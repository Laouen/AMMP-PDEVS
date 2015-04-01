#ifndef BOOST_SIMULATION_PDEVS_REACTION_H
#define BOOST_SIMULATION_PDEVS_REACTION_H
#include <string>
#include <vector>
#include <list>
#include <map>
#include <utility> // pair
#include <algorithm> // min, max, lowe_bound
#include <memory> // shared_ptr
#include <random> // rand
#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include
#include "../data-structures/types.hpp" // SetOfMolecules, Task, Rstate, Way, TaskQueue



using namespace boost::simulation::pdevs;
using namespace std;


template<class TIME, class MSG>
class reaction : public atomic<TIME, MSG>
{

private:
  // enzyme information
  string          _name;
  bool            _reversible;
  TIME            _rate;
  SetOfMolecules  _products_sctry;
  SetOfMolecules  _reactants_sctry;
  int             _amount;
  TIME            _interval_time;
  // elements bound
  SetOfMolecules  _reactants;
  SetOfMolecules  _products;
  TaskQueue<TIME> _tasks;


public:

  explicit reaction(
    const string&           other_name,
    const bool&             other_reversible,
    const TIME&             other_rate,
    const SetOfMolecules&   other_products_sctry,
    const SetOfMolecules&   other_reactants_sctry,
    const int               other_amount,
    const TIME&             other_interval_time
  ) noexcept :
  _name(other_name),
  _reversible(other_reversible),
  _rate(other_rate),
  _products_sctry(other_products_sctry),
  _reactants_sctry(other_reactants_sctry),
  _amount(other_amount),
  _interval_time(other_interval_time) {

    for (SetOfMolecules::const_iterator it = _reactants_sctry.cbegin(); it != _reactants_sctry.cend(); ++it) {
      _reactants[it->first] = 0;
    }

    for (SetOfMolecules::const_iterator it = _products_sctry.cbegin(); it != _products_sctry.cend(); ++it) {
      _products[it->first] = 0;
    }
  }

  void internal() noexcept {
    int amount_leaving;
    Task<TIME> to_reject;
    bool already_selected_them  = false;
    TIME current_time           = _tasks.front().time_left;

    while(_tasks.front().time_left <= current_time) {

      if ((_tasks.front().task_kind == SELECTING) && !already_selected_them) {
        already_selected_them = true;

        for (SetOfMolecules::iterator it = _reactants.begin(); it != _reactants.end(); ++it) {

          amount_leaving  = rand() % (it->second + 1);
          it->second      -= amount_leaving;
          addRejected(to_reject.rejected, it->first, amount_leaving);
        }

        for (SetOfMolecules::iterator it = _products.begin(); it != _products.end(); ++it) {

          amount_leaving  = rand() % (it->second + 1);
          it->second      -= amount_leaving;
          addRejected(to_reject.rejected, it->first, amount_leaving);
        }

        if (to_reject.rejected.size() > 0) {

          to_reject.time_left = TIME(0);
          to_reject.task_kind = REJECTING;
          this->innsertTask(to_reject);
        }
      } else if (_tasks.front().task_kind == REACTING) {

        _amount += _tasks.front().reaction.second;
      }

      _tasks.pop_front();
    }
  }

  TIME advance() const noexcept {

    return _tasks.front().time_left;
  }

  vector<MSG> out() const noexcept {
    MSG new_message;
    SetOfMolecules* curr_sctry;

    vector<MSG> result = {};
    TIME current_time  = _tasks.front().time_left;

    for (typename Task<TIME>::const_iterator it = _tasks.cbegin(); it != _tasks.cend(); ++it) {
      
      if (it->time_left > current_time) break;

      if (it->task_kind == REJECTING) {
        
        for (SetOfMolecules::const_iterator jt = it->rejected.cbegin(); jt != it->rejected.cend(); ++jt) {
          new_message.clear();
          new_message.specie = jt->first;
          new_message.amount = jt->second;
          result.push_back(new_message); 
        }
      } else if (it->task_kind == REACTING) {

        if (it->reaction.first == RTP)  curr_sctry = &_products_sctry;
        else                            curr_sctry = &_products_sctry;

        for (SetOfMolecules::const_iterator jt = curr_sctry->cbegin(); jt != curr_sctry->cend(); ++jt) {
          new_message.clear();
          new_message.specie = jt->first;
          new_message.amount = jt->second * it->reaction.second;
          result.push_back(new_message); 
        }
      }
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {
    int reactant_ready, product_ready, free_of_reactants, free_of_products, intersection;
    Task<TIME> to_reject, rtp, ptr;

    // inserting new metaboolits and rejecting the surplus
    this->bindMetabolitsAndDropSurplus(mb, to_reject);

    // looking for new reactions
    reactant_ready  = this->totalReadyFor(RTP);
    product_ready   = this->totalReadyFor(PTR);
    intersection    = rand() % min(reactant_ready, product_ready); // tengo que mejorar este random
    reactant_ready  -= intersection;
    product_ready   -= intersection;

    if (to_reject.rejected.size() > 0) {

      to_reject.time_left = TIME(0);
      to_reject.task_kind = REJECTING;
      this->innsertTask(to_reject);
    }

    if (reactant_ready > 0) {
      rtp.time_left = _rate;
      rtp.task_kind = REACTING;
      rtp.reaction  = make_pair(RTP, reactant_ready);
      this->innsertTask(ptr);
    }

    if (product_ready > 0) {
      ptr.time_left = _rate;
      ptr.task_kind = REACTING;
      ptr.reaction  = make_pair(PTR, reactant_ready);
      this->innsertTask(ptr);
    }

    this->deleteUsedMetabolics(reactant_ready, product_ready);
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
    internal();
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  void bindMetabolitsAndDropSurplus(const vector<MSG>& mb, Task<TIME>& tr) {
    int free_space, metabolites_taken_r, metabolites_taken_p;
    bool is_reactant, is_product;

    for (typename vector<MSG>::cont_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      
      metabolites_taken_r = metabolites_taken_p = 0;
      is_reactant = _reactants_sctry.find(it->specie) != _reactants_sctry.end();
      is_product  = _products_sctry.find(it->specie) != _products_sctry.end(); 

      // binding the allowed amount of metabolits
      if (is_reactant){

        free_space          = (_amount * _reactants_sctry.at(it->specie)) - _reactants.at(it->specie);
        metabolites_taken_r = min(it->amount, free_space);

        if (is_product) {
          metabolites_taken_r = rand() % (metabolites_taken_r + 1);
        }

        _reactants.at(it->specie) += metabolites_taken_r;
      } 

      if (is_product && _reversible) {
      
        free_space          = (_amount * _products_sctry.at(it->specie)) - _products.at(it->specie);
        metabolites_taken_p = min(it->amount, free_space);

        if (is_reactant) {
          metabolites_taken_p = min(metabolites_taken_p, it->amount - metabolites_taken_r);
        }

        _products.at(it->specie) += metabolites_taken_p;
      }
      
      // placing the surplus metabolites in the rejected list
      addRejected(tr.rejected, it->specie, it->amount - (metabolites_taken_p + metabolites_taken_r));
    }
  }

  void addRejected(SetOfMolecules& t, string n, int a){

    if (a > 0) {
      if (t.find(n) != t.end()) {
        t.at(n) += a;
      } else {
        t[n] = a;
      }
    }
  }

  // usando la stoichiometria restar a los reactantes y productos la cantidades utilizadas
  // _stoichiometry, _reactants, _products
  void deleteUsedMetabolics(int r, int p) {

    for (SetOfMolecules::iterator it = _reactants.begin(); it != _reactants.end(); ++it) {
      it->second -= r * _reactants_sctry.at(it->first); 
    } 

    for (SetOfMolecules::iterator it = _products.begin(); it != _products.end(); ++it) {
      it->second -= r * _products_sctry.at(it->first); 
    }

    _amount -= r + p;
  }

  // inserta la tarea t a la cola de tareas de manera de mantener el orden de menor a mayor time_left.
  void innsertTask(const Task<TIME>& t) {
    typename TaskQueue<TIME>::iterator it = lower_bound(_tasks.begin(), _tasks.end(), t);
    _tasks.insert(it, t);
  }

  // usando la stoichiometria y el espacio libre (segundo parametro) devuelve la cantidad total de elemento
  // listo a reaccionar.
  int totalReadyFor(Way d) const {
    int fr;

    SetOfMolecules *curr_metabolics, *curr_sctry;   
    if (d == RTP) {
      curr_metabolics = &_reactants;
      curr_sctry      = &_reactants_sctry;
      fr              = this->freeOf(_products);
    } else {
      curr_metabolics = &_products;
      curr_sctry      = &_products_sctry;
      fr              = this->freeOf(_reactants);
    }

    int m = 0;
    for (SetOfMolecules::const_iterator it = curr_metabolics->cbegin(); it != curr_metabolics->cend(); ++it)
      m = min(m, (it->second / curr_sctry->at(it->first)));

    return min(m,fr);
  }

  //usando la stoichiometry y el _amount calcula cuanto es la cantidad de enzymas que estan libres de
  // ese set de elementos. mira la maximo elemento que aparece y cuantas enzymas este ocupa.
  int freeOf(Way d) const {
    
    SetOfMolecules *curr_metabolics, *curr_sctry;   
    if (d == RTP) {
      curr_metabolics = &_reactants;
      curr_sctry      = &_reactants_sctry;
    } else {
      curr_metabolics = &_products;
      curr_sctry      = &_products_sctry;
    }

    int m = 0;
    for (SetOfMolecules::const_iterator it = curr_metabolics->cbegin(); it != curr_metabolics->cend(); ++it)
      m = max(m, (it->second / curr_sctry->at(it->first)) + 1); 

    return m;
  }


  ostream& show(ostream& os, const SetOfMolecules& to) const {
  
    os << "[";

    SetOfMolecules::const_iterator it = to.cbegin();
    while ( it != to.cend()) {
      os << "(" << it->first << "," << it->second << ")";
      ++it;
      if (it != to.cend()) os << ",";
    }
    os << "]" << endl;
    return os;
  }

};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
