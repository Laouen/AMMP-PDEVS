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
#include <limits> // numeric_limits

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // SetOfMolecules, Task, Rstate, Way, TaskQueue
#include "../data-structures/randomNumbers.hpp" // RealRandom

using namespace boost::simulation::pdevs;
using namespace std;


template<class TIME, class MSG>
class reaction : public atomic<TIME, MSG>
{

private:
  // enzyme information
  string                  _name;
  bool                    _reversible;
  TIME                    _rate;
  SetOfMolecules          _reactants_sctry;
  SetOfMolecules          _products_sctry;
  Integer                 _amount;
  TIME                    _interval_time;
  // elements bound
  SetOfMolecules          _reactants;
  SetOfMolecules          _products;
  TaskQueue<TIME>         _tasks;
  // used for uniform random number
  IntegerRandom<Integer>  _distribution;


public:

  explicit reaction(
    const string&             other_name,
    const bool&               other_reversible,
    const TIME&               other_rate,
    const SetOfMolecules&     other_reactants_sctry,
    const SetOfMolecules&     other_products_sctry,
    const Integer             other_amount,
    const TIME&               other_interval_time
  ) noexcept :
  _name(other_name),
  _reversible(other_reversible),
  _rate(other_rate),
  _reactants_sctry(other_reactants_sctry),
  _products_sctry(other_products_sctry),
  _amount(other_amount),
  _interval_time(other_interval_time) {

    random_device rd;
    _distribution.seed(rd());

    for (SetOfMolecules::const_iterator it = _reactants_sctry.cbegin(); it != _reactants_sctry.cend(); ++it) {
      _reactants[it->first] = 0;
    }

    for (SetOfMolecules::const_iterator it = _products_sctry.cbegin(); it != _products_sctry.cend(); ++it) {
      _products[it->first] = 0;
    }
  }

  void internal() noexcept {

    Task<TIME> to_reject;
    bool already_selected = false;

    // Updating time left
    this->updateTaskTimeLefts(_tasks.front().time_left);

    // Processing all the tasks with time left == 0 (happening now)
    for (typename TaskQueue<TIME>::iterator it = _tasks.begin(); it->time_left == 0; it = _tasks.erase(it)) {

      if ((it->task_kind == RState::SELECTING) && !already_selected) {

        this->selectFrom(_reactants, to_reject);
        
        this->selectFrom(_products, to_reject);

        this->lookForNewReactions();

        already_selected = true;
      
      } else if (it->task_kind == RState::REACTING) {

        _amount += it->reaction.second;
      } 
    }

    // add leaving metabolites to the tasks
    if (to_reject.rejected.size() > 0) {

      to_reject.time_left = TIME(0);
      to_reject.task_kind = RState::REJECTING;
      this->innsertTask(to_reject);
    }

    // if there is more metabolites set a new selection tasks in interval time
    this->setNextSelection();
  }

  TIME advance() const noexcept {
    if (_tasks.size() > 0) return _tasks.front().time_left;
    else                   return atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {
    MSG new_message;
    const SetOfMolecules* curr_sctry;

    vector<MSG> result = {};
    TIME current_time  = _tasks.front().time_left;

    for (typename TaskQueue<TIME>::const_iterator it = _tasks.cbegin(); it->time_left == current_time; ++it) {

      if (it->task_kind == RState::REJECTING) {
        
        for (SetOfMolecules::const_iterator jt = it->rejected.cbegin(); jt != it->rejected.cend(); ++jt) {
          new_message.clear();
          new_message.specie = jt->first;
          new_message.amount = jt->second;
          result.push_back(new_message); 
        }
      } else if (it->task_kind == RState::REACTING) {

        if (it->reaction.first == Way::RTP)      curr_sctry = &_products_sctry;
        else if (it->reaction.first == Way::PTR) curr_sctry = &_reactants_sctry;

        for (SetOfMolecules::const_iterator jt = curr_sctry->cbegin(); jt != curr_sctry->cend(); ++jt) {
          new_message.clear();
          new_message.specie = jt->first;
          new_message.amount = it->reaction.second * jt->second;
          result.push_back(new_message); 
        }
      }
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    Task<TIME> to_reject, new_selection;
    
    // Updating time left
    this->updateTaskTimeLefts(t);

    // inserting new metaboolits and rejecting the surplus
    this->bindMetabolitsAndDropSurplus(mb, to_reject);

    // looking for new reactions
    this->lookForNewReactions();

    // add rejecting surplus to the tasks
    if (to_reject.rejected.size() > 0) {

      to_reject.time_left = TIME(0);
      to_reject.task_kind = RState::REJECTING;
      this->innsertTask(to_reject);
    }

    // if there is more metabolites set a new selection tasks in interval time
    this->setNextSelection();
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

    internal();
    external(mb, TIME(0));
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  // It take the needed number of each specie in mb and rejected (by calling addRejected) the not needed part.
  void bindMetabolitsAndDropSurplus(const vector<MSG>& mb, Task<TIME>& tr) {
    Integer free_space, metabolites_taken_r, metabolites_taken_p, amount_for_r, amount_for_p;
    bool is_reactant, is_product;

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      
      metabolites_taken_r = 0;
      metabolites_taken_p = 0;
      is_reactant         = _reactants.find(it->specie) != _reactants.end();
      is_product          = _products.find(it->specie) != _products.end();

      // dividing the metabolit between reactants and products
      if (is_reactant && is_product && _reversible) {
        amount_for_r = _distribution.drawNumber(0, it->amount);
        amount_for_p = it->amount - amount_for_r;
      } else {
        amount_for_r = it->amount;
        amount_for_p = it->amount;
      }

      // binding the allowed amount of metabolits
      if (is_reactant){

        free_space                =   (_amount * _reactants_sctry.at(it->specie)) - _reactants.at(it->specie);
        metabolites_taken_r       =   min(amount_for_r, free_space);
        _reactants.at(it->specie) +=  metabolites_taken_r;
      } 

      if (is_product && _reversible) {
      
        free_space                =   (_amount * _products_sctry.at(it->specie)) - _products.at(it->specie);
        metabolites_taken_p       =   min(amount_for_p, free_space);
        _products.at(it->specie)  +=  metabolites_taken_p;
      }
      
      // placing the surplus metabolites in the rejected list
      addRejected(tr.rejected, it->specie, it->amount - (metabolites_taken_p + metabolites_taken_r));
    }
  }

  // Decrease the time left of all the current tasks in _tasks by the parameter t.
  void updateTaskTimeLefts(TIME t){

    for (typename TaskQueue<TIME>::iterator it = _tasks.begin(); it != _tasks.end(); ++it) {
      it->time_left -= t;
    }
  }

  void lookForNewReactions() {
    Integer reactant_ready, product_ready, intersection_range, intersection;
    Task<TIME> rtp, ptr;

    reactant_ready      = this->totalReadyFor(Way::RTP);
    product_ready       = this->totalReadyFor(Way::PTR);
    intersection_range  = min(reactant_ready, product_ready);
    intersection        =  _distribution.drawNumber(0, intersection_range);
    reactant_ready      -= intersection;
    product_ready       -= intersection;

    this->deleteUsedMetabolics(reactant_ready, product_ready);

    if (reactant_ready > 0) {
      rtp.time_left = _rate;
      rtp.task_kind = RState::REACTING;
      rtp.reaction  = make_pair(Way::RTP, reactant_ready);
      this->innsertTask(rtp);
    }

    if (product_ready > 0) {
      ptr.time_left = _rate;
      ptr.task_kind = RState::REACTING;
      ptr.reaction  = make_pair(Way::PTR, product_ready);
      this->innsertTask(ptr);
    }
  }

  void setNextSelection() {
    Task<TIME> new_selection;
    
    if ( (!isClean(_reactants) || !isClean(_products)) && !this->thereIsNextSelection() ) {

      new_selection.time_left = _interval_time;
      new_selection.task_kind = RState::SELECTING;
      this->innsertTask(new_selection);
    }
  }

  void selectFrom(SetOfMolecules& m, Task<TIME>& t) {
    Integer amount_leaving;

    for (SetOfMolecules::iterator it = m.begin(); it != m.end(); ++it) {

      amount_leaving  = _distribution.drawNumber(0, it->second);
      it->second      -= amount_leaving;
      addRejected(t.rejected, it->first, amount_leaving);
    }
  }

  // Add the rejected amount of the species specified in n by inserting/increasing the amount a in the set Of Molecule (garbage set) t.
  void addRejected(SetOfMolecules& t, string n, Integer a){

    if (a > 0) {
      if (t.find(n) != t.end()) {
        t.at(n) += a;
      } else {
        t[n] = a;
      }
    }
  }

  // It decrease the number of reactants and products removing the specified number by the parameters r and p.
  void deleteUsedMetabolics(Integer r, Integer p) {

    for (SetOfMolecules::iterator it = _reactants.begin(); it != _reactants.end(); ++it) {
      it->second -= r * _reactants_sctry.at(it->first); 
    } 

    for (SetOfMolecules::iterator it = _products.begin(); it != _products.end(); ++it) {
      it->second -= p * _products_sctry.at(it->first); 
    }

    _amount -= r + p;
  }

  // Insert a copy of the parameter t in the TaskQueue in an ordered way.
  void innsertTask(const Task<TIME>& t) {

    typename TaskQueue<TIME>::iterator it = lower_bound(_tasks.begin(), _tasks.end(), t);
    _tasks.insert(it, t);
  }

  // It return the total number of enzymes that are ready to react in the way specified by the parameter d.
  Integer totalReadyFor(Way d) const {
    Integer fr;

    const SetOfMolecules *curr_metabolics; 
    const SetOfMolecules *curr_sctry;

    if (d == Way::RTP) {
      curr_metabolics = &_reactants;
      curr_sctry      = &_reactants_sctry;
      fr              = this->freeFor(Way::RTP);
    } else {
      curr_metabolics = &_products;
      curr_sctry      = &_products_sctry;
      fr              = this->freeFor(Way::PTR);
    }

    Integer m = numeric_limits<Integer>::max();
    for (SetOfMolecules::const_iterator it = curr_metabolics->cbegin(); it != curr_metabolics->cend(); ++it)
      m = min(m, (Integer)floor(it->second / curr_sctry->at(it->first)));

    return min(m,fr);
  }

  //usando la stoichiometry y el _amount calcula cuanto es la cantidad de enzymas que estan libres de
  // ese set de elementos. mira la maximo elemento que aparece y cuantas enzymas este ocupa.
  Integer freeFor(Way d) const {
    
    const SetOfMolecules *curr_metabolics;
    const SetOfMolecules *curr_sctry;

    if (d == Way::RTP) {
      curr_metabolics = &_products;
      curr_sctry      = &_products_sctry;
    } else {
      curr_metabolics = &_reactants;
      curr_sctry      = &_reactants_sctry;
    }

    Integer m = 0;
    for (SetOfMolecules::const_iterator it = curr_metabolics->cbegin(); it != curr_metabolics->cend(); ++it)
      m = max( m, (Integer)ceil(it->second / curr_sctry->at(it->first)) ); 

    return _amount - m;
  }

  bool isClean(const SetOfMolecules& t) {

    bool result = true;
    for (SetOfMolecules::const_iterator it = t.cbegin(); it != t.cend(); ++it) {
      if (it->second != 0) {
        result = false;
        break;
      }
    }

    return result;
  }

  bool thereIsNextSelection() {
    bool result = false;

    for (typename TaskQueue<TIME>::iterator it = _tasks.begin(); it != _tasks.end(); ++it) {
      if ((it->task_kind == RState::SELECTING) && (it->time_left <= _interval_time)) {
        result = true;
        break;
      }
    }

    return result;
  }

  /*********************************************/
  /************** Testing functions ************/
  /*********************************************/

  ostream& show(ostream& os, const SetOfMolecules& to) {

    os << "[";

    SetOfMolecules::const_iterator it = to.cbegin();
    while ( it != to.cend()) {
      os << "(" << it->first << "," << it->second << ")";
      ++it;
      if (it != to.cend()) os << ",";
    }
    os << "]";
    return os;
  }

  ostream& show(ostream& os, const Task<TIME>& to) {

    string kind, w;
    if (to.task_kind == RState::SELECTING)        kind = "RState::SELECTING";
    else if (to.task_kind == RState::REJECTING)   kind = "RState::REJECTING";
    else if (to.task_kind == RState::REACTING)    kind = "RState::REACTING";

    os << "Task Kind: " << kind << endl;
    os << "Time left: " << to.time_left << endl;

    if(to.task_kind == RState::REJECTING) {
      show(os, to.rejected);
    } else if (to.task_kind == RState::REACTING) {

      if (to.reaction.first == Way::RTP)       w = "Way::RTP";
      else if (to.reaction.first == Way::PTR)  w = "Way::PTR";

      os << "way: " << w << " amount: " << to.reaction.second;
    }

    return os;
  }

  ostream& show(ostream& os, const TaskQueue<TIME>& to) {

    os << "Current Tasks in the queue: ";
    for (typename TaskQueue<TIME>::const_iterator it = to.cbegin(); it != to.cend(); ++it) {
      os << endl << endl;
      show(cout, *it);
    }
    return os;
  }

  ostream& show(ostream& os) {

    os << "name: "          << _name                            << endl;
    os << "rate: "          << _rate                            << endl;
    os << "reversible: "    << (_reversible ? "true" : "false") << endl;
    os << "interval time: " << _interval_time                   << endl;
    os << "free enzymes: "  << _amount                          << endl;
    
    os << "react sctry: ";
    show(os, _reactants_sctry);
    os << endl;
    
    os << "prod sctry: ";
    show(os, _products_sctry);
    os << endl;
    
    os << "reactants: ";
    show(os, _reactants);
    os << endl;
    
    os << "products: ";
    show(os, _products);
    os << endl;
    
    os << "schedulled tasks: " << endl;
    show(os, _tasks);
    os << endl;
    
    return os;
  }

};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
