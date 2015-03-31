#ifndef BOOST_SIMULATION_PDEVS_REACTION_H
#define BOOST_SIMULATION_PDEVS_REACTION_H
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <map>
#include <random>

// boost simalator include
#include <boost/simulation/pdevs/atomic.hpp>

using namespace boost::simulation::pdevs;
using namespace std;

/******************************************/
/********** Type definations **************/
/******************************************/

enum RState { REJECTING, REACTING, SELECTING, IDLE };
enum Way { RTP, PTR };

using SetOfMolecules  = map<string, int>;
using Stctry          = map<string, pair<string, int> >;

template<class TIME>
struct Task {
  TIME            time_left;
  RState          task_kind;
  SetOfMolecules  rejected;  
  pair<Way, int>  reaction;
};

template<class TIME>
using Heap = vector< Task<TIME> >;

/******************************************/
/******** End type definations ************/
/******************************************/

template<class TIME, class MSG>
class reaction : public atomic<TIME, MSG>
{

private:
  // enzyme information
  string          _name;
  bool            _reversible;
  TIME            _rate;
  Stctry          _stoichiometry;
  int             _amount;
  // elements bound
  SetOfMolecules  _reactants;
  SetOfMolecules  _products;
  Heap<TIME>      _tasks;
  // time information
  TIME            _interval_time;


public:

  explicit reaction(
    const string&             other_name,
    const bool&               other_reversible,
    const TIME&               other_rate,
    const stoichiometryDef&   other_stoichiometry,
    const int                 other_amount,
    const TIME&               other_interval_time
  ) noexcept :
  _name(other_name),
  _reversible(other_reversible),
  _rate(other_rate),
  _stoichiometry(other_stoichiometry),
  _amount(other_amount),
  _interval_time(other_interval_time) {

    for (stoichiometryDef::const_iterator it = _stoichiometry.cbegin(); it != _stoichiometry.cend(); ++it) {
      if (it->second.first == "reactant") 
        _reactants[it->first] = 0;
      else if (it->second.first == "product") 
        _products[it->first]  = 0;
    }
  }

  void internal() noexcept {

    if (_s.rejecting_species){

      _species_rejected.clear();
      _s.rejecting_species = false;
    } else if (_s.reaction_in_process != -1) {

      resetToZero(_reactants);
      resetToZero(_products);
      _s.reaction_in_process = -1;
      _next_internal    = _interval_time;
    } else {

      removeLeavingSpecies(_reactants, _species_rejected);
      removeLeavingSpecies(_products, _species_rejected);
      _s.rejecting_species = (_species_rejected.size() > 0);
      _next_internal      = _interval_time;
    }
  }

  TIME advance() const noexcept {
    return _tasks.front().time_left;
  }

  vector<MSG> out() const noexcept {
    
    vector<MSG>         result, current_species;
    TIME                next_reaction;
    ReactionData<TIME>  current_reaction_data;

    switch(_state) { 
      case REJECTING:
        result = gather(_rejected_species);
        break;
      case REACTING:
        next_reaction = _reactions_in_queue.front().first
        while(_reactions_in_queue.front().first == next_reaction) {
          current_reaction_data = _reactions_in_queue.front();
          current_species       = createProduct(current_reaction_data.direction, current_reaction_data.enzyme_amount);
          concat(result, current_species);
        }
        break;
      default:
        result.clear();
        break;
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {
    int free_space, metabolite_taken, reactant_ready, product_ready, free_of_reactants, free_of_products;
    Task<TIME> to_reject, rtp, ptr;

    // inserting new metaboolits and rejecting the surplus
    for (vector<MSG>::cont_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
        
      // binding the allowed amount of metabolits
      free_space       = (_amount * _stoichiometry.at(it->specie).second) - _reactants.at(it->specie);
      metabolite_taken = min(it->amount, free_space);
      this->insertMetabolit(it->specie, metabolites_taken);
      
      // placing the surplus metabolites in the rejected list
      addRejected(to_reject.rejected, it->specie, it->amount - metabolite_taken);
    }

    // looking for new reactions
    free_of_products    = this->freeOf(_products);
    free_of_reactants   = this->freeOf(_reactants);
    reactant_ready      = this->totalReadyFor(_reactants, free_of_products);
    product_ready       = this->totalReadyFor(_products, free_of_reactants);
    intersection        = rand() % min(reactant_ready, product_ready); // tengo que mejorar este random
    reactant_ready      -= intersection;
    product_ready       -= intersection;

    to_reject.time_left = TIME(0);
    to_reject.task_kind = REJECTING;
    this->innsertTask(to_reject);

    if (reactant_ready > 0) {
      rtp.time_left = _rate;
      rtp.task_kind = REACTING;
      rtp.reaction  = make_pair(RTP, reactant_ready);
      this->innsertTask(ptr);
    }

    if (product_ready > 0) {
      prt.time_left = _rate;
      prt.task_kind = REACTING;
      prt.reaction  = make_pair(PTR, reactant_ready);
      this->innsertTask(ptr);
    }

    this->deleteUsedMetabolics(reactant_ready, product_ready);
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

    if ((_s.reaction_in_process != -1) or _s.rejecting_species) {
      internal();
      external(mb, t);
    } else {
      external(mb, t);
    }
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  bool randomBool() {
    return ((rand() % 100) <= 50);
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

  }

  // inserta la tarea t a la cola de tareas de manera de mantener el orden de menor a mayor time_left.
  void innsertTask(const Task& t) {

  }

  // usando la stoichiometria y el espacio libre (segundo parametro) devuelve la cantidad total de elemento
  // listo a reaccionar.
  int totalReadyFor(const SetOfMolecules& elements, int free_space) const {

    return;
  }

  //usando la stoichiometry y el _amount calcula cuanto es la cantidad de enzymas que estan libres de
  // ese set de elementos. mira la maximo elemento que aparece y cuantas enzymas este ocupa.
  int freeOf(const SetOfMolecules& element) const {

  }

  // mirando si la especia y el amount agrega esa cantidad de elementos a los reactantes o productos segun corresponda
  // o sea, donde pertenesca la especie.
  void insertMetabolit(string specie, int amount) {

  }








  void removeLeavingSpecies(SetOfMolecules& molecules, SetOfMolecules& residues) {
    int totalAmount, totalToDrop;   

    for (SetOfMolecules::iterator it = molecules.begin(); it != molecules.end(); ++it) {
      totalToDrop = 0;
      totalAmount = it->second;
      for (int i = 0; i < totalAmount; ++i) {       
        if (_randomFunction()) ++totalToDrop;
      }
      it->second -= totalToDrop;

      if (residues.find(it->first) != residues.end()) residues.at(it->first) += totalToDrop;
      else residues[it->first] = totalToDrop;
    }
  }

  int readyToReact(const stoichiometryDef& stcmtry, const SetOfMolecules& rectnts, const SetOfMolecules& prdts, const int amount) const {
    int   result, ready_reactant_amount, ready_product_amount, interception;
    bool  reactantsFull, reactantsEmpty, productsFull, productsEmpty;

    if (_reversible) {
      ready_reactant_amount = readyReactant(rectnts, stcmtry, amount);
      ready_product_amount  = readyProduct(rectnts, stcmtry, amount);
      interception          = ready_product_amount - (amount - ready_reactant_amount);
      ready_reactant_amount -= interception;
      ready_product_amount  -= interception;
    }
    
    if (_reversible){

      reactantsFull   = isFull(rectnts, stcmtry);
      reactantsEmpty  = isEmpty(rectnts);
      productsFull    = isFull(prdts, stcmtry);
      productsEmpty   = isEmpty(prdts);
      
      if (productsEmpty && reactantsFull) {
        result = 0;
      } else if (reactantsEmpty && productsFull) {
        result = 1;
      } else {
        result = -1;
      }
    } else if (isFull(rectnts, stcmtry)) {
      result = 0;
    } else {
      result = -1;
    }

    return result;
  }

  bool isFull(const SetOfMolecules& molecules, const stoichiometryDef& stcmtry) const {
    bool result = true;

    for (SetOfMolecules::const_iterator it = molecules.cbegin(); it != molecules.cend(); ++it) {
      if((stcmtry.at(it->first)).second > it->second) {
        result = false;
        break;
      }
    }
    return result;
  }

  bool isEmpty(const SetOfMolecules& molecules) const {
    bool result = true;

    for (SetOfMolecules::const_iterator it = molecules.cbegin(); it != molecules.cend(); ++it) {
      if(it->second > 0) {
        result = false;
        break;
      }
    }
    return result;
  }

  void resetToZero(SetOfMolecules& molecules) const {
    for (SetOfMolecules::iterator it = molecules.begin(); it != molecules.end(); ++it) {
      it->second = 0;
    }
  }

  int bindNeededPart(SetOfMolecules& molecules, const stoichiometryDef& stcmtry, const MSG& element) {
    int free_space;
    int result = element.amount;


    if(molecules.find(element.specie) != molecules.end()) {
          
      free_space = stcmtry.at(element.specie).second - molecules.at(element.specie);
      if (element.amount <= free_space){
        molecules.at(element.specie) += element.amount;
        result = 0; 
      } else {
        molecules.at(element.specie) += free_space; 
        result = element.amount - free_space;
      }
    }
    return result; 
  }

  vector<MSG> gather(const SetOfMolecules& molecules) const {
    MSG current_message;
    vector<MSG> result;

    for (SetOfMolecules::const_iterator it = molecules.cbegin(); it != molecules.cend(); ++it) {
      current_message.clear();
      current_message.specie = it->first;
      current_message.amount = it->second;
      result.push_back(current_message);
    }

    return result;
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
