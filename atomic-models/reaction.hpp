#ifndef BOOST_SIMULATION_PDEVS_REACTION_H
#define BOOST_SIMULATION_PDEVS_REACTION_H
#include <string>
#include <utility>
#include <vector>
#include <map>

// boost simalator include
#include <boost/simulation/pdevs/atomic.hpp>

using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/


/******** vector types ************/
typedef map<string, int> setOfMolecules;
typedef map<string, pair<string, int> > stoichiometryDef;



template<class TIME, class MSG>
class reaction : public atomic<TIME, MSG>
{

private:
  // enzyme information
  string            _name;
  bool              _reversible;
  TIME              _rate;
  stoichiometryDef  _stoichiometry;
  // elements bound
  setOfMolecules    _reactants;
  setOfMolecules    _products;
  setOfMolecules    _species_rejected;
  // function to control which elements stay and which leave
  bool              (*_randomFunction)();
  // to know what it must to do next
  int               _is_reacting;
  bool              _rejecting_species;
  // time information
  TIME              _next_internal;
  TIME              _interval_time;


public:

  explicit reaction(
    const string& other_name,
    const bool& other_reversible,
    const TIME& other_rate,
    const TIME& other_interval_time,
    const stoichiometryDef& other_stoichiometry,
    decltype(_stoichiometry) other_randomFunction
  ) noexcept :
  _name(other_name),
  _reversible(other_reversible),
  _rate(other_rate),
  _reactants(),
  _products(),
  _stoichiometry(other_stoichiometry),
  _randomFunction(other_randomFunction),
  _is_reacting(-1),
  _rejecting_species(false),
  _next_internal(other_interval_time),
  _interval_time(other_interval_time) {

    _species_rejected.clear();

    for (stoichiometryDef::const_iterator it = _stoichiometry.cbegin(); it != _stoichiometry.cend(); ++it) {
      if (it->second.first == "reactant") 
        _reactants[it->first] = 0;
      else if (it->second.first == "product") 
        _products[it->first]  = 0;
    }
  }

  void internal() noexcept { 
    if (_is_reacting != -1) {

      resetToZero(_reactants);
      resetToZero(_products);
      _is_reacting    = -1;
      _next_internal  = _interval_time;
    } else if (_rejecting_species){

      _species_rejected.clear();
    }
  }

  TIME advance() const noexcept {
    if (_is_reacting != -1)
      return _rate;
    else if (_rejecting_species)
      return TIME(0);
    else
      return _next_internal;
  }

  vector<MSG> out() const noexcept {
    
    if (_is_reacting == 0) {
      return gather(_products);
    } else if (_is_reacting == 1) {
      return gather(_reactants);
    } else if (_rejecting_species) {
      return gather(_species_rejected);
    }
  }

 // falta manejar las moleculas que se son rechazadas
  void external(const vector<MSG>& mb, const TIME& t) noexcept {
    int rejected_amount;
    MSG rejected_specie;
    if (_is_reacting == -1) {

      deleteLeavingSpecies(_reactants);
      deleteLeavingSpecies(_products);

      for (vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {

        rejected_amount = min(bindNeededPart(_reactants, _stoichiometry, it), bindNeededPart(_products, _stoichiometry, it));
        
        if (rejected_amount > 0) {
          _species_rejected.push_back(make_pair(it->specie, rejected_amount));
        }
      }
    } else {

      for (vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {

        _species_rejected.push_back(make_pair(it->specie, it->amount));
      }
    }

    _is_reacting        = readyToReact(_stoichiometry, _reactants, _products);
    _rejecting_species  = (_species_rejected.size() > 0);
    _next_internal      = _interval_time - t;
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

    internal();
    external(mb, t);
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  void deleteLeavingSpecies(setOfMolecules& molecules) {
    int totalAmount;   

    for (setOfMolecules::iterator it = molecules.begin(); it != molecules.end(); ++it) {
      
      totalAmount = it->second;
      for (int i = 0; i < totalAmount; ++i) {       
        if (_randomFunction()) --(it->second);
      }
    }
  }

  int readyToReact(const stoichiometryDef& stcmtry, const setOfMolecules& rectnts, const setOfMolecules& prdts) {
    int   result;
    bool  reactantsFull, reactantsEmpty, productsFull, productsEmpty;
    
    if (_reversible){

      reactantsFull   = isFull(rectnts, stcmtry);
      reactantsEmpty  = isEmpty(rectnts, stcmtry);
      productsFull    = isFull(prdts, stcmtry);
      productsEmpty   = isEmpty(prdts, stcmtry);
      
      if (productsEmpty && reactantsFull) {
        result = 0;
      } else if (reactantsEmpty && productsFull) {
        result = 1
      } else {
        result = -1
      }
    } else if (isFull(rectnts, stcmtry)) {
      result = 0;
    } else {
      result = -1;
    }

    return result;
  }

  bool isFull(const setOfMolecules& molecules, const stoichiometryDef& stcmtry) {
    bool result = true;

    for (setOfMolecules::const_iterator it = molecules.cbegin(); it != molecules.cend(); ++it) {
      if(stcmtry.at(it->first).second.second > it->second) {
        result = false;
        break;
      }
    }
    return result;
  }

  void resetToZero(setOfMolecules& molecules) {
    for (setOfMolecules::iterator it = molecules.begin(); it != molecules.end(); ++it) {
      it->second = 0;
    }
  }

  int bindNeededPart(setOfMolecules& molecules, const stoichiometryDef& stcmtry, const MSG* element) {
    int free_space;
    int result = element->amount;


    if(molecules.find(element->specie) != molecules.end()) {
          
      free_space = stcmtry.at(element->specie).second.second - molecules.at(element->specie);
      if (element->amount <= free_space){
        molecules.at(element->specie) += element->amount;
        result = 0; 
      } else {
        molecules.at(element->specie) += free_space; 
        result = element->amount - free_space;
      }
    }
    return result; 
  }

  vector<MSG> gather(const setOfMolecules& molecules) {
    MSG current_message;
    vector<MSG> result;

    for (setOfMolecules::const_iterator it = molecules.cbegin(); it != molecules.cend(); ++it) {
      current_message.clear();
      current_message.specie = it->first;
      current_message.amount = it->second;
      result.push_back(current_message);
    }

    return result;
  }

};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
