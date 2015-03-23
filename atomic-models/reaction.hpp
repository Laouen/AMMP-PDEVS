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
  string            _name;
  bool              _reversible;
  TIME              _rate;
  setOfMolecules    _reactants;
  setOfMolecules    _products;
  stoichiometryDef  _stoichiometry;
  bool              (*_randomFunction)();
  bool              _is_reacting;
  TIME              _next_internal;
  TIME              _interval_time;


public:

  explicit reaction(string other_name, bool other_reversible, TIME other_rate, TIME other_interval_time stoichiometryDef other_stoichiometry, decltype(_stoichiometry) other_randomFunction) noexcept :
  _name(other_name),
  _reversible(other_reversible),
  _rate(other_rate),
  _reactants(),
  _products(),
  _stoichiometry(other_stoichiometry),
  _randomFunction(other_randomFunction),
  _isReacting(false),
  _next_internal(other_interval_time),
  _interval_time(other_interval_time) {

    for (stoichiometryDef::const_iterator it = _stoichiometry.cbegin(); it != _stoichiometry.cend(); ++it) {
      if (it->second.first == "reactant") 
        _reactants[it->first] = 0;
      else if (it->second.first == "product") 
        _products[it->first]  = 0;
    }
  }

  void internal() noexcept { 

  }

  TIME advance() const noexcept {
    if (_isReacting)
      return _rate;
    else
      return _next_internal;
  }

  vector<MSG> out() const noexcept {

  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {
    if (!_isReacting){

      deleteLeavingSpecies(_reactants);
      deleteLeavingSpecies(_products);

      for (vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {

        if(_reactants[it->specie] != _reactants.end()) {
          
          _reactants[it->specie] = max(_stoichiometry[it->specie].second, (_reactants[it->specie] + it->amount));
        
        } else if(_products[it->specie] != _products.end()) {

          _products[it->specie] = max(_stoichiometry[it->specie].second, (_products[it->specie] + it->amount));

        }
      }

      _is_reacting = readyForReact(_stoichiometry, _reactants, _products);
    }

    _next_internal = _interval_time - t;
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  void deleteLeavingSpecies(setOfMolecules) {
    cout << "function deleteLeavingSpecies not implemented" << endl;
  }

  bool readyForReact(stoichiometryDef stcmtry, setOfMolecules rectnts, setOfMolecules prdts) {
    
    cout << "function readyForReact not implemented" << endl;
    return true;
  }

};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
