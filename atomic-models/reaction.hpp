#ifndef BOOST_SIMULATION_PDEVS_REACTION_H
#define BOOST_SIMULATION_PDEVS_REACTION_H
#include <string>
#include <utility>
#include <vector>
#include <boost/simulation/pdevs/atomic.hpp>
#include "../data-structures/reaction_input.hpp"

using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** vector types ************/
typedef vector<molecule> vectorOfMolecules;


template<class TIME, class MSG>
class reaction : public atomic<TIME, MSG>
{

private:
  string _name;
  TIME _next;
  TIME _rate;
  vectorOfMolecules _reactants;
  vectorOfMolecules _products;
  vectorOfMolecules _enzymes;
  double _counter;

public:

  explicit reaction(reaction_input data) noexcept :
  _name(data.name),
  _next(atomic<TIME, MSG>::infinity),
  _rate(data.rate),
  _reactants(data.reactants),
  _products(data.products),
  _enzymes(data.enzymes),
  _counter(0) {}

  void internal() noexcept { 

    _counter  = 0;
    _next     = atomic<TIME, MSG>::infinity;
  }

  TIME advance() const noexcept {

    return _next;
  }

  vector<MSG> out() const noexcept {

    vector<MSG> result;

    for (vectorOfMolecules::const_iterator i = _products.cbegin(); i != _products.cend(); ++i){
      result.push_back( boost::any(make_pair(i->first, _counter)) );
    }
    return result; 
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    molecule current_molecule;

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      current_molecule = boost::any_cast<molecule>(*it);

      for (vectorOfMolecules::iterator i = _reactants.begin(); i != _reactants.end(); ++i){
        if (current_molecule.first == i->first ) {
          i->second += current_molecule.second;
        }
      }

      for (vectorOfMolecules::iterator i = _enzymes.begin(); i != _enzymes.end(); ++i){
        if (current_molecule.first == i->first ) {
          i->second = 1;
        }
      }
    }

    // if there is molecule enough a reaction start taking 0 time to finish.
    bool enough_molecule  = enoughMolecule(_reactants);
    bool reactant_present = allEnzymesPresent(_enzymes);

    if (enough_molecule && reactant_present) {

      _counter = minoritySpecies(_reactants);

      for (vectorOfMolecules::iterator i = _reactants.begin(); i != _reactants.end(); ++i){
        i->second -= _counter;
      }

      _next = _rate;
    } else {

      _next = atomic<TIME, MSG>::infinity;
    }
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {

      external(mb, t);
      internal();
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  bool enoughMolecule(vectorOfMolecules sp) {

    bool result = true;
    for (vectorOfMolecules::iterator i = sp.begin(); i != sp.end(); ++i){
      if (i->second <= 0) {

        result = false;
        break;
      }
    }
    return result;
  }

  bool allEnzymesPresent(vectorOfMolecules rt) {

    bool result = true;
    for (vectorOfMolecules::iterator i = rt.begin(); i != rt.end(); ++i){
      if (i->second == 0) {
        
        result = false;
        break;
      }
    }
    return result;
  }

  double minoritySpecies(vectorOfMolecules sp) {

    int result = sp.front().second;
    for (vectorOfMolecules::iterator i = sp.begin(); i != sp.end(); ++i){
      if (i->second < result)
        result = i->second;
    }
    return result;
  }
};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
