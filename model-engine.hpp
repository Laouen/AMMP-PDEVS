#ifndef BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H
#define BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H

#include <iostream>
#include <string>
#include <utility>
#include <map>

#include "parser/parser.hpp"
#include "data-structures/types.hpp"

using namespace std;

template<class TIME>
using vm_t  = vector< shared_ptr< model<TIME> > >;
template<class TIME>
using vmp_t = vector< pair< shared_ptr< model<TIME> >, shared_ptr< model<TIME> > > >;
template<class TIME>
using mm_t  = map<string, shared_ptr< model<TIME> > >;
template<class TIME,class MSG>
using vcm_t = vector< shared_ptr< flattened_coupled<TIME, MSG> > >;
template<class TIME,class MSG>
using cmm_t = map<string, shared_ptr< flattened_coupled<TIME, MSG> > >;

template<class TIME,class MSG>
class ModelEngine {

public:
	bool   _model_ready, _comment_mode;
	long double _cell_weight;
	string _e, _c, _p, _b;
	Parser_t _doc;

  // Parser information
  map<string, string>               _compartments;
  map<string, map<string, string> > _species;
  map<string, enzyme_parameter_t >  _reactions;
  enzyme_parameter_t                _biomass_info;

	// Addresses
	shared_ptr< map<string, Address_t> > _species_addresses;
  shared_ptr< map<string, Address_t> > _enzyme_addresses;

	// Models ptr
	shared_ptr<flattened_coupled<TIME, MSG>>   _biomass_model;
	map< string, cmm_t<TIME, MSG>> 	           _enzyme_models;
	cmm_t<TIME, MSG> 					                 _compartment_models;
	cmm_t<TIME, MSG>                           _enzyme_set_models;
	vcm_t<TIME, MSG>                           _organelle_models;
	shared_ptr<flattened_coupled<TIME, MSG>>   _extra_cellular_model;
	shared_ptr<flattened_coupled<TIME, MSG>>   _periplasm_model;
	shared_ptr<flattened_coupled<TIME, MSG>>   _cytoplasm_model;
	shared_ptr<flattened_coupled<TIME, MSG>>   _cell_model;
	shared_ptr< model<TIME> >                  _input;
	shared_ptr< model<TIME> >                  _show_request;

  // Constructors
	explicit ModelEngine(const long double cw, const char* fl, const string e, const string c, const string p, const string b, const Integer_t n) 
  : _model_ready(false), _comment_mode(false), _cell_weight(cw), _doc(fl, b, cw, n), _e(e), _c(c), _p(p), _b(b) {

    // Parsing the SBML file
    _doc.loadFile();
    _compartments  = _doc.getCompartments();
    _species       = _doc.getSpeciesByCompartment();
    _reactions     = _doc.getReactions();
    _biomass_info  = _doc.getBiomass();

    // Initializing the shared ptr for the addresses
    _species_addresses = make_shared< map<string, Address_t> >();
    _enzyme_addresses  = make_shared< map<string, Address_t> >();
  }

  void createEnzymeAddresses() {
    if (this->_comment_mode) cout << "[Model engine] Creating the enzyme addresses." << endl;

    Address_t na;

    for (map<string, map<string, string> >::const_iterator i = _species.cbegin(); i != _species.end(); ++i) {
      
      na.clear();
      na.push_back(i->first);
      na.push_back(i->first + "_s");
      for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
        _species_addresses->insert({j->first, na});
      }
    }
  }

  void addCompartment(string c) {
    if (_comment_mode) cout << "[Model engine] Adding new compartment: " << c << endl;
    _enzyme_models.insert({c + "_i", {}});
    _enzyme_models.insert({c + "_m", {}});
  }

  vector<string> getOrganelleIds() const {
    vector<string> result;

    for (map<string, string>::const_iterator i = _compartments.begin(); i != _compartments.end(); ++i) {
      if (isNotSpecial(i->first)) result.push_back(i->first);
    }

    return result;
  }

  void addEnzymeModel(string id, TIME it, TIME r, Integer_t a) {
    string place, sub_place;
    enzyme_parameter_t params = _reactions[id];

    place = getPlace(params);

    auto reactio = make_atomic_ptr< 
      reaction<TIME, MSG>, 
      const string, 
      const shared_ptr< map<string, Address_t> >, 
      const bool, 
      const TIME, 
      const SetOfMolecules_t&, 
      const SetOfMolecules_t&, 
      const Integer_t, 
      const TIME >(
        id, 
        _species_addresses, 
        params.reversible, 
        r, 
        params.reactants_sctry, 
        params.products_sctry, 
        a, 
        it
      );
  }


  /******************* helpers *************************/

  bool isNotSpecial(string id) const {
    return (id != _e) && (id != _p) && (id != _c);
  }

  string compartmentsOf(string id) {

    // TODO 
  }

  string getPlace(const enzyme_parameter_t& e) {

  string result;
  int amount_compartments;
  string curr_space;
  map<string, int> compartments;

  for (map<string, string>::const_iterator i = _compartments.cbegin(); i != _compartments.cend(); ++i) {
    compartments[i->first] = 0;
  }

  for (SetOfMolecules_t::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    compartments[compartmentsOf(jt->first)] += 1;
  }

  for (SetOfMolecules_t::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {  
    compartments[compartmentsOf(jt->first)] += 1;
  }


  ac = 0;
  for (map<string, int>::iterator i = compartments.begin(); i != compartments.end(); ++i) {   
    if (i->second > 0) ++ac;
  }

  switch(ac) {
    case 1:

      result = compartmentOfReactants(compartments) + "_i";
      break;
    case 2:

      if (thereAreOrganelleInvolved(compartments, sp)) {
        result = organelleOfReactants(compartments, sp) + "_m";
      
      } else if (compartments[_e] > 0){

        result = _p + "_um";
      } else if (compartments[_c] > 0) {

        result = _p + "_lm";
      } else if ((compartments[_e] > 0) && (compartments[_c] > 0)) {

        result = _p + "_tm";
      }
      break;
    case 3:

      result = _p + "_tm";
      break;
  }

  return result;
}
};

#endif // BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H