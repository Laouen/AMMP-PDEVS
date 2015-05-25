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

    // initializing the special compartment extra cellular, periplasm and cytoplasm
    _enzyme_models.insert({_e + "_i", {}});
    _enzyme_models.insert({_c + "_i", {}});
    _enzyme_models.insert({_p + "_i", {}});
    _enzyme_models.insert({_p + "_lm", {}});
    _enzyme_models.insert({_p + "_um", {}});
    _enzyme_models.insert({_p + "_tm", {}});
  }

  void createSpeciesAddresses() {
    if (this->_comment_mode) cout << "[Model engine] Creating the specie addresses." << endl;

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

  void addEnzymeModel(const string& id, TIME it, TIME r, Integer_t a) {
    assert(_reactions.find(id) != _reactions.end());
    
    if (_comment_mode) cout << "[Model engine] Adding new enzyme model: " << id << endl;

    string place, sub_place;
    enzyme_parameter_t params = _reactions[id];

    // reaction atomic model
    auto ereaction = make_atomic_ptr< 
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
    // filter atomic model
    auto efilter = make_atomic_ptr< filter<TIME, MSG>, const string>(id);
    
    // enzyme coupled model
    place = getPlace(params.reactants_sctry, params.products_sctry);
    _enzyme_models.at(place)[id] = make_shared<
      flattened_coupled<TIME, MSG>,
      vector<shared_ptr<model<TIME> > >,
      vector<shared_ptr<model<TIME> > >,
      vector<pair<shared_ptr<model<TIME>>, shared_ptr<model<TIME> > > >,
      vector<shared_ptr<model<TIME> > > 
    >(
      {efilter, ereaction}, 
      {efilter}, 
      {{efilter, ereaction}}, 
      {ereaction}
    );
  }

  void createEnzymeAddresses() {
    if (_comment_mode) cout << "[Model engine] Creating the enzyme addresses." << endl;
    Address_t na;
    string place;
    //string sub_place;

    _enzyme_addresses->clear();
    
    for (typename map< string, cmm_t<TIME, MSG>>::const_iterator i = _enzyme_models.cbegin(); i != _enzyme_models.cend(); ++i) {
      
      place       = i->first.substr(0, i->first.find_last_of("_"));
      //sub_place   = i->first.substr(i->first.find_last_of("_") + 1);
      
      na.clear();
      na.push_back(place);
      na.push_back(i->first);
      //if (sub_place == "i") na.push_back(place + "_s");

      for (typename cmm_t<TIME, MSG>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {

        na.push_back(j->first);
        _enzyme_addresses->insert({j->first, na});
        na.pop_back();
      }
    }
  }


  /******************* helpers *************************/

  bool isNotSpecial(const string& id) const {
    return (id != _e) && (id != _p) && (id != _c);
  }

  bool belong(string id, const map<string, string>& s) const {
    bool result = false;

    for (map<string, string>::const_iterator i = s.begin(); i != s.end(); ++i) {
      if (id == i->first) {
        result = true;
        break;
      }
    }

    return result;
  }

  string compartmentsOf(string id) const {
    string result;
    bool comp_found = false;

    for (map<string, map<string, string>>::const_iterator i = _species.begin(); i != _species.end(); ++i) {
      if (belong(id, i->second)) {
        result = i->first;
        comp_found = true;
        break;
      }
    }

    if (!comp_found) assert(comp_found && "some stoichiometry species does not belong to any compartments.");

    return result;
  }

  string getPlace(const SetOfMolecules_t& r, const SetOfMolecules_t& p) const {

    string result;
    int ac = 0;
    map<string, int> comps;

    for (map<string, string>::const_iterator i = _compartments.cbegin(); i != _compartments.cend(); ++i) {
      comps.insert({i->first, 0});
    }

    for (SetOfMolecules_t::const_iterator jt = r.cbegin(); jt != r.cend(); ++jt) {
      comps[this->compartmentsOf(jt->first)] += 1;
    }

    for (SetOfMolecules_t::const_iterator jt = p.cbegin(); jt != p.cend(); ++jt) {  
      comps[this->compartmentsOf(jt->first)] += 1;
    }

    for (map<string, int>::iterator i = comps.begin(); i != comps.end(); ++i) {   
      if (i->second > 0) ++ac;
    }

    switch(ac) {
      case 1:

        result = compOfReactant(comps) + "_i";
        break;
      case 2:

        if (isOrganelleTransport(comps)) {
          result = organelleOfReactant(comps) + "_m";
        
        } else if (comps[_e] > 0){

          result = _p + "_um";
        } else if (comps[_c] > 0) {

          result = _p + "_lm";
        } else if ((comps[_e] > 0) && (comps[_c] > 0)) {

          result = _p + "_tm";
        } else {

          assert(false && "The specie stoichiometry belong to a wrong conbination of 2 compartments.");
        }
        break;
      case 3:
        if ((comps[_e] > 0) && (comps[_p] > 0) && (comps[_c] > 0))
          result = _p + "_tm";
        else
          assert(false && "The specie stoichiometry belong to a wrong conbination of 3 compartments.");
        break;
    }

    return result;
  }

  string compOfReactant(const map<string, int>& c) const {

    string result;

    for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
      
      if (i->second > 0) {
        result = i->first;
        break;
      }
    }

    return result;
  }

  bool isOrganelleTransport(const map<string, int>& c) const {

    bool result = false;
    bool not_special;

    for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
      
      not_special = (i->first != _e) && (i->first != _p) && (i->first != _c);
      
      if (not_special && (i->second > 0)) {
        result = true;
        break;
      }
    }

    return result;
  }

  string organelleOfReactant(const map<string, int>& c) const {

    string result;
    bool not_special;

    for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
      
      not_special = (i->first != _e) && (i->first != _p) && (i->first != _c);
      
      if (not_special && (i->second > 0)) {
        result = i->first;
        break;
      }
    }

    return result;
  }

};

#endif // BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H