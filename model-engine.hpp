#ifndef BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H
#define BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H

#include <iostream>
#include <string>
#include <utility>
#include <map>

#include "parser/parser.hpp"
#include "data-structures/types.hpp"

using namespace std;

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
	// shared_ptr<flattened_coupled<TIME, MSG>> 		_biomass_model;
	// map< string, coupledModelsMap_t<TIME, MSG>> 	_enzyme_models;
	// coupledModelsMap_t<TIME, MSG> 					      _compartment_models;
	// coupledModelsMap_t<TIME, MSG>                _enzyme_set_models;
	// vectorOfCoupledModels_t<TIME, MSG>           _organelle_models;
	// shared_ptr<flattened_coupled<TIME, MSG>> 		_extra_cellular_model;
	// shared_ptr<flattened_coupled<TIME, MSG>> 		_periplasm_model;
	// shared_ptr<flattened_coupled<TIME, MSG>> 		_cytoplasm_model;
	// shared_ptr<flattened_coupled<TIME, MSG>> 		_cell_model;
	// shared_ptr< model<TIME> >						        _input;
	// shared_ptr< model<TIME> >						        _show_request;

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
};

#endif // BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H