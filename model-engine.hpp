#ifndef BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H
#define BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H

#include <iostream>
#include <string>
#include <utility>
#include <map>

#include "parser/parser.hpp"
#include "data-structures/types.hpp"

// Boost simalator include
#include <boost/simulation.hpp>

// atomic model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"
#include "atomic-models/space.hpp"
#include "atomic-models/biomass.hpp"

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
	string _e, _c, _p, _biomass_ID;
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
  : _model_ready(false), _comment_mode(false), _cell_weight(cw), _doc(fl, b, cw, n), _e(e), _c(c), _p(p), _biomass_ID(b) {

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

  void addCompartmentForEnzyme(const string& c) {
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

  void addCompartmentModel(const string& c, TIME it, TIME br, double v, double f) {
    if (_comment_mode) cout << "[Model engine] adding compartment: " + c << endl;

    map<string, enzyme_info_t> eis;
    enzyme_info_t ei;

    for (map<string, enzyme_parameter_t>::const_iterator i = _reactions.cbegin(); i != _reactions.end(); ++i) {

      if (belong(c, getCompartments(i->second))) {

        ei.location   = _enzyme_addresses->at(i->first);
        ei.reactants  = getReactants(i->second);
        eis.insert({i->first, ei});
      }
    }

    // creating filter and space atomic models
    auto cfilter = make_atomic_ptr< filter<TIME, MSG>, const string>(c + "_s");
    auto sm = make_atomic_ptr< 
      space<TIME, MSG>, 
      const string, 
      const TIME,
      const TIME,
      const Address_t&,
      const map<string, metabolite_info_t>&, 
      const map<string, enzyme_info_t>&, 
      const double, 
      const double >(
        c, 
        it,
        br,
        {"biomass", _biomass_ID}, // biomass address
        {}, // metabolites start as an empty map 
        eis, 
        v, 
        f
      );

    // creating compartment coupled model;
    _compartment_models[c] = make_shared<
      flattened_coupled<TIME, MSG>,
      vector<shared_ptr<model<TIME> > >,
      vector<shared_ptr<model<TIME> > >,
      vector<pair<shared_ptr<model<TIME>>, shared_ptr<model<TIME>> > >,
      vector<shared_ptr<model<TIME> > > 
    >(
      {cfilter, sm}, 
      {cfilter}, 
      {{cfilter, sm}}, 
      {sm}
    );
  }

  void createEnzymeSetModels() {
    if (_comment_mode) cout << "[Model engine] creating enzyme set models." << endl;
    vector<shared_ptr<model<TIME>>> models, eic, eoc;
    vmp_t<TIME> ic;

    for (typename map<string, cmm_t<TIME, MSG>>::const_iterator i = _enzyme_models.begin(); i != _enzyme_models.end(); ++i) {    
      
      models.clear();
      eic.clear();
      eoc.clear();
      ic.clear();

      auto esfilter = make_atomic_ptr< filter<TIME, MSG>, const string>(i->first);
      models.push_back(esfilter);
      eic.push_back(esfilter);
      
      for (typename cmm_t<TIME, MSG>::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
        models.push_back(j->second);
        eoc.push_back(j->second);
        ic.push_back({esfilter, j->second});
      }

      shared_ptr<flattened_coupled<TIME, MSG>> esm(new flattened_coupled<TIME, MSG>(models, eic, ic, eoc));
      _enzyme_set_models.insert({i->first, esm});
    }
  }

  void createCytoplasmModel() {
    if (_comment_mode) cout << "[Model engine] creating cytoplasm model." << endl;
    auto cytoplasm_filter     = make_atomic_ptr< filter<TIME, MSG>, const string>(_c);
    auto cytoplasm_space      = _compartment_models.at(_c);
    auto cytoplasm_inner      = _enzyme_set_models.at(_c + "_i");

    shared_ptr<flattened_coupled<TIME, MSG>> cm(new flattened_coupled<TIME, MSG>(
      {cytoplasm_filter, cytoplasm_space, cytoplasm_inner}, 
      {cytoplasm_filter}, 
      {
        {cytoplasm_filter, cytoplasm_space}, 
        {cytoplasm_space, cytoplasm_inner}, 
        {cytoplasm_inner, cytoplasm_space}
      }, 
      {cytoplasm_space}
    ));

    _cytoplasm_model = cm;
  }
  
  void createPeriplasmModel() {
    if (_comment_mode) cout << "[Model engine] creating periplasm model." << endl;

    // periplasm filter
    auto periplasm_filter = make_atomic_ptr< filter<TIME, MSG>, const string>(_p);
    // periplasm request filters
    auto out_filter = make_atomic_ptr< filter<TIME, MSG>, const string>(_p + "_or");
    auto bio_filter = make_atomic_ptr< filter<TIME, MSG>, const string>(_p + "_br");
    // periplasm output filters
    auto pof = make_atomic_ptr< filter<TIME, MSG>, const string>("output");
    auto pbf = make_atomic_ptr< filter<TIME, MSG>, const string>("biomass");
    
    // making a trivial coupled model of pof in order to solve the bug with coupled to atomic
    shared_ptr<flattened_coupled<TIME, MSG>> pocf(new flattened_coupled<TIME, MSG>(
      {pof}, {pof}, {}, {pof}
    ));

    // making a trivial coupled model of pbf in order to solve the bug with coupled to atomic
    shared_ptr<flattened_coupled<TIME, MSG>> pbcf(new flattened_coupled<TIME, MSG>(
      {pbf}, {pbf}, {}, {pbf}
    ));

    auto periplasm_space      = _compartment_models.at(_p);
    auto trans_membrane       = _enzyme_set_models.at(_p + "_tm");
    auto outer_membrane       = _enzyme_set_models.at(_p + "_um");
    auto inner_membrane       = _enzyme_set_models.at(_p + "_lm");
    auto periplasm_inner      = _enzyme_set_models.at(_p + "_i");

    shared_ptr<flattened_coupled<TIME, MSG>> periplasm_model(new flattened_coupled<TIME, MSG>(
      {periplasm_filter, out_filter, bio_filter, pocf, pbcf, periplasm_space, trans_membrane, outer_membrane, inner_membrane, periplasm_inner}, 
      {periplasm_filter, out_filter, bio_filter}, 
      {
        {out_filter, periplasm_space},
        {bio_filter, periplasm_space},
        {periplasm_filter, trans_membrane}, 
        {periplasm_filter, outer_membrane}, 
        {periplasm_filter, inner_membrane}, 
        {trans_membrane, periplasm_space}, 
        {outer_membrane, periplasm_space}, 
        {inner_membrane, periplasm_space}, 
        {periplasm_space, trans_membrane}, 
        {periplasm_space, outer_membrane}, 
        {periplasm_space, inner_membrane}, 
        {periplasm_space, periplasm_inner}, 
        {periplasm_inner, periplasm_space},
        {periplasm_space, pocf},
        {periplasm_space, pbcf}
      }, 
      {trans_membrane, outer_membrane, inner_membrane, pocf, pbcf}
    ));
  }
  
  void createExtraCellularModel() {
    if (_comment_mode) cout << "[Model engine] creating extra cellular model." << endl;

    auto extra_cellular_filter = make_atomic_ptr< filter<TIME, MSG>, const string>(_e);
    auto extra_cellular_space  = _compartment_models.at(_e);
    auto extra_cellular_inner  = _enzyme_set_models.at(_e + "_i");

    shared_ptr<flattened_coupled<TIME, MSG>> ecm(new flattened_coupled<TIME, MSG>(
      {extra_cellular_filter, extra_cellular_space, extra_cellular_inner}, 
      {extra_cellular_filter}, 
      {
        {extra_cellular_filter, extra_cellular_space},
        {extra_cellular_space, extra_cellular_inner},
        {extra_cellular_inner, extra_cellular_space}
      }, 
      {extra_cellular_space}
    ));

    _extra_cellular_model = ecm;
  }


  /******************* helpers *************************/
  vector<string> getReactants(const enzyme_parameter_t& e) const {
    vector<string> result;

    for (SetOfMolecules_t::const_iterator i = e.reactants_sctry.cbegin(); i != e.reactants_sctry.cend(); ++i) {
      result.push_back(i->first);
    }

    if (e.reversible) {
      for (SetOfMolecules_t::const_iterator i = e.products_sctry.cbegin(); i != e.products_sctry.cend(); ++i) {
        result.push_back(i->first);
      }
    }
    return result;
  }

  vector<string> getCompartments(const enzyme_parameter_t& e) {

    vector<string> result;
    map<string, int> comps;
    string cs;

    for (map<string, string>::const_iterator i = _compartments.cbegin(); i != _compartments.cend(); ++i) {
      comps.insert({i->first, 0});
    }

    for (SetOfMolecules_t::const_iterator i = e.reactants_sctry.cbegin(); i != e.reactants_sctry.cend(); ++i) {
      
      comps[this->compartmentsOf(i->first)] += 1;
    }

    if (e.reversible) {

      for (SetOfMolecules_t::const_iterator i = e.products_sctry.cbegin(); i != e.products_sctry.cend(); ++i) {
        
        comps[this->compartmentsOf(i->first)] += 1;
      }
    }

    for (map<string, int>::iterator i = comps.begin(); i != comps.end(); ++i) {   
      
      if (i->second > 0) result.push_back(i->first);
    }

    return result;
  }

  bool belong(const string& c, const vector<string>& comps) {
    bool result = false;

    for (vector<string>::const_iterator i = comps.begin(); i != comps.end(); ++i) {
      if (c == *i) {
        result = true;
        break;
      }
    }

    return result;
  }

  bool isNotSpecial(const string& id) const {
    
    return (id != _e) && (id != _p) && (id != _c);
  }

  string compartmentsOf(string id) const {
    string result;
    bool comp_found = false;

    for (map<string, map<string, string>>::const_iterator i = _species.begin(); i != _species.end(); ++i) {
      if (i->second.find(id) != i->second.cend()) {
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