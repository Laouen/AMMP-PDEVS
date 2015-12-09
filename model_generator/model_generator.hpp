#ifndef BOOST_SIMULATION_PDEVS_MODEL_GENERATOR_H
#define BOOST_SIMULATION_PDEVS_MODEL_GENERATOR_H

#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <vector>
#include <assert.h>

#include "../parser/parser.hpp"
#include "../data-structures/types.hpp"

// Boost simalator include
#include <boost/simulation.hpp>

#define COMMENTS false

// atomic model includes
#include "../atomic-models/filter.hpp"
#include "../atomic-models/reaction.hpp"
#include "../atomic-models/space.hpp"
#include "../atomic-models/biomass.hpp"

template<class TIME,class MSG>
class ModelGenerator {

private:
	Parser_t parser;
  TIME _it, _rt, _r, _br, _bit, _current_time; 
  shared_ptr< map<string, Address_t> > _species_addresses;
  string _biomass_ID, _p, _e, _c;
  map<string, double> _compartment_volums;
  
  // Models
  map< string, vm_t<TIME>> _reaction_models;
  mm_t<TIME> _reaction_set_models;
  mm_t<TIME> _space_models;
  mm_t<TIME> _bulk_solutions;
  shared_ptr<model<TIME>> _biomass_model;
  shared_ptr<model<TIME>> _periplasm_model;
  shared_ptr<model<TIME>> _cell_model;
  bool _comment_mode;

public:
	ModelGenerator(const char *filename, TIME other_it, TIME other_rt, TIME other_r, TIME other_br, TIME other_bit, TIME other_ct, bool other_comment_mode) 
  : parser(filename), _it(other_it), _rt(other_rt), _r(other_r), _br(other_br), _bit(other_bit), _current_time(other_ct), _species_addresses(new map<string, Address_t>), _comment_mode(other_comment_mode) {
    comment("Init ...");

    shared_ptr<map<string, Integer_t>> amounts(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> konSTPs(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> konPTSs(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> koffSTPs(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> koffPTSs(new map<string, Integer_t>());

    // TODO this is a temporal parametrization
    vector<string> reactID = parser.getReactionIDs();
    for (vector<string>::iterator i = reactID.begin(); i != reactID.end(); ++i) {
      amounts->insert({*i, 100});
      konSTPs->insert({*i, 1});
      konPTSs->insert({*i, 1});
      koffSTPs->insert({*i, 0});
      koffPTSs->insert({*i, 0});
    }

    _p = "p"; 
    _e = "e"; 
    _c = "c"; 
    _biomass_ID = "R_Ec_biomass_iJO1366_WT_53p95M";

    // data initialization
    parser.loadExternParameters(100, amounts, konSTPs, konPTSs, koffSTPs, koffPTSs, 1, 1, _p, _e, _c, _biomass_ID);
    map<string, string> places = parser.getCompartments();
    parser.getSpecieByCompartments();
    parser.getReactions();
    parser.getEnzymes();

    map<string, string> compartments = parser.getCompartments();
    for (map<string, string>::iterator comp = compartments.begin(); comp != compartments.end(); ++comp) {
      _compartment_volums.insert({comp->first, 0.000000000000000000000000001});
    }
    comment("End init.");
	}

  Parser_t& getParser() {
    return parser;
  }

  /*****************************************************/
  /*************** REACTION MODELS *********************/
  /*****************************************************/
  
  /**
   * it return the reaction models grouped by the reaction set that they belong
   *
   */
  map< string, vm_t<TIME>>& getReactionModels() {
    comment("Get reaction models  ...");
    if (_reaction_models.empty()) {   
      comment("Creating reaction models  ...");

      if (_species_addresses->empty()) this->createSpeciesAddresses();

      pair<string, string> place;
      string reaction_set, id, comp;
      map<string, reaction_info_t> reactions = parser.getReactions();
      
      for (map<string, reaction_info_t>::iterator r = reactions.begin(); r != reactions.end(); ++r) {

        place = parser.getCompAndSubComp(r->second.substrate_sctry, r->second.products_sctry);
        reaction_set = place.first + "_" + place.second;
        id = r->first;

        auto reaction_model = make_atomic_ptr< 
        reaction<TIME, MSG>, 
        const string,
        const shared_ptr< map<string, Address_t> >,
        const bool,
        const TIME,
        const map<string, SetOfMolecules_t>&,
        const map<string, SetOfMolecules_t>&,
        const map<string, Integer_t>,
        const map<string, Integer_t>,
        const double,
        const double,
        const TIME,
        const TIME>(
          id, 
          _species_addresses,
          r->second.reversible, 
          _r,
          splitStoichiometryByCmpartments(r->second.substrate_sctry),
          splitStoichiometryByCmpartments(r->second.products_sctry),
          getStoichiometryCompartments(r->second.substrate_sctry),
          getStoichiometryCompartments(r->second.products_sctry),
          r->second.koffSTP,
          r->second.koffPTS,
          _it,
          _rt
        );

        // filter atomic model
        auto rFilter = make_atomic_ptr< filter<TIME, MSG>, const string>(id);
        
        if (_reaction_models.find(reaction_set) == _reaction_models.end()){
          vm_t<TIME> empty_entry;
          _reaction_models.insert({reaction_set, empty_entry});
        }

        // enzyme flattened_coupled model
        _reaction_models.at(reaction_set).push_back(make_shared<
          flattened_coupled<TIME, MSG>,
          vector<shared_ptr<model<TIME> > >,
          vector<shared_ptr<model<TIME> > >,
          vector<pair<shared_ptr<model<TIME>>, shared_ptr<model<TIME> > > >,
          vector<shared_ptr<model<TIME> > > 
        >(
          {rFilter, reaction_model}, 
          {rFilter}, 
          {{rFilter, reaction_model}}, 
          {reaction_model}
        ));
      }
      comment("End creating reaction models.");
    }
    comment("End get reaction models.");
    return _reaction_models;
  }

  void createSpeciesAddresses() {
    comment("Getting specie addresses ...");

    Address_t na;
    map<string, map<string, string> > species = parser.getSpecieByCompartments();

    for (map<string, map<string, string> >::const_iterator i = species.cbegin(); i != species.cend(); ++i) {
      
      na.clear();
      na.push_back(i->first);
      //na.push_back(i->first + "_space");
      for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
        na.push_back(j->first);
        _species_addresses->insert({j->first, na});
        na.pop_back();
      }
    }

    comment("End getting specie addresses.");
  }

  void comment(string msg) {
    
    if (_comment_mode) cout << "[Model generator] " << msg << endl;
  }

  map<string, Integer_t> getStoichiometryCompartments(const SetOfMolecules_t& sctry) {
    map<string, Integer_t> sctry_comps;
    string comp;
    for (SetOfMolecules_t::const_iterator i = sctry.cbegin(); i != sctry.cend(); ++i) {
      comp = parser.specieComp(i->first);
      if (sctry_comps.find(comp) == sctry_comps.end())
        sctry_comps.insert({comp, 0});
    }
    return sctry_comps;
  }

  map<string, SetOfMolecules_t> splitStoichiometryByCmpartments(const SetOfMolecules_t& sctry) {
    map<string, SetOfMolecules_t> sctry_by_comps;
    string comp;
    for (SetOfMolecules_t::const_iterator i = sctry.cbegin(); i != sctry.cend(); ++i) {
      comp = parser.specieComp(i->first);
      if (sctry_by_comps.find(comp) == sctry_by_comps.end()){
        SetOfMolecules_t new_set;
        sctry_by_comps.insert({comp, new_set});
      } 
      sctry_by_comps.at(comp).insert({i->first, i->second});
    }
    return sctry_by_comps;
  }

  /*****************************************************/
  /************** ENZYME SET MODELS ********************/
  /*****************************************************/

  mm_t<TIME>& getEnzymeSetModels() {
    comment("Getting enzyme set models ...");
    if (_reaction_set_models.empty()) {
      comment("Creating enzyme set models ...");
      
      vm_t<TIME> models, eic, eoc;
      vmp_t<TIME> ic; 

      if (_reaction_models.empty()) this->getReactionModels();

      for (typename map<string, vm_t<TIME>>::const_iterator set = _reaction_models.cbegin(); set != _reaction_models.cend(); ++set) {
        models.clear(); eic.clear(); eoc.clear(); ic.clear();

        auto es_filter = make_atomic_ptr< filter<TIME, MSG>, const string>(set->first);
        models.push_back(es_filter);
        eic.push_back(es_filter);
        for(typename vm_t<TIME>::const_iterator react = set->second.begin(); react != set->second.end(); ++react) {
          models.push_back(*react);
          eoc.push_back(*react);
          ic.push_back({es_filter, *react});
        }

        _reaction_set_models.insert({set->first, make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc)});
      }

      comment("End creating enzyme set models.");
    }
    comment("End getting enzyme set models.");
    return _reaction_set_models;
  }

  /*****************************************************/
  /***************** SPACE MODELS **********************/
  /*****************************************************/

  mm_t<TIME>& getSpaceModels() {
    comment("Getting space models ...");
    if (_space_models.empty()) {
      comment("Creating space models ...");
      string id;
      
      if(_reaction_set_models.empty()) this->getEnzymeSetModels();

      map<string, map<string, enzyme_t>> enzymes_by_comps = parser.getEnzymesByCompartments();

      for (map<string, map<string, enzyme_t>>::iterator comp = enzymes_by_comps.begin(); comp != enzymes_by_comps.end(); ++comp) {
        id = comp->first;
        Address_t biomass_address = parser.getBiomass().location;

        auto space_model = make_atomic_ptr< 
        space<TIME, MSG>,
        const string,
        const TIME,
        const TIME,
        const TIME,
        const Address_t&,
        const SetOfMolecules_t  ,
        const map<string, enzyme_t>&,
        const double>(
          id, 
          _it,
          _br,
          _current_time,
          biomass_address,
          parser.getCompartmentMetabolites(id),
          comp->second,
          _compartment_volums.at(id)
        );
        auto space_model_coupled = makeCoupledModel(space_model);

        auto space_filter = make_atomic_ptr<filter<TIME, MSG>, const string>(id);
        auto space_filter_coupled = makeCoupledModel(space_filter);

        vm_t<TIME> models = {space_model_coupled, space_filter_coupled}; 
        vm_t<TIME> eic = {space_filter_coupled}; 
        vmp_t<TIME> ic = {{space_filter_coupled, space_model_coupled}};
        vm_t<TIME> eoc = {space_model_coupled};

        _space_models.insert({id, make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc)});
      }

      comment("End creating space models.");
    }
    comment("End getting space models.");
    return _space_models;
  }

  /*****************************************************/
  /**************** BIOMASS MODELS *********************/
  /*****************************************************/

  Address_t getCompartmentAddressesFromBiomass() {
    Address_t result;

    map<string, string> compartments = parser.getCompartments();

    for (map<string, string>::const_iterator comp = compartments.cbegin(); comp != compartments.cend(); ++comp) {
      result.push_back(comp->first);
      result.push_back(comp->first + "_" + _biomass_ID);
    }

    return result;
  }

  map<string, bool> getCollaborationsMap() {
    map<string, bool> result;
    map<string, string> comps = parser.getCompartments();
    for (map<string,string>::iterator c = comps.begin(); c != comps.end(); ++c) {
      result.insert({c->first, false});
    }
    return result;
  }

  shared_ptr<model<TIME>>& getBiomassModel() {
    comment("Getting biomass model ...");

    if (!_biomass_model) {
      comment("Creating biomass model ...");


      reaction_info_t biomass_info = parser.getBiomass();

      auto biomass_model = make_atomic_ptr< 
        biomass<TIME, MSG>,
        const string,
        const shared_ptr< map<string, Address_t>>,
        const SetOfMolecules_t&,
        const SetOfMolecules_t&,
        const Address_t,
        const TIME,
        const TIME>(
          _biomass_ID, 
          _species_addresses,
          biomass_info.substrate_sctry,
          biomass_info.products_sctry,
          this->getCompartmentAddressesFromBiomass(),
          _bit,
          _br,
          this->getCollaborationsMap()
        );
      auto biomass_model_coupled = makeCoupledModel(biomass_model);

      auto biomass_filter = make_atomic_ptr<filter<TIME, MSG>, const string>(parser.getBiomass().location.front());
      auto biomass_filter_coupled = makeCoupledModel(biomass_filter);
      
      vm_t<TIME> models = {biomass_model_coupled, biomass_filter_coupled}; 
      vm_t<TIME> eic = {biomass_filter_coupled}; 
      vmp_t<TIME> ic = {{biomass_filter_coupled, biomass_model_coupled}};
      vm_t<TIME> eoc = {biomass_model_coupled};

      _biomass_model = make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc);
      comment("End creating biomass model.");
    }
    comment("End getting biomass model ...");
    return _biomass_model;
  }

  shared_ptr<model<TIME>> makeCoupledModel(shared_ptr<model<TIME>> atomic_model) {
    vm_t<TIME> models(1, atomic_model); 
    vm_t<TIME> eic(1, atomic_model); 
    vm_t<TIME> eoc(1, atomic_model);
    vmp_t<TIME> ic;
    return make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc);
  }

  /*****************************************************/
  /************** BULKSOLUTION MODELS ******************/
  /*****************************************************/

  shared_ptr<model<TIME>>& getBulkSolutionModel(string comp) {
    comment("Getting " + comp + " model ...");

    if (_bulk_solutions.find(comp) == _bulk_solutions.end()) {
      comment("Creating " + comp + " model ...");
      
      shared_ptr<model<TIME>> bulk_inner = this->getEnzymeSetModels().at(comp + "_inner");
      shared_ptr<model<TIME>> bulk_space = this->getSpaceModels().at(comp);
      auto bulk_filter = make_atomic_ptr<filter<TIME, MSG>, const string>(comp);
      auto biomass_filter = make_atomic_ptr<filter<TIME, MSG>, const string>(comp + "_" + _biomass_ID);
      auto bulk_filter_coupled = makeCoupledModel(bulk_filter);
      auto biomass_filter_coupled = makeCoupledModel(biomass_filter);

      vm_t<TIME> models = {bulk_inner, bulk_space, bulk_filter_coupled, biomass_filter_coupled}; 
      vm_t<TIME> eic = {bulk_filter_coupled, biomass_filter_coupled}; 
      vm_t<TIME> eoc = {bulk_space};
      vmp_t<TIME> ic = {
        {bulk_filter_coupled, bulk_space}, 
        {biomass_filter_coupled, bulk_space}, 
        {bulk_space, bulk_inner}, 
        {bulk_inner, bulk_space}
      };
      _bulk_solutions.insert({comp, make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc)});

      comment("End creating " + comp + " model.");
    }
    comment("End getting " + comp + " model ...");
    return _bulk_solutions.at(comp);
  }

  /*****************************************************/
  /*************** PERIPLASM MODELS ********************/
  /*****************************************************/

  shared_ptr<model<TIME>>& getPeriplasmModel() {
    comment("Getting periplasm model ...");

    if (!_periplasm_model) {
      comment("Creating periplasm model ...");
      
      shared_ptr<model<TIME>> periplasm_inner = this->getEnzymeSetModels().at("p_inner");
      shared_ptr<model<TIME>> periplasm_outer_membrane = this->getEnzymeSetModels().at("p_outer_membrane");
      shared_ptr<model<TIME>> periplasm_inner_membrane = this->getEnzymeSetModels().at("p_inner_membrane");
      shared_ptr<model<TIME>> periplasm_trans_membrane = this->getEnzymeSetModels().at("p_trans_membrane");
      shared_ptr<model<TIME>> periplasm_space = this->getSpaceModels().at("p");
      
      auto periplasm_filter = make_atomic_ptr<filter<TIME, MSG>, const string>("p");
      auto biomass_filter = make_atomic_ptr<filter<TIME, MSG>, const string>("p_" + _biomass_ID);
      auto show_filter = make_atomic_ptr<filter<TIME, MSG>, const string>("p_show_request");

      auto periplasm_filter_coupled = makeCoupledModel(periplasm_filter);
      auto biomass_filter_coupled = makeCoupledModel(biomass_filter);
      auto show_filter_coupled = makeCoupledModel(show_filter);

      vm_t<TIME> models = {periplasm_inner, periplasm_space, periplasm_outer_membrane, periplasm_inner_membrane, periplasm_trans_membrane, periplasm_filter_coupled, biomass_filter_coupled, show_filter_coupled}; 
      vm_t<TIME> eic = {periplasm_filter_coupled, biomass_filter_coupled, show_filter_coupled}; 
      vm_t<TIME> eoc = {periplasm_space, periplasm_outer_membrane, periplasm_inner_membrane, periplasm_trans_membrane};
      vmp_t<TIME> ic = {
        {biomass_filter_coupled, periplasm_space}, 
        {show_filter_coupled, periplasm_space}, 
        {periplasm_filter_coupled, periplasm_outer_membrane}, 
        {periplasm_filter_coupled, periplasm_inner_membrane}, 
        {periplasm_filter_coupled, periplasm_trans_membrane},
        {periplasm_outer_membrane, periplasm_space}, 
        {periplasm_inner_membrane, periplasm_space}, 
        {periplasm_trans_membrane, periplasm_space},
        {periplasm_space, periplasm_outer_membrane}, 
        {periplasm_space, periplasm_inner_membrane}, 
        {periplasm_space, periplasm_trans_membrane},
        {periplasm_space, periplasm_inner}, 
        {periplasm_inner, periplasm_space}
      };
      _periplasm_model = make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc);

      comment("End creating periplasm model.");
    }
    comment("End getting periplasm model ...");
    return _periplasm_model;
  }

  /*****************************************************/
  /*************** PERIPLASM MODELS ********************/
  /*****************************************************/

  // TODO: Add organelles
  shared_ptr<model<TIME>>& getCellModel() {
    comment("Getting cell model ...");

    if (!_cell_model) {
      comment("Creating cell model ...");
      
      shared_ptr<model<TIME>> extracellular = this->getBulkSolutionModel(_e);
      shared_ptr<model<TIME>> cytoplasm = this->getBulkSolutionModel(_c);
      shared_ptr<model<TIME>> periplasm = this->getPeriplasmModel();
      shared_ptr<model<TIME>> biomass = this->getBiomassModel();

      vm_t<TIME> models = {extracellular, cytoplasm, periplasm, biomass}; 
      vm_t<TIME> eic = {extracellular, periplasm, cytoplasm}; 
      vm_t<TIME> eoc = {};
      vmp_t<TIME> ic = {
        {extracellular, periplasm},
        {periplasm, extracellular},
        {periplasm, cytoplasm},
        {cytoplasm, periplasm},
        {biomass, extracellular},
        {biomass, periplasm},
        {biomass, cytoplasm},
        {extracellular, biomass},
        {periplasm, biomass},
        {cytoplasm, biomass}
      };
      _cell_model = make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc);

      comment("End creating cell model.");
    }
    comment("End getting cell model ...");
    return _cell_model;
  }

};

#endif // BOOST_SIMULATION_PDEVS_MODEL_GENERATOR_H