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

// atomic model includes
#include "../atomic-models/filter.hpp"
#include "../atomic-models/reaction.hpp"
#include "../atomic-models/space.hpp"
#include "../atomic-models/biomass.hpp"

template<class TIME,class MSG>
class ModelGenerator {

private:
	Parser_t parser;
  TIME _it, _rt, _r, _br, _current_time; 
  shared_ptr< map<string, Address_t> > _species_addresses;
  string _biomass_ID;
  map<string, double> _compartment_volums;
  
  // Models
  map< string, vm_t<TIME>> _reaction_models;
  mm_t<TIME> _reaction_set_models;
  mm_t<TIME> _space_models;
  bool _comment_mode;

public:
	ModelGenerator(const char *filename, TIME other_it, TIME other_rt, TIME other_r, TIME other_br, TIME other_ct, bool other_comment_mode) 
  : parser(filename), _it(other_it), _rt(other_rt), _r(other_r), _br(other_br), _current_time(other_ct), _species_addresses(new map<string, Address_t>), _comment_mode(other_comment_mode) {
    comment("Init ...");

    shared_ptr<map<string, Integer_t>> amounts(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> konSTPs(new map<string, Integer_t>());
    shared_ptr<map<string, Integer_t>> konPTSs(new map<string, Integer_t>());

    // TODO this is a temporal parametrization
    vector<string> reactID = parser.getReactionIDs();
    for (vector<string>::iterator i = reactID.begin(); i != reactID.end(); ++i) {
      amounts->insert({*i, 100});
      konSTPs->insert({*i, 1});
      konPTSs->insert({*i, 1});
    }

    _biomass_ID = "R_Ec_biomass_iJO1366_WT_53p95M";

    // data initialization
    parser.loadExternParameters(100, amounts, konSTPs, konPTSs, 1, 1, "p", "e", "c", _biomass_ID);
    map<string, string> places = parser.getCompartments();
    parser.getSpecieByCompartments();
    parser.getReactions();
    parser.getEnzymes();

    map<string, string> compartments = parser.getCompartments();
    for (map<string, string>::iterator comp = compartments.begin(); comp != compartments.end(); ++comp) {
      _compartment_volums.insert({comp->first, 0});
    }
    comment("End init.");
	}

  /*****************************************************/
  /*************** REACTION MODELS *********************/
  /*****************************************************/
  
  /**
   * it return the reaction models grouped by the enzyme set that they belong
   *
   */
  map< string, vm_t<TIME>>& getReactionModels() {
    comment("Get reaction models  ...");
    if (_reaction_models.empty()) {   
      comment("Creating reaction models  ...");

      if (_species_addresses->empty()) this->createSpeciesAddresses();

      pair<string, string> place;
      string enzyme_set, id, comp;
      map<string, reaction_info_t> reactions = parser.getReactions();
      
      for (map<string, reaction_info_t>::iterator r = reactions.begin(); r != reactions.end(); ++r) {

        place = parser.getCompAndSubComp(r->second.substrate_sctry, r->second.products_sctry);
        enzyme_set = place.first + "_" + place.second;
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
          r->second.konSTP,
          r->second.konPTS,
          _it,
          _rt
        );

        // filter atomic model
        auto rFilter = make_atomic_ptr< filter<TIME, MSG>, const string>(id);
        
        if (_reaction_models.find(enzyme_set) == _reaction_models.end()){
          vm_t<TIME> empty_entry;
          _reaction_models.insert({enzyme_set, empty_entry});
        }

        // enzyme flattened_coupled model
        _reaction_models.at(enzyme_set).push_back(make_shared<
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
      na.push_back(i->first + "_s");
      for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
        _species_addresses->insert({j->first, na});
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
        Address_t biomass_address(1, _biomass_ID);

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

        // space flattened_coupled model
        vm_t<TIME> models(1, space_model); 
        vm_t<TIME> eic(1, space_model); 
        vm_t<TIME> eoc(1, space_model);
        vmp_t<TIME> ic; 
        _space_models.insert({id, make_shared<flattened_coupled<TIME, MSG>>(models, eic, ic, eoc)});
      }

      comment("End creating space models.");
    }
    comment("End getting space models.");
    return _space_models;
  }

};

#endif // BOOST_SIMULATION_PDEVS_MODEL_GENERATOR_H