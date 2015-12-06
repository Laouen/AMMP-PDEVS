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
  map< string, cmm_t<TIME, MSG>> _reaction_models;
  TIME _it;
  TIME _rt;
  TIME _r;
  shared_ptr< map<string, Address_t> > _species_addresses;
  bool _comment_mode;

public:
	ModelGenerator(const char *filename, TIME other_it, TIME other_rt, TIME other_r, bool other_comment_mode) 
  : parser(filename), _it(other_it), _rt(other_rt), _r(other_r), _comment_mode(other_comment_mode), _species_addresses(new map<string, Address_t>) {
    comment("[Model generator] Init ...");

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

    // data initialization
    parser.loadExternParameters(100, amounts, konSTPs, konPTSs, 1, 1, "p", "e", "c", "R_Ec_biomass_iJO1366_WT_53p95M");
    map<string, string> places = parser.getCompartments();
    parser.getSpecieByCompartments();
    parser.getReactions();
    parser.getEnzymes();

    // settin reaction model places
    for (map<string, string>::iterator p = places.begin(); p != places.end(); ++p) {
      cmm_t<TIME, MSG> empty_cmm;
      _reaction_models.insert({p->first, empty_cmm});
    }

    this->createSpeciesAddresses();
	}

  /*****************************************************/
  /*************** REACTION MODELS *********************/
  /*****************************************************/
  
  void createReactionModels() {
    comment("[Model generator] Creating reaction models  ...");
    string place, id, comp;
    map<string, Integer_t> substrate_comps;
    map<string, Integer_t> product_comps;
    map<string, SetOfMolecules_t> substrate_sctry_by_comps;
    map<string, SetOfMolecules_t> products_sctry_by_comps;
    map<string, reaction_info_t> reactions = parser.getReactions();
    
    for (map<string, reaction_info_t>::iterator r = reactions.begin(); r != reactions.end(); ++r) {

      substrate_comps.clear();
      product_comps.clear();
      place = parser.getCompAndSubComp(r->second.substrate_sctry, r->second.products_sctry).first;
      id = r->first;

      this->getStoichiometryCompartments(r->second.substrate_sctry, substrate_comps);
      this->getStoichiometryCompartments(r->second.products_sctry, product_comps);
      this->separateStoichiometryByCmpartments(r->second.substrate_sctry, substrate_sctry_by_comps);
      this->separateStoichiometryByCmpartments(r->second.products_sctry, products_sctry_by_comps);

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
        substrate_sctry_by_comps,
        products_sctry_by_comps,
        substrate_comps,
        product_comps,
        r->second.konSTP,
        r->second.konPTS,
        _it,
        _rt
      );

      // filter atomic model
      auto rFilter = make_atomic_ptr< filter<TIME, MSG>, const string>(id);
      
      // enzyme flattened_coupled model
      _reaction_models.at(place)[id] = make_shared<
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
      );
    }
  }

  void createSpeciesAddresses() {
    comment("[Model generator] Getting specie addresses ...");

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
  }

  void comment(string msg) {
    
    if (_comment_mode) cout << msg << endl;
  }

  void getStoichiometryCompartments(const SetOfMolecules_t& sctry, map<string, Integer_t>& sctry_comps) {
    string comp;
    for (SetOfMolecules_t::const_iterator i = sctry.cbegin(); i != sctry.cend(); ++i) {
      comp = parser.specieComp(i->first);
      if (sctry_comps.find(comp) == sctry_comps.end())
        sctry_comps.insert({comp, 0});
    }
  }

  void separateStoichiometryByCmpartments(const SetOfMolecules_t& sctry, map<string, SetOfMolecules_t> sctry_by_comps) {
    string comp;
    for (SetOfMolecules_t::const_iterator i = sctry.cbegin(); i != sctry.cend(); ++i) {
      comp = parser.specieComp(i->first);
      if (sctry_by_comps.find(comp) == sctry_by_comps.end()){
        SetOfMolecules_t new_set;
        sctry_by_comps.insert({comp, new_set});
      } 
      sctry_by_comps.at(comp).insert({i->first, i->second});
    }
  }

  map< string, cmm_t<TIME, MSG>> getReactionModels() {
    return _reaction_models;
  }

};

#endif // BOOST_SIMULATION_PDEVS_MODEL_GENERATOR_H