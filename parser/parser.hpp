#ifndef BOOST_SIMULATION_PDEVS_PARSER_H
#define BOOST_SIMULATION_PDEVS_PARSER_H

// STL
#include <list>
#include <map>
#include <utility> /* pair */
#include <memory> /* shared_ptr */
#include <string> /* atof */

// XML parser include
#include "../tinyXML/tinyxml.h"

// Data structures
#include "../data-structures/types.hpp"
#include "../data-structures/unit_definition.hpp"

using namespace std;


class Parser_t {

private:
	TiXmlDocument               _document;
  map<string, TiXmlElement*>  _models;
  string                      _biomass_ID;
  string                      _p, _e, _c;

  // data
  list<UnitDefinition_t>            _units;
  map<string, string>               _comps;
  map<string, map<string, string>>  _speciesByComps;
  map<string, reaction_info_t>      _reactions;
  bool                              _loaded;

  // extern parameters
  shared_ptr<map<string, Integer_t>> _amounts;
  shared_ptr<map<string, Integer_t>> _konSTPs;
  shared_ptr<map<string, Integer_t>> _konPTSs;
  Integer_t                          _cell_weight;
  long double                        _normalization;

  // private methods
  Integer_t getBiomassStoichiometryFrom(double);
  Integer_t getStoichiometryFrom(double);
  string specieComp(string);
  pair<string, string> getCompAndSubComp(const SetOfMolecules_t&, const SetOfMolecules_t&);
  Address_t getReactionAddress(const SetOfMolecules_t&, const SetOfMolecules_t&, string);

public:
  // Constructors
  Parser_t() 
  : _loaded(false) {};
	Parser_t(const char *filename, const string other_biomass_ID)
  : _document(filename), _biomass_ID(other_biomass_ID), _loaded(false), _cell_weight(0), _normalization(0) {};
	
  // methods to load the XML files in the class
  bool loadFile(const char *filename);
  bool loadFile();
  void loadExternParameters(
    shared_ptr<map<string, Integer_t>> other_amounts,
    shared_ptr<map<string, Integer_t>> other_konSTPs,
    shared_ptr<map<string, Integer_t>> other_konPTSs,
    long double other_cw,
    Integer_t other_n,
    string other_p,
    string other_e,
    string other_c
  );

  // methods to get the information from the current SBML file.
  list<UnitDefinition_t>& getUnitDefinitions();
  map<string, string>& getCompartments();
  map<string, map<string, string>>& getSpecieByCompartments();
  map<string, reaction_info_t>& getReactions();
  vector<string> getReactionIDs();
  reaction_info_t getBiomass();
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H