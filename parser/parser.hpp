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

  // data
  list<UnitDefinition_t>            _units;
  map<string, string>               _comps;
  map<string, map<string, string>>  _speciesByComps;

  // private methods
  Integer_t getBiomassStoichiometryFrom(long double, long double, Integer_t);
  Integer_t getStoichiometryFrom(long double);
  string specieComp(string);
  string getCompAndSubComp(const SetOfMolecules_t&, const SetOfMolecules_t&);
  Address_t getReactionAddress(const SetOfMolecules_t&, const SetOfMolecules_t&, string);

public:
  // Constructors
  Parser_t() = default;
	Parser_t(const char *filename, const string other_biomass_ID)
  : _document(filename), _biomass_ID(other_biomass_ID) {};
	
  // methods to load the XML files in the class
  bool loadFile(const char *filename);
  bool loadFile();

  // methods to get the information from the current SBML file.
  list<UnitDefinition_t>& getUnitDefinitions();
  map<string, string>& getCompartments();
  map<string, map<string, string>>& getSpecieByCompartments();
  map<string, reaction_info_t> getReactions();
  reaction_info_t getBiomass();
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H