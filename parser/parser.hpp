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
  long double                 _cell_weight;
  Integer_t                   _norm_number;

public:
  // Constructors
  Parser_t() = default;
	Parser_t(const char *filename, const string other_biomass_ID, const long double other_cell_weight, const Integer_t other_norm_number)
  : _document(filename), _biomass_ID(other_biomass_ID), _cell_weight(other_cell_weight), _norm_number(other_norm_number) {};
	
  // methods to load the XML files in the class
  bool loadFile(const char *filename);
  bool loadFile();

  // geting information from the current XML file loaded.
  list<UnitDefinition_t> getUnitDefinitions();
  map<string, string> getCompartments();
  map<string, map<string, string> > getSpecies();
  string getSpecieCompartment(string);
  map<string, enzyme_info_t > getEnzymes();
  enzyme_info_t getBiomass();
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H