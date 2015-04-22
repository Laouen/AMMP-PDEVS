#ifndef BOOST_SIMULATION_PDEVS_PARSER_H
#define BOOST_SIMULATION_PDEVS_PARSER_H

// STL
#include <list>
#include <map>
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
  map<string, TiXmlElement*> _models;

public:
  // Constructors
  Parser_t() = default;
	Parser_t(const char *filename)
  : _document(filename) {};
	
  // methods to load the XML files in the class
  bool loadFile(const char *filename);
  bool loadFile();

  // geting information from the current XML file loaded.
  list<UnitDefinition> getUnitDefinitions();
  map<string, string> getCompartments();
  map<string, map<string, string> > getSpeciesByCompartment();
  map<string, enzyme_parameter_t > getReactions();
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H