#ifndef BOOST_SIMULATION_PDEVS_PARSER_H
#define BOOST_SIMULATION_PDEVS_PARSER_H

// STL
#include <list>
#include <map>
#include <memory> // shared_ptr

// XML parser include
#include "../tinyXML/tinyxml.h"

// Data structures
//#include "../data-structures/types.hpp"
#include "../data-structures/unit_definition.hpp"

using namespace std;

template<class TIME>
struct enzyme_parameter {
  string&                           name;
  bool&                             reversible;
  TIME&                             rate;
  map<string, unsigned long long>&  reactants_sctry;
  map<string, unsigned long long>&  products_sctry;
  unsigned long long                amount;
  TIME&                             interval_time;
};


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
  map<string, enzyme_parameter<double> > getReactions(double, unsigned long long, double);
};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H