#ifndef BOOST_SIMULATION_PDEVS_PARSER_H
#define BOOST_SIMULATION_PDEVS_PARSER_H

// XML parser include
#include "tinyXML/tinyxml.h"

using namespace std;

class Parser_t {

private:
	TiXmlDocument _document;

public:
	Parser_t() = default;
	
  bool LoadFile(string filename);
  list<TiXmlElement*> getUnitDefinitions();
  list<TiXmlElement*> getCompartments();
  list<TiXmlElement*> getSpecies();
  list<TiXmlElement*> getReactions();


};

#endif // BOOST_SIMULATION_PDEVS_PARSER_H