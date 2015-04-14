#include "parser.hpp"

using namespace std;

bool Parser_t::loadFile(const char *fileName) {

	return _document.LoadFile(fileName);
}

bool Parser_t::loadFile() {

	return _document.LoadFile();
}

list<UnitDefinition> Parser_t::getUnitDefinitions() {

  string kind, name;
  double exponent, multiplier;
  int scale;
  list<UnitDefinition> result;
  TiXmlElement *listOfUnitDefinitions;
  list<Unit> list_of_units;

  name                = "listOfUnitDefinitions"; 
  TiXmlElement *root  = _document.FirstChildElement();
  TiXmlElement *model = root->FirstChildElement();
  for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if (it->Value() == name) {
      listOfUnitDefinitions = it;
      break;
    }
  }

  for (TiXmlElement *current_unit_definition = listOfUnitDefinitions->FirstChildElement(); current_unit_definition != NULL; current_unit_definition = current_unit_definition->NextSiblingElement()) {  
    for (TiXmlElement *list_of_sub_units = current_unit_definition->FirstChildElement(); list_of_sub_units != NULL; list_of_sub_units = list_of_sub_units->NextSiblingElement()) {

      list_of_units.clear();
      for(TiXmlElement * current_sub_unit = list_of_sub_units->FirstChildElement(); current_sub_unit != NULL; current_sub_unit = current_sub_unit->NextSiblingElement()){
        kind = current_sub_unit->Attribute("kind");

        if (current_sub_unit->Attribute("multiplier") != NULL) multiplier = stod(current_sub_unit->Attribute("multiplier"));
        else multiplier = 1;

        if (current_sub_unit->Attribute("scale") != NULL) scale = stod(current_sub_unit->Attribute("scale"));
        else scale = 0;

        if (current_sub_unit->Attribute("exponent") != NULL) exponent = stod(current_sub_unit->Attribute("exponent"));
        else exponent = 1;

        list_of_units.push_back(Unit(kind, exponent, scale, multiplier));
      }

      result.push_back(UnitDefinition(list_of_units, current_unit_definition->Attribute("id")));
    }
  }

  return result;
}

map<string, string> Parser_t::getCompartments() {
  
  map<string, string> result;

  string name;
  TiXmlElement *listOfCompartments;

  name                = "listOfCompartments"; 
  TiXmlElement *root  = _document.FirstChildElement();
  TiXmlElement *model = root->FirstChildElement();
  for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if (it->Value() == name) {
      listOfCompartments = it;
      break;
    }
  }

  for (TiXmlElement *it = listOfCompartments->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    
    result[it->Attribute("id")] = it->Attribute("name");
  }

  return result; 
}

map<string, map<string, string> > Parser_t::getSpeciesByCompartment() {
  map<string, map<string, string> > result;

  string listName, specieID;
  TiXmlElement *listOfSpecies;

  listName                = "listOfSpecies"; 
  TiXmlElement *root  = _document.FirstChildElement();
  TiXmlElement *model = root->FirstChildElement();
  for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if (it->Value() == listName) {
      listOfSpecies = it;
      break;
    }
  }

  for (TiXmlElement *it = listOfSpecies->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    
    specieID = it->Attribute("id");
    specieID.resize(specieID.length()-2);
    result[it->Attribute("compartment")][specieID] = it->Attribute("name");
  }

  return result; 
}

template<class TIME>
list< enzyme_parameter<TIME> > Parser_t::getReactions() {
  
  list< enzyme_parameter<TIME> > result;

  string name;
  TiXmlElement *listOfReactions;

  name                = "listOfReactions"; 
  TiXmlElement *root  = _document.FirstChildElement();
  TiXmlElement *model = root->FirstChildElement();
  for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if (it->Value() == name) {
      listOfReactions = it;
      break;
    }
  }

  // FALTA TERMINAR

  return result;
}
