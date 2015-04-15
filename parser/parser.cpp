#include "parser.hpp"

using namespace std;

bool Parser_t::loadFile(const char *fileName) {

  bool result = _document.LoadFile(fileName);

  if (result) {
    string name;
    TiXmlElement *root  = _document.FirstChildElement();
    TiXmlElement *model = root->FirstChildElement();
    for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      _models[it->Value()] = it;
    }
  }

	return result;
}

bool Parser_t::loadFile() {

  bool result = _document.LoadFile();

  if (result) {
    string name;
    TiXmlElement *root  = _document.FirstChildElement();
    TiXmlElement *model = root->FirstChildElement();
    for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      _models[it->Value()] = it;
    }
  }

  return result;
}

list<UnitDefinition> Parser_t::getUnitDefinitions() {

  string kind, name;
  double exponent, multiplier;
  int scale;
  list<UnitDefinition> result;
  list<Unit> list_of_units;

  for (TiXmlElement *current_unit_definition = _models["listOfUnitDefinitions"]->FirstChildElement(); current_unit_definition != NULL; current_unit_definition = current_unit_definition->NextSiblingElement()) {  
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

  for (TiXmlElement *it = _models["listOfCompartments"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    result[it->Attribute("id")] = it->Attribute("name");
  }

  return result; 
}

map<string, map<string, string> > Parser_t::getSpeciesByCompartment() {
  map<string, map<string, string> > result;
  string specieID;

  for (TiXmlElement *it = _models["listOfSpecies"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    
    specieID = it->Attribute("id");
    specieID.resize(specieID.length()-2);
    result[it->Attribute("compartment")][specieID] = it->Attribute("name");
  }

  return result; 
}

map<string, enzyme_parameter<double> > Parser_t::getReactions(double rate, unsigned long long amount, double interval_time) {
  
  map<string, enzyme_parameter<double> > result;
  string specieID;
  enzyme_parameter<double> new_enzyme;

  for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    
    result[it->Attribute("id")] = new_enzyme;
    enzyme_parameter current    = result[it->Attribute("id")];

    current.name          = it->Attribute("name");
    current.reversible    = (it->Attribute("reversible") != NULL) ? it->Attribute("reversible") : true;
    current.rate          = rate;
    current.amount        = amount;
    current.interval_time = interval_time;

    for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
      
      if ((string)it->Value() == "listOfReactants") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          specieID = lt->Attribute();
          specieID.resize(specieID.length()-2);
          current.reactants_sctry[specieID] = getStoichiometryFrom(lt->Attribute());
        } 
      } else if ((string)it->Value() == "listOfProducts") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          specieID = lt->Attribute();
          specieID.resize(specieID.length()-2);
          current.products_sctry[specieID] = getStoichiometryFrom(lt->Attribute());
        } 
      }
    }
  }


  return result;
}

// Get better this function
long long integer getStoichiometryFrom(double amount) {

  return (long long integer)amount;
}
