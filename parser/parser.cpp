#include "parser.hpp"

using namespace std;

/***************************************************************/
/******************** HELPER FUNCTIONS *************************/
/***************************************************************/

// is a temporary implemetation
unsigned long long getStoichiometryFrom(double amount) {

  return (unsigned long long)amount;
}

// is a temporary implemetation
unsigned long long getBiomassStoichiometryFrom(double amount) {

  return (unsigned long long)10000;
}

/***************************************************************/
/******************* END HELPER FUNCTIONS **********************/
/***************************************************************/

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

list<UnitDefinition_t> Parser_t::getUnitDefinitions() {

  string kind, name;
  double exponent, multiplier;
  int scale;
  list<UnitDefinition_t> result;
  list<Unit_t> list_of_units;

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

        list_of_units.push_back(Unit_t(kind, exponent, scale, multiplier));
      }

      result.push_back(UnitDefinition_t(list_of_units, current_unit_definition->Attribute("id")));
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

    result[it->Attribute("compartment")][it->Attribute("id")] = it->Attribute("name");
  }

  return result; 
}

map<string, enzyme_parameter_t > Parser_t::getReactions() {
  
  map<string, enzyme_parameter_t > result;
  string specieID, stoichiometry;
  enzyme_parameter_t new_enzyme;

  for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {

    if (it->Attribute("id") == _biomass_ID) continue;

    result[it->Attribute("id")] = new_enzyme;
    enzyme_parameter_t *current = &result[it->Attribute("id")];

    current->name = it->Attribute("name");

    if ( it->Attribute("reversible") == NULL )                current->reversible = true;
    else if ((string)it->Attribute("reversible") == "false")  current->reversible = false;
    else                                                      current->reversible = true;

    for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
      
      if ((string)jt->Value() == "listOfReactants") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID      = lt->Attribute("species");
          stoichiometry = (lt->Attribute("stoichiometry") == NULL) ? "1" : lt->Attribute("stoichiometry");
          current->reactants_sctry[specieID] = getStoichiometryFrom( stod(stoichiometry) );
        } 
      } else if ((string)jt->Value() == "listOfProducts") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID      = lt->Attribute("species");
          stoichiometry = (lt->Attribute("stoichiometry") == NULL) ? "1" : lt->Attribute("stoichiometry");
          current->products_sctry[specieID] = getStoichiometryFrom( stod(stoichiometry) );
        } 
      }
    }
  }


  return result;
}

enzyme_parameter_t Parser_t::getBiomass() {
  
  string specieID, stoichiometry;
  enzyme_parameter_t result;

  for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {

    if (it->Attribute("id") != _biomass_ID) continue;
    
    result.name = it->Attribute("name");

    if ( it->Attribute("reversible") == NULL )                result.reversible = true;
    else if ((string)it->Attribute("reversible") == "false")  result.reversible = false;
    else                                                      result.reversible = true;

    for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
      
      if ((string)jt->Value() == "listOfReactants") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID      = lt->Attribute("species");
          stoichiometry = (lt->Attribute("stoichiometry") == NULL) ? "1" : lt->Attribute("stoichiometry");
          result.reactants_sctry[specieID] = getBiomassStoichiometryFrom( stod(stoichiometry) );
        } 
      } else if ((string)jt->Value() == "listOfProducts") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID      = lt->Attribute("species");
          stoichiometry = (lt->Attribute("stoichiometry") == NULL) ? "1" : lt->Attribute("stoichiometry");
          result.products_sctry[specieID] = getBiomassStoichiometryFrom( stod(stoichiometry) );
        } 
      }
    }

    break;
  }


  return result;
}
