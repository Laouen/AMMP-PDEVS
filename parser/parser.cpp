#include "parser.hpp"

using namespace std;

/***************************************************************/
/******************** HELPER FUNCTIONS *************************/
/***************************************************************/

long double L   = 6.0221413e+23;
long double MOL = 1e-6;

// TODO this implementation should modify all the stoichimetry
Integer_t Parser_t::getStoichiometryFrom(long double amount) {

  return (unsigned long long)amount;
}

// is a temporary implemetation
Integer_t Parser_t::getBiomassStoichiometryFrom(long double amount, long double cw, Integer_t nm) {

  return (Integer_t) ((amount * MOL * L * cw) / nm);
}

string Parser_t::specieComp(string id) {
  assert(this->_loaded);

  string res;
  int a = 0;
  map<string, map<string, string>> s = this->getSpecieByCompartments();

  for (map<string, map<string, string>>::const_iterator i = s.begin(); i != s.end(); ++i) {
    if (i->second.find(id) != i->second.cend()) {
      res = i->first;
      ++a;
    }
  }

  if (a == 0) assert(false && "This specie does not belong to any compartment.");
  if (a > 1) assert(false && "More than one compartment for the specie " + id + ".");

  return res;
}

pair<string, string> Parser_t::getCompAndSubComp(const SetOfMolecules_t& st, const SetOfMolecules_t& pt) {
  assert(this->_loaded);

  pair<string, string> res;
  map<string, bool> comps;
  string c;

  for (SetOfMolecules_t::const_iterator jt = st.cbegin(); jt != st.cend(); ++jt) {
    c = this->specieComp(jt->first);
    if (comps.find(c) == comps.end()) comps.insert({c, true});
  }

  for (SetOfMolecules_t::const_iterator jt = pt.cbegin(); jt != pt.cend(); ++jt) {  
    c = this->specieComp(jt->first);
    if (comps.find(c) == comps.end()) comps.insert({c, true});
  }

  switch(ac.size()) {
  case 1:

    res.first = comps.begin()->first;
    res.second = "inner";
    break;
  case 2:

    res.first = _p;
    if (comps[_e] > 0){

      res = "outer_membrane";
    } else if (comps[_c] > 0) {

      res = "inner_membrane";
    } else if ((comps[_e] > 0) && (comps[_c] > 0)) {

      res = "trans_membrane";
    } else {

      assert(false && "The specie stoichiometry belong to a wrong conbination of 2 compartments.");
    }
    break;
  case 3:
    if ((comps[_e] > 0) && (comps[_p] > 0) && (comps[_c] > 0))
      res.first = _p;
      res.second = "trans_membrane";
    else
      assert(false && "The specie stoichiometry belong to a wrong conbination of 3 compartments.");
    break;
  }

  return res;
}

Address_t Parser_t::getReactionAddress(const SetOfMolecules_t& st, const SetOfMolecules_t& pt, string e) {
  assert(this->_loaded);

  Address_t to;
  pair<string, string> places = this->getCompAndSubComp(st, pt);

  to.push_back(places.first);
  to.push_back(places.second);
  to.push_back(e);

  return to;
}

/***************************************************************/
/******************* END HELPER FUNCTIONS **********************/
/***************************************************************/

bool Parser_t::loadFile(const char *fileName) {

  this->_loaded = _document.LoadFile(fileName);

  if (this->_loaded) {
    TiXmlElement *root  = _document.FirstChildElement();
    TiXmlElement *model = root->FirstChildElement();
    for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      _models[it->Value()] = it;
    }
  }

	return this->_loaded;
}

bool Parser_t::loadFile() {

  this->_loaded = _document.LoadFile();

  if (this->_loaded) {
    string name;
    TiXmlElement *root  = _document.FirstChildElement();
    TiXmlElement *model = root->FirstChildElement();
    for (TiXmlElement *it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      _models[it->Value()] = it;
    }
  }

  return this->_loaded;
}

list<UnitDefinition_t>& Parser_t::getUnitDefinitions() {
  assert(this->_loaded);

  if (this->_units.empty()) {

    string kind, name;
    double exponent, multiplier;
    int scale;
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

        this->_units.push_back(UnitDefinition_t(list_of_units, current_unit_definition->Attribute("id")));
      }
    }
  }

  return this->_units;
}

map<string, string>& Parser_t::getCompartments() {
  assert(this->_loaded);

  if (this->_comps.empty()){
    for (TiXmlElement *it = _models["listOfCompartments"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      this->_comps.insert({it->Attribute("id"), it->Attribute("name")});
    }
  }
  return this->_comps; 
}

map<string, map<string, string>>& Parser_t::getSpecieByCompartments(string s) {
  assert(this->_loaded);

  if (this->_speciesByComps.empty()) {
    string c, n, id;

    for (TiXmlElement *it = _models["listOfSpecies"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      c = it->Attribute("compartment");
      id = it->Attribute("id");
      n = it->Attribute("name");

      if (this->_speciesByComps.find(c) != this->_speciesByComps.end()){
        this->_speciesByComps.at(c).insert({id, n});
      } else {
        this->_speciesByComps.insert({c, {id, n}});
      }
    }
  }

  return this->_speciesByComps; 
}

map<string, reaction_info_t>& Parser_t::getReactions() {
  assert(this->_loaded);

  if (this->_reactions.empty()) {
    string specieID;
    Integer_t sctry_value;
    reaction_info_t react;

    for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      if (it->Attribute("id") == _biomass_ID) continue;

      react.clear();

      // Setting reversibilty
      if ((string)it->Attribute("reversible") == "false")  react.reversible = false;
      else react.reversible = true;

      // setting stoichiometries
      for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
        if ((string)jt->Value() == "listOfReactants") {

          for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
            
            specieID      = lt->Attribute("species");
            sctry_value = (lt->Attribute("stoichiometry") == NULL) ? 1 : stod(lt->Attribute("stoichiometry"));
            react.substrate_sctry.at(specieID) = getStoichiometryFrom(sctry_value);
          } 
        } else if ((string)jt->Value() == "listOfProducts") {

          for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
            
            specieID      = lt->Attribute("species");
            sctry_value = (lt->Attribute("stoichiometry") == NULL) ? 1 : stod(lt->Attribute("stoichiometry"));
            react.substrate_sctry.at(specieID) = getStoichiometryFrom(sctry_value);
          } 
        }
      }

      // setting location, amoun and both Knos
      react.location = this->getReactionAddress(react.substrate_sctry, react.products_sctry, it->Attribute("id"));
      react.amount = _amounts.at(it->Attribute("id"));
      react.konSTP = _knoSTPs.at(it->Attribute("id"));
      react.konPTS = _knoPTSs.at(it->Attribute("id"));

      // insert the new reaction to the res;
      this->_reactions.insert({it->Attribute("id"), react});
    }
  }


  return this->_reactions;
}

reaction_info_t Parser_t::getBiomass() {
  assert(this->_loaded);

  string specieID;
  Integer_t sctry_value;
  reaction_info_t react;

  for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {

    if (it->Attribute("id") != _biomass_ID) continue;
    
    // Setting reversibilty
    if ((string)it->Attribute("reversible") == "false")  react.reversible = false;
    else react.reversible = true;

    // setting stoichiometries
    for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
      if ((string)jt->Value() == "listOfReactants") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID    = lt->Attribute("species");
          sctry_value = (lt->Attribute("stoichiometry") == NULL) ? 1 : stod(lt->Attribute("stoichiometry"));
          react.substrate_sctry.at(specieID) = this->getBiomassStoichiometryFrom(sctry_value);
        } 
      } else if ((string)jt->Value() == "listOfProducts") {

        for (TiXmlElement *lt = jt->FirstChildElement(); lt != NULL; lt = lt->NextSiblingElement()) {
          
          specieID    = lt->Attribute("species");
          sctry_value = (lt->Attribute("stoichiometry") == NULL) ? 1 : stod(lt->Attribute("stoichiometry"));
          react.substrate_sctry.at(specieID) = this->getBiomassStoichiometryFrom(sctry_value);
        } 
      }
    }

    // setting location
    react.location = {this->_biomass_ID};

    break;
  }

  return react;
}
