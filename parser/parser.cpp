#include "parser.hpp"
#include <boost/algorithm/string.hpp>

using namespace std;

/***************************************************************/
/******************** HELPER FUNCTIONS *************************/
/***************************************************************/

// TODO this implementation should modify all the stoichimetry
Integer_t Parser_t::getStoichiometryFrom(double amount) {

  return (Integer_t)amount;
}

// is a temporary implemetation
Integer_t Parser_t::getBiomassStoichiometryFrom(double amount) {
  assert((this->_cell_weight != 0) && (this->_normalization != 0));
  
  return (Integer_t) ((amount * MOL * L * this->_cell_weight) / this->_normalization);
}

string Parser_t::specieComp(string id) {
  assert(this->_loaded);

  string res;
  int a = 0;
  map<string, map<string, string>> s = this->getSpecieByCompartments();

  for (map<string, map<string, string>>::const_iterator i = s.begin(); i != s.end(); ++i) {
    if (i->second.find(id) != i->second.end()) {
      res = i->first;
      a += 1;
    }
  }

  assert((a == 1));
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

  switch(comps.size()) {
  case 1:

    res.first = comps.begin()->first;
    res.second = "inner";
    break;
  case 2:

    res.first = this->_p;
    if (comps[this->_e] > 0){

      res.second = "outer_membrane";
    } else if (comps[this->_c] > 0) {

      res.second = "inner_membrane";
    } else if ((comps[this->_e] > 0) && (comps[this->_c] > 0)) {

      res.second = "trans_membrane";
    } else {

      assert(false && "The specie stoichiometry belong to a wrong conbination of 2 compartments.");
    }
    break;
  case 3:
    if ((comps[this->_e] > 0) && (comps[this->_p] > 0) && (comps[this->_c] > 0)) {
      res.first = this->_p;
      res.second = "trans_membrane";
    } else {
      cout << _e << " " << _p << " " << _c << endl;
      assert(false && "The specie stoichiometry belong to a wrong conbination of 3 compartments.");
    }
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

void Parser_t::setReactionSctry(TiXmlElement * p, SetOfMolecules_t& sctry) {
  string specieID;
  double sctry_value;

  for (TiXmlElement *it = p->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    specieID    = it->Attribute("species");
    sctry_value = (it->Attribute("stoichiometry") == NULL) ? 1 : stod(it->Attribute("stoichiometry"));
    
    assert(sctry.find(specieID) == sctry.end());    
    sctry.insert({specieID, getStoichiometryFrom(sctry_value)});
  } 
}

string Parser_t::getGAStringParameter(TiXmlElement *r) {

  TiXmlElement *body;
  TiXmlElement *p;
  string result;

  for (TiXmlElement *it = r->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if ((string)it->Value() != "notes") continue;
    body   = it->FirstChildElement();
    p      = body->FirstChildElement();
    result = p->GetText ();
  }

  return result;
}

vector<string> Parser_t::getEnzymesHandlerIDs(TiXmlElement *r) {

  string ga = getGAStringParameter(r);
  vector<string> tokens, result;
  boost::split(tokens, ga, boost::is_any_of(" "));

  for (vector<string>::iterator t = tokens.begin(); t != tokens.end(); ++t) {
    
    if ((*t) == "or") continue;
    if ((*t) == "and") continue;
    if ((*t) == "GENE_ASSOCIATION:") continue;

    if ((*t)[0] == '(') (*t).erase(0,1);
    if ((*t)[(*t).size()-1] == ')') (*t).erase((*t).size()-1,1);

    result.push_back(*t);
    cout << *t << endl;
  }

  return result;
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

/***************************************************************/
/******************* END HELPER FUNCTIONS **********************/
/***************************************************************/

Parser_t::Parser_t(const char *filename)
: _document(filename), _loaded(false) {

  this->loadFile();
};


void Parser_t::loadExternParameters(
  shared_ptr<map<string, Integer_t>> other_amounts,
  shared_ptr<map<string, Integer_t>> other_konSTPs,
  shared_ptr<map<string, Integer_t>> other_konPTSs,
  long double other_cw,
  Integer_t other_n,
  string other_p,
  string other_e,
  string other_c,
  string other_bID
) {
  
  _amounts = other_amounts;
  _konSTPs = other_konSTPs;
  _konPTSs = other_konPTSs;
  _cell_weight = other_cw;
  _normalization = other_n;
  _p = other_p;
  _e = other_e;
  _c = other_c;
  _biomass_ID = other_bID;
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

map<string, map<string, string>>& Parser_t::getSpecieByCompartments() {
  assert(this->_loaded);

  if (this->_speciesByComps.empty()) {
    string c, n, id;
    map<string, string> new_value;

    for (TiXmlElement *it = _models["listOfSpecies"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      c = it->Attribute("compartment");
      id = it->Attribute("id");
      n = it->Attribute("name");

      if (this->_speciesByComps.find(c) != this->_speciesByComps.end()){
        this->_speciesByComps.at(c).insert({id, n});
      } else {
        new_value = {{id, n}};
        this->_speciesByComps.insert({c, new_value});
      }
    }
  }

  return this->_speciesByComps; 
}

map<string, reaction_info_t>& Parser_t::getReactions() {
  assert(this->_loaded);
  assert(this->_konSTPs && this->_konPTSs);

  if (this->_reactions.empty()) {

    string reactID;
    reaction_info_t react;

    for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
      
      reactID = it->Attribute("id"); 
      
      if (reactID == _biomass_ID) continue;

      react.clear();

      // Setting reversibilty
      if ((it->Attribute("reversible") == NULL) || ((string)it->Attribute("reversible") == "false")) react.reversible = true;
      else react.reversible = false;

      // setting stoichiometries
      for (TiXmlElement *jt = it->FirstChildElement(); jt != NULL; jt = jt->NextSiblingElement()) {
        if ((string)jt->Value() == "listOfReactants") {
          setReactionSctry(jt, react.substrate_sctry);
        } else if ((string)jt->Value() == "listOfProducts") {
          setReactionSctry(jt, react.products_sctry);
        }
      }

      // setting location, id, Kons and Pons
      react.id = reactID;
      react.location = this->getReactionAddress(react.substrate_sctry, react.products_sctry, it->Attribute("id"));
      react.konSTP = this->_konSTPs->at(it->Attribute("id"));
      react.konPTS = this->_konPTSs->at(it->Attribute("id"));

      // insert the new reaction to the reactions map;
      this->_reactions.insert({reactID, react});
    }
  }

  return this->_reactions;
}

vector<string> Parser_t::getReactionIDs() {

  vector<string> res;  
  for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {
    if (it->Attribute("id") != _biomass_ID) res.push_back(it->Attribute("id"));
  }
  return res;
}

reaction_info_t Parser_t::getBiomass() {
  assert(this->_loaded);

  string specieID;
  double sctry_value;
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

map<string, enzyme_t>& Parser_t::getEnzymes() {
  assert(this->_loaded);
  assert(this->_amounts);
  assert(!this->_reactions.empty());

  string reactID;
  enzyme_t enz;
  vector<string> enzyme_hendlers;

  if (this->_enzymes.empty()) {
    for (TiXmlElement *it = _models["listOfReactions"]->FirstChildElement(); it != NULL; it = it->NextSiblingElement()) {

      reactID = it->Attribute("id");
      if (reactID == _biomass_ID) continue;

      enzyme_hendlers = getEnzymesHandlerIDs(it);
      /*for (vector<string>::iterator enzymeID = enzyme_hendlers.begin(); enzymeID != enzyme_hendlers.end(); ++enzymeID) {
        if (this->_enzymes.find(*enzymeID) == this->_enzymes.end()) {

          enz.clear();
          enz.id      = *enzymeID;
          enz.amount  = this->_amounts->at(enz.id);
          enz.handled_reactions.insert({reactID, this->_reactions.at(reactID)});
          this->_enzymes.insert({enz.id, enz});
        } else {

          this->_enzymes.at(*enzymeID).handled_reactions.insert({reactID, this->_reactions.at(reactID)});
        }
      } */
    }
  }

  return this->_enzymes;
}
