#include "parser.hpp"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

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

// For the moment this funtion is not considering organelles
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
  to.push_back(places.first + "_" + places.second);
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

// GA = Gene association
string Parser_t::getGAValue(TiXmlElement *r) {

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

  string ga = getGAValue(r);
  vector<string> tokens;
  boost::split(tokens, ga, boost::is_any_of(" "));

  return proccesExpresion(tokens, 1, tokens.size()-1);;
}

vector<string> Parser_t::proccesExpresion(vector<string>& tokens, int start, int end) {
  int i, j;
  vector<vector<string>> processed_tokens;
  vector<string> subExpresion, result; 

  i = start;

  while(i <= end) {
    if (tokens[i].front() != '(') {
      processed_tokens.push_back({tokens[i]});
      ++i;
    } else {
      j = getSubExpresion(tokens, i);
      tokens[i].erase(0,1);
      tokens[j].erase(tokens[j].size()-1);
      processed_tokens.push_back(proccesExpresion(tokens, i, j));
      i = j+1;
    }
  }

  while(processed_tokens.size() > 1) {
    subExpresion.clear();
    if (processed_tokens[1].front() == "or") {
      for (i = 0; i < processed_tokens[0].size(); ++i) {
        subExpresion.push_back(processed_tokens[0][i]);
      }
      for (i = 0; i < processed_tokens[2].size(); ++i) {
        subExpresion.push_back(processed_tokens[2][i]);
      }
    } else if (processed_tokens[1].front() == "and") {
      for (i = 0; i < processed_tokens[0].size(); ++i) {
        for (j = 0; j < processed_tokens[2].size(); ++j) {
          subExpresion.push_back(processed_tokens[0][i] + "-" + processed_tokens[2][j]);
        }
      }
    }

    processed_tokens.erase(processed_tokens.begin(),processed_tokens.begin()+3);
    processed_tokens.insert(processed_tokens.begin(),subExpresion);
  }

  if (!processed_tokens.empty()) result = processed_tokens.front();
  return result;
}

int Parser_t::getSubExpresion(const vector<string>& tokens, int start) {
  assert(tokens[start].front() == '(');
  int i, to_close, end;

  to_close = 0;
  end = start-1;
  do {
    ++end;
    for (i=0; i<tokens[end].length(); ++i) {
      if (tokens[end].at(i) == '(') ++to_close;
      else break;
    }
    for (i=tokens[end].length()-1; i>=0; --i) {
      if (tokens[end].at(i) == ')') --to_close;
      else break;
    }
  } while(to_close > 0);

  return end;
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

vector<string> Parser_t::getEnzymeCompartments(const enzyme_t& en) {
  vector<string> compartments;
  for (map<string, reaction_info_t>::const_iterator r = en.handled_reactions.cbegin(); r != en.handled_reactions.cend(); ++r) {
    for (SetOfMolecules_t::const_iterator substrate = r->second.substrate_sctry.begin(); substrate != r->second.substrate_sctry.end(); ++substrate) {
      compartments.push_back(this->specieComp(substrate->first));
    }

    if (r->second.reversible) {
      for (SetOfMolecules_t::const_iterator product = r->second.products_sctry.begin(); product != r->second.products_sctry.end(); ++product) {
        compartments.push_back(this->specieComp(product->first));
      }
    }
  }

  return compartments;
}

/***************************************************************/
/******************* END HELPER FUNCTIONS **********************/
/***************************************************************/

Parser_t::Parser_t(const char *filename)
: _document(filename), _loaded(false), _enzyme_ID_counter(0) {

  this->loadFile();
};

void Parser_t::loadExternParameters(
  Integer_t other_default_amount,
  shared_ptr<map<string, Integer_t>> other_amounts,
  shared_ptr<map<string, Integer_t>> other_konSTPs,
  shared_ptr<map<string, Integer_t>> other_konPTSs,
  shared_ptr<map<string, Integer_t>> other_koffSTPs,
  shared_ptr<map<string, Integer_t>> other_koffPTSs,
  long double other_cw,
  Integer_t other_n,
  string other_p,
  string other_e,
  string other_c,
  string other_bID
) {
  
  _default_amount = other_default_amount;
  _amounts = other_amounts;
  _konSTPs = other_konSTPs;
  _konPTSs = other_konPTSs;
  _koffSTPs = other_koffSTPs;
  _koffPTSs = other_koffPTSs;
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
      if ((it->Attribute("reversible") == NULL) || ((string)it->Attribute("reversible") == "false")) react.reversible = false;
      else react.reversible = true;

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
      react.koffSTP = this->_koffSTPs->at(it->Attribute("id"));
      react.koffPTS = this->_koffPTSs->at(it->Attribute("id"));

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
        setReactionSctry(jt, react.substrate_sctry);
      } else if ((string)jt->Value() == "listOfProducts") {
        setReactionSctry(jt, react.products_sctry);
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
      if (!enzyme_hendlers.empty()) { 
        for (vector<string>::iterator enzymeID = enzyme_hendlers.begin(); enzymeID != enzyme_hendlers.end(); ++enzymeID) {
          if (this->_enzymes.find(*enzymeID) == this->_enzymes.end()) {

            enz.clear();
            enz.id = *enzymeID;
            if (this->_amounts->find(enz.id) != this->_amounts->end()) 
              enz.amount = this->_amounts->at(enz.id);
            else 
              enz.amount = this->_default_amount;
            enz.handled_reactions.insert({reactID, this->_reactions.at(reactID)});
            this->_enzymes.insert({enz.id, enz});
          } else {

            this->_enzymes.at(*enzymeID).handled_reactions.insert({reactID, this->_reactions.at(reactID)});
          }
        }
      } else {
        enz.clear();
        enz.id = "not_handled_" + _enzyme_ID_counter;
        enz.amount = this->_default_amount;
        enz.handled_reactions.insert({reactID, this->_reactions.at(reactID)});
        this->_enzymes.insert({enz.id, enz});
        ++_enzyme_ID_counter;
      }
    }
  }

  return this->_enzymes;
}

map<string, map<string, enzyme_t>> Parser_t::getEnzymesByCompartments() {
  map<string, map<string, enzyme_t>> result;
  vector<string> compartments;

  map<string, enzyme_t> enzymes = this->getEnzymes();

  for (map<string, enzyme_t>::iterator en = enzymes.begin(); en != enzymes.end(); ++en) {
    compartments = this->getEnzymeCompartments(en->second);
    for (vector<string>::iterator comp = compartments.begin(); comp != compartments.end(); ++comp) {
      if (result.find(*comp) == result.end()) {
        map<string, enzyme_t> empty_entry;
        result.insert({*comp, empty_entry});
      }
      result.at(*comp).insert({en->first, en->second});
    }
  }

  return result;
}

SetOfMolecules_t Parser_t::getCompartmentMetabolites(string comp) {
  SetOfMolecules_t result;
  map<string, string> species = this->getSpecieByCompartments().at(comp); 
  for (map<string, string>::iterator s = species.begin(); s != species.end(); ++s) {
    result.insert({s->first, 0});
  } 
  return result;
}