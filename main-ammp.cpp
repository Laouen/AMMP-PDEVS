// STD includes
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <memory>

// Boost simalator include
#include <boost/simulation.hpp>

// atomic model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"
#include "atomic-models/space.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/types.hpp"
#include "data-structures/randomNumbers.hpp"

// tinyXML parser
#include "parser/parser.hpp"

#define TIXML_USE_STL

using namespace std;
using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;


/***************************************/
/********* Type definations ************/
/***************************************/

typedef double Time;
typedef chrono::high_resolution_clock hclock;
typedef vector< shared_ptr< model<Time> > > vectorOfModels;
typedef vector< pair< shared_ptr< model<Time> >, shared_ptr< model<Time> > > > vectorOfModelPairs;
typedef map<string, shared_ptr< model<Time> > > modelsMap;



/***************************************/
/******** End type definations *********/
/***************************************/

/***************************************/
/********** Helper functions ***********/
/***************************************/

string compartmentOfReactants(const map<string, int>& c) {

  string result;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    if (i->second > 0) {
      result = i->first;
      break;
    }
  }

  return result;
}

string organelleOfReactants(const map<string, int>& c, string* sp) {

  string result;
  bool is_not_special;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    is_not_special = (i->first != sp[0]) && (i->first != sp[1]) && (i->first != sp[2]);
    
    if (is_not_special && (i->second > 0)) {
      result = i->first;
      break;
    }
  }

  return result;
}

bool thereAreOrganelleInvolved(const map<string, int>& c, string* sp) {

  bool result = false;
  bool is_not_special;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    is_not_special = (i->first != sp[0]) && (i->first != sp[1]) && (i->first != sp[2]);
    
    if (is_not_special && (i->second > 0)) {
      result = true;
      break;
    }
  }

  return result;
}

string getPlace(const enzyme_parameter_t& e, const map<string, map<string, string> >& s, const map<string, string>& c, string* sp) {

  string result;
  int amount_compartments;
  string curr_space;
  map<string, int> compartments;

  for (map<string, string>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    compartments[i->first] = 0;
  }

  for (SetOfMolecules::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }

  for (SetOfMolecules::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }


  amount_compartments = 0;
  for (map<string, int>::iterator i = compartments.begin(); i != compartments.end(); ++i) {   
    if (i->second > 0) ++amount_compartments;
  }

  switch(amount_compartments) {
    case 1:

      result = compartmentOfReactants(compartments) + "_i";
      break;
    case 2:

      if (thereAreOrganelleInvolved(compartments, sp)) {
        result = organelleOfReactants(compartments, sp) + "_m";
      
      } else if (compartments[sp[0]] > 0){

        result = sp[1] + "_um";
      } else if (compartments[sp[2]] > 0) {

        result = sp[1] + "_lm";
      }
      break;
    case 3:

      result = sp[1] + "_tm";
      break;
  }

  return result;
}

/***************************************/
/******** End helper functions *********/
/***************************************/

int main(int argc, char* argv[]) {


  /**************************************************************************************************************/
  /************************************* Parsing the SBML file **************************************************/
  /**************************************************************************************************************/

  if (argc <= 1){
      cout << "an SBML file is required." << endl;
      exit(1);
  }

  Parser_t input_doc(argv[1]);
  input_doc.loadFile();

  map<string, string>               compartements     = input_doc.getCompartments();
  map<string, map<string, string> > species           = input_doc.getSpeciesByCompartment();
  map<string, enzyme_parameter_t >  reactions         = input_doc.getReactions();
  string                            special_places[3]  = {"e", "p", "c"};

  cout << "creating enzymes atomic models and ordering them by places" << endl;
  map< string, modelsMap >  reaction_models;
  Time                      interval_time, rate;
  Integer                   amount;
  bool                      isCompartment;
  string                    place, sub_place;

  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    isCompartment = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (isCompartment) {
      reaction_models[it->first + "_i"] = {};
      reaction_models[it->first + "_m"] = {};
    }
  }

  reaction_models[special_places[0] + "_i"] = {};
  reaction_models[special_places[2] + "_i"] = {};
  reaction_models[special_places[1] + "_i"] = {};
  reaction_models[special_places[1] + "_lm"] = {};
  reaction_models[special_places[1] + "_um"] = {};
  reaction_models[special_places[1] + "_tm"] = {};

  cout << reaction_models.size() << endl;


  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    interval_time   = 0.001;
    rate            = 0.001;
    amount          = 1;

    place = getPlace(it->second, species, compartements, special_places);

    reaction_models.at(place)[it->first] = make_atomic_ptr< reaction<Time, Message>, string, bool, Time, SetOfMolecules&, SetOfMolecules&, Integer, Time >(it->first, it->second.reversible, rate, it->second.reactants_sctry, it->second.products_sctry, amount, interval_time);
  }

  cout << "creating addresses to send metabolites to enzymes" << endl;
  map<string, Address> reaction_addresses;
  Address new_address;

  for (map<string, modelsMap>::const_iterator it = reaction_models.cbegin(); it != reaction_models.cend(); ++it) {
    

    place       = it->first.substr(0, it->first.find_last_of("_"));
    sub_place   = it->first.substr(it->first.find_last_of("_") + 1);

    
    new_address.clear();
    new_address.push_back(place);
    new_address.push_back(it->first);
    if (sub_place == "i") new_address.push_back(place + "_s");

    for (modelsMap::const_iterator jt = it->second.cbegin(); jt != it->second.cend(); ++jt) {

      new_address.push_back(jt->first);
      reaction_addresses[jt->first] = new_address;
      new_address.pop_back();
    }
  }

  for (map<string, Address>::iterator i = reaction_addresses.begin(); i != reaction_addresses.end(); ++i) {
    cout << i->first << ": " << i->second << endl;
  }

  cout << reaction_addresses.size() << endl;


  /**************************************************************************************************************/
  /*********************************** End creating unit definitions ********************************************/
  /**************************************************************************************************************/

  /**************************************************************************************************************/
  /*********************************** End parsing the SBML file ************************************************/
  /**************************************************************************************************************/



  /**************************************************************************************************************/
  /************************************** Making the models *****************************************************/
  /**************************************************************************************************************/



  /**************************************************************************************************************/
  /************************************ End making the models ***************************************************/
  /**************************************************************************************************************/

  return 0;
}