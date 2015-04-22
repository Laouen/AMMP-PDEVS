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



/***************************************/
/******** End type definations *********/
/***************************************/

/***************************************/
/********** Helper functions ***********/
/***************************************/

Address getLocation(const enzyme_parameter_t& e, const map<string, map<string, string> >& s) {

  return {};
}




void addLocation(map<string, enzyme_info_t>& enzymes, string new_location) {

  for (map<string, enzyme_info_t>::iterator it = enzymes.begin(); it != enzymes.end(); ++it){
    it->second.location.push_back(new_location);
  }
}

bool thereIsReactantsFrom(const vector<string>& reactants, const map<string, string>& species) {
  bool result = false;
  map<string, string>::const_iterator endOfMap = species.cend();

  for (vector<string>::const_iterator it = reactants.cbegin(); it != reactants.cend(); ++it) {
    if(species.find(*it) != endOfMap) {
      result = true;
      break;
    }
  }

  return result;
}

map<string, enzyme_info_t > getEnzymes(const map<string, enzyme_parameter_t >& reactions, const map<string, string>& species) {

  map<string, enzyme_info_t > result;
  string reactionID;
  enzyme_info_t new_enzyme;
  map<string, string>::const_iterator endOfMap = species.cend(); 

  for (map<string, enzyme_parameter_t >::const_iterator it = reactions.cbegin(); it != reactions.cend(); ++it) {

    reactionID          = it->first;
    new_enzyme.location = { it->first };
    new_enzyme.reactants.clear();

    for (SetOfMolecules::const_iterator jt = it->second.reactants_sctry.cbegin(); jt != it->second.reactants_sctry.cend(); ++jt) {
      new_enzyme.reactants.push_back(jt->first);
    }

    if (it->second.reversible) {
      for (SetOfMolecules::const_iterator jt = it->second.products_sctry.cbegin(); jt != it->second.products_sctry.cend(); ++jt) {
        new_enzyme.reactants.push_back(jt->first);
      }
    }

    if(thereIsReactantsFrom(new_enzyme.reactants, species))
      result[reactionID] = new_enzyme;
  }
  
  return result;
}

/***************************************/
/******** END Helper functions *********/
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

  map<string, string>               compartements = input_doc.getCompartments();
  map<string, map<string, string> > species       = input_doc.getSpeciesByCompartment();
  map<string, enzyme_parameter_t >  reactions     = input_doc.getReactions();

  cout << "creating enzymes atomic models" << endl;
  vectorOfModels              reaction_models;
  Time                        interval_time, rate;
  Integer                     amount;
  Address                     location;

  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    //cout << "interval time for: " << ct->second << endl;
    //cin >> interval_time;
    //cout << "volume time for: " << ct->second << endl;
    //cin >> volume;
    //cout << "factor time for: " << ct->second << endl;
    //cin >> factor;
    //cout << "enzymes initial amount: " << endl;
    //cin >> amount;
    interval_time   = 0.001;
    rate            = 0.001;
    amount          = 1;
    location        = getLocation(it->second, species);

    reaction_models.push_back( make_atomic_ptr< reaction<Time, Message>, string, Address, bool, Time, SetOfMolecules&, SetOfMolecules&, Integer, Time >(it->first, location, it->second.reversible, rate, it->second.reactants_sctry, it->second.products_sctry, amount, interval_time) );
  }


  
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