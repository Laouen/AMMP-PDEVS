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

  map<string, string>               compartements   = input_doc.getCompartments();
  map<string, map<string, string> > species         = input_doc.getSpeciesByCompartment();
  map<string, enzyme_parameter_t >  reactions       = input_doc.getReactions();
  string                            e_ID            = "e";
  string                            p_ID            = "p";
  string                            c_ID            = "c";

  cout << "creating enzymes atomic models" << endl;
  vectorOfModels              reaction_models;
  Time                        interval_time, rate;
  Integer                     amount;
  Address                     space_location;

  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    interval_time   = 0.001;
    rate            = 0.001;
    amount          = 1;

    reaction_models.push_back( make_atomic_ptr< reaction<Time, Message>, string, bool, Time, SetOfMolecules&, SetOfMolecules&, Integer, Time >(it->first, it->second.reversible, rate, it->second.reactants_sctry, it->second.products_sctry, amount, interval_time) );
  }

/*
  int cantComps;
  string k;
  map<string, int> comps; 

  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    comps["c"] = 0;
    comps["p"] = 0;
    comps["e"] = 0;

    for (SetOfMolecules::iterator jt = it->second.reactants_sctry.begin(); jt != it->second.reactants_sctry.end(); ++jt) {
      k = jt->first.back();
      if (comps.find(k) == comps.end()) cout << jt->first << endl;
      comps[k] += 1;
    }

    for (SetOfMolecules::iterator jt = it->second.products_sctry.begin(); jt != it->second.products_sctry.end(); ++jt) {
      k = jt->first.back();
      if (comps.find(k) == comps.end()) cout << jt->first << endl;
      comps[k] += 1;
    }
    cantComps = 0;
    if (comps["c"] > 0) ++cantComps;
    if (comps["p"] > 0) ++cantComps;
    if (comps["e"] > 0) ++cantComps;

    if (cantComps > 2) {
      cout << it->first << endl;
    }


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