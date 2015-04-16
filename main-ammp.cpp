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

map<string, enzyme_parameter > getMetabolites(const map<string, enzyme_parameter >& reactions) {

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
  map<string, enzyme_parameter >    reactions     = input_doc.getReactions();
  /*
  << "creating space atomic models"
  Time interval_time;
  map<string, metabolite_info_t> enzymes;
  map<string, enzyme_info_t> metabolites;
  double volume, factor, amount;
  vectorOfModels spaces = {};

  for (map<string, string>::iterator ct = compartements.begin(); ct != compartements.end(); ++ct) {
    
    cout << "interval time for: " << i->second << endl;
    cin >> interval_time;
    cout << "volume time for: " << i->second << endl;
    cin >> volume;
    cout << "factor time for: " << i->second << endl;
    cin >> factor;
    cout << "species initial amount: " << endl;
    cin >> amount;

    enzymes     = getEnzymes(reactions);

    spaces.push_back( make_atomic_ptr< space<Time, Message>, Time, map<string, metabolite_info_t>, map<string, enzyme_info_t>, double, double >(interval_time, metabolites, enzymes, volume, factor) );
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