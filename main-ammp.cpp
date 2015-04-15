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


int main(int argc, char* argv[]) {


  /**************************************************************************************************************/
  /************************************* Parsing the SBML file **************************************************/
  /**************************************************************************************************************/

  if (argc <= 1){
      cout << "The SBML file is required." << endl;
      exit(1);
  }

  Parser_t input_doc(argv[1]);
  input_doc.loadFile();

  input_doc.getReactions(0.001, 2, 0.01);
  /*
  list<UnitDefinition> units = input_doc.getUnitDefinitions();
  for (list<UnitDefinition>::iterator it = units.begin(); it != units.end(); ++it) {
    cout << it->unitName() << endl;
  }

  cout << endl;
  map<string, string> compartements = input_doc.getCompartments();
  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    cout << "id: " << it->first << " - name: " << it->second << endl;
  }

  cout << endl;
  map<string, map<string, string> > species = input_doc.getSpeciesByCompartment();

  for (map<string, map<string, string> >::iterator it = species.begin(); it != species.end(); ++it) {
    cout << " compartement: " << it->first << endl << endl;
   for (map<string, string>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
     cout << "id: " << jt->first << " - name: " << jt->second << endl;
   }
  }
  */

  
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