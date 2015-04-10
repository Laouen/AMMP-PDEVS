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

  TiXmlDocument doc;

  TiXmlElement *root, *lists, *listOfReactions;
  map<string, TiXmlElement *> model_lists;

  if (argc <= 1){
      cout << "The SBML file is required." << endl;
      exit(1);
  }

  if ( !doc.LoadFile(argv[1]) ){
      cout << "Fail loading SBML file." << endl;
      exit(1);
  }

  root  = doc.FirstChildElement();
  lists = root->FirstChildElement();

  for(TiXmlElement *it = lists->FirstChildElement(); it != NULL; it = it->NextSiblingElement()){
    model_lists[it->Value()] = it;
    cout << it->Value() << endl;
  }

  /**************************************************************************************************************/
  /********************************* Creating unit definitions **************************************************/
  /**************************************************************************************************************/

  // declaring variables
  string kind, factors_in_string;
  double exponent, multiplier;
  int scale;
  list<Unit> list_of_units;
  list<UnitDefinition> list_of_units_definitions;

  for (TiXmlElement * current_unit_definition = model_lists["listOfUnitDefinitions"]->FirstChildElement(); current_unit_definition != NULL; current_unit_definition = current_unit_definition->NextSiblingElement()){
    for(TiXmlElement * current_list_of_units = current_unit_definition->FirstChildElement(); current_list_of_units != NULL; current_list_of_units = current_list_of_units->NextSiblingElement()){

      list_of_units.clear();
      for(TiXmlElement * current_unit = current_list_of_units->FirstChildElement(); current_unit != NULL; current_unit = current_unit->NextSiblingElement()){
        kind = current_unit->Attribute("kind");

        if (current_unit->Attribute("multiplier") != NULL) multiplier = stod(current_unit->Attribute("multiplier"));
        else multiplier = 1;

        if (current_unit->Attribute("scale") != NULL) scale = stod(current_unit->Attribute("scale"));
        else scale = 0;

        if (current_unit->Attribute("exponent") != NULL) exponent = stod(current_unit->Attribute("exponent"));
        else exponent = 1;

        list_of_units.push_back(Unit(kind, exponent, scale, multiplier));
      }

      list_of_units_definitions.push_back(UnitDefinition(list_of_units, current_unit_definition->Attribute("id")));
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