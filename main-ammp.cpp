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

// boost simalator include
#include <boost/simulation.hpp>

// XML parser include
#include "tinyXML/tinyxml.h"

// model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/message.hpp"

#define TIXML_USE_STL

using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;
using namespace std;


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


//int main(int argc, char* argv[]) {
//
//
//  /**************************************************************************************************************/
//  /************************************* Parsing the SBML file **************************************************/
//  /**************************************************************************************************************/
//
//  TiXmlDocument doc;
//
//  TiXmlElement * root, * model, * listOfReactions;
//  map<string, TiXmlElement *> model_lists;
//
//  if (argc <= 1){
//      cout << "one argument must to be defined" << endl;
//      exit(1);
//  }
//
//  if ( !doc.LoadFile(argv[1]) ){
//      cout << "fail loading" << endl;
//      exit(1);
//  }
//
//  root  = doc.FirstChildElement();
//  model = root->FirstChildElement();
//
//  for(TiXmlElement * it = model->FirstChildElement(); it != NULL; it = it->NextSiblingElement()){
//    model_lists[it->Value()] = it;
//  }
//
//  /**************************************************************************************************************/
//  /********************************* Creating unit definitions **************************************************/
//  /**************************************************************************************************************/
//
//  // declaring variables
//  string kind, factors_in_string;
//  double exponent, multiplier;
//  int scale;
//  list<Unit> list_of_units;
//  list<UnitDefinition> list_of_units_definitions;
//
//  for (TiXmlElement * current_unit_definition = model_lists["listOfUnitDefinitions"]->FirstChildElement(); current_unit_definition != NULL; current_unit_definition = current_unit_definition->NextSiblingElement()){
//    for(TiXmlElement * current_list_of_units = current_unit_definition->FirstChildElement(); current_list_of_units != NULL; current_list_of_units = current_list_of_units->NextSiblingElement()){
//
//      list_of_units.clear();
//      for(TiXmlElement * current_unit = current_list_of_units->FirstChildElement(); current_unit != NULL; current_unit = current_unit->NextSiblingElement()){
//        kind = current_unit->Attribute("kind");
//
//        if (current_unit->Attribute("multiplier") != NULL) multiplier = stod(current_unit->Attribute("multiplier"));
//        else multiplier = 1;
//
//        if (current_unit->Attribute("scale") != NULL) scale = stod(current_unit->Attribute("scale"));
//        else scale = 0;
//
//        if (current_unit->Attribute("exponent") != NULL) exponent = stod(current_unit->Attribute("exponent"));
//        else exponent = 1;
//
//        list_of_units.push_back(Unit(kind, exponent, scale, multiplier));
//      }
//
//      list_of_units_definitions.push_back(UnitDefinition(list_of_units, current_unit_definition->Attribute("id")));
//    }
//  }
//
//  /**************************************************************************************************************/
//  /*********************************** End creating unit definitions ********************************************/
//  /**************************************************************************************************************/
//
//  /**************************************************************************************************************/
//  /*********************************** End parsing the SBML file ************************************************/
//  /**************************************************************************************************************/
//
//
//
//  /**************************************************************************************************************/
//  /************************************** Making the models *****************************************************/
//  /**************************************************************************************************************/
//
//
//
//  /**************************************************************************************************************/
//  /************************************ End making the models ***************************************************/
//  /**************************************************************************************************************/
//
//  return 0;
//}


/**************************************************************************************************************/
/***************************************** Testing models *****************************************************/
/**************************************************************************************************************/

typedef map<string, pair<string, int> > stoichiometryDef;

bool randomDesition() {
  return ((rand() % 100) <= 50);
}

int main () {

  srand(time(NULL));

 /*
  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  piss->str("1 {organelle,org} {enzyme set,inner} {enzyme,enzyme 10} | ADP 2");
  
  auto pf = make_atomic_ptr<external_events<Time, Message, Time, string >, shared_ptr<istringstream>, Time>(piss, Time(0),
    [](const string& s, Time& t_next, Message& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
    int delimiter;
    string send_to;
    string collector;
    string thrash;
    stringstream ss;
    Message msg_out;

    ss.str(s);
    ss >> t_next;

    ss >> collector;

    while (collector != "|") {     
      send_to = "";

      do {

        if (send_to != "") send_to += " ";
        send_to += collector;
        ss >> collector;

      } while (send_to.at(send_to.size()-1) != '}');
      
      send_to.erase(send_to.begin());
      send_to.erase(send_to.end()-1);
      delimiter = send_to.find(",");


      msg_out.sendTo(send_to.substr(0, delimiter), send_to.substr(delimiter+1, -1));
    }

    ss >> msg_out.specie;
    ss >> msg_out.amount;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });


  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time, Message> > root( new flattened_coupled<Time, Message>{{pf, cell}, {}, {{pf, cell}}, {cell}});

  cout << "Preparing runner" << endl;
  Time initial_time{0};
  runner<Time, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << "specie: " << m.specie << endl << "amount: " << m.amount << endl; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;

  /**************************************************************************************************************/
  /********************************************* End Models *****************************************************/
  /**************************************************************************************************************/

  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




