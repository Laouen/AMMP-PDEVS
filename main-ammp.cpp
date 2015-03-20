#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <boost/simulation.hpp>
#include "atomic-models/reaction.hpp"
#include "atomic-models/controler.hpp"
#include "data-structures/reaction_input.hpp"
#include "data-structures/unit_definition.hpp"
#include "data-structures/message.hpp"
#include "atomic-models/filter.hpp"
#include "tinyXML/tinyxml.h"

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

int main () {

  /**************************************************************************************************************/
  /**************************************** Testing filters *****************************************************/
  /**************************************************************************************************************/

  /*********************************************** TEST 1 *******************************************************/
  /*
  cout << "Creating the atomic models for the filters:" << endl;
  cout << "enzyme set: Membrane -> cytoplasm: Cytoplasm -> compartment: Cytoplasm space" << endl;
  auto membrane_filter = make_atomic_ptr< filter<Time>, string, string >("enzyme set", "Membrane");
  auto cytoplasm_filter = make_atomic_ptr< filter<Time>, string, string >("cytoplasm", "Cytoplasm");
  auto cytoplasm_space_filter = make_atomic_ptr< filter<Time>, string, string >("compartment", "Cytoplasm space");

  cout << "Creating the atomic models for the filters that should not let the message pass:" << endl;
  cout << "4 enzyme set type, 4 cytoplasm type, 4 compartment type, 4 enzyme types, 4 organelle types" << endl;
  auto mf1 = make_atomic_ptr< filter<Time>, string, string >("enzyme set", "Inner");
  auto mf2 = make_atomic_ptr< filter<Time>, string, string >("enzyme set", "membrane");
  auto mf3 = make_atomic_ptr< filter<Time>, string, string >("enzyme set", "Cytoplasm space");
  auto mf4 = make_atomic_ptr< filter<Time>, string, string >("enzyme set", "");
  auto cf1 = make_atomic_ptr< filter<Time>, string, string >("cytoplasm", "cyto");
  auto cf2 = make_atomic_ptr< filter<Time>, string, string >("cytoplasm", "cytoplasm");
  auto cf3 = make_atomic_ptr< filter<Time>, string, string >("cytoplasm", "Membrane");
  auto cf4 = make_atomic_ptr< filter<Time>, string, string >("cytoplasm", "");
  auto csf1 = make_atomic_ptr< filter<Time>, string, string >("compartment", "Cytoplasm-space");
  auto csf2 = make_atomic_ptr< filter<Time>, string, string >("compartment", "cytoplasm space");
  auto csf3 = make_atomic_ptr< filter<Time>, string, string >("compartment", "Cytoplasm");
  auto csf4 = make_atomic_ptr< filter<Time>, string, string >("compartment", "");
  auto e1 = make_atomic_ptr< filter<Time>, string, string >("enzyme", "Membrane");
  auto e2 = make_atomic_ptr< filter<Time>, string, string >("enzyme", "Cytoplasm");
  auto e3 = make_atomic_ptr< filter<Time>, string, string >("enzyme", "Cytoplasm space");
  auto e4 = make_atomic_ptr< filter<Time>, string, string >("enzyme", "enzym");
  auto o1 = make_atomic_ptr< filter<Time>, string, string >("organelle", "Membrane");
  auto o2 = make_atomic_ptr< filter<Time>, string, string >("organelle", "Cytoplasm");
  auto o3 = make_atomic_ptr< filter<Time>, string, string >("organelle", "Cytoplasm space");
  auto o4 = make_atomic_ptr< filter<Time>, string, string >("organelle", "organ");

  /*********************************************** TEST 2 *******************************************************/

  cout << "Creating the atomic models for the filters:" << endl;
  cout << "enzyme set: Membrane -> cytoplasm: Cytoplasm -> compartment: Cytoplasm space" << endl;
  auto organelle_filter = make_atomic_ptr< filter<Time, Message>, string, string >("organelle", "org");
  auto inner_filter = make_atomic_ptr< filter<Time, Message>, string, string >("enzyme set", "inner");

  
  vectorOfModels models = {organelle_filter, inner_filter};

  for (int i = 0; i < 100; ++i){
    auto new_enzyme_filter = make_atomic_ptr< filter<Time, Message>, string, string >("enzyme", "enzyme " + to_string(i));
    models.push_back(new_enzyme_filter);
  }

  /************************************ COUPLED MODEL FOR BOTH TESTS **********************************************/

  cout << "Coupling the models into the filter_test_model" << endl;
  vectorOfModels      eic = {organelle_filter};
  vectorOfModels      eoc = models;
  vectorOfModelPairs  ic;

  for (int i = 1; i < models.size(); ++i) {
    ic.push_back(make_pair(organelle_filter, models[i]));
  }

  for (int i = 2; i < models.size(); ++i) {
    ic.push_back(make_pair(inner_filter, models[i]));
    for (int j = 2; j < models.size(); ++j) {
      if (i != j) {
        ic.push_back(make_pair(models[i], models[j]));
      }
    }
  }

  shared_ptr< flattened_coupled<Time, Message> > filter_test_model( new flattened_coupled<Time, Message>{models, eic, ic, eoc});

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  piss->str("1 {organelle,org} {enzyme set,inner} {enzyme,enzyme 10} | for-enzyme-10 2");
  
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
  shared_ptr< flattened_coupled<Time, Message> > root( new flattened_coupled<Time, Message>{{pf, filter_test_model}, {}, {{pf, filter_test_model}}, {filter_test_model}});

  cout << "Preparing runner" << endl;
  Time initial_time{0};
  runner<Time, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << "specie: " << m.specie << endl << "amount: " << m.amount << endl; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;

  /**************************************************************************************************************/
  /************************************** End testing filters ***************************************************/
  /**************************************************************************************************************/

  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




