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
#include "data-structures/types.hpp"

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

  srand(time(NULL));
  SetOfMolecules react_stoichiometry;
  SetOfMolecules prod_stoichiometry;
  int action, enzymes;
  double rate, interval_time, current_time;
  bool reversible, add;
  bool non_stop = true;
  map<string, shared_ptr< reaction<Time, Message> > > reactions;
  map<string, Message> messages;
  vector<Message> messages_to_send, output;
  string name, reaction_name, trash;
  Message msg;

  // initial reaction
  react_stoichiometry["A"] = 1;
  react_stoichiometry["B"] = 2;
  prod_stoichiometry["C"] = 3;
  prod_stoichiometry["D"] = 1;
  reactions["r1"] = make_shared< reaction<Time, Message> > (reaction<Time, Message>("r1", true, 0.3, react_stoichiometry, prod_stoichiometry, 1, 0.5) );

  msg.specie = "A";
  msg.amount = 8;
  messages["m1"] = msg;

  msg.specie = "B";
  msg.amount = 15;
  messages["m2"] = msg;

  msg.specie = "C";
  msg.amount = 23;
  messages["m3"] = msg;

  msg.specie = "D";
  msg.amount = 5;
  messages["m4"] = msg;

  while(non_stop) {
    system("clear");
    cout << "chois next action:" << endl;
    cout << "1 \t create new reaction" << endl;
    cout << "2 \t create new message" << endl;
    cout << "3 \t show messages" << endl;
    cout << "4 \t show reactions" << endl;
    cout << "5 \t make external" << endl;
    cout << "6 \t make internal" << endl;
    cout << "7 \t make out" << endl;
    cout << "8 \t make advance" << endl;
    cout << "9 \t make confluence" << endl;
    cout << "other \t exit" << endl;
    cin >> action;

    switch(action) {
      case 1:
        system("clear");
        react_stoichiometry.clear();
        prod_stoichiometry.clear();

        cout << "name: ";
        cin >> reaction_name;
        
        cout << "make react_stoichiometry:" << endl;
        cout << "add specie: ";
        cin >> add;
        while(add) {
          cout << "new specie name: ";
          cin >> name;
          cout << "amount: ";
          cin >> react_stoichiometry[name];
          cout << "add specie: ";
          cin >> add;
        }

        cout << "make prod_stoichiometry:" << endl;
        cout << "add specie: ";
        cin >> add;
        while(add) {
          cout << "new specie name: ";
          cin >> name;
          cout << "amount: ";
          cin >> prod_stoichiometry[name];
          cout << "add specie: ";
          cin >> add;
        }

        cout << "reversible: ";
        cin >> reversible;

        cout << "rate time: ";
        cin >> rate;

        cout << "interval time: ";
        cin >> interval_time;

        cout << "total enzymes: ";
        cin >> enzymes;

        reactions[reaction_name] = make_shared< reaction<Time, Message> > (reaction<Time, Message>(reaction_name, reversible, rate, react_stoichiometry, prod_stoichiometry, enzymes, interval_time) );
        break;

      case 2:
        system("clear");
        cout << "message name: ";
        cin >> name;

        cout << "specie: ";
        cin >> msg.specie;

        cout << "amount: ";
        cin >> msg.amount;

        messages[name] = msg;
        break;

      case 3:
        system("clear");
        for (map<string, Message>::iterator i = messages.begin(); i != messages.end(); ++i) {
          cout << i->first << ":" << endl;
          cout << i->second << endl;
        }

        cout << "come back to menu: ";
        cin >> trash;
        break;
      
      case 4:
        system("clear");
        for (map<string, shared_ptr< reaction<Time, Message> > >::iterator i = reactions.begin(); i != reactions.end(); ++i) {
          cout << i->first << ":" << endl;
          i->second->show(cout);
          cout << endl;
        }

        cout << "come back to menu: ";
        cin >> trash;
        break;

      case 5:
        system("clear");
        messages_to_send.clear();

        cout << "reaction name: ";
        cin >> reaction_name;

        cout << "add messages: ";
        cin >> add;

        while(add) {

          cout << "message name: ";
          cin >> name;
          messages_to_send.push_back(messages[name]);

          cout << "add messages: ";
          cin >> add;
        }

        cout << "time elapsed: ";
        cin >> current_time;

        reactions[reaction_name]->external(messages_to_send, current_time);
        break;

      case 6:
        system("clear");
        
        cout << "reaction name: ";
        cin >> reaction_name;

        reactions[reaction_name]->internal();
        break;

      case 7:
        system("clear");
        
        cout << "reaction name: ";
        cin >> reaction_name;

        output = reactions[reaction_name]->out();

        for (vector<Message>::iterator i = output.begin(); i != output.end(); ++i) {
          cout << *i << endl;
        }

        cout << "come back to menu: ";
        cin >> trash;
        break;

      case 8:
        system("clear");
        
        cout << "reaction name: ";
        cin >> reaction_name;

        cout << "next internal: ";
        cout << reactions[reaction_name]->advance();
        cout << endl;

        cout << "come back to menu: ";
        cin >> trash;
        break;

      case 9:
        system("clear");
        messages_to_send.clear();

        cout << "reaction name: ";
        cin >> reaction_name;

        cout << "add messages: ";
        cin >> add;

        while(add) {

          cout << "message name: ";
          cin >> name;
          messages_to_send.push_back(messages[name]);

          cout << "add messages: ";
          cin >> add;
        }

        cout << "time elapsed: ";
        cin >> current_time;

        reactions[reaction_name]->confluence(messages_to_send, current_time);
        break;

      default:
        system("clear");
        non_stop = false;
    }
  }
  
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




