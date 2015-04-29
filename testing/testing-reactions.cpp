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

typedef double Time_t;
typedef chrono::high_resolution_clock hclock_t;
typedef vector< shared_ptr< model<Time_t> > > vectorOfModels_t;
typedef vector< pair< shared_ptr< model<Time_t> >, shared_ptr< model<Time_t> > > > vectorOfModelPairs_t;



/***************************************/
/******** End type definations *********/
/***************************************/


/**************************************************************************************************************/
/***************************************** Testing models *****************************************************/
/**************************************************************************************************************/


int main () {

  srand(time(NULL));
  SetOfMolecules_t react_stoichiometry;
  SetOfMolecules_t prod_stoichiometry;
  int action, enzymes;
  double rate, interval_time, current_time;
  bool reversible, add;
  bool non_stop = true;
  map<string, shared_ptr< reaction<Time_t, Message_t> > > reactions;
  map<string, Message_t> messages;
  vector<Message_t> messages_to_send, output;
  string name, reaction_name, trash;
  Message_t msg;

  // initial reaction
  react_stoichiometry["A"] = 1;
  react_stoichiometry["B"] = 2;
  prod_stoichiometry["C"] = 3;
  prod_stoichiometry["D"] = 1;
  reactions["r1"] = make_shared< reaction<Time_t, Message_t> > (reaction<Time_t, Message_t>("r1", true, 0.3, react_stoichiometry, prod_stoichiometry, 1, 0.5) );

  // initial menssages
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

  // menu
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

        reactions[reaction_name] = make_shared< reaction<Time_t, Message_t> > (reaction<Time_t, Message_t>(reaction_name, reversible, rate, react_stoichiometry, prod_stoichiometry, enzymes, interval_time) );
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
        for (map<string, Message_t>::iterator i = messages.begin(); i != messages.end(); ++i) {
          cout << i->first << ":" << endl;
          cout << i->second << endl;
        }

        cout << "come back to menu: ";
        cin >> trash;
        break;
      
      case 4:
        system("clear");
        for (map<string, shared_ptr< reaction<Time_t, Message_t> > >::iterator i = reactions.begin(); i != reactions.end(); ++i) {
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

        for (vector<Message_t>::iterator i = output.begin(); i != output.end(); ++i) {
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

  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




