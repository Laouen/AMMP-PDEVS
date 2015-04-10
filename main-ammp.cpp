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
#include "atomic-models/space.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/types.hpp"
#include "data-structures/randomNumbers.hpp"

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


/**************************************************************************************************************/
/***************************************** Testing models *****************************************************/
/**************************************************************************************************************/


int main () {


  try {

    cout << "Creating the space models" << endl;
    map<string, metabolite_info_t> metabolites;
    map<string, enzyme_info_t> enzymes;
    auto cytoplasm = make_atomic_ptr<space<Time, Message>, Time, map<string, metabolite_info_t>, map<string, enzyme_info_t>, double, double >(0.2, metabolites, enzymes, 250.0, 1.0);

    //shared_ptr<space<Time, Message>> c = dynamic_pointer_cast<space<Time, Message>>(cytoplasm);
    //c->weightedRandomBool(10000);

    cout << "Creating the model to insert the input from stream" << endl;
    auto piss = make_shared<istringstream>();
    string input = "";
    int a = 1;
    for (double i = 0.001; i < 0.109; i += 0.001) {
      input += to_string(i) + " | A 1 \n ";
      if ((a % 2) == 0 ) input += to_string(i) + " | B 1 \n ";
      if ((a % 3) == 0 ) input += to_string(i) + " | C 1 \n ";
      if ((a % 4) == 0 ) input += to_string(i) + " | D 1 \n ";
      if ((a % 5) == 0 ) input += to_string(i) + " | E 1 \n ";
      if ((a % 6) == 0 ) input += to_string(i) + " | F 1 \n ";
      if ((a % 7) == 0 ) input += to_string(i) + " | G 1 \n ";
      if ((a % 8) == 0 ) input += to_string(i) + " | H 1 \n ";
      if ((a % 9) == 0 ) input += to_string(i) + " | I 1 \n ";
      if ((a % 10) == 0 ) input += to_string(i) + " | J 1 \n ";
      if ((a % 11) == 0 ) input += to_string(i) + " | K 1 \n ";
      if ((a % 12) == 0 ) input += to_string(i) + " | L 1 \n ";
      ++a;
    }
    input.pop_back();
    input.pop_back();
    input.pop_back();
    input.pop_back();
    piss->str(input);
    
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


        msg_out.to.push_back("ID");
      }

      ss >> msg_out.specie;
      ss >> msg_out.amount;

      m_next = msg_out;
      ss >> thrash;
      if ( 0 != thrash.size()) throw exception();
    });


    cout << "Coupling the input to the model" << endl;
    shared_ptr< flattened_coupled<Time, Message> > root( new flattened_coupled<Time, Message>{{pf}, {}, {}, {pf}});

    cout << "Preparing runner" << endl;
    Time initial_time{0};
    runner<Time, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

    cout << "Starting simulation until passivate" << endl;

    auto start = hclock::now(); //to measure simulation execution time

    r.runUntilPassivate();

    auto elapsed = chrono::duration_cast< chrono::duration< Time, ratio<1> > > (hclock::now() - start).count();

    cout << "Simulation took:" << elapsed << "sec" << endl;

    /**************************************************************************************************************/
    /********************************************* End Models *****************************************************/
    /**************************************************************************************************************/
  } catch(exception& e) {
    
    cout << "exception found:" << endl;
    cout << e.what() << endl;
  }
  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




