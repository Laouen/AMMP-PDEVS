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

typedef double Time_t;
typedef chrono::high_resolution_clock hclock_t;
typedef vector< shared_ptr< model<Time_t> > > vectorOfModels_t;
typedef vector< pair< shared_ptr< model<Time_t> >, shared_ptr< model<Time_t> > > > vectorOfModelPairs_t;



/***************************************/
/******** End type definations *********/
/***************************************/

int main () {

  /**************************************************************************************************************/
  /**************************************** Testing filters *****************************************************/
  /**************************************************************************************************************/

  /*********************************************** TEST 1 *******************************************************/
  
  cout << "Creating the atomic models for the filters:" << endl;
  cout << "enzyme set: Membrane -> cytoplasm: Cytoplasm -> compartment: Cytoplasm space" << endl;
  auto membrane_filter = make_atomic_ptr< filter<Time_t>, string, string >("enzyme set", "Membrane");
  auto cytoplasm_filter = make_atomic_ptr< filter<Time_t>, string, string >("cytoplasm", "Cytoplasm");
  auto cytoplasm_space_filter = make_atomic_ptr< filter<Time_t>, string, string >("compartment", "Cytoplasm space");

  cout << "Creating the atomic models for the filters that should not let the message pass:" << endl;
  cout << "4 enzyme set type, 4 cytoplasm type, 4 compartment type, 4 enzyme types, 4 organelle types" << endl;
  auto mf1 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme set", "Inner");
  auto mf2 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme set", "membrane");
  auto mf3 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme set", "Cytoplasm space");
  auto mf4 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme set", "");
  auto cf1 = make_atomic_ptr< filter<Time_t>, string, string >("cytoplasm", "cyto");
  auto cf2 = make_atomic_ptr< filter<Time_t>, string, string >("cytoplasm", "cytoplasm");
  auto cf3 = make_atomic_ptr< filter<Time_t>, string, string >("cytoplasm", "Membrane");
  auto cf4 = make_atomic_ptr< filter<Time_t>, string, string >("cytoplasm", "");
  auto csf1 = make_atomic_ptr< filter<Time_t>, string, string >("compartment", "Cytoplasm-space");
  auto csf2 = make_atomic_ptr< filter<Time_t>, string, string >("compartment", "cytoplasm space");
  auto csf3 = make_atomic_ptr< filter<Time_t>, string, string >("compartment", "Cytoplasm");
  auto csf4 = make_atomic_ptr< filter<Time_t>, string, string >("compartment", "");
  auto e1 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme", "Membrane");
  auto e2 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme", "Cytoplasm");
  auto e3 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme", "Cytoplasm space");
  auto e4 = make_atomic_ptr< filter<Time_t>, string, string >("enzyme", "enzym");
  auto o1 = make_atomic_ptr< filter<Time_t>, string, string >("organelle", "Membrane");
  auto o2 = make_atomic_ptr< filter<Time_t>, string, string >("organelle", "Cytoplasm");
  auto o3 = make_atomic_ptr< filter<Time_t>, string, string >("organelle", "Cytoplasm space");
  auto o4 = make_atomic_ptr< filter<Time_t>, string, string >("organelle", "organ");

  /*********************************************** TEST 2 *******************************************************/

  cout << "Creating the atomic models for the filters:" << endl;
  cout << "enzyme set: Membrane -> cytoplasm: Cytoplasm -> compartment: Cytoplasm space" << endl;
  auto organelle_filter = make_atomic_ptr< filter<Time_t, Message_t>, string, string >("organelle", "org");
  auto inner_filter = make_atomic_ptr< filter<Time_t, Message_t>, string, string >("enzyme set", "inner");

  
  vectorOfModels_t models = {organelle_filter, inner_filter};

  for (int i = 0; i < 100; ++i){
    auto new_enzyme_filter = make_atomic_ptr< filter<Time_t, Message_t>, string, string >("enzyme", "enzyme " + to_string(i));
    models.push_back(new_enzyme_filter);
  }

  /************************************ COUPLED MODEL FOR BOTH TESTS **********************************************/

  cout << "Coupling the models into the filter_test_model" << endl;
  vectorOfModels_t      eic = {organelle_filter};
  vectorOfModels_t      eoc = models;
  vectorOfModelPairs_t  ic;

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

  shared_ptr< flattened_coupled<Time_t, Message_t> > filter_test_model( new flattened_coupled<Time_t, Message_t>{models, eic, ic, eoc});

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  piss->str("1 {organelle,org} {enzyme set,inner} {enzyme,enzyme 10} | for-enzyme-10 2");
  
  auto pf = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
    int delimiter;
    string send_to;
    string collector;
    string thrash;
    stringstream ss;
    Message_t msg_out;

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
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{{pf, filter_test_model}, {}, {{pf, filter_test_model}}, {filter_test_model}});

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "specie: " << m.specie << endl << "amount: " << m.amount << endl; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;

  /**************************************************************************************************************/
  /************************************** End testing filters ***************************************************/
  /**************************************************************************************************************/

  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




