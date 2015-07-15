// STD includes
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <memory>

// data structures
#include "data-structures/britime.hpp"
#include "data-structures/types.hpp"

// Boost simalator include
#include <boost/simulation.hpp>

// model engine
#include "model-engine.hpp"

// include vendors
#include "vendors/pdevs_tools.hpp"

// TEMPORAL FOR TEST
#include "parser/parser.hpp"

#define TIXML_USE_STL

using namespace std;
using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;

using Time_t   = BRITime;
using hclock_t = chrono::high_resolution_clock;


int main(int argc, char* argv[]) {

  bool comment_mode = true;

  if (argc <= 1){
    cout << "An SBML file is required." << endl;
    exit(1);
  }

  return 0;
}
  /*
  long double cell_weight   = 280 * 1e-15;  
  Integer_t norm_number     = 1;  
  string e                  = "e";
  string c                  = "c";
  string p                  = "p";
  string biomass_ID         = "R_Ec_biomass_iJO1366_WT_53p95M";

  ModelEngine<Time_t, Message_t> m(cell_weight, argv[1], e, c, p, biomass_ID, norm_number);
  m._comment_mode = comment_mode;
  
  m.createSpeciesAddresses();

  // adding the organelles compartments
  for (map<string, string> ::const_iterator i = m._compartments.begin(); i != m._compartments.end(); ++i) {
    if (m.isNotSpecial(i->first)) m.addCompartmentForEnzyme(i->first);
  }

  // creating the enzymes models
  m._comment_mode = false;
  for (map<string, enzyme_parameter_t >::const_iterator i = m._reactions.begin(); i != m._reactions.end(); ++i) {
      m.addEnzymeModel(i->first, BRITime(1,10), BRITime(1,100), Integer_t(100));
  }
  m._comment_mode = comment_mode;

  m.createEnzymeAddresses();

  m.createCytoplasmModel(BRITime(1,100), BRITime(1,100), 250, 1);

  m.createExtraCellularModel(BRITime(1,100), BRITime(1,100), 0, 1);
  
  m.createPeriplasmModel(BRITime(1,100), BRITime(1,100), 15, 1);

  m.createBiomassModel(BRITime(1,10), BRITime(1,100));

  // creating the organelles
  for (map<string, string> ::const_iterator i = m._compartments.begin(); i != m._compartments.end(); ++i) {
    if (m.isNotSpecial(i->first)) m.addOrganelleModel(i->first, BRITime(1,100), BRITime(1,100), 0, 1);
  }

  m.createCellModel();


  /*****************************************************************************************************/
  /************************************** Runing Simulation ********************************************/
  /*****************************************************************************************************/
  /*
  if (comment_mode) cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string partial_input = "";
  string input = "";
  string cp;
  for (Time_t t(1,100); t < Time_t(1,50); t += Time_t(1,100)) {

    for (map<string, map<string, string>>::const_iterator i = m._species.begin(); i != m._species.end(); ++i){
      if (i->first == "p") cp = "p_init";
      else cp = i->first;
      partial_input = t.toString() + " " + cp + " " + i->first + "_s | ";
      
      for (map<string,string>::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
        input += partial_input + j->first + " 100000 \n ";
      }
    }
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);

  auto pf = make_atomic_ptr<input_stream<Time_t, Message_t, Time_t, Message_t >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    int delimiter;
    string collector;
    string thrash;
    stringstream ss;
    Message_t msg_out;

    ss.str(s);
    ss >> t_next;
    ss >> collector;

    while (collector != "|") {     
  
      msg_out.to.push_back(collector);
      ss >>collector;
    }

    ss >> msg_out.specie;
    ss >> msg_out.amount;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });


  if (comment_mode) cout << "Creating the model to show the space state from stream" << endl;
  piss = make_shared<istringstream>();
  input = "";

  for (Time_t i(1,1000); i <= Time_t(300); i += Time_t(1,1000)) {

    input += i.toString() + " \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);

  auto so = make_atomic_ptr<input_stream<Time_t, Message_t, Time_t, Message_t >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    string thrash;
    stringstream ss;
    Message_t msg_out;

    ss.str(s);
    ss >> t_next;

    msg_out.to            = {"e", "e_s", "c", "c_s", "p_or", "p_s"};
    msg_out.show_request  = true;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });

  if (comment_mode) cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{
    {pf, so, m._cell_model}, 
    {}, 
    {
      {pf, m._cell_model},
      {so, m._cell_model}
    }, 
    {m._cell_model}
  });

  //pdevs_tools::pdevs_coupling_diagram<Time_t, Message_t> pd{*root};
  //string puml_result = pd.get_plant_uml();
  //cout << puml_result;
  //return 0;
  
  if (comment_mode) cout << "Preparing runner" << endl;
  runner<Time_t, Message_t> r(root, Time_t(0), cout, [](ostream& os, Message_t m){  os << m; });

  if (comment_mode) cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntil(Time_t(300000, 1));

  auto elapsed = chrono::duration_cast< chrono::duration< double, ratio<1> > > (hclock_t::now() - start).count();

  if (comment_mode) cout << "Simulation took:" << elapsed << "sec" << endl;

  return 0;
}
*/