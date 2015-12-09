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
#include "model_generator/model_generator.hpp"

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

  if (argc <= 1){
    cout << "An SBML file is required." << endl;
    exit(1);
  }

  ModelGenerator<Time_t,Message_t> mg(argv[1], Time_t(1,1), Time_t(1,1000), Time_t(1,100), Time_t(1,1000), Time_t(1,1), Time_t(0,1), true);
  shared_ptr<model<Time_t>>& cell_model = mg.getCellModel();

  /*****************************************************************************************************/
  /************************************** Runing Simulation ********************************************/
  /*****************************************************************************************************/

  //TODO(lao) emprolijar y modularizar este codigo. parametrizar los valores como input
  bool comment_mode = true;

  Parser_t p = mg.getParser();
  map<string, map<string, string>> species =  p.getSpecieByCompartments();

  if (comment_mode) cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "1/10000 e e_input | A_e 100 \n 2/1 e e_input | A_e 200 ";
  //string input = "";
  /*
  string partial_input = "";
  for (Time_t t(1,100); t < Time_t(1,50); t += Time_t(1,100)) {

    for (map<string, map<string, string>>::const_iterator i = species.begin(); i != species.end(); ++i){
      partial_input = t.toString() + " " + i->first + " " + i->first + "_s | ";
      
      for (map<string,string>::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
        input += partial_input + j->first + " 100 \n ";
      }
    }
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  */
  piss->str(input);

  auto pf = make_atomic_ptr<input_stream<Time_t, Message_t, Time_t, Message_t >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    int delimiter, amount;
    string collector, specie;
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

    ss >> specie;
    ss >> amount;

    msg_out.metabolites.insert({specie, amount});
    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });

  auto coupled_pf = mg.makeCoupledModel(pf);


  if (comment_mode) cout << "Creating the model to show the space state from stream" << endl;
  piss = make_shared<istringstream>();
  input = "";

  for (Time_t i(1,1); i <= Time_t(30000); i += Time_t(1,1)) {

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

    msg_out.to = {"e", "c", "p", "p_show_request"};
    msg_out.show_request = true;
    msg_out.biomass_request = false;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });

  auto coupled_so = mg.makeCoupledModel(so);

  if (comment_mode) cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{
    {coupled_pf, coupled_so, cell_model}, 
    {}, 
    {
      {coupled_pf, cell_model},
      {coupled_so, cell_model}
    }, 
    {}
  });

  //pdevs_tools::pdevs_coupling_diagram<Time_t, Message_t> pd{*root};
  //string puml_result = pd.get_plant_uml();
  //cout << puml_result;
  //return 0;
  
  if (comment_mode) cout << "Preparing runner" << endl;
  runner<Time_t, Message_t> r(root, Time_t(0), cout, [](ostream& os, Message_t m){  os << m; });

  Time_t until(30000,1);
  if (comment_mode) cout << "Starting simulation until time " << until << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntil(until);

  auto elapsed = chrono::duration_cast< chrono::duration< double, ratio<1> > > (hclock_t::now() - start).count();

  if (comment_mode) cout << "Simulation took:" << elapsed << "sec" << endl;

  return 0;
}