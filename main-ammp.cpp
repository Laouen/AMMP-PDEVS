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

typedef chrono::high_resolution_clock hclock;
typedef initializer_list< shared_ptr< model<double> > > modelsList;
typedef initializer_list< pair< shared_ptr< model<double> >, shared_ptr< model<double> > > > pairOfModelsList;


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

  cout << "Creating the atomic models for the filters:" << endl;
  cout << "enzyme set: Membrane -> cytoplasm: Cytoplasm ->compartment: Cytoplasm space" << endl;
  auto membrane_filter = make_atomic_ptr< filter<double>, string, string >("enzyme set", "Membrane");
  auto cytoplasm_filter = make_atomic_ptr< filter<double>, string, string >("cytoplasm", "Cytoplasm");
  auto cytoplasm_space_filter = make_atomic_ptr< filter<double>, string, string >("compartment", "Cytoplasm space");

  cout << "Creating the atomic models for the filters that should not let the message pass:" << endl;
  cout << "4 enzyme set type, 4 cytoplasm type, 4 compartment type, 4 enzyme types, 4 organelle types" << endl;
  auto mf1 = make_atomic_ptr< filter<double>, string, string >("enzyme set", "Inner");
  auto mf2 = make_atomic_ptr< filter<double>, string, string >("enzyme set", "membrane");
  auto mf3 = make_atomic_ptr< filter<double>, string, string >("enzyme set", "Cytoplasm space");
  auto mf4 = make_atomic_ptr< filter<double>, string, string >("enzyme set", "");
  auto cf1 = make_atomic_ptr< filter<double>, string, string >("cytoplasm", "cyto");
  auto cf2 = make_atomic_ptr< filter<double>, string, string >("cytoplasm", "cytoplasm");
  auto cf3 = make_atomic_ptr< filter<double>, string, string >("cytoplasm", "Membrane");
  auto cf4 = make_atomic_ptr< filter<double>, string, string >("cytoplasm", "");
  auto csf1 = make_atomic_ptr< filter<double>, string, string >("compartment", "Cytoplasm-space");
  auto csf2 = make_atomic_ptr< filter<double>, string, string >("compartment", "cytoplasm space");
  auto csf3 = make_atomic_ptr< filter<double>, string, string >("compartment", "Cytoplasm");
  auto csf4 = make_atomic_ptr< filter<double>, string, string >("compartment", "");
  auto e1 = make_atomic_ptr< filter<double>, string, string >("enzyme", "Membrane");
  auto e2 = make_atomic_ptr< filter<double>, string, string >("enzyme", "Cytoplasm");
  auto e3 = make_atomic_ptr< filter<double>, string, string >("enzyme", "Cytoplasm space");
  auto e4 = make_atomic_ptr< filter<double>, string, string >("enzyme", "enzym");
  auto o1 = make_atomic_ptr< filter<double>, string, string >("organelle", "Membrane");
  auto o2 = make_atomic_ptr< filter<double>, string, string >("organelle", "Cytoplasm");
  auto o3 = make_atomic_ptr< filter<double>, string, string >("organelle", "Cytoplasm space");
  auto o4 = make_atomic_ptr< filter<double>, string, string >("organelle", "organ");

  cout << "Coupling the models into the filter_test_model" << endl;
  modelsList models   = {membrane_filter, cytoplasm_filter, cytoplasm_space_filter, mf1, mf2, mf3, mf4, cf1, cf2, cf3, cf4, csf1, csf2, csf3, csf4, e1, e2, e3, e4, o1, o2, o3, o4};
  modelsList eic      = {membrane_filter, e1, mf2, cf3, csf4, o1, e2, mf3, cf4, csf1, csf2, csf3};
  pairOfModelsList ic = {{membrane_filter, cytoplasm_filter}, {cytoplasm_filter, cytoplasm_space_filter} , {membrane_filter, e1}, {membrane_filter, e2}, {membrane_filter, e3}, {membrane_filter, e4}, {membrane_filter, o1}, {membrane_filter, o2}, {membrane_filter, o3}, {membrane_filter, o4}, {cytoplasm_space_filter, mf1}, {cytoplasm_space_filter, mf2}, {cytoplasm_space_filter, mf3}, {cytoplasm_space_filter, mf4}, {cytoplasm_filter, cf1}, {cytoplasm_filter, cf2}, {cytoplasm_filter, cf3}, {cytoplasm_filter, cf4}};
  modelsList eoc      = {membrane_filter, cytoplasm_filter, cytoplasm_space_filter, mf1, mf2, mf3, mf4, cf1, cf2, cf3, cf4, csf1, csf2, csf3, csf4, e1, e2, e3, e4, o1, o2, o3, o4};

  shared_ptr< flattened_coupled<double, Message> > filter_test_model( new flattened_coupled<double, Message>{models, eic, ic, eoc});

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  piss->str("1 {enzyme set,Membrane} {cytoplasm,Cytoplasm} {compartment,Cytoplasm space} | Must-pass-all-the-filters 2 \n 2 {cytoplasm,Cytoplasm} | Musn't-pass-all-the-filters 3 \n 3 {compartment,Cytoplasm space} | Musn't-pass-all-the-filters 4 \n 4 {enzyme set,Membrane} {cytoplasm,Cytoplasm} {compartment,} | Must-pass-the-first-two-filters 5 \n 5 {enzyme set,Membrane} {cytoplasm,wron name} {compartment,Cytoplasm space} | Must-pass-only-the-first-filters 6 \n 5 {enzyme set,Membrane} {cytoplasm,wron name} {compartment,Cytoplasm space} | Musnt-pass-the-filters 6");
  
  auto pf = make_atomic_ptr<external_events<double, Message, double, string >, shared_ptr<istringstream>, double>(piss, double(0),
    [](const string& s, double& t_next, Message& m_next)->void{ 

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
  shared_ptr< flattened_coupled<double, Message> > root( new flattened_coupled<double, Message>{{pf, filter_test_model}, {}, {{pf, filter_test_model}}, {filter_test_model}});

  cout << "Preparing runner" << endl;
  double initial_time{0};
  runner<double, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << m; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< double, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;

  /**************************************************************************************************************/
  /************************************** End testing filters ***************************************************/
  /**************************************************************************************************************/

  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/


/*
int main(int argc, char* argv[]){

  if (argc < 2){
    cout << "show usage" << endl;
    exit(1);
  }

  cout << "Preparing reactions parameters" << endl;
  reaction_input step1_input("step1", 0, {"atpc", "glucosec"}, {"g6pc", "h", "adp"}, {"ifhec"});
  reaction_input step2_input("step2", 0, {"g6pc"}, {"f6pc"}, {"ifpgisomerase"});
  reaction_input step3_input("step3", 0, {"atpc", "f6pc"}, {"fructose_16_bisphosphate", "adp"}, {"ifpfk"});
  reaction_input step4_input("step4", 0, {"fructose_16_bisphosphate"}, {"dhp"}, {"aldolase"});
  reaction_input step5_input("step5", 0, {"fructose_16_bisphosphate"}, {"gdp"}, {"aldolase"});
  reaction_input step6_input("step6", 0, {"gdp", "nadc", "pc"}, {"nadc", "h", "_13bpg"}, {"g3pd"});
  reaction_input step7_input("step7", 0, {"_13bpg", "adp"}, {"_3_phosphoglycerate", "atpc"}, {"pgk"});
  reaction_input step8_input("step8", 0, {"_3_phosphoglycerate"}, {"_2_phosphoglycerate"}, {"pgm"});
  reaction_input step9_input("step9", 0, {"_2_phosphoglycerate"}, {"pepc", "h2o"}, {"enolase"});
  reaction_input step10_input("step10", 0, {"pepc", "adp"}, {"pyruvate", "atp"}, {"pyruvate_kinase"});
  reaction_input step4to5_input("step4to5", 0, {"dhp"}, {"gdp"}, {"isomerase"});

  cout << "Creating the atomic models for the 3 steps" << endl;
  auto step1 = make_atomic_ptr< reaction<double, boost::any>, reaction_input >(step1_input);
  auto step2 = make_atomic_ptr< reaction<double, boost::any>, reaction_input >(step2_input);
  auto step3 = make_atomic_ptr< reaction<double, boost::any>, reaction_input >(step3_input);

  cout << "Creating the atomic models for the 3 filters" << endl;
  auto filter1 = make_atomic_ptr< filter<double, boost::any>, string >("step1");
  auto filter2 = make_atomic_ptr< filter<double, boost::any>, string >("step2");
  auto filter3 = make_atomic_ptr< filter<double, boost::any>, string >("step3");

  cout << "Creating the atomic model for the controler" << endl;
  auto glyco_ctrl = make_atomic_ptr< controler<double, boost::any>, vector<reaction_input> >({step1_input, step2_input, step3_input});


  cout << "Coupling the models into the glyco" << endl;
  modelsList models   = {step1, step2, step3, glyco_ctrl, filter1, filter2, filter3};
  modelsList eic      = {glyco_ctrl};
  pairOfModelsList ic = {{step1, glyco_ctrl}, {step2, glyco_ctrl}, {step3, glyco_ctrl}, {glyco_ctrl, filter1}, {glyco_ctrl, filter2}, {glyco_ctrl, filter3}, {filter1, step1}, {filter2, step2}, {filter3, step3}};
  modelsList eoc      = {glyco_ctrl};

  shared_ptr< flattened_coupled<double, boost::any> > glyco( new flattened_coupled<double, boost::any>{models, eic, ic, eoc});

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  piss->str("0 glucosec 3 \n 2 atpc 6 \n 5 ifhec 1 \n 6 ifpfk 1 \n 7 ifpgisomerase 1");
  auto pf = make_atomic_ptr<external_events<double, boost::any, double, pair<string, double> >, shared_ptr<istringstream>, double>(piss, double(0),
    [](const string& s, double& t_next, boost::any& m_next)->void{ //parsing function
      //intermediary vars for casting
      double tmp_next;
      double amoutn_next_out; 
      string element_next_out;
      stringstream ss;
      ss.str(s);
      ss >> tmp_next;
      t_next = tmp_next;
      ss >> element_next_out;
      ss >> amoutn_next_out;
      m_next = static_cast<boost::any>(make_pair(element_next_out, amoutn_next_out));
      string thrash;
      ss >> thrash;
      if ( 0 != thrash.size()) throw exception();
    });

  cout << "Coupling the glyco to the input" << endl;
  shared_ptr< flattened_coupled<double, boost::any> > root( new flattened_coupled<double, boost::any>{{pf, glyco}, {}, {{pf, glyco}}, {glyco}});

  cout << "Preparing runner" << endl;
  double initial_time{0};
  runner<double, boost::any> r(root, initial_time, cout, 
    [](ostream& os, boost::any m){ 
        
        pair<string, molecule> msg = boost::any_cast< pair<string, molecule> >(m);
        os << "send To: " << msg.first << " - Molecule: " << msg.second.first << " - Amount: " << msg.second.second;
        //molecule msg = boost::any_cast< pair<string, double> >(m);
        //os << "Molecule: " << msg.first << " - amount: " << msg.second;
    });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< double, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;
  return 0;

}
*/


