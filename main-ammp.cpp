#include <iostream>
#include <map>
#include <chrono>
#include <algorithm>
#include <utility>
#include <vector>
#include <boost/simulation.hpp>
#include "atomic-models/reaction.hpp"
#include "atomic-models/controler.hpp"
#include "data-structures/reaction_input.hpp"
#include "atomic-models/filter.hpp"
#include "tinyXML/tinyxml.h"

#define TIXML_USE_STL

using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;
using namespace std;


/*************************************
*********type definations*************
*************************************/

/********** In Out types ************/
typedef chrono::high_resolution_clock hclock;
typedef initializer_list< shared_ptr< model<double> > > modelList;
typedef initializer_list< pair< shared_ptr< model<double> >, shared_ptr< model<double> > > > modelPairList;

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** Vector types ************/
typedef vector<molecule> vectorOfMolecules;

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
  modelList models  = {step1, step2, step3, glyco_ctrl, filter1, filter2, filter3};
  modelList eic     = {glyco_ctrl};
  modelPairList ic  = {{step1, glyco_ctrl}, {step2, glyco_ctrl}, {step3, glyco_ctrl}, {glyco_ctrl, filter1}, {glyco_ctrl, filter2}, {glyco_ctrl, filter3}, {filter1, step1}, {filter2, step2}, {filter3, step3}};
  modelList eoc     = {glyco_ctrl};

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

//esta es una prueba de como parsear el imput, me falta interpretar los mat. voy a crear una function que haga eso

int main(int argc, char* argv[]) {

    int count;
    string elemName;
    TiXmlDocument doc;


    TiXmlElement * root, * model, * list, * listOfReactions;
    if (argc <= 1){
        cout << "one argument must to be defined" << endl;
        exit(1);
    }


    if ( !doc.LoadFile(argv[1]) ){
        cout << "fail loading" << endl;
        exit(1);
    }

    root = doc.FirstChildElement();
    model = root->FirstChildElement();

    for(list = model->FirstChildElement(); list != NULL; list = list->NextSiblingElement()){
        elemName = list->Value();
        if (elemName == "listOfReactions"){
            listOfReactions = list;
        }
    }
    list = listOfReactions->FirstChildElement();

    count = 0;
    for(list = listOfReactions->FirstChildElement(); list != NULL; list = list->NextSiblingElement()){
        ++count;
    }

    cout << "Cantidad de reacciones: " << count << endl;
  return 0;
}


