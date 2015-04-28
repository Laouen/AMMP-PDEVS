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
#include <memory>

// Boost simalator include
#include <boost/simulation.hpp>

// atomic model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"
#include "atomic-models/space.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/types.hpp"
#include "data-structures/randomNumbers.hpp"

// tinyXML parser
#include "parser/parser.hpp"

#define TIXML_USE_STL

using namespace std;
using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;


/***************************************/
/********* Type definations ************/
/***************************************/

typedef double Time;
typedef chrono::high_resolution_clock hclock;
typedef vector< shared_ptr< model<Time> > > vectorOfModels;
typedef vector< pair< shared_ptr< model<Time> >, shared_ptr< model<Time> > > > vectorOfModelPairs;
typedef map<string, shared_ptr< model<Time> > > modelsMap;
typedef map<string, shared_ptr< flattened_coupled<Time, Message> > > coupledModelsMap;



/***************************************/
/******** End type definations *********/
/***************************************/

/***************************************/
/********** Helper functions ***********/
/***************************************/

vector<string> getReactants(const enzyme_parameter_t& e) {
  vector<string> result;

  for (SetOfMolecules::const_iterator it = e.reactants_sctry.cbegin(); it != e.reactants_sctry.cend(); ++it) {
    result.push_back(it->first);
  }

  if (e.reversible) {
    for (SetOfMolecules::const_iterator it = e.products_sctry.cbegin(); it != e.products_sctry.cend(); ++it) {
      result.push_back(it->first);
    }
  }
  return result;
}

string compartmentOfReactants(const map<string, int>& c) {

  string result;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    if (i->second > 0) {
      result = i->first;
      break;
    }
  }

  return result;
}

string organelleOfReactants(const map<string, int>& c, string* sp) {

  string result;
  bool is_not_special;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    is_not_special = (i->first != sp[0]) && (i->first != sp[1]) && (i->first != sp[2]);
    
    if (is_not_special && (i->second > 0)) {
      result = i->first;
      break;
    }
  }

  return result;
}

bool thereAreOrganelleInvolved(const map<string, int>& c, string* sp) {

  bool result = false;
  bool is_not_special;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    is_not_special = (i->first != sp[0]) && (i->first != sp[1]) && (i->first != sp[2]);
    
    if (is_not_special && (i->second > 0)) {
      result = true;
      break;
    }
  }

  return result;
}

string getPlace(const enzyme_parameter_t& e, const map<string, map<string, string> >& s, const map<string, string>& c, string* sp) {

  string result;
  int amount_compartments;
  string curr_space;
  map<string, int> compartments;

  for (map<string, string>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    compartments[i->first] = 0;
  }

  for (SetOfMolecules::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }

  for (SetOfMolecules::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }


  amount_compartments = 0;
  for (map<string, int>::iterator i = compartments.begin(); i != compartments.end(); ++i) {   
    if (i->second > 0) ++amount_compartments;
  }

  switch(amount_compartments) {
    case 1:

      result = compartmentOfReactants(compartments) + "_i";
      break;
    case 2:

      if (thereAreOrganelleInvolved(compartments, sp)) {
        result = organelleOfReactants(compartments, sp) + "_m";
      
      } else if (compartments[sp[0]] > 0){

        result = sp[1] + "_um";
      } else if (compartments[sp[2]] > 0) {

        result = sp[1] + "_lm";
      }
      break;
    case 3:

      result = sp[1] + "_tm";
      break;
  }

  return result;
}

vector<string> getCompartments(const enzyme_parameter_t& e, const map<string, map<string, string> >& s, const map<string, string>& c, string* sp) {

  vector<string> result;
  string curr_space;
  map<string, int> compartments;

  for (map<string, string>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    compartments[i->first] = 0;
  }

  for (SetOfMolecules::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }

  if (e.reversible) {

    for (SetOfMolecules::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {
      
      for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
        if (i->second.find(jt->first) != i->second.end()) {
          
          curr_space = i->first;
          break;
        }
      }
      compartments[curr_space] += 1;
    }
  }

  for (map<string, int>::iterator i = compartments.begin(); i != compartments.end(); ++i) {   
    
    if (i->second > 0) result.push_back(i->first);
  }

  return result;
}

/***************************************/
/******** End helper functions *********/
/***************************************/

int main(int argc, char* argv[]) {


  /**************************************************************************************************************/
  /************************************* Parsing the SBML file **************************************************/
  /**************************************************************************************************************/

  if (argc <= 1){
      cout << "An SBML file is required." << endl;
      exit(1);
  }

  Parser_t input_doc(argv[1]);
  input_doc.loadFile();

  map<string, string>               compartements     = input_doc.getCompartments();
  map<string, map<string, string> > species           = input_doc.getSpeciesByCompartment();
  map<string, enzyme_parameter_t >  reactions         = input_doc.getReactions();
  string                            special_places[3]  = {"e", "p", "c"};
  string                            place, sub_place;

  cout << "Creating space addresses." << endl;
  shared_ptr< map<string, Address> > species_addresses = make_shared< map<string, Address> >();
  Address new_address;

  for (map<string, map<string, string> >::const_iterator i = species.cbegin(); i != species.end(); ++i) {
    
    new_address.clear();
    new_address.push_back(i->first);
    for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
      
      new_address.push_back(i->first + "_s");
      species_addresses->emplace(j->first, new_address);
      new_address.pop_back();
    }
  }

  cout << "Creating enzymes atomic models and ordering them by places." << endl;
  map< string, modelsMap >  reaction_models;
  Time                      interval_time, rate;
  Integer                   amount;
  bool                      isCompartment;

  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    isCompartment = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (isCompartment) {
      reaction_models[it->first + "_i"] = {};
      reaction_models[it->first + "_m"] = {};
    }
  }

  reaction_models[special_places[0] + "_i"] = {};
  reaction_models[special_places[2] + "_i"] = {};
  reaction_models[special_places[1] + "_i"] = {};
  reaction_models[special_places[1] + "_lm"] = {};
  reaction_models[special_places[1] + "_um"] = {};
  reaction_models[special_places[1] + "_tm"] = {};



  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    interval_time   = 0.1;
    rate            = 0.01;
    amount          = 3;

    place = getPlace(it->second, species, compartements, special_places);

    reaction_models.at(place)[it->first] = make_atomic_ptr< reaction<Time, Message>, string, shared_ptr< map<string, Address> >, bool, Time, SetOfMolecules&, SetOfMolecules&, Integer, Time >(it->first, species_addresses, it->second.reversible, rate, it->second.reactants_sctry, it->second.products_sctry, amount, interval_time);
  }

  cout << "Creating enzyme coupled models with filters." << endl;
  map< string, coupledModelsMap >  enzyme_models;

  for (map<string, modelsMap>::iterator it = reaction_models.begin(); it != reaction_models.end(); ++it) {    
    for (modelsMap::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {

      auto new_filter = make_atomic_ptr< filter<Time, Message>, string>(jt->first);
      enzyme_models[it->first][jt->first] = make_shared<
              flattened_coupled<Time, Message>,
              vector<shared_ptr<model<Time> > >,
              vector<shared_ptr<model<Time> > >,
              vector<pair<shared_ptr<model<Time>>, shared_ptr<model<Time> > > >,
              vector<shared_ptr<model<Time> > > 
            >(
              {new_filter, jt->second}, 
              {new_filter}, 
              {{new_filter, jt->second}}, 
              {jt->second}
            );

    }
  }
  

  cout << "Creating enzyme addresses." << endl;
  map<string, Address> reaction_addresses;

  for (map<string, modelsMap>::const_iterator it = reaction_models.cbegin(); it != reaction_models.cend(); ++it) {
    

    place       = it->first.substr(0, it->first.find_last_of("_"));
    sub_place   = it->first.substr(it->first.find_last_of("_") + 1);

    
    new_address.clear();
    new_address.push_back(place);
    new_address.push_back(it->first);
    if (sub_place == "i") new_address.push_back(place + "_s");

    for (modelsMap::const_iterator jt = it->second.cbegin(); jt != it->second.cend(); ++jt) {

      new_address.push_back(jt->first);
      reaction_addresses[jt->first] = new_address;
      new_address.pop_back();
    }
  }


  cout << "Getting enzyme information for spaces creation." << endl;
  map<string, map<string, enzyme_info_t> > enzyme_informations;
  enzyme_info_t new_enzyme;
  vector<string> compartments_who_use_the_enzyme;

  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    isCompartment = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (isCompartment) {
      enzyme_informations[it->first] = {};
    }
  }

  enzyme_informations[special_places[0]] = {};
  enzyme_informations[special_places[2]] = {};
  enzyme_informations[special_places[1]] = {};


  for (map<string, enzyme_parameter_t >::const_iterator it = reactions.cbegin(); it != reactions.end(); ++it) {
    
    new_enzyme.location   = reaction_addresses[it->first];
    new_enzyme.reactants  = getReactants(it->second); 

    compartments_who_use_the_enzyme = getCompartments(it->second, species, compartements, special_places);

    for (vector<string>::iterator i = compartments_who_use_the_enzyme.begin(); i != compartments_who_use_the_enzyme.end(); ++i) {
      
      enzyme_informations.at(*i)[it->first] = new_enzyme;
    }
  }


  cout << "Creating space atomic models." << endl;
  modelsMap space_models;
  double volume, factor;
  map<string, metabolite_info_t> metabolites;


  for (map<string, string>::const_iterator it = compartements.cbegin(); it != compartements.cend(); ++it) {
    
    interval_time = 0.05;
    volume        = 5;
    factor        = 1;

    space_models[it->first] = make_atomic_ptr< space<Time, Message>, string, Time, map<string, metabolite_info_t>&, map<string, enzyme_info_t>&, double, double >(it->first, interval_time, metabolites, enzyme_informations[it->first], volume, factor);
  }



  /*****************************************************************************************************/
  /**************************** Testing enzyme coupled model *******************************************/
  /*****************************************************************************************************/


  cout << "Testing spaces atomic models" << endl;
  //shared_ptr< space<Time, Message> > s = dynamic_pointer_cast< space<Time, Message> >( space_models["c"] );
  auto enzyme = enzyme_models.at("p_lm").at("R_F6Pt6_2pp");

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.001; i < 0.009; i += 0.001) {

    input += to_string(i) + " R_F6Pt6_2pp | M_pi_c 1 \n ";
    input += to_string(i) + " R_F6Pt6_2pp | M_f6p_p 1 \n ";
    //input += to_string(i) + " | M_hdcea_c 1 \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

  auto pf = make_atomic_ptr<external_events<Time, Message, Time, string >, shared_ptr<istringstream>, Time>(piss, Time(0),
    [](const string& s, Time& t_next, Message& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
    int delimiter;
    string collector;
    string thrash;
    stringstream ss;
    Message msg_out;

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

  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time, Message> > root( new flattened_coupled<Time, Message>{{pf, enzyme}, {}, {{pf, enzyme}}, {enzyme}});

  cout << "Preparing runner" << endl;
  Time initial_time{0};
  runner<Time, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;


  /*****************************************************************************************************/
  /************************** Testing space atomic model ***********************************************/
  /*****************************************************************************************************/



/*
  cout << "Testing spaces atomic models" << endl;
  //shared_ptr< space<Time, Message> > s = dynamic_pointer_cast< space<Time, Message> >( space_models["c"] );
  auto cytoplasm = space_models["c"];

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.001; i < 0.109; i += 0.001) {

    input += to_string(i) + " | M_adphep_DASH_DD_c 1 \n ";
    input += to_string(i) + " | M_acon_DASH_T_c 1 \n ";
    //input += to_string(i) + " | M_hdcea_c 1 \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

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
  shared_ptr< flattened_coupled<Time, Message> > root( new flattened_coupled<Time, Message>{{pf, cytoplasm}, {}, {{pf, cytoplasm}}, {cytoplasm}});

  cout << "Preparing runner" << endl;
  Time initial_time{0};
  runner<Time, Message> r(root, initial_time, cout, [](ostream& os, Message m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time, ratio<1> > > (hclock::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;



  /*****************************************************************************************************/
  /*********************** Testing reaction atomic model ***********************************************/
  /*****************************************************************************************************/

/*
  shared_ptr< reaction<Time, Message> > m = dynamic_pointer_cast< reaction<Time, Message> >( reaction_models["c_i"].at("R_2AGPEAT120") );

  Message m1;
  Message m2;
  Message m3;

  m1.specie = "M_2agpe120_c";
  m1.amount = 10;
  m2.specie = "M_atp_c";
  m2.amount = 2;
  m3.specie = "M_ddca_c";
  m3.amount = 5;

  cout << endl;
  m->show(cout);
  cout << endl;
  m->external({m1, m2, m3}, 0);
  m->show(cout);
  cout << endl;
  vector<Message> ms = m->out();
  for(vector<Message>::iterator i = ms.begin(); i != ms.end(); ++i) {
    cout << *i << endl;
  }
  m->internal();
  m->advance();
  m->show(cout);
  cout << endl;
  ms = m->out();
  for(vector<Message>::iterator i = ms.begin(); i != ms.end(); ++i) {
    cout << *i << endl;
  }
  m->internal();
  m->advance();
  m->show(cout);
  cout << endl;


  /**************************************************************************************************************/
  /*********************************** End creating unit definitions ********************************************/
  /**************************************************************************************************************/

  /**************************************************************************************************************/
  /*********************************** End parsing the SBML file ************************************************/
  /**************************************************************************************************************/



  /**************************************************************************************************************/
  /************************************** Making the models *****************************************************/
  /**************************************************************************************************************/



  /**************************************************************************************************************/
  /************************************ End making the models ***************************************************/
  /**************************************************************************************************************/

  return 0;
}