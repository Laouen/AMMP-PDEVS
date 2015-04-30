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

using Time_t                  = double;
using hclock_t                = chrono::high_resolution_clock;
using vectorOfModels_t        = vector< shared_ptr< model<Time_t> > >;
using vectorOfModelPairs_t    = vector< pair< shared_ptr< model<Time_t> >, shared_ptr< model<Time_t> > > >;
using modelsMap_t             = map<string, shared_ptr< model<Time_t> > >;
using vectorOfCoupledModels_t = vector< shared_ptr< flattened_coupled<Time_t, Message_t> > >;
using coupledModelsMap_t      = map<string, shared_ptr< flattened_coupled<Time_t, Message_t> > >;



/***************************************/
/******** End type definations *********/
/***************************************/

/***************************************/
/********** Helper functions ***********/
/***************************************/

vector<string> getReactants(const enzyme_parameter_t& e) {
  vector<string> result;

  for (SetOfMolecules_t::const_iterator it = e.reactants_sctry.cbegin(); it != e.reactants_sctry.cend(); ++it) {
    result.push_back(it->first);
  }

  if (e.reversible) {
    for (SetOfMolecules_t::const_iterator it = e.products_sctry.cbegin(); it != e.products_sctry.cend(); ++it) {
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

  for (SetOfMolecules_t::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }

  for (SetOfMolecules_t::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {
    
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

  for (SetOfMolecules_t::const_iterator jt = e.reactants_sctry.cbegin(); jt != e.reactants_sctry.cend(); ++jt) {
    
    for (map<string, map<string, string> >::const_iterator i = s.cbegin(); i != s.cend(); ++i) {
      if (i->second.find(jt->first) != i->second.end()) {
        
        curr_space = i->first;
        break;
      }
    }
    compartments[curr_space] += 1;
  }

  if (e.reversible) {

    for (SetOfMolecules_t::const_iterator jt = e.products_sctry.cbegin(); jt != e.products_sctry.cend(); ++jt) {
      
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

  cout << "Creating space addresses for use in the reaction atomic models." << endl;
  shared_ptr< map<string, Address_t> > species_addresses = make_shared< map<string, Address_t> >();
  Address_t new_address;

  for (map<string, map<string, string> >::const_iterator i = species.cbegin(); i != species.end(); ++i) {
    
    new_address.clear();
    new_address.push_back(i->first);
    for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
      
      new_address.push_back(i->first + "_s");
      species_addresses->emplace(j->first, new_address);
      new_address.pop_back();
    }
  }

  cout << "Creating reaction atomic models and ordering them by places." << endl;
  map< string, modelsMap_t >  reaction_models;
  Time_t                      interval_time, rate;
  Integer_t                   amount;
  bool                      is_not_special;

  // setting all the existent organelle places in the map
  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    is_not_special = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (is_not_special) {
      reaction_models[it->first + "_i"] = {};
      reaction_models[it->first + "_m"] = {};
    }
  }

  // setting the special places extra cellular, periplasm and cytoplasm places in the map
  reaction_models[special_places[0] + "_i"]   = {};
  reaction_models[special_places[2] + "_i"]   = {};
  reaction_models[special_places[1] + "_i"]   = {};
  reaction_models[special_places[1] + "_lm"]  = {};
  reaction_models[special_places[1] + "_um"]  = {};
  reaction_models[special_places[1] + "_tm"]  = {};



  // creating the reactions atomic models
  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    interval_time   = 0.0001;
    rate            = 0.001;
    amount          = 3;

    place = getPlace(it->second, species, compartements, special_places);

    reaction_models.at(place)[it->first] = make_atomic_ptr< 
      reaction<Time_t, Message_t>, 
      const string, 
      const shared_ptr< map<string, Address_t> >, 
      const bool, 
      const Time_t, 
      const SetOfMolecules_t&, 
      const SetOfMolecules_t&, 
      const Integer_t, 
      const Time_t >(
        it->first, 
        species_addresses, 
        it->second.reversible, 
        rate, 
        it->second.reactants_sctry, 
        it->second.products_sctry, 
        amount, interval_time
      );
  }

  cout << "Creating enzyme coupled models with filters." << endl;
  map< string, coupledModelsMap_t >  enzyme_models;

  for (map<string, modelsMap_t>::iterator it = reaction_models.begin(); it != reaction_models.end(); ++it) {    
    for (modelsMap_t::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {

      auto enzyme_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(jt->first);
      enzyme_models[it->first][jt->first] = make_shared<
        flattened_coupled<Time_t, Message_t>,
        vector<shared_ptr<model<Time_t> > >,
        vector<shared_ptr<model<Time_t> > >,
        vector<pair<shared_ptr<model<Time_t>>, shared_ptr<model<Time_t> > > >,
        vector<shared_ptr<model<Time_t> > > 
      >(
        {enzyme_filter, jt->second}, 
        {enzyme_filter}, 
        {{enzyme_filter, jt->second}}, 
        {jt->second}
      );

    }
  }
  

  cout << "Creating enzyme addresses." << endl;
  map<string, Address_t> enzyme_addresses;

  for (map<string, coupledModelsMap_t>::const_iterator it = enzyme_models.cbegin(); it != enzyme_models.cend(); ++it) {
    
    place       = it->first.substr(0, it->first.find_last_of("_"));
    sub_place   = it->first.substr(it->first.find_last_of("_") + 1);
    
    new_address.clear();
    new_address.push_back(place);
    new_address.push_back(it->first);
    if (sub_place == "i") new_address.push_back(place + "_s");

    for (coupledModelsMap_t::const_iterator jt = it->second.cbegin(); jt != it->second.cend(); ++jt) {

      new_address.push_back(jt->first);
      enzyme_addresses[jt->first] = new_address;
      new_address.pop_back();
    }
  }


  cout << "Getting enzyme information for spaces creation." << endl;
  map<string, map<string, enzyme_info_t> > enzyme_informations;
  enzyme_info_t new_enzyme;
  vector<string> compartments_who_use_the_enzyme;

  // setting the organelle places in the map
  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    is_not_special = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (is_not_special) enzyme_informations[it->first] = {};
  }

  // setting the special places extra cellular, periplasm and cytoplasm places in the map
  enzyme_informations[special_places[0]] = {};
  enzyme_informations[special_places[2]] = {};
  enzyme_informations[special_places[1]] = {};

  // saving the enzyme informations in the map
  for (map<string, enzyme_parameter_t >::const_iterator it = reactions.cbegin(); it != reactions.end(); ++it) {
    
    new_enzyme.location   = enzyme_addresses[it->first];
    new_enzyme.reactants  = getReactants(it->second); 

    compartments_who_use_the_enzyme = getCompartments(it->second, species, compartements, special_places);

    for (vector<string>::iterator i = compartments_who_use_the_enzyme.begin(); i != compartments_who_use_the_enzyme.end(); ++i) {
      
      enzyme_informations.at(*i)[it->first] = new_enzyme;
    }
  }


  cout << "Creating space atomic models." << endl;
  modelsMap_t space_models;
  double volume, factor;
  map<string, metabolite_info_t> metabolites;


  for (map<string, string>::const_iterator it = compartements.cbegin(); it != compartements.cend(); ++it) {
    
    interval_time = 0.5;
    volume        = 5;
    factor        = 1;

    space_models[it->first] = make_atomic_ptr< 
      space<Time_t, Message_t>, 
      const string, 
      const Time_t,
      const map<string, metabolite_info_t>&, 
      const map<string, enzyme_info_t>&, 
      const double, 
      const double >(
        it->first, 
        interval_time, 
        metabolites, 
        enzyme_informations[it->first], 
        volume, 
        factor
      );
  }


  cout << "Creating compartment coupled models" << endl;
  coupledModelsMap_t compartment_models;
   
  for (modelsMap_t::iterator it = space_models.begin(); it != space_models.end(); ++it) {

    auto compartment_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(it->first + "_s");
    compartment_models[it->first] = make_shared<
      flattened_coupled<Time_t, Message_t>,
      vector<shared_ptr<model<Time_t> > >,
      vector<shared_ptr<model<Time_t> > >,
      vector<pair<shared_ptr<model<Time_t>>, shared_ptr<model<Time_t>> > >,
      vector<shared_ptr<model<Time_t> > > 
    >(
      {compartment_filter, it->second}, 
      {compartment_filter}, 
      {{compartment_filter, it->second}}, 
      {it->second}
    );
  }


  cout << "Creating enzyme set coupled models with filters." << endl;
  coupledModelsMap_t  enzyme_set_models;
  vector<shared_ptr<model<Time_t> > > models, eic, eoc;
  vectorOfModelPairs_t ic;

  for (map<string, coupledModelsMap_t>::iterator it = enzyme_models.begin(); it != enzyme_models.end(); ++it) {    
    
    models.clear();
    eic.clear();
    eoc.clear();
    ic.clear();

    auto enzyme_set_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(it->first);
    models.push_back(enzyme_set_filter);
    eic.push_back(enzyme_set_filter);
    
    for (coupledModelsMap_t::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
      models.push_back(jt->second);
      eoc.push_back(jt->second);
      ic.push_back({enzyme_set_filter, jt->second});
    }

    shared_ptr<flattened_coupled<Time_t, Message_t>> new_enzyme_set(new flattened_coupled<Time_t, Message_t>(models, eic, ic, eoc));
    enzyme_set_models[it->first] = new_enzyme_set;
  }


  cout << "Creating cytoplasm bulk solution coupled model." << endl;
  auto cytoplasm_filter     = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[2]);
  auto cytoplasm_or_filter  = make_atomic_ptr< filter<Time_t, Message_t>, const string>("show_state");
  auto cytoplasm_os_filter  = make_atomic_ptr< filter<Time_t, Message_t>, const string>("sending_output");
  auto cytoplasm_space      = compartment_models.at(special_places[2]);
  auto cytoplasm_inner      = enzyme_set_models.at(special_places[2] + "_i");
  
  shared_ptr<flattened_coupled<Time_t, Message_t>> cytoplasm_cos_filter(new flattened_coupled<Time_t, Message_t>(
    {cytoplasm_os_filter}, 
    {cytoplasm_os_filter}, 
    {}, 
    {cytoplasm_os_filter}
  ));

  shared_ptr<flattened_coupled<Time_t, Message_t>> cytoplasm_model(new flattened_coupled<Time_t, Message_t>(
    {cytoplasm_filter, cytoplasm_or_filter, cytoplasm_cos_filter, cytoplasm_space, cytoplasm_inner}, 
    {cytoplasm_filter, cytoplasm_or_filter}, 
    {
      {cytoplasm_or_filter, cytoplasm_space}, 
      {cytoplasm_space, cytoplasm_cos_filter}, 
      {cytoplasm_filter, cytoplasm_space}, 
      {cytoplasm_space, cytoplasm_inner}, 
      {cytoplasm_inner, cytoplasm_space}
    }, 
    {cytoplasm_cos_filter}
  ));

  cout << "Creating extra cellular bulk solution coupled model." << endl;
  auto extra_cellular_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[0]);
  auto extra_cellular_space  = compartment_models.at(special_places[0]);
  auto extra_cellular_inner  = enzyme_set_models.at(special_places[0] + "_i");
  shared_ptr<flattened_coupled<Time_t, Message_t>> extra_cellular_model(new flattened_coupled<Time_t, Message_t>(
    {extra_cellular_filter, extra_cellular_space, extra_cellular_inner}, 
    {extra_cellular_filter}, 
    {{extra_cellular_filter, extra_cellular_space}, {extra_cellular_space, extra_cellular_inner}, {extra_cellular_inner, extra_cellular_space}}, 
    {extra_cellular_space}
  ));

  cout << "Creating periplasm coupled model." << endl;
  auto periplasm_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[1]);
  auto periplasm_space  = compartment_models.at(special_places[1]);
  auto trans_membrane   = enzyme_set_models.at(special_places[1] + "_tm");
  auto outer_membrane   = enzyme_set_models.at(special_places[1] + "_um");
  auto inner_membrane   = enzyme_set_models.at(special_places[1] + "_lm");
  auto periplasm_inner  = enzyme_set_models.at(special_places[1] + "_i");
  shared_ptr<flattened_coupled<Time_t, Message_t>> periplasm_model(new flattened_coupled<Time_t, Message_t>(
    {periplasm_filter, periplasm_space, trans_membrane, outer_membrane, inner_membrane, periplasm_inner}, 
    {extra_cellular_filter}, 
    {{extra_cellular_filter, trans_membrane}, {extra_cellular_filter, outer_membrane}, {extra_cellular_filter, inner_membrane}, {trans_membrane, periplasm_space}, {outer_membrane, periplasm_space}, {inner_membrane, periplasm_space}, {periplasm_space, trans_membrane}, {periplasm_space, outer_membrane}, {periplasm_space, inner_membrane}, {periplasm_space, periplasm_inner}, {periplasm_inner, periplasm_space}}, 
    {trans_membrane, outer_membrane, inner_membrane}
  ));

  cout << "Creating organelles coupled models." << endl;
  vectorOfCoupledModels_t organelle_models;

  for (coupledModelsMap_t::iterator it = compartment_models.begin(); it != compartment_models.end(); ++it) {

    is_not_special = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]);
    
    if(is_not_special) {
      
      auto organelle_filter   = make_atomic_ptr< filter<Time_t, Message_t>, const string>(it->first);
      auto organelle_space    = compartment_models.at(it->first);
      auto organelle_membrane = enzyme_set_models.at(it->first + "_m");
      auto organelle_inner    = enzyme_set_models.at(it->first + "_i");
      organelle_models.push_back( 
        make_shared<
          flattened_coupled<Time_t, Message_t>,
          vector<shared_ptr<model<Time_t> > >,
          vector<shared_ptr<model<Time_t> > >,
          vector<pair<shared_ptr<model<Time_t>>, shared_ptr<model<Time_t>> > >,
          vector<shared_ptr<model<Time_t> > > 
        >(
          {organelle_filter, organelle_membrane, organelle_space, organelle_inner}, 
          {organelle_filter}, 
          {{organelle_filter, organelle_membrane}, {organelle_membrane, organelle_space}, {organelle_space, organelle_inner}, {organelle_inner, organelle_space}, {organelle_space, organelle_membrane}}, 
          {organelle_membrane}
        )
      );
    }
  }

  cout << "Creating the cell coupled model." << endl;

  vectorOfModels_t cell_models = {extra_cellular_model, periplasm_model, cytoplasm_model};
  vectorOfModels_t cell_eic = {extra_cellular_model, cytoplasm_model};
  vectorOfModels_t cell_eoc = {extra_cellular_model, periplasm_model, cytoplasm_model};
  vectorOfModelPairs_t cell_ic{ {extra_cellular_model, periplasm_model}, {periplasm_model, cytoplasm_model}, {cytoplasm_model, periplasm_model}, {periplasm_model, extra_cellular_model} };
  
  for (vectorOfCoupledModels_t::iterator it = organelle_models.begin(); it != organelle_models.end(); ++it) {
    cell_ic.push_back({cytoplasm_model, *it});
    cell_ic.push_back({*it, cytoplasm_model});
  }

  shared_ptr<flattened_coupled<Time_t, Message_t>> cell_model(new flattened_coupled<Time_t, Message_t>(
    cell_models,
    cell_eic,
    cell_ic,
    cell_eoc
  ));

  /*****************************************************************************************************/
  /****************************** Testing cytoplasm coupled model *************************************/
  /*****************************************************************************************************/

  cout << "Testing cytoplasm coupled model with filter" << endl;

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.001; i < 0.002; i += 0.001) {

    input += to_string(i) + "c c_s | A_c 1 \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

  auto pf = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
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


  cout << "Creating the model to show the space state from stream" << endl;
  piss = make_shared<istringstream>();
  input = "";

  for (double i = 0.0; i < 500.0; i += 1.0) {

    input += to_string(i) + " \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

  auto so = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
    string thrash;
    stringstream ss;
    Message_t msg_out;

    ss.str(s);
    ss >> t_next;

    msg_out.to      = {"show_state", "c_s"};
    msg_out.specie  = "";
    msg_out.amount  = 0;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });


  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{
    {pf, so, cytoplasm_model}, 
    {}, 
    {{pf, cytoplasm_model}, {so, cytoplasm_model}}, 
    {cytoplasm_model}
  });

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;


  /*****************************************************************************************************/
  /****************************** Testing enzyme set coupled model *************************************/
  /*****************************************************************************************************/

/*
  cout << "Testing enzyme set coupled models with filters" << endl;
  //shared_ptr< space<Time_t, Message_t> > s = dynamic_pointer_cast< space<Time_t, Message_t> >( space_models["c"] );
  auto enzyme_set = enzyme_set_models.at("p_i");

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.001; i < 0.009; i += 0.001) {

    input += to_string(i) + " p_i R_AGM4PCPpp | M_anhgm4p_p 1 \n ";
    input += to_string(i) + " p_i R_AGM4PCPpp | M_h_p 1 \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

  auto pf = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
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

  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{{pf, enzyme_set}, {}, {{pf, enzyme_set}}, {enzyme_set}});

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;


  /*****************************************************************************************************/
  /*************************** Testing compartment coupled model ***************************************/
  /*****************************************************************************************************/

/*
  cout << "Testing compartment coupled models" << endl;
  //shared_ptr< space<Time_t, Message_t> > s = dynamic_pointer_cast< space<Time_t, Message_t> >( space_models["c"] );
  auto cytoplasm = compartment_models["c"];

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.001; i < 0.109; i += 0.001) {

    input += to_string(i) + "c_s | M_adphep_DASH_DD_c 1 \n ";
    input += to_string(i) + "c_s | M_acon_DASH_T_c 1 \n ";
    //input += to_string(i) + " | M_hdcea_c 1 \n ";
  }
  input.pop_back();
  input.pop_back();
  input.pop_back();
  piss->str(input);
  
  //cout << input << endl;

  auto pf = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
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

  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{{pf, cytoplasm}, {}, {{pf, cytoplasm}}, {cytoplasm}});

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;



  /*****************************************************************************************************/
  /**************************** Testing enzyme coupled model *******************************************/
  /*****************************************************************************************************/
/*

  cout << "Testing enzyme coupled models" << endl;
  //shared_ptr< space<Time_t, Message_t> > s = dynamic_pointer_cast< space<Time_t, Message_t> >( space_models["c"] );
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

  auto pf = make_atomic_ptr<external_events<Time_t, Message_t, Time_t, string >, shared_ptr<istringstream>, Time_t>(piss, Time_t(0),
    [](const string& s, Time_t& t_next, Message_t& m_next)->void{ 

    // Parsing function
    // Intermediary vars for casting
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

  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{{pf, enzyme}, {}, {{pf, enzyme}}, {enzyme}});

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;


  /*****************************************************************************************************/
  /************************** Testing space atomic model ***********************************************/
  /*****************************************************************************************************/



/*
  cout << "Testing spaces atomic models" << endl;
  //shared_ptr< space<Time_t, Message_t> > s = dynamic_pointer_cast< space<Time_t, Message_t> >( space_models["c"] );
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


      msg_out.to.push_back("ID");
    }

    ss >> msg_out.specie;
    ss >> msg_out.amount;

    m_next = msg_out;
    ss >> thrash;
    if ( 0 != thrash.size()) throw exception();
  });

  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{{pf, cytoplasm}, {}, {{pf, cytoplasm}}, {cytoplasm}});

  cout << "Preparing runner" << endl;
  Time_t initial_time{0};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << "To: " << m.to << endl << "specie: " << m.specie << endl << "amount: " << m.amount; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntilPassivate();

  auto elapsed = chrono::duration_cast< chrono::duration< Time_t, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;



  /*****************************************************************************************************/
  /*********************** Testing reaction atomic model ***********************************************/
  /*****************************************************************************************************/

/*
  shared_ptr< reaction<Time_t, Message_t> > m = dynamic_pointer_cast< reaction<Time_t, Message_t> >( reaction_models["c_i"].at("R_2AGPEAT120") );

  Message_t m1;
  Message_t m2;
  Message_t m3;

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
  vector<Message_t> ms = m->out();
  for(vector<Message_t>::iterator i = ms.begin(); i != ms.end(); ++i) {
    cout << *i << endl;
  }
  m->internal();
  m->advance();
  m->show(cout);
  cout << endl;
  ms = m->out();
  for(vector<Message_t>::iterator i = ms.begin(); i != ms.end(); ++i) {
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