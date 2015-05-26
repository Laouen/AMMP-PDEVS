// STD includes
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <utility> /* pair */
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <memory>

//vendors
#include "vendors/britime.hpp"

// Boost simalator include
#include <boost/simulation.hpp>

// atomic model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"
#include "atomic-models/space.hpp"
#include "atomic-models/biomass.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/types.hpp"
#include "data-structures/randomNumbers.hpp"

// tinyXML parser
#include "parser/parser.hpp"

// model engine
#include "model-engine.hpp"

#define TIXML_USE_STL

using namespace std;
using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;


/***************************************/
/********* Type definations ************/
/***************************************/

using Time_t                  = BRITime;
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

vector<string> getReactants_m(const enzyme_parameter_t& e) {
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

string compartmentOfReactants_m(const map<string, int>& c) {

  string result;

  for (map<string, int>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    
    if (i->second > 0) {
      result = i->first;
      break;
    }
  }

  return result;
}

string compartmentOfReactants_m(const map<string, int>& c, string* sp) {

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

bool thereAreOrganelleInvolved_m(const map<string, int>& c, string* sp) {

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

string getPlace_m(const enzyme_parameter_t& e, const map<string, map<string, string> >& s, const map<string, string>& c, string* sp) {

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

      result = compartmentOfReactants_m(compartments) + "_i";
      break;
    case 2:

      if (thereAreOrganelleInvolved_m(compartments, sp)) {
        result = compartmentOfReactants_m(compartments, sp) + "_m";
      
      } else if (compartments[sp[0]] > 0){

        result = sp[1] + "_um";
      } else if (compartments[sp[2]] > 0) {

        result = sp[1] + "_lm";
      } else if ((compartments[sp[0]] > 0) && (compartments[sp[2]] > 0)) {

        result = sp[1] + "_tm";
      }
      break;
    case 3:

      result = sp[1] + "_tm";
      break;
  }

  return result;
}

vector<string> getCompartments_m(const enzyme_parameter_t& e, const map<string, map<string, string> >& s, const map<string, string>& c, string* sp) {

  vector<string> result;
  string curr_space;
  map<string, int> compartments;

  for (map<string, string>::const_iterator i = c.cbegin(); i != c.cend(); ++i) {
    compartments.insert({i->first, 0});
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

  if (argc <= 1){
    cout << "An SBML file is required." << endl;
    exit(1);
  }

  long double cell_weight   = 280 * 1e-15;  
  Integer_t norm_number     = 1;  
  string e                  = "e";
  string c                  = "c";
  string p                  = "p";
  string biomass_ID         = "R_Ec_biomass_iJO1366_WT_53p95M";

  ModelEngine<Time_t, Message_t> m(cell_weight, argv[1], e, c, p, biomass_ID, norm_number);
  m._comment_mode = true;
  
  m.createSpeciesAddresses();

  // adding the organelles compartments
  for (map<string, string> ::const_iterator i = m._compartments.begin(); i != m._compartments.end(); ++i) {
    if (m.isNotSpecial(i->first)) m.addCompartmentForEnzyme(i->first);
  }

  // creating the enzymes models
  m._comment_mode = false;
  for (map<string, enzyme_parameter_t >::const_iterator i = m._reactions.begin(); i != m._reactions.end(); ++i) {
      m.addEnzymeModel(i->first, BRITime(1), BRITime(1,100), Integer_t(40));
  }
  m._comment_mode = true;

  m.createEnzymeAddresses();

  // create compartments models
  for (map<string, string>::const_iterator i = m._compartments.begin(); i != m._compartments.end(); ++i) {
    m.addCompartmentModel(i->first, BRITime(1,100), BRITime(1,100), 0, 1);
  }

  m.createEnzymeSetModels();

  return 0;

  //long double cell_weight   = 280 * 1e-15;  
  //Integer_t norm_number     = 1;  
  string special_places[3]  = {"e", "p", "c"};
  //string biomass_ID         = "R_Ec_biomass_iJO1366_WT_53p95M";

  /**************************************************************************************************************/
  /************************************* Parsing the SBML file **************************************************/
  /**************************************************************************************************************/

  if (argc <= 1){
    cout << "An SBML file is required." << endl;
    exit(1);
  }

  Parser_t input_doc(argv[1], biomass_ID, cell_weight, norm_number);
  input_doc.loadFile();

  map<string, string>               compartements       = input_doc.getCompartments();
  map<string, map<string, string> > species             = input_doc.getSpeciesByCompartment();
  map<string, enzyme_parameter_t >  reactions           = input_doc.getReactions();
  enzyme_parameter_t                biomass_info        = input_doc.getBiomass();



  /**************************************************************************************************************/
  /*********************************** End parsing the SBML file ************************************************/
  /**************************************************************************************************************/




  /**************************************************************************************************************/
  /******************************* Creating the map of species addresses ****************************************/
  /**************************************************************************************************************/

  /*
  Brieff: This map store for each specie, the filter keys to get in the corresponding space where it belong.
  */




  cout << "Creating species addresses for use in the reaction atomic models and biomass atomic model." << endl;
  shared_ptr< map<string, Address_t> > species_addresses = make_shared< map<string, Address_t> >();
  Address_t new_address;

  for (map<string, map<string, string> >::const_iterator i = species.cbegin(); i != species.end(); ++i) {
    
    new_address.clear();
    new_address.push_back(i->first);
    for (map<string, string>::const_iterator j = i->second.cbegin(); j != i->second.cend(); ++j) {
      
      new_address.push_back(i->first + "_s");
      species_addresses->insert({j->first, new_address});
      new_address.pop_back();
    }
  }


  /**************************************************************************************************************/
  /********************************* Setting the reaction and enzyme map ****************************************/
  /**************************************************************************************************************/

  /* 
  Brief: setting the places in the maps reaction_models and enzyme_models.
  */



  cout << "Setting the reaction and enzyme model map and enzyme model map." << endl;
  map< string, modelsMap_t >          reaction_models;
  map< string, coupledModelsMap_t >   enzyme_models;
  bool                                is_not_special;

  // setting all the existent organelle places in the maps
  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    is_not_special = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (is_not_special) {
      reaction_models.insert({it->first + "_i", {}});
      reaction_models.insert({it->first + "_m", {}});
      enzyme_models.insert({it->first + "_i", {}});
      enzyme_models.insert({it->first + "_m", {}});
    }
  }

  // setting the special places extra cellular, periplasm and cytoplasm places in the map
  reaction_models.insert({special_places[0] + "_i", {}});
  reaction_models.insert({special_places[2] + "_i", {}});
  reaction_models.insert({special_places[1] + "_i", {}});
  reaction_models.insert({special_places[1] + "_lm", {}});
  reaction_models.insert({special_places[1] + "_um", {}});
  reaction_models.insert({special_places[1] + "_tm", {}});
  enzyme_models.insert({special_places[0] + "_i", {}});
  enzyme_models.insert({special_places[2] + "_i", {}});
  enzyme_models.insert({special_places[1] + "_i", {}});
  enzyme_models.insert({special_places[1] + "_lm", {}});
  enzyme_models.insert({special_places[1] + "_um", {}});
  enzyme_models.insert({special_places[1] + "_tm", {}});


  /**************************************************************************************************************/
  /******************************* End Setting the reaction and enzyme map **************************************/
  /**************************************************************************************************************/



  /**************************************************************************************************************/
  /********************* Creating the reaction atomic models and enzyme coupled model ***************************/
  /**************************************************************************************************************/

  /* 
  brief: all the reaction atomic model are placed in the correct place where they belong to.
  After that, the enzymes are created by coupling the reaction atomic model with its corresponding filter. 
  */


  cout << "Creating reaction atomic models and ordering them by places." << endl;
  Time_t    interval_time, rate;
  Integer_t amount;
  string    place, sub_place;

  // creating the reactions atomic models
  for (map<string, enzyme_parameter_t >::iterator it = reactions.begin(); it != reactions.end(); ++it) {
    
    interval_time   = BRITime(1,1);
    rate            = BRITime(1,100);
    amount          = 40;

    place = getPlace_m(it->second, species, compartements, special_places);

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

  for (map<string, modelsMap_t>::iterator it = reaction_models.begin(); it != reaction_models.end(); ++it) {    
    for (modelsMap_t::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {

      auto enzyme_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(jt->first);
      enzyme_models.at(it->first)[jt->first] = make_shared<
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

  /**************************************************************************************************************/
  /******************* End reating the reaction atomic models and enzyme coupled model **************************/
  /**************************************************************************************************************/


  /**************************************************************************************************************/
  /****************************** Getting the enzyme coupled models addresses ***********************************/
  /**************************************************************************************************************/

  /* 
  Biref: enzyme_addresses store for each ezyme, the filter keys to get to them.
  */
  

  cout << "Getting enzyme addresses." << endl;
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
      enzyme_addresses.insert({jt->first, new_address});
      new_address.pop_back();
    }
  }

  /**************************************************************************************************************/
  /**************************** End getting the enzyme coupled models addresses *********************************/
  /**************************************************************************************************************/


  /**************************************************************************************************************/
  /*************************** Getting the space parameters for the atomic model ********************************/
  /**************************************************************************************************************/

  /* 
  Brief: The atomic models space need some parameters, this part is responsible to collect those parameters.
  */

  cout << "Getting enzyme information for spaces creation." << endl;
  map<string, map<string, enzyme_info_t> > enzyme_informations;
  enzyme_info_t new_enzyme_information;
  vector<string> compartments_who_use_the_enzyme;

  // setting the organelle places in the map
  for (map<string, string>::iterator it = compartements.begin(); it != compartements.end(); ++it) {
    
    is_not_special = (it->first != special_places[0]) && (it->first != special_places[1]) && (it->first != special_places[2]) ;
    if (is_not_special) enzyme_informations.insert({it->first, {}});
  }

  // setting the special places extra cellular, periplasm and cytoplasm places in the map
  enzyme_informations.insert({special_places[0], {}});
  enzyme_informations.insert({special_places[2], {}});
  enzyme_informations.insert({special_places[1], {}});

  // saving the enzyme informations in the map
  for (map<string, enzyme_parameter_t >::const_iterator it = reactions.cbegin(); it != reactions.end(); ++it) {

    new_enzyme_information.location   = enzyme_addresses[it->first];
    new_enzyme_information.reactants  = getReactants_m(it->second); 
    compartments_who_use_the_enzyme   = getCompartments_m(it->second, species, compartements, special_places);

    for (vector<string>::iterator i = compartments_who_use_the_enzyme.begin(); i != compartments_who_use_the_enzyme.end(); ++i) {
      
      enzyme_informations.at(*i)[it->first] = new_enzyme_information;
    }
  }

  /**************************************************************************************************************/
  /*************************** End getting the space parameters for the atomic model ****************************/
  /**************************************************************************************************************/

  /**************************************************************************************************************/
  /************************* Creating the space atomic and compartment coupled models ***************************/
  /**************************************************************************************************************/

  /* 
  Brief: all the space atomic model are created using the parameters.
  After that, the compartmen are created by coupling the space atomic model with its corresponding filter. 
  */

  cout << "Creating space atomic models." << endl;
  modelsMap_t space_models;
  double volume, factor;
  map<string, metabolite_info_t> metabolites;
  Time_t biomass_rate;
  Address_t biomass_address;
  for (map<string, string>::const_iterator it = compartements.cbegin(); it != compartements.cend(); ++it) {
    
    interval_time   = BRITime(1,100);
    biomass_rate    = BRITime(1,100);
    biomass_address = {"biomass", biomass_ID};
    volume          = 0;
    factor          = 1;

    space_models[it->first] = make_atomic_ptr< 
      space<Time_t, Message_t>, 
      const string, 
      const Time_t,
      const Time_t,
      const Address_t&,
      const map<string, metabolite_info_t>&, 
      const map<string, enzyme_info_t>&, 
      const double, 
      const double >(
        it->first, 
        interval_time,
        biomass_rate,
        biomass_address,
        metabolites, 
        enzyme_informations.at(it->first), 
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
      {
        {compartment_filter, it->second}
      }, 
      {it->second}
    );
  }


  /**************************************************************************************************************/
  /*********************** End Creating the space atomic and compartment coupled models *************************/
  /**************************************************************************************************************/


  /**************************************************************************************************************/
  /********************************** Creaatin Enzyme set coupled model *****************************************/
  /**************************************************************************************************************/



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
    enzyme_set_models.insert({it->first, new_enzyme_set});
  }


  /**************************************************************************************************************/
  /******************************** End creaatin Enzyme set coupled model ***************************************/
  /**************************************************************************************************************/

  /**************************************************************************************************************/
  /**************************************** coupling everything *************************************************/
  /**************************************************************************************************************/


  cout << "Creating cytoplasm coupled model." << endl;
  auto cytoplasm_filter     = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[2]);
  auto cytoplasm_space      = compartment_models.at(special_places[2]);
  auto cytoplasm_inner      = enzyme_set_models.at(special_places[2] + "_i");

  shared_ptr<flattened_coupled<Time_t, Message_t>> cytoplasm_model(new flattened_coupled<Time_t, Message_t>(
    {cytoplasm_filter, cytoplasm_space, cytoplasm_inner}, 
    {cytoplasm_filter}, 
    {
      {cytoplasm_filter, cytoplasm_space}, 
      {cytoplasm_space, cytoplasm_inner}, 
      {cytoplasm_inner, cytoplasm_space}
    }, 
    {cytoplasm_space}
  ));

  cout << "Creating extra cellular coupled model." << endl;
  auto extra_cellular_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[0]);
  auto extra_cellular_space  = compartment_models.at(special_places[0]);
  auto extra_cellular_inner  = enzyme_set_models.at(special_places[0] + "_i");

  shared_ptr<flattened_coupled<Time_t, Message_t>> extra_cellular_model(new flattened_coupled<Time_t, Message_t>(
    {extra_cellular_filter, extra_cellular_space, extra_cellular_inner}, 
    {extra_cellular_filter}, 
    {
      {extra_cellular_filter, extra_cellular_space},
      {extra_cellular_space, extra_cellular_inner},
      {extra_cellular_inner, extra_cellular_space}
    }, 
    {extra_cellular_space}
  ));

  cout << "Creating periplasm coupled model." << endl;
  auto periplasm_filter        = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[1]);
  auto periplasm_or_filter     = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[1] + "_or");
  auto periplasm_br_filter     = make_atomic_ptr< filter<Time_t, Message_t>, const string>(special_places[1] + "_br");

  auto periplasm_output_filter  = make_atomic_ptr< filter<Time_t, Message_t>, const string>("output");
  auto periplasm_biomass_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>("biomass");
  
  // making a trivial coupled model of periplasm_output_filter in order to solve the bug with coupled to atomic
  shared_ptr<flattened_coupled<Time_t, Message_t>> periplasm_output_coupled_filter(new flattened_coupled<Time_t, Message_t>(
    {periplasm_output_filter}, 
    {periplasm_output_filter}, 
    {}, 
    {periplasm_output_filter}
  ));

  // making a trivial coupled model of periplasm_biomass_filter in order to solve the bug with coupled to atomic
  shared_ptr<flattened_coupled<Time_t, Message_t>> periplasm_biomass_coupled_filter(new flattened_coupled<Time_t, Message_t>(
    {periplasm_biomass_filter}, 
    {periplasm_biomass_filter}, 
    {}, 
    {periplasm_biomass_filter}
  ));

  auto periplasm_space      = compartment_models.at(special_places[1]);
  auto trans_membrane       = enzyme_set_models.at(special_places[1] + "_tm");
  auto outer_membrane       = enzyme_set_models.at(special_places[1] + "_um");
  auto inner_membrane       = enzyme_set_models.at(special_places[1] + "_lm");
  auto periplasm_inner      = enzyme_set_models.at(special_places[1] + "_i");

  shared_ptr<flattened_coupled<Time_t, Message_t>> periplasm_model(new flattened_coupled<Time_t, Message_t>(
    {periplasm_filter, periplasm_or_filter, periplasm_br_filter, periplasm_output_coupled_filter, periplasm_biomass_coupled_filter, periplasm_space, trans_membrane, outer_membrane, inner_membrane, periplasm_inner}, 
    {periplasm_filter, periplasm_or_filter, periplasm_br_filter}, 
    {
      {periplasm_or_filter, periplasm_space},
      {periplasm_br_filter, periplasm_space},
      {periplasm_filter, trans_membrane}, 
      {periplasm_filter, outer_membrane}, 
      {periplasm_filter, inner_membrane}, 
      {trans_membrane, periplasm_space}, 
      {outer_membrane, periplasm_space}, 
      {inner_membrane, periplasm_space}, 
      {periplasm_space, trans_membrane}, 
      {periplasm_space, outer_membrane}, 
      {periplasm_space, inner_membrane}, 
      {periplasm_space, periplasm_inner}, 
      {periplasm_inner, periplasm_space},
      {periplasm_space, periplasm_output_coupled_filter},
      {periplasm_space, periplasm_biomass_coupled_filter}
    }, 
    {trans_membrane, outer_membrane, inner_membrane, periplasm_output_coupled_filter, periplasm_biomass_coupled_filter}
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
          {
            {organelle_filter, organelle_membrane},
            {organelle_membrane, organelle_space},
            {organelle_space, organelle_inner},
            {organelle_inner, organelle_space},
            {organelle_space, organelle_membrane}
          }, 
          {organelle_membrane}
        )
      );
    }
  }

  cout << "Creating the biomass reaction model" << endl;
  auto biomass_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>(biomass_ID);
  auto biomass_atomic = make_atomic_ptr< 
    biomass<Time_t, Message_t>, 
    const string,
    const shared_ptr< map<string, Address_t> >,
    const SetOfMolecules_t&,
    const SetOfMolecules_t&,
    const Address_t,
    const Time_t,
    const Time_t >(
      biomass_ID,
      species_addresses,
      biomass_info.reactants_sctry,
      biomass_info.products_sctry,
      {"e", "e_s", "p_br", "p_s", "c", "c_s"}, // put all the adresses
      BRITime(1,10), // interval time
      BRITime(1,100) // rate_time
    );

  shared_ptr<flattened_coupled<Time_t, Message_t>> biomass_model(new flattened_coupled<Time_t, Message_t>(
    {biomass_filter, biomass_atomic}, 
    {biomass_filter}, 
    {
      {biomass_filter, biomass_atomic}
    }, 
    {biomass_atomic}
  ));

  cout << "Creating the cell coupled model." << endl;
  auto output_filter = make_atomic_ptr< filter<Time_t, Message_t>, const string>("output");
  
  // making a trivial coupled model of output_filter in order to solve the bug with coupled to atomic
  shared_ptr<flattened_coupled<Time_t, Message_t>> output_coupled_filter(new flattened_coupled<Time_t, Message_t>(
    {output_filter}, 
    {output_filter}, 
    {}, 
    {output_filter}
  ));

  vectorOfModels_t cell_models  = {extra_cellular_model, periplasm_model, cytoplasm_model, /*biomass_model,*/ output_coupled_filter};
  vectorOfModels_t cell_eic     = {extra_cellular_model, periplasm_model, cytoplasm_model};
  vectorOfModels_t cell_eoc     = {output_coupled_filter};
  vectorOfModelPairs_t cell_ic  = {
    {extra_cellular_model, periplasm_model},
    {periplasm_model, cytoplasm_model},
    {cytoplasm_model, periplasm_model},
    {periplasm_model, extra_cellular_model},
    // biomass
    /*{biomass_model, extra_cellular_model},
    {biomass_model, periplasm_model},
    {biomass_model, cytoplasm_model},
    {extra_cellular_model, biomass_model},
    {periplasm_model, biomass_model},
    {cytoplasm_model, biomass_model},*/
    // outputs
    {extra_cellular_model, output_coupled_filter},
    {cytoplasm_model, output_coupled_filter},
    {periplasm_model, output_coupled_filter}
  };
  
  for (vectorOfCoupledModels_t::iterator it = organelle_models.begin(); it != organelle_models.end(); ++it) {
    cell_ic.push_back({cytoplasm_model, *it});
    cell_ic.push_back({*it, cytoplasm_model});
  }

  shared_ptr<flattened_coupled<Time_t, Message_t>> cell_model(new flattened_coupled<Time_t, Message_t>(cell_models, cell_eic, cell_ic, cell_eoc));


  /*****************************************************************************************************/
  /************************************** Runing Simulation ********************************************/
  /*****************************************************************************************************/

  cout << "Testing cytoplasm coupled model with filter" << endl;

  cout << "Creating the model to insert the input from stream" << endl;
  auto piss = make_shared<istringstream>();
  string input = "";

  for (double i = 0.01; i < 0.02; i += 0.1) {

    input += "1 " + to_string((int)(i*10000)) + " " + "c c_s | A_c 1 \n ";
    //input += "1 " + to_string((int)(i*10000)) + " " + "e e_s | A_e 1 \n ";
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


  cout << "Creating the model to show the space state from stream" << endl;
  piss = make_shared<istringstream>();
  input = "";

  for (double i = 0.1; i <= 0.1; i += 0.1) {

    input += "1 " + to_string((int)(i*100)) + " \n ";
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


  cout << "Coupling the input to the model" << endl;
  shared_ptr< flattened_coupled<Time_t, Message_t> > root( new flattened_coupled<Time_t, Message_t>{
    {pf, so, cell_model}, 
    {}, 
    {
      {pf, cell_model},
      {so, cell_model}
    }, 
    {cell_model}
  });

  cout << "Preparing runner" << endl;
  Time_t initial_time{0, 1};
  runner<Time_t, Message_t> r(root, initial_time, cout, [](ostream& os, Message_t m){  os << m; });

  cout << "Starting simulation until passivate" << endl;

  auto start = hclock_t::now(); //to measure simulation execution time

  r.runUntil(Time_t(3, 10));

  auto elapsed = chrono::duration_cast< chrono::duration< double, ratio<1> > > (hclock_t::now() - start).count();

  cout << "Simulation took:" << elapsed << "sec" << endl;

  return 0;
}
