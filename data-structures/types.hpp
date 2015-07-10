#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>

#include <boost/simulation.hpp>


using namespace std;
using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;

/******************************************/
/********** Enums and renames *************/
/******************************************/

enum class SState_t { SHOWING = 0, SELECTING_FOR_BIOMAS = 1, SELECTING_FOR_REACTION = 2, SENDING_BIOMAS = 3, SENDING_REACTIONS = 4};
enum class BState_t { ENOUGH = 0, NOT_ENOUGH = 1, NOTHING = 2, START = 3 };
enum class Way_t { STP, PTS };

using Integer_t = unsigned long long;
using SetOfMolecules_t  = map<string, Integer_t>;



/******************************************/
/******** End enums and renames ***********/
/******************************************/

/******************************************/
/************** models type ***************/
/******************************************/

template<class TIME>
using vm_t  = vector<shared_ptr< model<TIME>>>;
template<class TIME>
using vmp_t = vector<pair< shared_ptr< model<TIME>>, shared_ptr< model<TIME>>>>;
template<class TIME>
using mm_t  = map<string, shared_ptr< model<TIME>>>;
template<class TIME,class MSG>
using vcm_t = vector< shared_ptr<flattened_coupled<TIME, MSG>>>;
template<class TIME,class MSG>
using cmm_t = map<string, shared_ptr<flattened_coupled<TIME, MSG>>>;

/******************************************/
/************ End models type *************/
/******************************************/

/******************************************/
/************ Parser types ****************/
/******************************************/


struct enzyme_parameter_t {
  
  string            name;
  bool              reversible;
  SetOfMolecules_t  reactants_sctry;
  SetOfMolecules_t  products_sctry;
};


/******************************************/
/*********** End parser types *************/
/******************************************/


/*******************************************/
/**************** RTask_t ******************/
/*******************************************/

template<class TIME>
struct RTask_t {
  TIME        time_left;
  Way_t       direction; 
  Integer_t   amount;

  RTask_t() {}

  RTask_t(const TIME& other_t, const Way_t& other_d, const Integer_t& other_a) 
  : time_left(other_t), direction(other_d), amount(other_a) {}

  RTask_t(const RTask_t<TIME>& other) 
  : time_left(other.time_left), direction(other.direction), amount(other.amount) {}

  inline bool operator<(const RTask_t<TIME>& o) const {

    return (time_left < o.time_left);
  }

  inline bool operator==(const RTask_t<TIME>& o)  const { 

    return (time_left == o.time_left) && (direction == o.direction) && (amount == o.amount);
  }

};

template<class TIME>
using RTaskQueue_t = list< RTask_t<TIME> >;

/*******************************************/
/************** End RTask_t ****************/
/*******************************************/









/*******************************************/
/**************** STask_t ******************/
/*******************************************/



template<class TIME, class MSG>
struct STask_t {
  TIME          time_left;
  SState_t      task_kind;
  vector<MSG>   msgs;

  STask_t() {}

  STask_t(const STask_t<TIME, MSG>& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    msgs   = other.msgs;
  }

  inline bool operator<(const STask_t<TIME, MSG>& o)  const {

    bool result;
    if (time_left != o.time_left) result = (time_left < o.time_left);
    else                          result = (task_kind < o.task_kind);

    return result;
  }

  inline bool operator==(const STask_t<TIME, MSG>& o)  const {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if ((task_kind == SState_t::SENDING_REACTIONS) || (task_kind == SState_t::SENDING_BIOMAS)) {
      result = result && (msgs == o.msgs);
    }  

    return result;
  }

};

template<class TIME, class MSG>
using STaskQueue_t = list< STask_t<TIME, MSG> >;


/*******************************************/
/************** End STask_t ****************/
/*******************************************/





/*******************************************/
/**************** Message_t ****************/
/*******************************************/

using Address_t = list<string>;

struct Message_t {
  
  Address_t to;
  string from;
  SetOfMolecules_t metabolites;
  Way_t react_direction;
  Integer_t react_amount;
  bool show_request;
  bool biomass_request;

  Message_t(const Address_t& other_to, string other_from, const SetOfMolecules_t& other_m, Way_t other_rd, Integer_t other_ra, bool other_sr, bool other_br)
  : to(other_to), from(other_from), metabolites(other_m), react_direction(other_rd), react_amount(other_ra), show_request(other_sr), biomass_request(other_br) {}

  Message_t()
  :to(), from(""), metabolites(), react_direction(), react_amount(), show_request(false), biomass_request(false) {}

  Message_t(const Message_t& other)
  : to(other.to), from(other.from), metabolites(other.metabolites), react_direction(other.react_direction), react_amount(other.react_amount), show_request(other.show_request), biomass_request(other.biomass_request) {}

  Message_t(Message_t& other)
  : to(other.to), from(other.from), metabolites(other.metabolites), react_direction(other.react_direction), react_amount(other.react_amount), show_request(other.show_request), biomass_request(other.biomass_request) {}

  void clear() {
    to.clear();
    from = "";
    metabolites.clear();
    show_request = false;
    biomass_request = false;
  }
};

ostream& operator<<(ostream& os, const Message_t& msg);
ostream& operator<<(ostream& os, const Address_t& to);
ostream& operator<<(ostream& os, const vector<string>& m);
ostream& operator<<(ostream& os, const SetOfMolecules_t& m);
ostream& operator<<(ostream& os, const SState_t& s);
ostream& operator<<(ostream& os, const BState_t& s);
ostream& operator<<(ostream& os, const RState_t& s);
ostream& operator<<(ostream& os, const Way_t& s);

/*******************************************/
/************** End Message_t **************/
/*******************************************/







/*******************************************/
/************* Data info type **************/
/*******************************************/

struct metabolite_info_t {
  
  Integer_t         amount;
  vector<Address_t> enzymes;

  metabolite_info_t()
  : amount(0), enzymes() {}

  metabolite_info_t(Integer_t other_amount, const vector<Address_t>& other_enzymes)
  : amount(other_amount), enzymes(other_enzymes) {}

  metabolite_info_t(const metabolite_info_t& other)
  : amount(other.amount), enzymes(other.enzymes) {}
};

// TODO this type must be removed after the new implementation.
struct enzyme_info_t {

  Address_t         location;
  Integer_t         amount;
  SetOfMolecules_t  reactants_sctry;
  SetOfMolecules_t  products_sctry;
  double            konSTP;
  double            konPTS;
  bool              reversible;

  enzyme_info_t()
  : location(), reactants() {}

  enzyme_info_t(const Address_t& other_location, const vector<string>& other_reactants, double other_konSTP, double other_konPTS, bool other_reversible)
  : location(other_location), reactants(other_reactants), konSTP(other_konSTP), konPTS(other_konPTS), reversible(other_reversible) {}

  enzyme_info_t(const enzyme_info_t& other)
  : location(other.location), reactants(other.reactants), konSTP(other.konSTP), konPTS(other.konPTS), reversible(other.reversible) {}
};

metabolite_info_t& metabolite_info_t::operator=(metabolite_info_t &&) = default;

/*******************************************/
/*********** End Data info type ************/
/*******************************************/






/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H