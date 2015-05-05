#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>


using namespace std;

/******************************************/
/********** Enums and renames *************/
/******************************************/

enum class RState_t { SELECTING = 0, REJECTING = 1, REACTING = 2 };
enum class SState_t { SHOWING = 0, SENDING_BIOMAS = 1, SENDING_REACTION = 2, SELECTING = 3, IDLE = 4 };
enum class BState_t { ENOUGH = 0, NOT_ENOUGH = 1, NOTHING = 2 };
enum class Way_t { RTP, PTR };

using Integer_t = unsigned long long;
using SetOfMolecules_t  = map<string, Integer_t>;



/******************************************/
/******** End enums and renames ***********/
/******************************************/







/*******************************************/
/**************** RTask_t ******************/
/*******************************************/

struct enzyme_parameter_t {
  
  string            name;
  bool              reversible;
  SetOfMolecules_t  reactants_sctry;
  SetOfMolecules_t  products_sctry;
};

template<class TIME>
struct RTask_t {
  TIME                    time_left;
  RState_t                task_kind;
  SetOfMolecules_t        rejected;  
  pair<Way_t, Integer_t>  reaction;

  RTask_t() {}

  RTask_t(const RTask_t& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    rejected  = other.rejected;
    reaction  = other.reaction;
  }

  inline bool operator<(const RTask_t<TIME>& o) {

    bool result;
    if (time_left != o.time_left) result = (time_left < o.time_left);
    else                          result = (task_kind < o.task_kind);

    return result;
  }

  inline bool operator==(const RTask_t<TIME>& o) {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if (task_kind == RState_t::REJECTING) {
      result = result && (rejected == o.rejected);
    }

    if (task_kind == RState_t::REACTING) {
      result = result && (reaction == o.reaction);
    }    

    return result;
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



template<class TIME>
struct STask_t {
  TIME                    time_left;
  SState_t                task_kind;
  SetOfMolecules_t        to_send;

  STask_t() {}

  STask_t(const STask_t& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    to_send   = other.to_send;
  }

  inline bool operator<(const STask_t<TIME>& o) {

    bool result;
    if (time_left != o.time_left) result = (time_left < o.time_left);
    else                          result = (task_kind < o.task_kind);

    return result;
  }

  inline bool operator==(const STask_t<TIME>& o) {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if ((task_kind == SState_t::SENDING_REACTION) || (task_kind == SState_t::SENDING_BIOMAS)) {
      result = result && (to_send == o.to_send);
    }  

    return result;
  }

};

template<class TIME>
using STaskQueue_t = list< STask_t<TIME> >;


/*******************************************/
/************** End STask_t ****************/
/*******************************************/








/*******************************************/
/**************** Message_t ****************/
/*******************************************/

using Address_t = list<string>;


struct Message_t {
  
  Address_t to;
  string specie;
  Integer_t amount;
  bool show_request;

  Message_t(const Address_t& other_to, const string& other_specie, const Integer_t& other_amount)
  : to(other_to), specie(other_specie), amount(other_amount), show_request(false) {}

  Message_t()
  :to(), specie(""), amount(0), show_request(false) {}

  Message_t(const Message_t& other)
  : to(other.to), specie(other.specie), amount(other.amount), show_request(other.show_request) {}

  Message_t(Message_t* other)
  : to(other->to), specie(other->specie), amount(other->amount), show_request(other->show_request) {}

  void clear() {
    to.clear();
    specie        = "";
    amount        = 0;
    show_request  = false;
  }
};

ostream& operator<<(ostream& os, Message_t msg);
ostream& operator<<(ostream& os, Address_t to);
ostream& operator<<(ostream& os, vector<string> m);
ostream& operator<<(ostream& os, SetOfMolecules_t m);

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

struct enzyme_info_t {

  Address_t          location;
  vector<string>  reactants;

  enzyme_info_t()
  : location(), reactants() {}

  enzyme_info_t(const Address_t& other_location, const vector<string>& other_reactants)
  : location(other_location), reactants(other_reactants) {}

  enzyme_info_t(const enzyme_info_t& other)
  : location(other.location), reactants(other.reactants) {}
};

/*******************************************/
/*********** End Data info type ************/
/*******************************************/






/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H