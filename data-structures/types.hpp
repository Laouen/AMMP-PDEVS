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

enum class RState_t { REJECTING, REACTING, SELECTING };
enum class SState_t { SENDING, SELECTING, IDLE };
enum class Way_t { RTP, PTR };

using Integer_t = unsigned long long;
using SetOfMolecules_t  = map<string, Integer_t>;



/******************************************/
/******** End enums and renames ***********/
/******************************************/

struct enzyme_parameter_t {
  
  string          name;
  bool            reversible;
  SetOfMolecules_t  reactants_sctry;
  SetOfMolecules_t  products_sctry;
};

template<class TIME>
struct Task_t {
  TIME                time_left;
  RState_t              task_kind;
  SetOfMolecules_t      rejected;  
  pair<Way_t, Integer_t>  reaction;

  Task_t() {}

  Task_t(const Task_t& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    rejected  = other.rejected;
    reaction  = other.reaction;
  }

  inline bool operator<(const Task_t<TIME>& o) {

    bool result;
    if (time_left != o.time_left) {
      
      result = (time_left < o.time_left);
    } else {

      if ((task_kind == RState_t::SELECTING) && (o.task_kind != RState_t::SELECTING)) result = true;
      else if ((task_kind == RState_t::REJECTING) && (o.task_kind == RState_t::REACTING)) result = true;
      else result = false;
    }

    return result;
  }

  inline bool operator==(const Task_t<TIME>& o) {

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
using TaskQueue = list< Task_t<TIME> >;


/*******************************************/
/**************** Message_t ******************/
/*******************************************/

using Address_t = list<string>;


struct Message_t {
  
  Address_t to;
  string specie;
  Integer_t amount;

  Message_t(const Address_t& other_to, const string& other_specie, const Integer_t& other_amount)
  : to(other_to), specie(other_specie), amount(other_amount) {}

  Message_t()
  :to(), specie(""), amount(0) {}

  Message_t(const Message_t& other)
  : to(other.to), specie(other.specie), amount(other.amount) {}

  Message_t(Message_t* other)
  : to(other->to), specie(other->specie), amount(other->amount) {}

  void clear() {
    to.clear();
    specie = "";
    amount = 0;
  }
};

ostream& operator<<(ostream& os, Message_t msg);
ostream& operator<<(ostream& os, Address_t to);
ostream& operator<<(ostream& os, vector<string> m);

/*******************************************/
/************** End Message_t ****************/
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


/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H