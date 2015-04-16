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

enum class RState { REJECTING, REACTING, SELECTING };
enum class SState { SENDING, SELECTING, IDLE };
enum class Way { RTP, PTR };

using Integer = unsigned long long;
using SetOfMolecules  = map<string, Integer>;



/******************************************/
/******** End enums and renames ***********/
/******************************************/

template<class TIME>
struct Task {
  TIME                time_left;
  RState              task_kind;
  SetOfMolecules      rejected;  
  pair<Way, Integer>  reaction;

  Task() {}

  Task(const Task& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    rejected  = other.rejected;
    reaction  = other.reaction;
  }

  inline bool operator<(const Task<TIME>& o) {

    bool result;
    if (time_left != o.time_left) {
      
      result = (time_left < o.time_left);
    } else {

      if ((task_kind == RState::SELECTING) && (o.task_kind != RState::SELECTING)) result = true;
      else if ((task_kind == RState::REJECTING) && (o.task_kind == RState::REACTING)) result = true;
      else result = false;
    }

    return result;
  }

  inline bool operator==(const Task<TIME>& o) {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if (task_kind == RState::REJECTING) {
      result = result && (rejected == o.rejected);
    }

    if (task_kind == RState::REACTING) {
      result = result && (reaction == o.reaction);
    }    

    return result;
  }

};

template<class TIME>
using TaskQueue = list< Task<TIME> >;


/*******************************************/
/**************** Message ******************/
/*******************************************/

using Address = list<string>;


struct Message {
  
  Address to;
  string specie;
  Integer amount;

  Message(Address other_to, string other_specie, Integer other_amount)
  : to(other_to), specie(other_specie), amount(other_amount) {}

  Message()
  :to(), specie(""), amount(0) {}

  Message(const Message& other)
  : to(other.to), specie(other.specie), amount(other.amount) {}

  Message(Message* other)
  : to(other->to), specie(other->specie), amount(other->amount) {}

  void clear() {
    to.clear();
    specie = "";
    amount = 0;
  }
};

ostream& operator<<(ostream& os, Message msg);
ostream& operator<<(ostream& os, Address to);

/*******************************************/
/************** End Message ****************/
/*******************************************/

struct metabolite_info_t {
  
  Integer         amount;
  vector<Address> enzymes;

  metabolite_info_t()
  : amount(0), enzymes() {}

  metabolite_info_t(Integer other_amount, const vector<Address>& other_enzymes)
  : amount(other_amount), enzymes(other_enzymes) {}

  metabolite_info_t(const metabolite_info_t& other)
  : amount(other.amount), enzymes(other.enzymes) {}
};

struct enzyme_info_t {

  Address          location;
  vector<string>  reactants;

  enzyme_info_t()
  : location(), reactants() {}

  enzyme_info_t(const Address& other_location, const vector<string>& other_reactants)
  : location(other_location), reactants(other_reactants) {}

  enzyme_info_t(const enzyme_info_t& other)
  : location(other.location), reactants(other.reactants) {}
};


/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H