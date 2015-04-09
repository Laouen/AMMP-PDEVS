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

long double e = 2.71828182845904523536028747135266249775724709369995L;

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

struct Adress {

  string _organelle;
  string _cytoplasm;
  string _ensyme_set;
  string _compartment;
  string _enzyme;

  Adress(string other_organelle, string other_cytoplasm, string other_ensyme_set, string other_compartment, string other_enzyme)
  : _organelle(other_organelle), _cytoplasm(other_cytoplasm), _ensyme_set(other_ensyme_set), _compartment(other_compartment), _enzyme(other_enzyme) {}

  Adress()
  : _organelle(""), _cytoplasm(""), _ensyme_set(""), _compartment(""), _enzyme("") {}

  Adress(const Adress& other)
  : _organelle(other._organelle), _cytoplasm(other._cytoplasm), _ensyme_set(other._ensyme_set), _compartment(other._compartment), _enzyme(other._enzyme) {}
  
  Adress(Adress* other)
  : _organelle(other->_organelle), _cytoplasm(other->_cytoplasm), _ensyme_set(other->_ensyme_set), _compartment(other->_compartment), _enzyme(other->_enzyme) {}



  string atModel(string model_type) const{
    
    string result = "";
    if(model_type == "organelle") {
      
      result = _organelle;
    } else if(model_type == "cytoplasm") {
      
      result = _cytoplasm;
    } else if(model_type == "enzyme set") {
      
      result = _ensyme_set;
    } else if(model_type == "compartment") {
      
      result = _compartment;
    } else if(model_type == "enzyme") {
      
      result = _enzyme;
    }

    return result;
  }

  void setModel(string model_type, string new_value) {

    if(model_type == "organelle") {
      
      _organelle = new_value;
    } else if(model_type == "cytoplasm") {
      
      _cytoplasm = new_value;
    } else if(model_type == "enzyme set") {
      
      _ensyme_set = new_value;
    } else if(model_type == "compartment") {
      
      _compartment = new_value;
    } else if(model_type == "enzyme") {
      
      _enzyme = new_value;
    }
  }

  void clear() {
    _organelle    = "";
    _cytoplasm    = "";
    _ensyme_set   = "";
    _compartment  = "";
    _enzyme     = "";
  }

};

ostream& operator<<(ostream& os, Adress to) {
  os << "Organelle: " << to._organelle << endl;
  os << "Cytoplasm: " << to._cytoplasm << endl;
  os << "Enzyme set: " << to._ensyme_set << endl;
  os << "Compartment: " << to._compartment << endl;
  os << "Enzyme: " << to._enzyme << endl;
  return os;
}

struct Message {
  
  Adress to;
  string specie;
  Integer amount;

  Message(Adress other_to, string other_specie, Integer other_amount)
  : to(other_to), specie(other_specie), amount(other_amount) {}

  Message()
  :to(), specie(""), amount(0) {}

  Message(const Message& other)
  : to(other.to), specie(other.specie), amount(other.amount) {}

  Message(Message* other)
  : to(other->to), specie(other->specie), amount(other->amount) {}

  bool sendTo(string model_type, string new_value) {
    to.setModel(model_type, new_value);

    return (to.atModel(model_type) == new_value);
  }

  void clear() {
    to.clear();
    specie = "";
    amount = 0;
  }
};

ostream& operator<<(ostream& os, Message msg) {
  os << "To: " << endl << msg.to << endl;
  os << "Specie: " << msg.specie << endl;
  os << "Amount: " << msg.amount << endl;
  return os;
}

/*******************************************/
/************** End Message ****************/
/*******************************************/


struct metabolite_info_t {
  
  Integer         amount;
  bool            to_send;
  vector<Adress>  enzymes;

  metabolite_info_t()
  : amount(0), to_send(false), enzymes() {}

  metabolite_info_t(string other_specie, Integer other_amount, bool other_to_send, vector<Adress> other_enzymes)
  : amount(other_amount), to_send(other_to_send), enzymes(other_enzymes) {}

  metabolite_info_t(const metabolite_info_t& other)
  : amount(other.amount), to_send(other.to_send), enzymes(other.enzymes) {}
};

struct enzyme_info_t {

  Adress          location;
  vector<string>  reactants;

  enzyme_info_t()
  : location(), reactants() {}

  enzyme_info_t(Adress other_location, vector<string> other_reactants)
  : location(other_location), reactants(other_reactants) {}

  enzyme_info_t(const enzyme_info_t& other)
  : location(other.location), reactants(other.reactants) {}
};


/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H