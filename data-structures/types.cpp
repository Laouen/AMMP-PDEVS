#include "types.hpp"

ostream& operator<<(ostream& os, Address_t to) {
  
  os << "[";
  Address_t::iterator i = to.begin();
  while(i != to.end()){
    os << *i;
    ++i;
    if (i != to.end()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, Message_t msg) {
  os << "To: " << msg.to << endl;
  os << "Specie: " << msg.specie << endl;
  os << "Amount: " << msg.amount << endl;
  os << "show request: " << (msg.show_request ? "true" : "false") << endl;
  os << "send biomass: " << (msg.biomass_request ? "true" : "false");
  return os;
}

ostream& operator<<(ostream& os, vector<string> m) {
  
  os << "[";
  vector<string>::iterator i = m.begin();
  while(i != m.end()){
    os << *i;
    ++i;
    if (i != m.end()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, SetOfMolecules_t m) {
  
  os << "[";
  SetOfMolecules_t::iterator i = m.begin();
  while(i != m.end()){
    os << i->second << "-" << i->first;
    ++i;
    if (i != m.end()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, SState_t s) {

  switch(s) {
    case SState_t::SHOWING:
      os << "SHOWING";
      break;
    case SState_t::SENDING_BIOMAS:
      os << "SENDING_BIOMAS";
      break;
    case SState_t::SENDING_REACTIONS:
      os << "SENDING_REACTIONS";
      break;
    case SState_t::SELECTING_FOR_BIOMAS:
      os << "SELECTING_FOR_BIOMAS";
      break;
    case SState_t::SELECTING_FOR_REACTION:
      os << "SELECTING_FOR_REACTION";
      break;
  }

  return os;
}

ostream& operator<<(ostream& os, BState_t s) {

  switch(s) {
    case BState_t::START:
      os << "START";
      break;
    case BState_t::NOTHING:
      os << "NOTHING";
      break;
    case BState_t::NOT_ENOUGH:
      os << "NOT_ENOUGH";
      break;
    case BState_t::ENOUGH:
      os << "ENOUGH";
      break;
  }

  return os;
}