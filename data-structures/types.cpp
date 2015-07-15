#include "types.hpp"

ostream& operator<<(ostream& os, const Address_t& to) {
  
  os << "[";
  Address_t::const_iterator i = to.cbegin();
  while(i != to.cend()){
    os << *i;
    ++i;
    if (i != to.cend()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, const Message_t& msg) {
  os << "To: " << msg.to << endl;
  os << "Metabolites: " << msg.metabolites << endl;
  os << "reaction way: " << msg.react_direction << endl;
  os << "show request: " << (msg.show_request ? "true" : "false") << endl;
  os << "send biomass: " << (msg.biomass_request ? "true" : "false");
  return os;
}

ostream& operator<<(ostream& os, const vector<string>& m) {
  
  os << "[";
  vector<string>::const_iterator i = m.cbegin();
  while(i != m.cend()){
    os << *i;
    ++i;
    if (i != m.cend()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, const SetOfMolecules_t& m) {
  
  os << "[";
  SetOfMolecules_t::const_iterator i = m.cbegin();
  while(i != m.cend()){
    os << i->second << "-" << i->first;
    ++i;
    if (i != m.cend()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, const SState_t& s) {

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

ostream& operator<<(ostream& os, const BState_t& s) {

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

ostream& operator<<(ostream& os, const Way_t& s) {

  switch(s) {
    case Way_t::STP:
      os << "STP";
      break;
    case Way_t::PTS:
      os << "PTS";
      break;
  }

  return os;
}