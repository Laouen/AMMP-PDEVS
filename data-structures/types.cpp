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
  os << "To: " << endl << msg.to << endl;
  os << "Specie: " << msg.specie << endl;
  os << "Amount: " << msg.amount << endl;
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