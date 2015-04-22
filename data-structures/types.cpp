#include "types.hpp"

ostream& operator<<(ostream& os, Address to) {
  
  os << "[";
  Address::iterator i = to.begin();
  while(i != to.end()){
    os << *i;
    ++i;
    if (i != to.end()) os << ", ";
  }
  os << "]";
  return os;
}

ostream& operator<<(ostream& os, Message msg) {
  os << "To: " << endl << msg.to << endl;
  os << "Specie: " << msg.specie << endl;
  os << "Amount: " << msg.amount << endl;
  return os;
}