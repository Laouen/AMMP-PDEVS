#ifndef MESSAGE_TYPES_H
#define MESSAGE_TYPES_H

#include <string>
#include <list>
#include <map>

#include "types.hpp" // MetaboliteAmounts, Address_t, Integer_t, Way_t

using namespace std;

namespace msg_event {
  struct Metabolites {

    Address_t to;
    string from;
    Way_t react_direction;
    Integer_t react_amount;
    MetaboliteAmounts metabolites;

    Metabolites() {
      this->clear();
    }

    void clear() {
      to.clear();
      from = "";
      react_amount = 0;
      metabolites.clear();
    }
  };
}

#endif