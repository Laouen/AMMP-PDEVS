#ifndef MESSAGE_TYPES_H
#define MESSAGE_TYPES_H

#include <string>
#include <list>
#include <map>

#include "../structures/types.hpp" // MetaboliteAmounts, Address_t, Integer, Way

using namespace std;

namespace msg_event {
  struct Metabolites {

    Address_t to;
    string from;
    Way react_direction;
    Integer react_amount;
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