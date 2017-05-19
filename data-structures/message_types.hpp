#ifndef MESSAGE_TYPES_H
#define MESSAGE_TYPES_H

#include <string>
#include <list>
#include <map>

#include "types.hpp" // SetOfMolecules_t, Address_t, Integer_t, Way_t

using namespace std;

namespace messages {

  struct Metabolites {

    Address_t to;
    string from;
    Way_t react_direction;
    Integer_t react_amount;
    SetOfMolecules_t metabolites;

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