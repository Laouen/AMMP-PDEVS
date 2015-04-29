#include "unit_definition.hpp"
#include <list>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

UnitDefinition_t::UnitDefinition_t(list<Unit_t> other_list_of_units, string other_name) 
: list_of_units(other_list_of_units), name(other_name) {

  general_multiplier      = 1;
  general_scale_exponent  = 0;

  for (list<Unit_t>::iterator it = list_of_units.begin(); it != list_of_units.end(); ++it){

    general_multiplier *= pow(it->multiplier, it->exponent);
    general_scale_exponent += it->scale * it->exponent;
  }
}

double UnitDefinition_t::convertFrom(double amount) const {

  return amount * (general_multiplier * pow(10, general_scale_exponent));
}

double UnitDefinition_t::convertTo(double amount) const {

  return amount / (general_multiplier * pow(10, general_scale_exponent));
}

string UnitDefinition_t::unitName() const {

  return name;
}