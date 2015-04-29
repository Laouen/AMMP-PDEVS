#ifndef AMMP_PDEVS_UNIT_DEFINITION_H
#define AMMP_PDEVS_UNIT_DEFINITION_H
#include <list>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

struct Unit_t {
  // unit type
  string  kind;
  // unit dimention
  int     scale;
  // conversion factors
  double  exponent;
  double  multiplier;

  Unit_t(string other_kind, double other_exponent, int other_scale, double other_multiplier)
  : kind(other_kind), exponent(other_exponent), scale(other_scale), multiplier(other_multiplier) {}

  double convertFrom(double amount) {

    return amount * pow(multiplier * pow(10, scale), exponent);
  }

  double convertTo(double amount) {

    return amount / pow(multiplier * pow(10, scale), exponent);
  }

};

class UnitDefinition_t {
  
public:

  // constructors
  UnitDefinition_t();
  UnitDefinition_t(list<Unit_t>, string);

  //methods
  double convertTo(double amount) const;
  double convertFrom(double amount) const;
  string unitName() const;

private:
  string name;
  list<Unit_t> list_of_units;
  double general_multiplier;
  double general_scale_exponent;
};

#endif // AMMP_PDEVS_UNIT_DEFINITION_H