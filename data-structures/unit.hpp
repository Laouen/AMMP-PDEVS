#ifndef AMMP_PDEVS_UNIT_DEFINITION_H
#define AMMP_PDEVS_UNIT_DEFINITION_H
#include <list>

using namespace std;

struct Unit {
  string  kind;
  double  exponent;
  int     scale;
  double  multiplier;

  Unit(string other_kind, double other_exponent, int other_scale, double other_multiplier)
  : kind(other_kind), exponent(other_exponent), scale(other_scale), multiplier(other_multiplier) {}


};

class UnitDefinition {
  
public:

  // constructors
  UnitDefinition();
  UnitDefinition(list<unit> other_list_of_units)
  : list_of_units(other_list_of_units) {}

  //methods
  double convert(double amount);

private:
  list<unit> list_of_units;
};

#endif // AMMP_PDEVS_UNIT_DEFINITION_H