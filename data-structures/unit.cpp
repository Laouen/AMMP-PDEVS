#include "unit.hpp"
#include <cmath>

using namespace std;

double UnitDefinition::convert(double amount) {
  double m = 1;
  double s = 0;

  for (list<unit>::iterator it = list_of_units.begin(); it != list_of_units.end(); ++it){

    m *= it->multiplier;
    s += it->scale * it->exponent;
  }

  return m * pow(10, s);
}