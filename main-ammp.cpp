// STD includes
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <algorithm>
#include <utility>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <limits>

// boost simalator include
#include <boost/simulation.hpp>

// XML parser include
#include "tinyXML/tinyxml.h"

// model includes
#include "atomic-models/filter.hpp"
#include "atomic-models/reaction.hpp"
#include "atomic-models/space.hpp"

// data structure includes
#include "data-structures/unit_definition.hpp"
#include "data-structures/types.hpp"
#include "data-structures/randomNumbers.hpp"

#define TIXML_USE_STL

using namespace boost::simulation;
using namespace boost::simulation::pdevs;
using namespace boost::simulation::pdevs::basic_models;
using namespace std;


/***************************************/
/********* Type definations ************/
/***************************************/

typedef double Time;
typedef chrono::high_resolution_clock hclock;
typedef vector< shared_ptr< model<Time> > > vectorOfModels;
typedef vector< pair< shared_ptr< model<Time> >, shared_ptr< model<Time> > > > vectorOfModelPairs;



/***************************************/
/******** End type definations *********/
/***************************************/


/**************************************************************************************************************/
/***************************************** Testing models *****************************************************/
/**************************************************************************************************************/

bool weightedRandomBool(Integer a, long double _volume, long double _factor, RealRandom<double>& _distribution) {

  //double proportion = _volume / (double)a;
  //double threshold  = (double)1.0 / pow( (double)e, (double)_factor*proportion );

  return (_distribution.drawNumber(0.0, 1.0) < 0.5);

}

int main () {

  cout.precision(numeric_limits< long double >::digits);

  map<string, metabolite_info_t> metabolites;
  map<string, enzyme_info_t> enzymes;
  long double volume = 1.0;
  long double factor = 1.0;
  Integer     amount = 10000;
  Integer     count  = Integer(1000000)*Integer(100);

  random_device rd;
  RealRandom<double> distribution(rd());
  Integer total = 0;

  for (Integer i = 0; i < count; ++i){
    if ((i % 1000000) == 0) cout << i << endl;
    if (weightedRandomBool(amount, volume, factor, distribution)) ++total;
  }

  cout << total << " " << count << endl;
  cout << (long double)( (double)total / (double)count ) << endl;

  space<Time, Message> new_space(0.001, metabolites, enzymes, 0.052, 1.0);
  
  return 0;
}


/**************************************************************************************************************/
/*************************************** End Testing models ***************************************************/
/**************************************************************************************************************/




