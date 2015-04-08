#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include "../data-structures/types.hpp"


using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** vector types ************/
typedef vector<molecule> vectorOfMolecules;

template<class TIME, class MSG>
class space : public atomic<TIME, MSG>
{
private:
  TIME _next_internal;
  TIME _interval_time;

public:

  explicit space(const vector<reaction_input> reactions) noexcept {
  }

  void internal() noexcept { 
  }

  TIME advance() const noexcept {
    
    return _next_internal;
  }

  vector<MSG> out() const noexcept {
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    internal();
    external(mb, TIME(0));
  }
};

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
