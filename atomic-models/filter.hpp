#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <string>
#include <utility>
#include <map>
#include <memory>

#include <boost/simulation/pdevs/atomic.hpp>

#include "../data-structures/types.hpp" /* Address_t */


using namespace boost::simulation::pdevs;
using namespace boost::simulation;
using namespace std;

template<class TIME, class MSG>
class filter : public pdevs::atomic<TIME, MSG>
{
private:

  string      _accepted_input;
  vector<MSG> _filtered_input;

  // constant values
  TIME ZERO;

public:

  explicit filter(const string other_accepted_input) noexcept :
  _accepted_input(other_accepted_input) {

    _filtered_input.clear();

    // Constant values;
    ZERO = TIME(0);
  }

  void internal() noexcept {
    //if (_accepted_input == "output") cout << "internal" << endl;
    _filtered_input.clear();
  }

  TIME advance() const noexcept {
    //if (_accepted_input == "output") cout << "advance " << endl;

    TIME result = (!_filtered_input.empty()) ? ZERO : pdevs::atomic<TIME, MSG>::infinity;
    //if (_accepted_input == "output")  cout << "advance return: " << result << endl;
    return result;
  }

  vector<MSG> out() const noexcept {
    //if (_accepted_input == "output") cout << "out" << endl;
    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    //if (_accepted_input == "output") cout << "external "  << endl;
    
    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      //if (_accepted_input == "output")  cout << "input " << it->to << endl;
      if (isAcceptedInpunt(_accepted_input, it->to)) {

        _filtered_input.push_back(*it);
      }
    }

    //if (_accepted_input == "output") cout << "end external" << endl;
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {
    //if (_accepted_input == "output") cout << "confluence "  << endl;
    internal();
    external(mb, ZERO);
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  bool isAcceptedInpunt(string a, Address_t to) {

    bool result = false;
    for (Address_t::iterator i = to.begin(); i != to.end(); ++i) {
      if( a == (*i)) {
        result = true;
        break;
      }
    }

    return result;
  }

};

#endif // BOOST_SIMULATION_PDEVS_FILTER_H
