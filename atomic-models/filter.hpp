#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include <memory>

#include "../data-structures/types.hpp" /* Address_t */


using namespace boost::simulation::pdevs;
using namespace std;

template<class TIME, class MSG>
class filter : public atomic<TIME, MSG>
{
private:

  string      _accepted_input;
  vector<MSG> _filtered_input;

public:

  explicit filter(const string other_accepted_input) noexcept :
  _accepted_input(other_accepted_input) {

    _filtered_input.clear();
  }

  void internal() noexcept {
    if (_accepted_input == "output")cout << "internal" << endl;
    _filtered_input.clear();
  }

  TIME advance() const noexcept {

    TIME result = (_filtered_input.size() > 0) ? TIME(0) : atomic<TIME, MSG>::infinity;
    if (_accepted_input == "output") cout << "advance " << endl;

    return result;
  }

  vector<MSG> out() const noexcept {
    
    if (_accepted_input == "output") cout << "out" << endl;
    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    if (_accepted_input == "output") cout << "external "  << endl;
    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
    
      if (isAcceptedInpunt(_accepted_input, it->to)) {

        _filtered_input.push_back(*it);
      }
    }
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    internal();
    external(mb, TIME(0));
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
