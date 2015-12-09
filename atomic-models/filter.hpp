#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <string>
#include <utility>
#include <map>
#include <memory>

#include <boost/simulation/pdevs/atomic.hpp>

#include "../data-structures/types.hpp" /* Address_t */

#define MINIMUM_TIME_FOR_FILTER TIME(1,1000);

using namespace boost::simulation::pdevs;
using namespace boost::simulation;
using namespace std;

template<class TIME, class MSG>
class filter : public pdevs::atomic<TIME, MSG>
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
    comment("internal init.");
    _filtered_input.clear();
    comment("internal end.");
  }

  TIME advance() const noexcept {
    comment("advance init.");
    TIME result;
    if (_filtered_input.empty()) result = pdevs::atomic<TIME, MSG>::infinity;
    else result = MINIMUM_TIME_FOR_FILTER;

    if (result <= TIME(0,1)) cout << _accepted_input << " " << result << endl;
    comment("advance time result " + result.toStringAsDouble());
    return result;
  }

  vector<MSG> out() const noexcept {
    comment("out init.");
    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    comment("external init.");
    
    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      if (isAcceptedInpunt(_accepted_input, it->to)) {
        _filtered_input.push_back(*it);
      }
    }

    comment("external end.");
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {
    comment("conluence init.");
    internal();
    external(mb, ZERO);
    comment("conluence end.");
  }

  /***************************************
  ********* helper functions *************
  ***************************************/
  void comment(string msg) const {
    if (COMMENTS) cout << "[filter " << _accepted_input << "] " << msg << endl;
  }

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
