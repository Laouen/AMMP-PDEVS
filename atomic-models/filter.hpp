#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include "../data-structures/message.hpp"


using namespace boost::simulation::pdevs;
using namespace std;

template<class TIME, class MSG>
class filter : public atomic<TIME, MSG>
{
private:

  string      _model_type;
  string      _acepted_input;
  vector<MSG> _filtered_input;

public:

  explicit filter(string other_model_type, string other_acepted_input) noexcept :
  _model_type(other_model_type),
  _acepted_input(other_acepted_input) {

    _filtered_input.clear();
  }

  void internal() noexcept {
    cout << _acepted_input << " internal." << endl;
    _filtered_input.clear();
  }

  TIME advance() const noexcept {
    cout << _acepted_input << " advance." << endl;
    return (_filtered_input.size() > 0) ? TIME(0) : atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {
    cout << _acepted_input << " out." << endl;
    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    cout << _acepted_input << " external." << endl;
    cout << "current state catn msg: " << _filtered_input.size() << endl;
    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      if ( (_acepted_input != "") && (it->to.atModel(_model_type) == _acepted_input) ) {
        MSG new_message(*it);
        new_message.specie += " " + _acepted_input;
        _filtered_input.push_back(new_message);
      }
    }
    cout << "new state catn msg: " << _filtered_input.size() << endl;
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {
    cout << _acepted_input << " confluence." << endl;
    internal();
    external(mb, t);
  }

};

#endif // BOOST_SIMULATION_PDEVS_FILTER_H
