#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include <memory>


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

    _filtered_input.clear();
  }

  TIME advance() const noexcept {

    return (_filtered_input.size() > 0) ? TIME(0) : atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {

    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      if (it->to.atModel(_model_type) == _acepted_input) {
        _filtered_input.push_back(*it);
      }
    }

  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    internal();
    external(mb, t);
  }

};

#endif // BOOST_SIMULATION_PDEVS_FILTER_H
