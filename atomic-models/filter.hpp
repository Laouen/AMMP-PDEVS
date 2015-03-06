#ifndef BOOST_SIMULATION_PDEVS_FILTER_H
#define BOOST_SIMULATION_PDEVS_FILTER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>


using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** vector types ************/
typedef vector<molecule> vectorOfMolecules;

typedef pair<string, molecule> tp;

template<class TIME, class MSG>
class filter : public atomic<TIME, MSG>
{
private:

  string      _aceptedInput;
  vector<MSG> _filteredInput;
  bool        _hasToSend;

public:

  explicit filter(string other_aceptedInput) noexcept :
  _aceptedInput(other_aceptedInput),
  _hasToSend(false) {
    _filteredInput.clear();
  }

  void internal() noexcept { 
    _hasToSend = false;
    _filteredInput.clear();
  }

  TIME advance() const noexcept {
    return _hasToSend?TIME(0):atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {

    return _filteredInput; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      
      tp input = boost::any_cast<tp>(*it); 
      if ( input.first == _aceptedInput){
        _filteredInput.push_back(boost::any(input.second));
      }
    }

    _hasToSend = (_filteredInput.size() > 0);
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
    internal();
  }

};

#endif // BOOST_SIMULATION_PDEVS_FILTER_H
