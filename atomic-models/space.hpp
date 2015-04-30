#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <string>
#include <utility>
#include <map>
#include <math.h>
#include <random>
#include <limits>
#include <memory>

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // metabolite_info_t, enzyme_info_t, SState_t, Integer_t
#include "../data-structures/randomNumbers.hpp" // IntegerRandom_t


using namespace boost::simulation::pdevs;
using namespace std;

long double e = 2.71828182845904523536028747135266249775724709369995L;

template<class TIME, class MSG>
class space : public atomic<TIME, MSG>
{
private:
  string                            _id;
  TIME                              _next_internal;
  TIME                              _interval_time;
  map<string, metabolite_info_t>    _metabolites;
  map<string, enzyme_info_t>        _enzymes;
  double                            _volume;
  double                            _factor;
  SState_t                          _s;
  vector<MSG>                       _output;
  // used for uniform random numbers
  RealRandom_t<double>              _real_random;
  IntegerRandom_t<Integer_t>        _integer_random;
  // used to show the states
  bool                              _show_state;

public:

  explicit space(
    const string                          other_id,
    const TIME                            other_interval_time,
    const map<string, metabolite_info_t>& other_metabolites,
    const map<string, enzyme_info_t>&     other_enzymes,
    const double                          other_volume,
    const double                          other_factor
    ) noexcept :
  _id(other_id),
  _interval_time(other_interval_time),
  _metabolites(other_metabolites),
  _enzymes(other_enzymes),
  _volume(other_volume),
  _factor(other_factor),
  _show_state(false) {

    random_device real_rd;
    _real_random.seed(real_rd());
    random_device integer_rd;
    _integer_random.seed(integer_rd());

    if (this->thereIsMetabolites()) {
      
      _s              = SState_t::SELECTING;
      _next_internal  = _interval_time;
    } else {
      
      _s              = SState_t::IDLE;
      _next_internal  = atomic<TIME, MSG>::infinity;
    }
  }

  void internal() noexcept {
    
    vector<Integer_t> distributed_reactants = {};
    MSG current_message;

    if (_show_state) {
    
    _show_state = false;
    } else {
      if (_s == SState_t::SELECTING) {

        for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {

          if (this->weightedRandomBool(it->second.amount)){

            distributed_reactants.clear();
            distributed_reactants.resize(it->second.enzymes.size());
            randomDistribution(distributed_reactants, it->second.amount);
            current_message.specie = it->first; 

            for (int i = 0; i < distributed_reactants.size(); ++i) {
              
              current_message.amount  = distributed_reactants[i];
              it->second.amount       -= distributed_reactants[i];
              current_message.to      = it->second.enzymes[i];

              _output.push_back(current_message);
            }
          }
        }

        _next_internal  = TIME(0);
        _s              = SState_t::SENDING;
      } else if (_s == SState_t::SENDING) {

        _output.clear();

        if (this->thereIsMetabolites()) {
          
          _s              = SState_t::SELECTING;
          _next_internal  = _interval_time;
        } else {
          
          _s              = SState_t::IDLE;
          _next_internal  = atomic<TIME, MSG>::infinity;
        }
      }
    }
    
  }

  TIME advance() const noexcept {
    
    TIME result = _show_state ? TIME(0) : _next_internal; 
    
    return result;
  }

  vector<MSG> out() const noexcept {
    
    vector<MSG> result;

    if(_show_state) {
      for (map<string, metabolite_info_t>::const_iterator it = _metabolites.cbegin(); it != _metabolites.cend(); ++it) {

        result.push_back( MSG({"sending_output"}, it->first, it->second.amount) );
      }
    } else {

      result = _output;
    }
    
    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      
      if (isShowRequest(it->to)) _show_state = true;
      else this->addToMetabolites(it->specie, it->amount);
    }

    if (_next_internal != atomic<TIME, MSG>::infinity) {
      _next_internal = _next_internal - t;
    } else if (this->thereIsMetabolites()) {
      
      _s              = SState_t::SELECTING;
      _next_internal  = _interval_time;

    }
    
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {
    
    internal();
    external(mb, TIME(0));
    
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  void addToMetabolites(string n, Integer_t a){
    metabolite_info_t new_metabolite;

    if (a > 0) {
      if (_metabolites.find(n) != _metabolites.end()) {

        _metabolites.at(n).amount += a;
      } else {
        
        new_metabolite.amount   = a;
        new_metabolite.enzymes  = lookForEnzymes(n);
        _metabolites.insert({n, new_metabolite});
      }
    }
  }

  vector<Address_t> lookForEnzymes(string n) const {
    vector<Address_t> result;

    for (map<string, enzyme_info_t>::const_iterator it = _enzymes.cbegin(); it != _enzymes.cend(); ++it) {
      if (belong(n, it->second.reactants)) {
        result.push_back(it->second.location);
      }
    }

    return result;
  }

  bool belong(const string& n, const vector<string>& ls) const {
    bool result = false;

    for (vector<string>::const_iterator it = ls.cbegin(); it != ls.cend(); ++it) {
      if (n == *it) {
        result = true;
        break;
      }
    }
    return result;
  }

  bool weightedRandomBool(Integer_t a) {

    double proportion;

    if (a > 0)  proportion = _volume / (double)a;
    else        proportion = numeric_limits<double>::infinity();
    
    double threshold  = (double)1.0 / pow( (double)e, (double)_factor*proportion );


    return (_real_random.drawNumber(0.0, 1.0) < threshold);
  }

  void randomDistribution(vector<Integer_t>& ds, Integer_t a) {
    Integer_t current_amount;
    Integer_t reactions_amount = ds.size();
    
    if (reactions_amount > 0) {

      for (int i = 0; i < reactions_amount; ++i) {
        ds[i] = 0;
      }

      for (Integer_t i = 0; i < a; ++i) {
        ds[_integer_random.drawNumber(0, reactions_amount - 1)] += 1;
      }
    }
  }

  bool thereIsMetabolites() const {

    bool result = false;
    for (map<string, metabolite_info_t>::const_iterator it = _metabolites.cbegin(); it != _metabolites.cend(); ++it) {
      if (it->second.amount > 0) {
        result = true;
        break;
      }
    }
    return result;
  }

  bool isShowRequest(const Address_t& a) const{

    bool result = false;
    string show_requested_string = "show_state";

    for (Address_t::const_iterator it = a.cbegin(); it != a.cend(); ++it) {
      if (*it == show_requested_string) { 
        result = true;
        break;
      }
    }

    return result;
  }
};
  

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
