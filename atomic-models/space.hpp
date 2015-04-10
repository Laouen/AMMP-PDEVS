#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <string>
#include <utility>
#include <map>
#include <math.h>
#include <random>

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // number e, metabolite_info_t, enzyme_info_t, SState, Integer
#include "../data-structures/randomNumbers.hpp" // IntegerRandom


using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/

template<class TIME, class MSG>
class space : public atomic<TIME, MSG>
{
private:
  TIME                              _next_internal;
  TIME                              _interval_time;
  map<string, metabolite_info_t>    _metabolites;
  map<string, enzyme_info_t>        _enzymes;
  double                            _volume;
  double                            _factor;
  SState                            _s;
  // used for uniform random numbers
  RealRandom<double>                _distribution;

public:

  explicit space(
    const TIME                            other_interval_time,
    const map<string, metabolite_info_t>& other_metabolites,
    const map<string, enzyme_info_t>&     other_enzymes,
    const double                          other_volume,
    const double                          other_factor
    ) noexcept :
  _interval_time(other_interval_time),
  _metabolites(other_metabolites),
  _enzymes(other_enzymes),
  _volume(other_volume),
  _factor(other_factor) {

    random_device rd;
    _distribution.seed(rd());

    if (this->thereIsMetabolites()) {
      
      _s              = SState::SELECTING;
      _next_internal  = _interval_time;
    } else {
      
      _s = SState::IDLE;
      _next_internal = atomic<TIME, MSG>::infinity;
    }
  }

  void internal() noexcept {

    if (_s == SState::SELECTING) {
      
      for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {
        it->second.to_send = this->weightedRandomBool(it->second.amount);
      }

      _next_internal  = TIME(0);
      _s              = SState::SENDING;
    } else if (_s == SState::SENDING) {

      for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {

        if (it->second.to_send) {
          it->second.amount  = 0;
          it->second.to_send = false;
        }
      }

      if (this->thereIsMetabolites()) {
        
        _s              = SState::SELECTING;
        _next_internal  = _interval_time;
      } else {
        
        _s              = SState::IDLE;
        _next_internal  = atomic<TIME, MSG>::infinity;
      }
    }
  }

  TIME advance() const noexcept {
    
    return _next_internal;
  }

  vector<MSG> out() const noexcept {

    vector<MSG> result                    = {};
    vector<Integer> distributed_reactants = {};
    MSG current_message;

    if (_s == SState::SENDING) {
      for (map<string, metabolite_info_t>::const_iterator it = _metabolites.cbegin(); it != _metabolites.cend(); ++it) {
        
        if (it->second.to_send) {
          distributed_reactants.clear();
          distributed_reactants.resize(it->second.enzymes.size());
          randomDistribution(distributed_reactants, it->second.amount);
          current_message.specie = it->first; 

          for (int i = 0; i < distributed_reactants.size(); ++it) {

            current_message.amount  = distributed_reactants[i];
            current_message.to      = it->second.enzymes[i];
            result.push_back(current_message);
          }
        } 
      }
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      
      this->addToMetabolites(it->specie, it->amount);
    }

    if (_next_internal != atomic<TIME, MSG>::infinity) {
      
      _next_internal = _interval_time - t;
    } else if (this->thereIsMetabolites()) {
      
      _s              = SState::SELECTING;
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

  void addToMetabolites(string n, Integer a){
    metabolite_info_t new_metabolite;

    if (a > 0) {
      if (_metabolites.find(n) != _metabolites.end()) {

        _metabolites.at(n).amount += a;
      } else {
        
        new_metabolite.amount   = a;
        new_metabolite.to_send  = false;
        new_metabolite.enzymes  = lookForEnzymes(n);
        _metabolites[n]         = new_metabolite;
      }
    }
  }

  vector<Adress> lookForEnzymes(string n) const {
    vector<Adress> result;

    for (map<string, enzyme_info_t>::const_iterator it = _enzymes.cbegin(); it != _enzymes.cend(); ++it) {
      if (belong(n, it->second.reactants)) {
        result.push_back(it->second.location);
      }
    }

    return result;
  }

  bool belong(string n, const vector<string>& ls) const {
    bool result = false;

    for (vector<string>::const_iterator it = ls.cbegin(); it != ls.cend(); ++it) {
      if (n == *it) {
        result = true;
        break;
      }
    }
    return result;
  }

  bool weightedRandomBool(Integer a) {

    double proportion = _volume / (double)a;
    double threshold  = (double)1.0 / pow( (double)e, (double)_factor*proportion );
    cout << "proportion: " << proportion << endl;
    cout << "threshold: " << threshold << endl;
    return (_distribution.drawNumber(0.0, 1.0) < threshold);
  }

  void randomDistribution(vector<Integer>& ds, Integer a) const {
    Integer current_amount;

    for (int i = 0; i < ds.size(); ++i) {
      ds[i] = 0;
    }

    for (Integer i = 0; i < a; ++i) {
      ds[rand() % ds.size()] += 1;
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
};
  

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
