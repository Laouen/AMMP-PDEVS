#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include <math.h>
#include <random>
#include "../data-structures/types.hpp" // number e, metabolite_info_t, enzyme_info_t, SState, Integer


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
  SState                            _s;
  // used for
  random_device                     _rd;
  mt19937                           _generator;
  uniform_real_distribution<double> _distribution;

public:

  explicit space(
    TIME                            other_next_internal,
    TIME                            other_interval_time,
    map<string, metabolite_info_t>  other_metabolites
    map<string, enzyme_info_t>      other_enzymes,
    double                          other_volume,
    ) noexcept :
  _next_internal(other_next_internal),
  _interval_time(other_interval_time),
  _metabolites(other_metabolites),
  _enzymes(other_enzymes),
  _volume(other_volume),
  _s(SState::SELECTING),
  _rd(),
  _generator(_rd),
  _distribution(0.0,1.0) {}

  void internal() noexcept {

    if (_s == SState::SELECTING) {
      
      for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {
        it->to_send = this->weightedRandomBool(it->amount);
      }

      _next_internal  = TIME(0);
      _s              = SState::SENDING;
    } else if (_s == SState::SENDING) {

      for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {

        if (it->to_send) {
          it->amount  = 0;
          it->to_send = false;
        }
      }

      _next_internal  = _interval_time;
      _s              = SState::SELECTING;
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
          distributed_reactants.resize(it->second.enzymes.size())
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

    for (vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {
      
      this->addToMetabolites(it->specie, it->amount);
    }

    _next_internal = _interval_time - t;
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
        result true;
        break;
      }
    }
    return result;
  }

  bool weightedRandomBool(Integer a) const {

    double proportion = _volume / (double)a;
    double threshold  = (double)1.0 / pow( (double)e, (double)_factor*proportion )
    
    return _distribution(_generator) < threshold ;
  }

  void randomDistribution(vector<Integer>& ds, Integer a) {
    Integer current_amount;

    for (int i = 0; i < ds.size(); ++i) {
      ds[i] = 0;
    }

    for (Integer i = 0; i < a; ++i) {
      ds[rand() % ds.size()] += 1;
    }
  }
};
  

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
