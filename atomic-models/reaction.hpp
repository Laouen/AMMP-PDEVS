#ifndef BOOST_SIMULATION_PDEVS_REACTION_H
#define BOOST_SIMULATION_PDEVS_REACTION_H
#include <string>
#include <vector>
#include <list>
#include <map>
#include <utility> // pair
#include <algorithm> // min, max, lowe_bound
#include <memory> // shared_ptr
#include <limits> // numeric_limits

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // SetOfMolecules_t, RTask_t, Way_t, RTaskQueue_t
#include "../data-structures/randomNumbers.hpp" // RealRandom

using namespace boost::simulation::pdevs;
using namespace boost::simulation;
using namespace std;

BRITime ZERO(0);

template<class TIME, class MSG>
class reaction : public pdevs::atomic<TIME, MSG>
{

private:
  // enzyme information
  string                                _id;
  shared_ptr< map<string, Address_t> >  _addresses;
  bool                                  _reversible;
  TIME                                  _rate;
  SetOfMolecules_t                      _substrate_sctry;
  SetOfMolecules_t                      _products_sctry;
  map<string, Integer_t>                _substrate_comps;
  map<string, Integer_t>                _product_comps;
  TIME                                  _it;
  RTaskQueue_t<TIME>                    _tasks;
  // used for uniform random number 
  IntegerRandom_t<Integer_t>            _distribution;

public:

  explicit reaction(
    const string                                other_id,
    const shared_ptr< map<string, Address_t> >  other_addresses,
    const bool                                  other_reversible,
    const TIME                                  other_rate,
    const SetOfMolecules_t&                     other_substrate_sctry,
    const SetOfMolecules_t&                     other_products_sctry,
    const map<string, Integer_t>                other_substrate_comps,
    const map<string, Integer_t>                other_product_comps,
    const TIME                                  other_it
  ) noexcept :
  _id(other_id),
  _addresses(other_addresses),
  _reversible(other_reversible),
  _rate(other_rate),
  _substrate_sctry(other_substrate_sctry),
  _products_sctry(other_products_sctry),
  _substrate_comps(other_substrate_comps),
  _product_comps(other_product_comps),  
  _it(other_it) {

    random_device rd;
    _distribution.seed(rd());

  }

  void internal() noexcept {

    // Updating time left
    this->updateTaskTimeLefts(_tasks.front().time_left);

    // removing all the task made in the las out function
    this->removeFinishedTasks();
  }

  TIME advance() const noexcept {
  
    if (!_tasks.empty()) return _tasks.front().time_left;
    else                 return pdevs::atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {
    
    MSG new_message;
    const SetOfMolecules_t* curr_sctry;

    vector<MSG> result = {};
    TIME current_time  = _tasks.front().time_left;

    for (typename RTaskQueue_t<TIME>::const_iterator it = _tasks.cbegin(); (it != _tasks.cend()) && (it->time_left == current_time); ++it) {

      if (it->direction == Way_t::STP) curr_sctry = &_products_sctry;
      else curr_sctry = &_substrate_sctry;

      for (SetOfMolecules_t::const_iterator jt = curr_sctry->cbegin(); jt != curr_sctry->cend(); ++jt) {
        new_message.clear();
        new_message.to = _addresses->at(jt->first);
        new_message.metabolites.insert({jt->first, it->amount*jt->second});
        result.push_back(new_message); 
      }
    }

    unifyMessages(result);
    
    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {
    // Updating time left
    this->updateTaskTimeLefts(t);

    // inserting new metabolites
    this->bindMetabolits(mb);

    // looking for new reactions
    this->lookForNewReactions();
  }

  virtual void confluence(const vector<MSG>& mb, const TIME& t) noexcept {
    
    internal();
    external(mb, ZERO);
    
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  // Decrease the time left of all the current tasks in _tasks by the parameter t.
  void updateTaskTimeLefts(TIME t){

    for (typename RTaskQueue_t<TIME>::iterator it = _tasks.begin(); it != _tasks.end(); ++it) {
      it->time_left -= t;
    }
  }

  void removeFinishedTasks() {

    while(!_tasks.empty() && (_tasks.front().time_left == ZERO)) {
      _tasks.pop_front();
    }
  }

  // It take the needed number of each specie in mb and rejected (by calling addRejected) the not needed part.
  void bindMetabolits(const vector<MSG>& mb) {

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {

      if (it->react_direction == Way_t::STP) {
        _substrate_comps.at(it->from) += it->react_amount;
      } else {
        _product_comps.at(it->from) += it->react_amount;
      }
    }
  }

  void lookForNewReactions() {
    
    Integer_t stp_ready = totalReadyFor(_substrate_comps); 
    Integer_t pts_ready = totalReadyFor(_product_comps);

    if (stp_ready > 0) {
      RTask_t<TIME> stp_task(_rate, Way_t::STP, stp_ready);
      this->insertTask(stp_task);
      this->removeMetabolites(_substrate_comps, stp_ready);
    }

    if (pts_ready > 0) {
      RTask_t<TIME> pts_task(_rate, Way_t::STP, pts_ready);
      this->insertTask(pts_task);
      this->removeMetabolites(_product_comps, pts_ready);
    }
  }

  void removeMetabolites(map<string, Integer_t>& comp, Integer_t a) {

    for (map<string, Integer_t>::iterator it = comp.begin(); it != comp.end(); ++it) {
      it->second -= a;  
    }
  }

  // TODO test this function specially
  Integer_t totalReadyFor(const map<string, Integer_t>& comp) {
    
    Integer_t result = numeric_limits<Integer_t>::max();
    for (map<string, Integer_t>::const_iterator it = comp.cbegin(); it != comp.cend(); ++it) {
      if (result > it->second) result = it->second;
    }

    return result;
  }

  void insertTask(const RTask_t<TIME>& t) {

    typename STaskQueue_t<TIME>::iterator it = lower_bound(_tasks.begin(), _tasks.end(), t);
    _tasks.insert(it, t);
  }

  void unifyMessages(vector<MSG>& m) const {

    map<Address_t, MSG> unMsgs;

    for (vector<MSG>::iterator it = m.begin(); it != m.end(); ++it) {
      insertMessage(unMsgs, *it);
    }

    m.clear();

    for (map<Address_t, MSG>::iterator it = unMsgs.begin(); it != unMsgs.end(); ++it) {
      m.push_back(it->second);
    }
  }

  void insertMessage(map<Address_t, MSG>& ms, MSG& m) const {

    if (ms.find(m.to) != ms.end()) {
      addMultipleMetabolites(ms.at(m.to).metabolites, m.metabolites);
    } else {
      ms.insert({m.to, m}); // TODO: change all the initializer_list because they don't work on windows
    }
  }

  void addMultipleMetabolites(SetOfMolecules_t& m, const SetOfMolecules_t& om) {
  
    for (SetOfMolecules_t::const_iterator it = om.cbegin(); it != om.cend(); ++it) {
      addMetabolite(m, it->first, it->second);
    }
  }

  void addMetabolite(SetOfMolecules_t& m, string n, Integer_t a){

    if (a > 0) {
      if (m.find(n) != m.end()) {
        m.at(n) += a;
      } else {
        m.insert({n, a}); // TODO: change all the initializer_list because they don't work on windows
      }
    }
  }
};

#endif // BOOST_SIMULATION_PDEVS_REACTION_H
