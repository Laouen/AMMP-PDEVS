/**
 * Copyright (c) 2017, Laouen Mayal Louan Belloli
 * Carleton University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PMGBBP_PDEVS_REACTION_ATOMIC_MODEL_HPP
#define PMGBBP_PDEVS_REACTION_ATOMIC_MODEL_HPP

#include <limits> // numeric_limits
#include <cstddef>
#include <string>
#include <fstream>

#include <vector>
#include <list>
#include <map>
#include <utility> // pair
#include <algorithm> // min, max, lowe_bound
#include <memory> // shared_ptr
#include <limits> 

#include <cadmium/modeling/ports.hpp>
#include <cadmium/modeling/message_bag.hpp>

#include "../structures/types.hpp" // MetaboliteAmounts, RTask_t, Way_t, RTaskQueue_t
#include "../libs/randomNumbers.hpp" // RealRandom

using namespace std;
using namespace cadmium;

template<typename MSG>
struct reaction_defs{
    // custom ports
    struct out : public out_port<MSG> {
    };
    struct in : public in_port<MSG> {
    };
};

template<class MSG, class TIME>
class reaction {

  using defs=reaction_defs<MSG>;

public:

  /**
   * 
   * @author Laouen Mayal Louan Belloli
   * @date 16 May 2017
   * 
   * @struct reaction::state_type reaction.hpp
   * 
   * @brief This struct stores all the variables related with the state of a 
   * reaction atomic model.
   * 
   */
  struct state_type{
      string                                id;
      shared_ptr<map<string, Address_t>>    addresses; // TODO: This must be replaced with the use of ports
      bool                                  reversible; // TODO: check where is this field used
      TIME                                  rate;
      map<string, MetaboliteAmounts>         substrate_sctry; // the stoichiometry is separated by compartments
      map<string, MetaboliteAmounts>         products_sctry; // the stoichiometry is separated by compartments
      map<string, Integer>                substrate_comps;
      map<string, Integer>                product_comps;
      double                                koff_STP;
      double                                koff_PTS;
      TIME                                  interval_time; // TODO: Check where is this field used
      TIME                                  reject_time;

      RTaskQueue_t<TIME, MSG>               tasks;
  };

  reaction() noexcept {
    // real_random is initialized with a random generator engine
    random_device real_rd; // Random generator engine
    real_random.seed(real_rd());
  }

  /**
   * @brief Default constructor
   * @details Construct a new instance of a reaction model using the state passed
   * as parameter to initialize the new instance.
   * 
   * @param initialized_state A reaction::state_type already initialized.
   */
  explicit reaction(const state_type& initialized_state) noexcept {

    this->state = initialized_state;

    // real_random is initialized with a random generator engine
    random_device real_rd; // Random generator engine
    real_random.seed(real_rd());
  }

  RealRandom<double> real_random; // used for uniform random number
  state_type state;

  // ports definition
  using input_ports=std::tuple<typename defs::in>;
  using output_ports=std::tuple<typename defs::out>;

  void internal_transition() {
    comment("reaction::internal function.");

    // Updating time left
    this->updateTaskTimeLefts(state.tasks.front().time_left);

    // Removing all the task made in the las out function
    this->removeFinishedTasks();
  }

  void external_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
    comment("reaction::external_transition function.");

    // Updating time left
    this->updateTaskTimeLefts(e);

    // Inserting new acepted metabolites 
    map<string, pair<int, int>> rejected = {}; // first = STP, second = PTS
    this->bindMetabolits(mbs, rejected);

    // Puting the rejected metabolites in msg to be send it.
    vector<MSG> rejected_task;
    collectInMessage(rejected, rejected_task);
    unifyMessages(rejected_task);

    // Adding a task for the rejected metabolites
    RTask_t<TIME, MSG> new_task(state.reject_time, rejected_task);
    insertTask(new_task);

    // looking for new reactions
    this->lookForNewReactions();
  }

  void confluence_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
    comment("reaction::confluence function");
    internal_transition();
    external_transition(TIME(0), mbs);
  }

  typename make_message_bags<output_ports>::type output() const {
      comment("reaction::output function.");

      map<string, MetaboliteAmounts>::const_iterator jt;
      typename RTaskQueue_t<TIME, MSG>::const_iterator it;
      typename vector<MSG>::iterator mt

      typename make_message_bags<output_ports>::type bags;
      MSG new_message;
      const map<string, MetaboliteAmounts>* curr_sctry;

      vector<MSG> result = {};
      TIME current_time = state.tasks.front().time_left;

      for (it = state.tasks.cbegin(); (it != state.tasks.cend()) && (it->time_left == current_time); ++it) {
          if (it->task_kind == RState_t::REACTING) {

              if (it->direction == Way_t::STP) curr_sctry = &state.products_sctry;
              else curr_sctry = &state.substrate_sctry;

              for (jt = curr_sctry->cbegin(); jt != curr_sctry->cend(); ++jt) {

                  for (MetaboliteAmounts::const_iterator mt = jt->second.cbegin(); mt != jt->second.cend(); ++mt) {
                        new_message.clear();
                        new_message.to = state.addresses->at(mt->first); // TODO: change this by ports routing
                        new_message.metabolites.insert({mt->first, it->amount*mt->second});
                        result.push_back(new_message);
                  }
              }
          } else {
              result.insert(result.end(), it->toSend.begin(), it->toSend.end());
          }
      }

      this->unifyMessages(result);
      for (mt = result.begin(); mt != result.end(); ++mt) {
            cadmium::get_messages<typename defs::out>(bags).emplace_back(*mt);
      }
    
    return bags;
  }

  TIME time_advance() const {
    comment("reaction::time_advance function.");
    
    TIME result;
    if (!state.tasks.empty()) result = state.tasks.front().time_left;
    else result = TIME("inf");
    
    return result;
  }

  friend std::ostringstream& operator<<(std::ostringstream& os, const typename reaction<MSG,TIME>::state_type& i) {}

  /***************************************
  ********* helper functions *************
  ***************************************/

  void comment(string msg) const {
    cout << "[reaction " << state.id << "] " << msg << endl;
  }

  // Decreases the time left of all the current tasks in state.tasks by the parameter t.
  void updateTaskTimeLefts(TIME t){

    for (typename RTaskQueue_t<TIME, MSG>::iterator it = state.tasks.begin(); it != state.tasks.end(); ++it) {
      it->time_left -= t;
    }
  }

  void removeFinishedTasks() {

    while(!state.tasks.empty() && (state.tasks.front().time_left <= TIME(0))) {
      state.tasks.pop_front();
    }
  }

  void bindMetabolits(typename make_message_bags<input_ports>::type mbs, map<string, pair<int, int>>& r) {

    for (const auto &x : get_messages<typename defs::in>(mbs)) {

      if (x.react_direction == Way_t::STP) {

        // TODO: replace this by a single step using uniform distribution to calculate the rejected and accepted amount
        for (int i = 0; i < x.react_amount; ++i) {
          if (aceptedMetabolites(state.koff_STP)) state.substrate_comps.at(x.from) += 1;
          else increaseRejected(r, x.from, Way_t::STP); // r.first = STP, r.second = PTS
        }
      } else {

        // TODO: replace this by a single step using uniform distribution to calculate the rejected and accepted amount
        for (int i = 0; i < x.react_amount; ++i) {
          if (aceptedMetabolites(state.koff_PTS)) state.product_comps.at(x.from) += 1;
          else increaseRejected(r, x.from, Way_t::PTS); // r.first = STP, r.second = PTS
        }
      }
    }
  }

  bool aceptedMetabolites(double k) {

    return (real_random.drawNumber(0.0, 1.0) > k);
  }

  void increaseRejected(map<string, pair<int, int>>& r, string f, Way_t w) {

    if (w == Way_t::STP) {

      if(r.find(f) != r.end()) {
        r.at(f).first += 1;
      } else {
        r.insert({f, make_pair(1,0)});
      }
    } else {

      if(r.find(f) != r.end()) {
        r.at(f).second += 1;
      } else {
        r.insert({f, make_pair(0,1)});
      }
    }
  }

  void collectInMessage(const map<string, pair<int, int>>& r, vector<MSG>& ts) const {

    MSG m;
    for (map<string, pair<int, int> >::const_iterator it = r.cbegin(); it != r.cend(); ++it) {

      if (it->second.first > 0) {
        assert(state.substrate_sctry.find(it->first) != state.substrate_sctry.end());
        
        for (MetaboliteAmounts::const_iterator jt = state.substrate_sctry.at(it->first).cbegin(); jt != state.substrate_sctry.at(it->first).cend(); ++jt) {
          
          m.clear();
          m.to = state.addresses->at(jt->first);
          m.metabolites.insert({jt->first, it->second.first*jt->second});
          ts.push_back(m); 
        }
      }

      if (it->second.second > 0) {
        assert(state.products_sctry.find(it->first) != state.products_sctry.end());
        
        for (MetaboliteAmounts::const_iterator jt = state.products_sctry.at(it->first).cbegin(); jt != state.products_sctry.at(it->first).cend(); ++jt) {
          
          m.clear();
          m.to = state.addresses->at(jt->first);
          m.metabolites.insert({jt->first, it->second.second*jt->second});
          ts.push_back(m); 
        }
      }
    }
  }

  void lookForNewReactions() {
    
    Integer stp_ready = totalReadyFor(state.substrate_comps);
    Integer pts_ready = totalReadyFor(state.product_comps);

    if (stp_ready > 0) {
      RTask_t<TIME, MSG> stp_task(state.rate, Way_t::STP, stp_ready);
      this->insertTask(stp_task);
      this->removeMetabolites(state.substrate_comps, stp_ready);
    }

    if (pts_ready > 0) {
      RTask_t<TIME, MSG> pts_task(state.rate, Way_t::STP, pts_ready);
      this->insertTask(pts_task);
      this->removeMetabolites(state.product_comps, pts_ready);
    }
  }

  void removeMetabolites(map<string, Integer>& comp, Integer a) {

    for (map<string, Integer>::iterator it = comp.begin(); it != comp.end(); ++it) {
      it->second -= a;  
    }
  }

  // TODO: test this function specially
  Integer totalReadyFor(const map<string, Integer>& comp) {

    Integer result = comp.cbegin()->second;
    for (map<string, Integer>::const_iterator it = comp.cbegin(); it != comp.cend(); ++it) {
      if (result > it->second) result = it->second;
    }

    return result;
  }

  void insertTask(const RTask_t<TIME, MSG>& t) {

    typename RTaskQueue_t<TIME, MSG>::iterator it = lower_bound(state.tasks.begin(), state.tasks.end(), t);
    state.tasks.insert(it, t);
  }

  void unifyMessages(vector<MSG>& m) const {

    map<Address_t, MSG> unMsgs;

    for (typename vector<MSG>::iterator it = m.begin(); it != m.end(); ++it) {
      insertMessage(unMsgs, *it);
    }

    m.clear();

    for (typename map<Address_t, MSG>::iterator it = unMsgs.begin(); it != unMsgs.end(); ++it) {
      m.push_back(it->second);
    }
  }

  void insertMessage(map<Address_t, MSG>& ms, MSG& m) const {

    if (ms.find(m.to) != ms.end()) {
      addMultipleMetabolites(ms.at(m.to).metabolites, m.metabolites);
    } else {
      ms.insert({m.to, m});
    }
  }

  void addMultipleMetabolites(MetaboliteAmounts& m, const MetaboliteAmounts& om) const {
  
    for (MetaboliteAmounts::const_iterator it = om.cbegin(); it != om.cend(); ++it) {
      addMetabolite(m, it->first, it->second);
    }
  }

  void addMetabolite(MetaboliteAmounts& m, string n, Integer a) const {

    if (a > 0) {
      if (m.find(n) != m.end()) {
        m.at(n) += a;
      } else {
        m.insert({n, a});
      }
    }
  }
};

#endif