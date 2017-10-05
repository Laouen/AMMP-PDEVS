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

#ifndef PMGBP_PDEVS_MODEL_REACTION_HPP
#define PMGBP_PDEVS_MODEL_REACTION_HPP

#include <limits> // numeric_limits
#include <cstddef>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <list>
#include <map>
#include <utility> // pair
#include <algorithm> // min, max, lowe_bound
#include <memory> // shared_ptr
#include <limits> 

#include <cadmium/modeling/ports.hpp>
#include <cadmium/modeling/message_bag.hpp>

#include <types.hpp> // MetaboliteAmounts, RTask_t, Way, RTaskQueue_t
#include <Random.hpp> // RealRandom
#include <TaskScheduler.hpp>
#include <Logger.hpp>
#include <TupleOperators.hpp>

#include "../structures/reaction.hpp"

namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;
using namespace pmgbp::types;
using namespace pmgbp::structs::reaction;

template<class PORTS, class TIME>
class reaction {
public:

    using Product=typename PORTS::output_type;

    using input_ports=typename PORTS::input_ports;
    using output_ports=typename PORTS::output_ports;

    using output_bags=typename make_message_bags<output_ports>::type;
    using input_bags=typename make_message_bags<input_ports>::type;

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
        string id;
        TIME rate;
        map<string, MetaboliteAmounts> substrate_sctry; // the stoichiometry is separated by compartments
        map<string, MetaboliteAmounts> products_sctry; // the stoichiometry is separated by compartments
        map<string, Integer> substrate_comps;
        map<string, Integer> product_comps;
        double koff_STP;
        double koff_PTS;
        TIME reject_time;

        RoutingTable<string> routing_table;
        TaskScheduler<TIME, output_bags> tasks;
    };

    RealRandom<double> real_random; // used for uniform random number
    state_type state;

    reaction() noexcept = default;

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
        this->logger.setModuleName("Reaction_" + this->state.id);
    }

    void internal_transition() {
        this->logger.info("Begin internal_transition");
        this->state.tasks.advance();
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin external_transition");

        this->state.tasks.update(e);

        // Inserting new accepted metabolites
        map<pair<string, Way>, int> rejected = {}; // first = STP, second = PTS
        this->bindMetabolites(mbs, rejected);

        // New task for the rejected metabolites to be send it.
        output_bags rejected_metabolites;
        sendBackRejected(rejected, rejected_metabolites);
        this->state.tasks.add(this->state.reject_time, rejected_metabolites);

        // looking for new reactions
        output_bags products;
        this->lookForNewReactions(products);
        this->state.tasks.add(this->state.rate, products);
        this->logger.info("End external_transition");
    }

    void confluence_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin confluence_transition");
        internal_transition();
        external_transition(TIME::zero(), mbs);
        this->logger.info("End confluence_transition");
    }

    output_bags output() const {
        this->logger.info("Begin output");

        output_bags bags;

        list<output_bags> outputs = this->state.tasks.next();
        for (const auto &current_bags: outputs) {
            pmgbp::tuple::merge(bags, current_bags);
        }

        this->logger.info("End confluence_transition");
        return bags;
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");

        TIME result = this->state.tasks.time_advance();

//        if (result < TIME::zero()) {
//            this->logger.error("Bad time: negative time: " + result);
//        }

        this->logger.info("End time_advance");
        return result;
    }

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename reaction<PORTS,TIME>::state_type& i) {}

private:

    Logger logger;

    /***************************************
    ********* helper functions *************
    ***************************************/

    void push_to_correct_port(string metabolite_id, output_bags& bags, const Product& m) const {
        int port_number = this->state.routing_table.at(metabolite_id);
        pmgbp::tuple::get<Product>(bags, port_number).emplace_back(m);
    }

    void bindMetabolites(typename make_message_bags<input_ports>::type mbs,
                         map<pair<string, Way>, int> &rejected) {

        for (const auto &x : get_messages<typename PORTS::in>(mbs)) {

            if (x.reaction_direction == Way::STP) {

                // TODO: replace this by a single step using uniform distribution to calculate the rejected and accepted amount
                for (int i = 0; i < x.reaction_amount; ++i) {
                    if (acceptedMetabolites(state.koff_STP)) state.substrate_comps.at(x.from) += 1;
                    else increaseRejected(rejected, x.from, Way::STP);
                }
            } else {

                // TODO: replace this by a single step using uniform distribution to calculate the rejected and accepted amount
                for (int i = 0; i < x.reaction_amount; ++i) {
                    if (acceptedMetabolites(state.koff_PTS)) state.product_comps.at(x.from) += 1;
                    else increaseRejected(rejected, x.from, Way::PTS);
                }
            }
        }
    }

    bool acceptedMetabolites(double k) {
        return (real_random.drawNumber(0.0, 1.0) > k);
    }

    void increaseRejected(map<pair<string, Way>, int>& rejected, const string& compartment, Way direction) {

        // insert rejected information compartment and rejected reaction direction
        pair<string, Way> key = make_pair(compartment, direction);
        if(rejected.find(key) != rejected.end()) {
            rejected.at(key) += 1;
        } else {
            rejected.insert({key, 1});
        }
    }

    void sendBackRejected(const map<pair<string, Way>, int> &rejected, output_bags &bags) const {

        Product message;
        for ( const auto &it : rejected) {

            if (it.first.second == Way::STP) {
                assert(state.substrate_sctry.find(it.first.first) != state.substrate_sctry.end());

                for (const auto &metabolite : state.substrate_sctry.at(it.first.first)) {

                    message.clear();
                    message.metabolites.insert({metabolite.first, it.second*metabolite.second});
                    this->push_to_correct_port(metabolite.first, bags, message);
                }
            } else {
                assert(state.products_sctry.find(it.first.first) != state.products_sctry.end());

                for (const auto &metabolite : state.products_sctry.at(it.first.first)) {

                    message.clear();
                    message.metabolites.insert({metabolite.first, it.second*metabolite.second});
                    this->push_to_correct_port(metabolite.first, bags, message);
                }
            }
        }
    }

    void lookForNewReactions(output_bags& bags) {

        Product message;
        Integer stp_ready = totalReadyFor(state.substrate_comps);
        Integer pts_ready = totalReadyFor(state.product_comps);

        if (stp_ready > 0) {

            this->removeMetabolites(state.substrate_comps, stp_ready);
            for (const auto &compartment_sctry : state.products_sctry) {

                for (const auto &metabolite : compartment_sctry.second) {
                    message.clear();
                    message.metabolites.insert({metabolite.first, stp_ready * metabolite.second});
                    this->push_to_correct_port(compartment_sctry.first, bags, message);
                }
            }
        }

        if (pts_ready > 0) {

            this->removeMetabolites(state.product_comps, pts_ready);
            for (const auto &compartment_sctry : state.substrate_sctry) {

                for (const auto &metabolite : compartment_sctry.second) {
                    message.clear();
                    message.metabolites.insert({metabolite.first, pts_ready * metabolite.second});
                    this->push_to_correct_port(compartment_sctry.first, bags, message);
                }
            }
        }
    }

    void removeMetabolites(map<string, Integer>& comp, Integer a) {

        for (auto &it : comp) {
            it.second -= a;
        }
    }

    // TODO: test this function specially
    Integer totalReadyFor(const map<string, Integer>& comp) {

        Integer result = comp.cbegin()->second;
        for (const auto &it : comp) {
            if (result > it.second) {
                result = it.second;
            }
        }

        return result;
    }

};
}
}

#endif //PMGBP_PDEVS_MODEL_REACTION_HPP