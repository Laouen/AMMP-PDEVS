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

#ifndef PMGBP_PDEVS_SPACE_HPP
#define PMGBP_PDEVS_SPACE_HPP

#include <limits>
#include <string>

#include <cadmium/modeling/ports.hpp>
#include <cadmium/modeling/message_bag.hpp>

#include "../libs/types.hpp" // reaction_info_t, SState_t, Integer_t
#include "../libs/randomNumbers.hpp" // RealRandom
#include "../libs/Logger.hpp"

using namespace std;
using namespace cadmium;

// TODO: Ports must be defined by the model generator
template <typename MSG>
struct space_defs{
    // custom ports
    struct out : public out_port<MSG> {
    };
    struct in : public in_port<MSG> {
    };
};

template <class MSG, class TIME>
class space {

    using defs=space_defs<MSG>;

public:
    /**
     * @author Laouen Mayal Louan Belloli
     *
     * @struct space::state_type space.hpp
     *
     * @brief Represents a valid P-DEVS atomic space internal model state
     */
    // TODO: Implement generic state class with state validations
    struct state_type {
        std::string                 id;
        TIME                        current_time;
        TIME                        internal_time;
        TIME                        biomass_request;
        Address_t                   biomass_address;
        SetOfMolecules_t            metabolites;
        map<std::string, enzyme_t>  enzymes;
        double                     volume;

        STaskQueue_t<TIME, MSG>     tasks;
    };

    RealRandom_t<double>       _real_random;
    IntegerRandom_t<Integer_t>  _integer_random;

    state_type                  _state;

    Logger logger;

    // ports definition
    using input_ports=std::tuple<typename defs::in>;
    using output_ports=std::tuple<typename defs::out>;

    /********* Space constructors *************/

    /**
     * @brief Default constructor
     */
    space() noexcept { this->initialize_random_engines(); }

    /**
     * @brief Constructs a new space instance using the internal state passed as parameter as the
     * initial model state.
     *
     * @param state_other - The model initial internal state.
     * @tparam state_other - space::state_type
     */
    explicit space(const state_type& state_other) noexcept {
        this->_state = state_other;
        this->initialize_random_engines();
        this->logger.setModuleName(this->_state.id);
    }

    /********* Space constructors *************/

    /********** P-DEVS functions **************/

    void internal_transition() {
        this->logger.info("Begin internal_transition");

        MSG cm;
        STask_t<TIME, MSG> sr; // sr = selected_reactants, sb = selected_biomass (TODO: check sb)

        bool reaction_selection_done = false;

        this->_state.current_time += this->_state.tasks.front().time_left;
        this->updateTaskTimeLefts(this->_state.tasks.front().time_left);

        // For all the tasks that are happening now.
        // Because The task time_lefts were updated, the current time is ZERO.
        typename STaskQueue_t<TIME, MSG>::iterator it;
        for (it = _tasks.begin(); !_tasks.empty() && (it->time_left == ZERO); it = _tasks.erase(it)) {

            if (it->task_kind != SState_t::SELECTING_FOR_REACTION) continue;

            if (!reaction_selection_done) {
                reaction_selection_done = true;

                // set a new task to send the selected metabolites.
                sr.task_kind  = SState_t::SENDING_REACTIONS;
                sr.time_left  = TIME_TO_SEND_FOR_REACTION;
                this->selectMetabolitesToReact(sr.msgs);
                this->unifyMessages(sr.msgs);
                if (!sr.msgs.empty()) this->insertTask(sr);
            }
        }

        // setting new selection
        this->setNextSelection();
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
        this->logger.info("Begin external_transition");

        bool select_biomass = false;
        bool show_metabolites = false;

        this->_state.current_time += t;
        this->updateTaskTimeLefts(t);

        for (const auto &x : get_messages<typename defs::in>(mbs)) {
            if (x.biomass_request) {
                select_biomass = true;
                continue;
            }

            this->addMultipleMetabolites(this->_state.metabolites, x.metabolites);
        }

        if (select_biomass) this->selectForBiomass();
        this->setNextSelection();

        this->logger.info("End external_transition");
    }

    void confluence_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
        this->logger.info("Begin confluence_transition");
        this->logger.info("End confluence_transition");
    }

    typename make_message_bags<output_ports>::type output() const {
        this->logger.info("Begin output");

        typename STaskQueue_t<TIME, MSG>::const_iterator it;
        typename vector<MSG>::iterator rt;
        vector<MSG> result;
        MSG b_msg;
        TIME current_time  = _tasks.front().time_left;

        // for all the tasks currently occurring. These tasks are processed now.
        for (it = _tasks.cbegin(); (it != _tasks.end()) && (it->time_left == current_time); ++it) {
            if (it->task_kind == SState_t::SELECTING_FOR_REACTION) continue;
            result.insert(result.end(), it->msgs.cbegin(), it->msgs.cend());
        }

        typename make_message_bags<output_ports>::type bags;
        for (rt = result.begin(); rt != result.end(); ++rt) {
            cadmium::get_messages<typename defs::out>(bags).emplace_back(*rt);
        }

        this->logger.info("End output");
        return bags;
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");

        TIME result;
        if (!this->_state.tasks.empty()) result = this->_state.tasks.front().time_left;
        else result = TIME("inf");

        if (result <= TIME(0)) {
            this->logger.error("Bad time: negative advance_time: " + result);
        }

        this->logger.info("End time_advance");
        return result;
    }

    /********** P-DEVS functions **************/

private:

    void initialize_random_engines() {

        // The random attributes must be initialized with a random generator variables
        random_device real_rd;
        this->_real_random.seed(real_rd());

        random_device integer_rd;
        this->_integer_random.seed(integer_rd());
    }

    void selectMetabolitesToReact(vector<MSG> &m) {
        MSG cm;
        double rv, total, partial;
        map<string, double> sons, pons;
        enzyme_t en;
        reaction_info_t re;
        vector<string> enzyme_IDs;

        // Iterators
        vector<string>::iterator it;
        SetOfMolecules_t::iterator st;
        map<string, double>::iterator i;

        // Enzyme are individually considered
        this->unfoldEnzymes(enzyme_IDs);
        // Enzymes are randomly iterated
        this->shuffleEnzymes(enzyme_IDs);

        for (it = enzyme_IDs.begin(); it != enzyme_IDs.end(); ++it) {
            partial = 0.0; total = 0.0, rv = 0.0;
            sons.clear(); pons.clear(); re.clear(); en.clear();
            en = this->_state.enzymes.at(*it);

            this->collectOns(en.handled_reactions, sons, pons);


            // sons + pons can't be greater than 1. If that happen, they are normalized
            // if sons + pons is smaller than 1, there is a chance that the enzyme does'nt react
            total = this->sumAll(sons) + this->sumAll(pons);
            if (total > 1) {
                this->normalize(sons, total);
                this->normalize(pons, total);
            }

            // The interval [0,1] is  divided in pieces:
            // {[0,son1), [son1, son1+son2),
            // ... ,
            // [son1 + ... + sonk, son1 + ... + sonk + pon1),
            // ... ,
            // [son1 + ... + sonk + pon1 + ... + ponk, 1)}
            // Depending to which intervals rv belongs, the enzyme triggers
            // the corresponding reaction or do nothing (last interval).

            rv = _real_random.drawNumber(0.0, 1.0);

            for (i = sons.begin(); i != sons.end(); ++i) {

                partial += i->second;
                if (rv < partial) {
                    // send message to trigger the reaction
                    re = en.handled_reactions.at(i->first);
                    cm.clear();
                    cm.to = re.location; // TODO: use new routing mechanism
                    cm.from = _id;
                    cm.react_direction = Way_t::STP;
                    cm.react_amount = 1;
                    m.push_back(cm);
                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (st = re.substrate_sctry.begin(); st != re.substrate_sctry.end(); ++st) {
                    if (_metabolites.find(st->first) != _metabolites.end()) {
                        assert(_metabolites.at(st->first) >= st->second);
                        _metabolites.at(st->first) -= st->second;
                    }
                }
                // once the reaction is set the enzyme was processed and it moves on to the
                // next enzyme
                continue;
            }

            // if re is empty, then no one of the STP reactions have ben triggered
            // and the search continues with the PTS reactions.
            for (i = pons.begin(); i != pons.end(); ++i) {

                partial += i->second;
                if (rv < partial) {
                    // send message to trigger the reaction
                    re = en.handled_reactions.at(i->first);
                    cm.clear();
                    cm.to = re.location;
                    cm.from = _id;
                    cm.react_direction = Way_t::PTS;
                    cm.react_amount = 1;
                    m.push_back(cm);
                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (st = re.products_sctry.begin(); st != re.products_sctry.end(); ++st) {
                    if (_metabolites.find(st->first) != _metabolites.end()) {
                        assert(_metabolites.at(st->first) >= st->second);
                        _metabolites.at(st->first) -= st->second;
                    }
                }
            }
        }
    }

    void unfoldEnzymes(vector<string>& ce) const {
        map<string, enzyme_t>::const_iterator it;

        for (it = this->_state.enzymes.cbegin(); it != this->_state.enzymes.cend(); ++it) {
            ce.insert(ce.end(), it->second.amount, it->second.id);
        }
    }

    // TODO test this function specially
    void shuffleEnzymes(vector<string>& ce) const {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(ce.begin(), ce.end(), g);
    }

    void collectOns(const map<string, reaction_info_t>& r,
                    map<string, double>& s,
                    map<string, double>& p) {

        map<string, reaction_info_t>::const_iterator it;
        double threshold;

        for (it = r.cbegin(); it != r.cend(); ++it) {

            // calculating the sons and pons
            if (this->thereAreEnaughFor(it->second.substrate_sctry)) {
                threshold = this->bindingThreshold(it->second.substrate_sctry, it->second.konSTP);
                s.insert({it->first, threshold});
            } else {
                s.insert({it->first, 0});
            }

            if (it->second.reversible && this->thereAreEnaughFor(it->second.products_sctry)) {
                threshold = this->bindingThreshold(it->second.products_sctry, it->second.konPTS);
                p.insert({it->first, threshold});
            } else {
                p.insert({it->first, 0});
            }
        }
    }

    // TODO test this function specially
    double bindingThreshold(const SetOfMolecules_t &sctry, double kon) const {
        // calculation of the concentrations [A][B][C]

        double concentration = 1.0;
        SetOfMolecules_t::const_iterator it;
        for (it = sctry.cbegin(); it != sctry.cend(); ++it) {
            if (this->_state.metabolites.find(it->first) != this->_state.metabolites.end()) {
                concentration *= this->._state.metabolites.at(it->first) / (L * _volume);
            }
        }

        if (concentration == 0.0)
            return 0.0;

        return exp(-(1.0 / (concentration * kon)));
    }
};


#endif //PMGBP_PDEVS_SPACE_HPP
