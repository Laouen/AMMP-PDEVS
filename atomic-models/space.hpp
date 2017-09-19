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
#include <assert.h>

#include <cadmium/modeling/ports.hpp>
#include <cadmium/modeling/message_bag.hpp>
#include <algorithm>

#include "../libs/types.hpp" // reaction_info_t, SpaceState, Integer_t
#include "../libs/randomNumbers.hpp" // RealRandom
#include "../libs/Logger.hpp"
#include "../libs/TaskScheduler.hpp"

using namespace std;
using namespace cadmium;

#define TIME_TO_SEND_FOR_REACTION TIME(0,0,0,1) // 1 millisecond

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
        TIME                        biomass_request_time;
        Address_t                   biomass_address;
        MetaboliteAmounts           metabolites;
        map<std::string, enzyme_t>  enzymes;
        double                     volume;

        TaskScheduler<TIME, SpaceTask<MSG>> tasks;
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

        this->_state.current_time += this->_state.tasks.time_advance();

        if (this->_state.tasks.is_in_next(SpaceTask(SpaceState::SELECTING_FOR_REACTION))) {

            // advance() must be called after the is_in_next() and before to add the
            // SENDING_REACTION task
            this->_state.tasks.advance();

            // set a new task to send the selected metabolites.
            // selected_reactants = selected_reactants
            SpaceTask<MSG> selected_reactants(SpaceState::SENDING_REACTIONS);
            this->selectMetabolitesToReact(selected_reactants.msgs);
            this->mergeMessages(selected_reactants.msgs);
            if (!selected_reactants.msgs.empty()) {
                this->_state.tasks.add(TIME_TO_SEND_FOR_REACTION, selected_reactants);
            }
        } else {

            this->_state.tasks.advance();
        }


        // setting new selection
        this->setNextSelection();
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
        this->logger.info("Begin external_transition");

        bool select_biomass = false;
        bool show_metabolites = false;

        this->_state.current_time += e;
        this->_state.tasks.update(e);

        for (const auto &x : get_messages<typename defs::in>(mbs)) {
            if (x.biomass_request) {
                select_biomass = true;
            } else {
                this->addMultipleMetabolites(this->_state.metabolites, x.metabolites);
            }
        }

        if (select_biomass) this->selectForBiomass();
        this->setNextSelection();

        this->logger.info("End external_transition");
    }

    void confluence_transition(TIME e, typename make_message_bags<input_ports>::type mbs) {
        this->logger.info("Begin confluence_transition");
        internal_transition();
        external_transition(mbs, TIME::zero());
        this->logger.info("End confluence_transition");
    }

    typename make_message_bags<output_ports>::type output() const {
        this->logger.info("Begin output");

        typename std::list<MSG>::const_iterator it;
        typename vector<MSG>::iterator rt;
        typename make_message_bags<output_ports>::type bags;

        std::list<MSG> current_tasks = this->_state.tasks.next();
        for (it = current_tasks.cbegin(); it != current_tasks.cend(); ++it) {
            if (it->task_kind == SpaceState::SELECTING_FOR_REACTION) continue;

            for (rt = it->msgs.cbegin(); rt != it->msgs.cend(); ++rt) {
                cadmium::get_messages<typename defs::out>(bags).emplace_back(*rt);
            }
        }

        this->logger.info("End output");
        return bags;
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");

        TIME result = this->_state.tasks.time_advance();

        if (result <= TIME::zero()) {
            this->logger.error("Bad time: negative time: " + result);
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
        MetaboliteAmounts::iterator st;
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
                    cm.from = this->_state.id;
                    cm.react_direction = Way_t::STP;
                    cm.react_amount = 1;
                    m.push_back(cm);
                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (st = re.substrate_sctry.begin(); st != re.substrate_sctry.end(); ++st) {
                    if (this->_state.metabolites.find(st->first) != this->_state.metabolites.end()) {
                        assert(this->_state.metabolites.at(st->first) >= st->second);
                        this->_state.metabolites.at(st->first) -= st->second;
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
                    cm.from = this->_state.id;
                    cm.react_direction = Way_t::PTS;
                    cm.react_amount = 1;
                    m.push_back(cm);
                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (st = re.products_sctry.begin(); st != re.products_sctry.end(); ++st) {
                    if (this->_state.metabolites.find(st->first) != this->_state.metabolites.end()) {
                        assert(this->_state.metabolites.at(st->first) >= st->second);
                        this->_state.metabolites.at(st->first) -= st->second;
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

    // TODO test this function specially and put this in the random number class
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
            if (this->thereAreEnoughFor(it->second.substrate_sctry)) {
                threshold = this->bindingThreshold(it->second.substrate_sctry, it->second.konSTP);
                s.insert({it->first, threshold});
            } else {
                s.insert({it->first, 0});
            }

            if (it->second.reversible && this->thereAreEnoughFor(it->second.products_sctry)) {
                threshold = this->bindingThreshold(it->second.products_sctry, it->second.konPTS);
                p.insert({it->first, threshold});
            } else {
                p.insert({it->first, 0});
            }
        }
    }

    // TODO test this function specially
    double bindingThreshold(const MetaboliteAmounts &sctry, double kon) const {
        // calculation of the concentrations [A][B][C]

        double concentration = 1.0;
        MetaboliteAmounts::const_iterator it;
        for (it = sctry.cbegin(); it != sctry.cend(); ++it) {
            if (this->_state.metabolites.find(it->first) != this->_state.metabolites.end()) {
                concentration *= this->._state.metabolites.at(it->first) / (L * this->_state.volume);
            }
        }

        if (concentration == 0.0)
            return 0.0;

        return exp(-(1.0 / (concentration * kon)));
    }

    bool thereAreEnoughFor(const MetaboliteAmounts &stcry) const {
        MetaboliteAmounts::const_iterator it;
        bool local_metabolite, not_enough;

        for (it = stcry.begin(); it != stcry.end(); ++it) {
            local_metabolite = this->_state.metabolites.find(it->first) != this->_state.metabolites.end();
            not_enough = this->_state.metabolites.at(it->first) < it->second;
            if(local_metabolite && not_enough) {
                return false;
            }
        }

        return true;
    }

    /**
     * Takes all the metabolites from om with an amount grater than 0 and add them to m.
     */
    void addMultipleMetabolites(MetaboliteAmounts& m, const MetaboliteAmounts& om) {
        MetaboliteAmounts::const_iterator it;

        for (it = om.cbegin(); it != om.cend(); ++it) {
            if (m.find(it->first) != m.end()) {
                m.at(it->first) += it->second;
            } else {
                m.insert({it->first, it->second});
            }
        }
    }

    /**
     * @brief Merges all message unifying those with the same receiver address
     * @param m The non grouped messages to Unify
     */
    void mergeMessages(vector<MSG> &m) const {

        typename vector<MSG>::iterator it;
        typename map<Address_t, MSG>::iterator mt;
        map<Address_t, MSG> unMsgs;

        for (it = m.begin(); it != m.end(); ++it) {
            this->insertMessageMerging(unMsgs, *it);
        }

        m.clear();

        for (mt = unMsgs.begin(); mt != unMsgs.end(); ++mt) {
            m.push_back(mt->second);
        }
    }

    void insertMessageMerging(map<Address_t, MSG> &ms, MSG &m) const {

        if (m.react_amount > 0) {
            if (ms.find(m.to) != ms.end()) {
                ms.at(m.to).react_amount += m.react_amount;
            } else {
                ms.insert({m.to, m});
            }
        }
    }

    /**
     * @brief Looks if there is metabolites to send and in this case, if the space have
     * not already programed a selection task to send metabolites, it will program one.
     */
    void setNextSelection() {

        if ( this->thereIsMetabolites() && !this->thereIsNextSelection() ) {
            SpaceTask<MSG> selection_task(SpaceState::SELECTING_FOR_REACTION);
            this->_state.tasks.add(this->_state.internal_time, selection_task);
        }
    }

    /**
     * @brief Tells if there is or not metabolites in the space.
     * @return True if there is metabolites in the space.
     */
    bool thereIsMetabolites() const {

        MetaboliteAmounts::const_iterator it;
        for (it = this->_state.metabolites.cbegin(); it != this->_state.metabolites.cend(); ++it) {
            if (it->second > 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Looks if there is a selection task already programed.
     * @return True if there is a selection task, otherwise false.
     */
    bool thereIsNextSelection() const {
        return this->_state.tasks.exists(SpaceState::SELECTING_FOR_REACTION);
    }

    void selectForBiomass() {

        // Look for metabolites to send
        MSG cm;
        cm.to = this->_state.biomass_address;
        cm.from = this->_state.id;
        addMultipleMetabolites(cm.metabolites, this->_state.metabolites);

        // Set a new task for output() to send the selected metabolites.
        SpaceTask<MSG> send_biomass(SpaceState::SENDING_BIOMASS);
        send_biomass.msgs.push_back(cm);
        this->_state.tasks.add(this->_state.biomass_request_time, send_biomass);

        // Once the metabolite are all send to biomass, there is no more metabolites in the space.
        this->removeAllMetabolites();

    }

    void removeAllMetabolites() {

        MetaboliteAmounts::iterator it;
        for (it = this->_state.metabolites.begin(); it != this->_state.metabolites.end(); ++it) {
            it->second = 0;
        }
    }

    double sumAll(const std::map<string, double>& ons) const {
        double result = 0;

        std::map<string, double>::const_iterator it;
        for (it = ons.cbegin(); it != ons.cend(); ++it) {
            result += it->second;
        }

        return result;
    }

    void normalize(std::map<string, double>& ons, double t) {
        std::map<string, double>::const_iterator it;
        for (it = ons.begin(); it != ons.end(); ++it) {
            it->second = it->second / t;
        }
    }
};


#endif //PMGBP_PDEVS_SPACE_HPP
