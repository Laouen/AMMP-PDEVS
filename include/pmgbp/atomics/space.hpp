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

#ifndef PMGBP_PDEVS_MODEL_SPACE_HPP
#define PMGBP_PDEVS_MODEL_SPACE_HPP

#include <limits>
#include <string>
#include <cassert>
#include <algorithm>
#include <cmath>

#include <cadmium/modeling/message_bag.hpp>

#include <tinyxml2.h>

#include <pmgbp/lib/Random.hpp> // RealRandom
#include <pmgbp/lib/Logger.hpp>
#include <pmgbp/lib/TaskScheduler.hpp>
#include <pmgbp/lib/TupleOperators.hpp>

#include <pmgbp/structures/types.hpp> // ReactionInfo, Integer, RoutingTable
#include <pmgbp/structures/space.hpp> // Status, Task

#define TIME_TO_SEND_FOR_REACTION TIME({0,0,0,1}) // 1 millisecond
namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;
using namespace pmgbp::types;
using namespace pmgbp::structs::space;

    /**
     * @author Laouen Mayal Louan Belloli
     *
     * @struct space::space space.hpp
     *
     * @brief Represents a valid P-DEVS atomic space model
     *
     * @typedef PORTS the ports struct ...
     * @typedef TIME The type of the time class
     *
     */
template<typename PORTS, class TIME>
class space {

public:

    using Reactant=typename PORTS::reactant_type;

    using input_ports=typename PORTS::input_ports;
    using output_ports=typename PORTS::output_ports;

    using output_bags=typename make_message_bags<output_ports>::type;
    using input_bags=typename make_message_bags<input_ports>::type;

    struct state_type {
        string id;
        TIME interval_time;
        MetaboliteAmounts metabolites;
        map<string, Enzyme> enzymes;
        RoutingTable<EnzymeAddress> routing_table;
        long double volume;
        
        //TODO: Take the volume from the .xml parameter
        //long double volume = 0.0000000000000000001;

        TaskScheduler<TIME, Task<output_ports>> tasks;
    };

    state_type state;

    /********* Space constructors *************/

    space() = default;

    /**
     * @brief Constructs a new space atomic model instance using the internal state passed
     * as parameter as the initial model state.
     *
     * @param state_other - The model initial internal state.
     * @tparam state_other - space::state_type
     */
    explicit space(const state_type &state_other) noexcept {
        this->state = state_other;
        this->logger.setModuleName("Space_" + this->state.id);

        // Initialize random generators
        this->initialize_random_engines();
    }

    /**
     * @brief Parser constructor
     * @details Construct a new space atomic model instance by opening and parsing the xml
     * file in the path xml_file.
     *
     * @param xml_file path where the xml file containing all the parameters is located.
     * @param id model id.
     */
    explicit space(const char* xml_file, const char* id) {
        this->state.id = id;
        logger.setModuleName("Space_" + this->state.id);
        logger.debug("Loading from XML");

        // Initialize random generators
        this->initialize_random_engines();

        tinyxml2::XMLDocument doc;
        tinyxml2::XMLError opened = doc.LoadFile(xml_file);
        assert(opened == tinyxml2::XML_SUCCESS);

        tinyxml2::XMLElement* root = doc.RootElement()
                ->FirstChildElement("spaces")
                ->FirstChildElement(id);


        // Load interval_time
        this->state.interval_time = TIME(root->FirstChildElement("intervalTime")->GetText());

        // Load metabolites
        tinyxml2::XMLElement* metabolites = root->FirstChildElement("metabolites");
        tinyxml2::XMLElement* metabolite = metabolites->FirstChildElement();
        while (metabolite != nullptr) {
            string specie = metabolite->Attribute("id");
            Integer amount = Integer(std::stoi(metabolite->Attribute("amount")));
            this->state.metabolites.insert({specie, amount});
            metabolite = metabolite->NextSiblingElement();
        }

        // Load enzymes
        string specie_id;
        Integer specie_amount;

        tinyxml2::XMLElement* enzymes = root->FirstChildElement("enzymes");
        for (tinyxml2::XMLElement* enzyme_entry = enzymes->FirstChildElement(); enzyme_entry != nullptr; enzyme_entry = enzyme_entry->NextSiblingElement()) {

            string enzyme_id = enzyme_entry->Value();
            logger.debug("Loading enzyme " + enzyme_id);

            // Load the enzyme address
            tinyxml2::XMLElement* address = enzyme_entry->FirstChildElement("address");
            EnzymeAddress enzyme_location(address->Attribute("cid"), address->Attribute("esn"));

            Integer enzyme_amount = Integer(std::stoi(enzyme_entry->Attribute("amount")));

            // Load handled reactions
            map<string, ReactionInfo> handled_reactions;
            tinyxml2::XMLElement* reactions = enzyme_entry->FirstChildElement("reactions");
            for (tinyxml2::XMLElement* handled_reaction = reactions->FirstChildElement(); handled_reaction != nullptr; handled_reaction = handled_reaction->NextSiblingElement()) {

                string reaction_id = handled_reaction->Attribute("id");
                logger.debug("Loading enzyme " + enzyme_id + " reaction " + reaction_id);

                tinyxml2::XMLElement* reaction_parameters = doc.RootElement()
                        ->FirstChildElement("reactions")
                        ->FirstChildElement(reaction_id.c_str());

                // Load the reaction stoichiometry
                tinyxml2::XMLElement* compartment_sctry = reaction_parameters
                        ->FirstChildElement("stoichiometryByCompartments")
                        ->FirstChildElement();

                while (compartment_sctry != nullptr && std::string(compartment_sctry->Attribute("cid")) != this->state.id) {
                    compartment_sctry = compartment_sctry->NextSiblingElement();
                }

                // If the reaction don't consume nor produce any metabolite from the compartment
                // the compartment must not handle the reaction
                if (compartment_sctry == nullptr) continue;

                double konSTP = std::stod(reaction_parameters->FirstChildElement("konSTP")->GetText());
                double konPTS = std::stod(reaction_parameters->FirstChildElement("konPTS")->GetText());
                double koffSTP = std::stod(reaction_parameters->FirstChildElement("koffSTP")->GetText());
                double koffPTS = std::stod(reaction_parameters->FirstChildElement("koffPTS")->GetText());
                bool reversible =reaction_parameters->FirstChildElement("reversible")->GetText() == "true";

                // Load the substrate stoichiometry
                MetaboliteAmounts substrate_sctry;
                tinyxml2::XMLElement* substrate = compartment_sctry->FirstChildElement("substrate");
                if (substrate != nullptr) {
                    for (tinyxml2::XMLElement* species = substrate->FirstChildElement(); species != nullptr; species = species->NextSiblingElement()) {
                        specie_id = species->Attribute("id");
                        specie_amount = Integer(std::stoi(species->Attribute("amount")));
                        substrate_sctry.insert({specie_id, specie_amount});
                    }
                }

                // Load the product stoichiometry
                MetaboliteAmounts products_sctry;
                tinyxml2::XMLElement* product = compartment_sctry->FirstChildElement("product");
                if (product != nullptr) {
                    for (tinyxml2::XMLElement* species = product->FirstChildElement(); species != nullptr; species = species->NextSiblingElement()) {
                        specie_id = species->Attribute("id");
                        specie_amount = Integer(std::stoi(species->Attribute("amount")));
                        products_sctry.insert({specie_id, specie_amount});
                    }
                }

                // Save the reaction information
                ReactionInfo reaction_information = ReactionInfo(reaction_id,
                                                                 substrate_sctry,
                                                                 products_sctry,
                                                                 konSTP,
                                                                 konPTS,
                                                                 koffSTP,
                                                                 koffPTS,
                                                                 reversible);

                handled_reactions.insert({reaction_id, reaction_information});
            }

            // If all the enzyme reactions are not related with the compartment, then, the compartment
            // mustn't handle the enzyme at all.
            if (!handled_reactions.empty()) {
                Enzyme enzyme(enzyme_id, enzyme_location, enzyme_amount, handled_reactions);
                this->state.enzymes.insert({enzyme_id + ":" + enzyme.location.str(), enzyme});
            }

            logger.debug("Loaded enzyme " + enzyme_id);
        }

        logger.debug("Loading routing table");
        // Load routing_table
        tinyxml2::XMLElement* routing_table;
        tinyxml2::XMLElement* entry;
        int port_number;

        routing_table = root->FirstChildElement("routingTable");
        entry = routing_table->FirstChildElement();
        while (entry != nullptr) {
            EnzymeAddress enzyme_address(entry->Attribute("cid"), entry->Attribute("esn"));
            port_number = std::stoi(entry->Attribute("port"));
            this->state.routing_table.insert(enzyme_address, port_number);
            entry = entry->NextSiblingElement();
        }
    }

    /********* Space constructors *************/

    /********** P-DEVS functions **************/

    void internal_transition() {
        this->logger.info("Begin internal_transition");

        if (this->state.tasks.is_in_next(Task<output_ports>(Status::SELECTING_FOR_REACTION))) {

            // advance() must be called after the is_in_next() and before
            // selecting new metabolites to react
            // Status::SENDING_REACTION task
            this->state.tasks.advance();

            // set a new task to send the selected metabolites.
            // selected_reactants = selected_reactants
            Task<output_ports> selected_reactants(Status::SENDING_REACTIONS);
            this->selectMetabolitesToReact(selected_reactants.message_bags);
            if (!pmgbp::tuple::empty(selected_reactants.message_bags)) {
                pmgbp::tuple::map(selected_reactants.message_bags, space::mergeMessages);
                this->state.tasks.add(TIME_TO_SEND_FOR_REACTION, selected_reactants);
            }
        } else {

            this->state.tasks.advance();
        }

        // setting new selection
        this->setNextSelection();
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin external_transition");

        std::ostringstream oss;
        oss << "External: ";
        pmgbp::tuple::print(oss, mbs);
        this->logger.debug(oss.str());

        bool show_metabolites = false;

        this->state.tasks.update(e);

        // Receive new metabolites
        for (const auto &x : get_messages<typename PORTS::in_0_product>(mbs)) {
            this->addMultipleMetabolites(this->state.metabolites, x.metabolites);
        }

        // Receive released enzymes
        for (const auto &x : get_messages<typename PORTS::in_0_information>(mbs)) {
            string eid = x.enzyme_id + ":" + x.location.str();
            this->state.enzymes.at(eid).amount += x.released_enzymes;
        }

        this->setNextSelection();

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

        list<Task<output_ports>> current_tasks = this->state.tasks.next();
        for (const auto &task : current_tasks) {
            if (task.kind == Status::SELECTING_FOR_REACTION) continue;
            pmgbp::tuple::merge(bags, task.message_bags);
        }

        std::ostringstream oss;
        oss << "Output: ";
        pmgbp::tuple::print(oss, bags);
        this->logger.debug(oss.str());

        this->logger.info("End output");
        return bags;
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");

        TIME result = this->state.tasks.time_advance();

        if (result == TIME::infinity()) {
            result = this->state.interval_time;
        }

        this->logger.info("End time_advance");
        return result;
    }

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename space<PORTS,TIME>::state_type& s) {
        bool separate = false;

        os << "{";
        os << "\"model_class\":\"space\",";
        os << "\"id\":\"" << s.id << "\",";
        os << "\"enzymes\": [";
        for(const auto& enzyme : s.enzymes) {
            if (separate) {
                os << ",";
            }
            separate = true;
            os << "{";
            os << "\"id\":\"" << enzyme.second.id << "\",";
            os << "\"amount\":" << enzyme.second.amount;
            os << "}";
        }
        os << "],";

        separate = false;
        os << "\"metabolites\": [";
        for(const auto& metabolite : s.metabolites) {
            if (separate) {
                os << ", ";
            }
            separate = true;
            os << "{";
            os << "\"id\":\"" << metabolite.first << "\",";
            os << "\"amount\":" << metabolite.second;
            os << "}";
        }

        os << "]";
        os << "}";

        return os;
    }

    /********** P-DEVS functions **************/

private:

    /*********** Private attributes **********/

    RealRandom<double> real_random;
    IntegerRandom<Integer> integer_random;
    Logger logger;

    /*********** Private attributes **********/

    void initialize_random_engines() {

        // The random attributes must be initialized with a random generator variables
        random_device real_rd;
        this->real_random.seed(real_rd());

        random_device integer_rd;
        this->integer_random.seed(integer_rd());
    }

    void push_to_correct_port(EnzymeAddress address, output_bags& bags, const Reactant& p) {
        int port_number = this->state.routing_table.at(address);
        pmgbp::tuple::get<Reactant>(bags, port_number).emplace_back(p);
    }

    void selectMetabolitesToReact(output_bags& bags) {

        Reactant reactant;
        double rv, total, partial;
        map<string, double> sons, pons;
        ReactionInfo re;
        vector<string> unfolded_enzymes;

        // Enzyme are individually considered
        this->unfoldEnzymes(unfolded_enzymes);

        // Enzymes are randomly iterated
        this->shuffleEnzymes(unfolded_enzymes);

        for (auto &eid : unfolded_enzymes) {

            partial;
            sons.clear();
            pons.clear();
            re.clear();

            Enzyme& enzyme = this->state.enzymes.at(eid);

            this->collectOns(enzyme.handled_reactions, sons, pons);

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

            rv = real_random.drawNumber(0.0, 1.0);

            partial = 0.0;
            for (auto &son : sons) {

                partial += son.second;
                if (rv < partial) {

                    // send message to trigger the reaction
                    re = enzyme.handled_reactions.at(son.first);
                    reactant.clear();
                    reactant.rid = re.id;
                    reactant.enzyme_id = enzyme.id;
                    reactant.from = this->state.id;
                    reactant.reaction_direction = Way::STP;
                    reactant.reaction_amount = 1;
                    this->push_to_correct_port(enzyme.location, bags, reactant);

                    // update enzyme amount
                    enzyme.amount--;

                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (auto &metabolite : re.substrate_sctry) {
                    if (this->state.metabolites.find(metabolite.first) !=
                        this->state.metabolites.end()) {
                        assert(this->state.metabolites.at(metabolite.first) >= metabolite.second);
                        this->state.metabolites.at(metabolite.first) -= metabolite.second;
                    }
                }
                // once the reaction is set the enzyme was processed and it moves on to the
                // next enzyme
                continue;
            }

            // TODO(Lao): By adding the reaction direction to the Pons and Sons I can merge both lists
            // and then just iterate over a sngle list instead of this.
            // if re is empty, then no one of the STP reactions have ben triggered
            // and the search continues with the PTS reactions.
            for (auto &pon : pons) {

                partial += pon.second;
                if (rv < partial) {

                    // send message to trigger the reaction
                    re = enzyme.handled_reactions.at(pon.first);
                    reactant.clear();
                    reactant.rid = re.id;
                    reactant.enzyme_id = enzyme.id;
                    reactant.from = this->state.id;
                    reactant.reaction_direction = Way::PTS;
                    reactant.reaction_amount = 1;
                    this->push_to_correct_port(enzyme.location, bags, reactant);

                    // update enzyme amount
                    enzyme.amount--;

                    break;
                }
            }

            // update the metabolite amount in the space
            if (!re.empty()) {
                for (auto &metabolite : re.products_sctry) {
                    if (this->state.metabolites.find(metabolite.first) !=
                        this->state.metabolites.end()) {
                        assert(this->state.metabolites.at(metabolite.first) >= metabolite.second);
                        this->state.metabolites.at(metabolite.first) -= metabolite.second;
                    }
                }
            }
        }
    }

    void unfoldEnzymes(vector<string> &ce) const {
        for (const auto &enzyme : this->state.enzymes) {
            ce.insert(ce.end(), enzyme.second.amount, enzyme.first);
        }
    }

    // TODO: test this function specially and put this in the random number class
    void shuffleEnzymes(vector<string> &ce) const {
        random_device rd;
        mt19937 g(rd());
        shuffle(ce.begin(), ce.end(), g);
    }

    void collectOns(const map<string, ReactionInfo>& reactions, map<string, double>& son, map<string, double>& pon) {

        long double threshold;

        for (const auto &reaction : reactions) {

            // calculating the sons and pons
            if (this->thereAreEnoughFor(reaction.second.substrate_sctry)) {
                threshold = this->bindingThreshold(reaction.second.substrate_sctry,
                                                   reaction.second.konSTP);
                son.insert({reaction.first, threshold});
            } else {
                son.insert({reaction.first, 0});
            }

            if (reaction.second.reversible && this->thereAreEnoughFor(reaction.second.products_sctry)) {
                threshold = this->bindingThreshold(reaction.second.products_sctry,
                                                   reaction.second.konPTS);
                pon.insert({reaction.first, threshold});
            } else {
                pon.insert({reaction.first, 0});
            }
        }
    }

    // TODO: test this function specially
    long double bindingThreshold(const MetaboliteAmounts &sctry, double kon) const {
        // concentrations calculation [A][B][C]

        double concentration = 1.0;
        for (const auto &metabolite : sctry) {
            if (this->state.metabolites.find(metabolite.first) != this->state.metabolites.end()) {
                Integer amount = this->state.metabolites.at(metabolite.first);
                concentration *= amount / (L * this->state.volume);
            }
        }

        if (concentration == 0.0)
            return 0.0;
        return 0.2;
        // exp(-(1.0 / (concentration * kon)));
    }

    bool thereAreEnoughFor(const MetaboliteAmounts &stcry) const {
        bool is_local, not_enough, there_is_enough;
        there_is_enough = false;

        for (const auto &metabolite : stcry) {
            is_local = this->state.metabolites.find(metabolite.first) != this->state.metabolites.end();
            not_enough = this->state.metabolites.at(metabolite.first) < metabolite.second;
            if (is_local && not_enough) {
                return false;
            }
            if (is_local) {
                there_is_enough = true;
            }
        }
        return there_is_enough;
    }

    /**
     * Takes all the metabolites from om with an amount grater than 0 and add them to m.
     */
    void addMultipleMetabolites(MetaboliteAmounts &m, const MetaboliteAmounts &om) {

        for (const auto &metabolite : om) {
            if (m.find(metabolite.first) != m.end()) {
                m.at(metabolite.first) += metabolite.second;
            } else {
                m.insert({metabolite.first, metabolite.second});
            }
        }
    }

    /**
     * @brief Merges all message unifying those with the same receiver address
     * @param messages The non grouped messages to Unify
     */
    static void mergeMessages(cadmium::bag<Reactant> &messages) {
        map<string, Reactant> merged_messages;

        for (auto &product : messages) {
            space::insertMessageMerging(merged_messages, product);
        }

        messages.clear();

        for (auto &merged_product : merged_messages) {
            messages.emplace_back(merged_product.second);
        }
    }

    static void insertMessageMerging(std::map<string, Reactant>& ms, Reactant &m) {

        if (m.reaction_amount > 0) {
            if (ms.find(m.rid) != ms.end()) {
                ms.at(m.rid).reaction_amount += m.reaction_amount;
            } else {
                ms.insert({m.rid, m});
            }
        }
    }

    /**
     * @brief Looks if there is metabolites to send and in this case, if the space have
     * not already programed a selection task to send metabolites, it will program one.
     */
    void setNextSelection() {

        if (this->thereIsMetabolites() && !this->thereIsNextSelection()) {
            Task<output_ports> selection_task(Status::SELECTING_FOR_REACTION);
            this->state.tasks.add(this->state.interval_time, selection_task);
        }
    }

    /**
     * @brief Tells if there is or not metabolites in the space.
     * @return True if there is metabolites in the space.
     */
    bool thereIsMetabolites() const {

        for (const auto &metabolite : this->state.metabolites) {
            if (metabolite.second > 0) {
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
        return this->state.tasks.exists(Task<output_ports >(Status::SELECTING_FOR_REACTION));
    }

    double sumAll(const map<string, double> &ons) const {
        double result = 0;

        for (const auto &on : ons) {
            result += on.second;
        }
        return result;
    }

    void normalize(map<string, double> &ons, double t) {
        for (auto &on : ons) {
            on.second /= t;
        }
    }
};

}
}

#endif //PMGBP_PDEVS_MODEL_SPACE_HPP
