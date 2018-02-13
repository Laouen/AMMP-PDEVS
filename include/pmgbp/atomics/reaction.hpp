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

#include <pmgbp/lib/Logger.hpp>
#include <pmgbp/lib/Random.hpp> // RealRandom
#include <pmgbp/lib/TaskScheduler.hpp>
#include <pmgbp/lib/TupleOperators.hpp>

#include <pmgbp/structures/types.hpp> // MetaboliteAmounts, RTask_t, Way, RTaskQueue_t
#include <pmgbp/structures/reaction.hpp>

#include <tinyxml2.h>


namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;
using namespace pmgbp::types;
using namespace pmgbp::structs::reaction;

template<class PORTS, class TIME>
class reaction_template {
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
        TIME reject_rate;

        RoutingTable<string> routing_table;
        TaskScheduler<TIME, output_bags> tasks;
    };

    state_type state;

    reaction_template() = default;

    /**
     * @brief Default constructor
     * @details Construct a new instance of a reaction model using the state passed
     * as parameter to initialize the new instance.
     *
     * @param state_other A reaction::state_type already initialized.
     */
    explicit reaction_template(const state_type& state_other) noexcept {
        this->state = state_other;
        this->logger.setModuleName("Reaction_" + this->state.id);

        // Initialize random generators
        this->initialize_random_engines();
    }

    /**
     * @brief Parser constructor
     * @details Construct a new reaction atomic model instance by opening and parsing the xml
     * file in the path xml_file.
     *
     * @param xml_file path where the xml file containing all the parameters is located.
     * @param id model id.
     */
    explicit reaction_template(const char* xml_file, const char* id) {
        this->state.id = id;
        this->logger.setModuleName("Reaction_" + this->state.id);

        // Initialize random generators
        this->initialize_random_engines();

        tinyxml2::XMLDocument doc;
        tinyxml2::XMLError opened = doc.LoadFile(xml_file);
        assert(opened == tinyxml2::XML_SUCCESS);

        tinyxml2::XMLElement* root = doc.RootElement()
                ->FirstChildElement("reactions")
                ->FirstChildElement(id);


        // Read simple state parameters
        this->state.rate = TIME(root->FirstChildElement("rate")->GetText());
        this->state.reject_rate = TIME(root->FirstChildElement("rejectRate")->GetText());
        this->state.koff_STP = std::stod(root->FirstChildElement("koffSTP")->GetText());
        this->state.koff_PTS = std::stod(root->FirstChildElement("koffPTS")->GetText());

        // Read stoichiometry by compartments
        tinyxml2::XMLElement* compartment;
        tinyxml2::XMLElement* stoichiometry_specie;
        MetaboliteAmounts  substrate_sctry, products_sctry;
        string specie_id, cid;
        int specie_amount;

        compartment = root
                ->FirstChildElement("stoichiometryByCompartments")
                ->FirstChildElement("compartment");
        while (compartment != nullptr) {

            cid = compartment->FirstChildElement("id")->GetText();

            substrate_sctry.clear();
            stoichiometry_specie = compartment->FirstChildElement("substrate");
            if (stoichiometry_specie != nullptr) {
                stoichiometry_specie = stoichiometry_specie->FirstChildElement();
            }
            while (stoichiometry_specie != nullptr) {
                specie_id = stoichiometry_specie->Attribute("id");
                specie_amount = std::stoi(stoichiometry_specie->Attribute("amount"));
                substrate_sctry.insert({specie_id, specie_amount});

                stoichiometry_specie = stoichiometry_specie->NextSiblingElement();
            }
            if (!substrate_sctry.empty()) {
                this->state.substrate_sctry.insert({cid, substrate_sctry});
                this->state.substrate_comps.insert({cid, 0});
            }

            products_sctry.clear();
            stoichiometry_specie = compartment->FirstChildElement("product");
            if (stoichiometry_specie != nullptr) {
                stoichiometry_specie = stoichiometry_specie->FirstChildElement();
            }
            while (stoichiometry_specie != nullptr) {
                specie_id = stoichiometry_specie->Attribute("id");
                specie_amount = std::stoi(stoichiometry_specie->Attribute("amount"));
                products_sctry.insert({specie_id, specie_amount});

                stoichiometry_specie = stoichiometry_specie->NextSiblingElement();
            }
            if (!products_sctry.empty()) {
                this->state.products_sctry.insert({cid, products_sctry});
                this->state.product_comps.insert({cid, 0});
            }

            compartment = compartment->NextSiblingElement();
        }

        // Read routing_table
        tinyxml2::XMLElement* routing_table;
        tinyxml2::XMLElement* entry;
        int port_number;
        string metabolite_id;

        routing_table = root->FirstChildElement("routingTable");
        entry = routing_table->FirstChildElement();
        while (entry != nullptr) {
            metabolite_id = entry->Attribute("metaboliteId");
            port_number = std::stoi(entry->Attribute("port"));
            this->state.routing_table.insert(metabolite_id, port_number);

            entry = entry->NextSiblingElement();
        }
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
        this->state.tasks.add(this->state.reject_rate, rejected_metabolites);

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

        this->logger.info("End time_advance");
        return result;
    }

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename reaction_template<PORTS,TIME>::state_type& s) {
        os << "{";
        os << "\"model_class\":\"reaction\",";
        os << "\"id\":\"" << s.id << "\",";
        os << "\"reactions_in_progress\": [";

        bool separate = false;
        for (const auto& task : s.tasks.queue()) {
            if (separate) {
                os << ",";
            }
            separate = true;

            os << "{";
            os << "\"time_left\":\"" << task.time_left << "\",";
            os << "\"reaction_amount\":" << task.task_elements.size();
            os << "}";
        }

        os << "]";
        os << "}";

        return os;
    }

private:

    /*********** Private attributes **********/

    RealRandom<double> real_random; // used for uniform random number
    Logger logger;

    /*********** Private attributes **********/

    /***************************************
    ********* helper functions *************
    ***************************************/

    void initialize_random_engines() {

        // real_random is initialized with a random generator engine
        random_device real_rd; // Random generator engine
        this->real_random.seed(real_rd());
    }

    void push_to_correct_port(string metabolite_id, output_bags& bags, const Product& m) const {
        int port_number = this->state.routing_table.at(metabolite_id);
        pmgbp::tuple::get<Product>(bags, port_number).emplace_back(m);
    }

    void bindMetabolites(typename make_message_bags<input_ports>::type mbs,
                         map<pair<string, Way>, int> &rejected) {

        for (const auto &x : get_messages<typename PORTS::in_0>(mbs)) {

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
        return (this->real_random.drawNumber(0.0, 1.0) > k);
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
                    this->push_to_correct_port(metabolite.first, bags, message);
                }
            }
        }

        if (pts_ready > 0) {

            this->removeMetabolites(state.product_comps, pts_ready);
            for (const auto &compartment_sctry : state.substrate_sctry) {

                for (const auto &metabolite : compartment_sctry.second) {
                    message.clear();
                    message.metabolites.insert({metabolite.first, pts_ready * metabolite.second});
                    this->push_to_correct_port(metabolite.first, bags, message);
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

struct reaction_ports {

    struct out_0: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_1: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_2: public cadmium::out_port<pmgbp::types::Product>{};

    struct in_0: public cadmium::in_port<pmgbp::types::Reactant>{};

    using output_type=pmgbp::types::Product;
    using input_type=pmgbp::types::Reactant;

    using output_ports=std::tuple<
        out_0,
        out_1,
        out_2
    >;
    using input_ports=std::tuple<in_0>;
};

template<typename TIME>
class reaction : public reaction_template<pmgbp::models::reaction_ports, TIME> {
public:
    reaction() = default;
    explicit reaction(const char* xml_file, const char* id) : reaction_template<pmgbp::models::reaction_ports, TIME>(xml_file, id) {}
};

}
}

#endif //PMGBP_PDEVS_MODEL_REACTION_HPP