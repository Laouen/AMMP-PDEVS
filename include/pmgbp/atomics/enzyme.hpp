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

#ifndef PMGBP_PDEVS_MODEL_ENZYME_HPP
#define PMGBP_PDEVS_MODEL_ENZYME_HPP

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
#include <pmgbp/structures/space.hpp> // EnzymeAddress
#include <pmgbp/structures/reaction.hpp>

#include <tinyxml2.h>

namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;
using namespace pmgbp::types;
using namespace pmgbp::structs::reaction;

struct enzyme_ports {

    struct out_0_product: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_1_product: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_2_product: public cadmium::out_port<pmgbp::types::Product>{};

    struct out_0_information: public cadmium::out_port<pmgbp::types::Information>{};
    struct out_1_information: public cadmium::out_port<pmgbp::types::Information>{};
    struct out_2_information: public cadmium::out_port<pmgbp::types::Information>{};


    struct in_0: public cadmium::in_port<pmgbp::types::Reactant>{};

    using product_type=pmgbp::types::Product;
    using information_type=pmgbp::types::Information;

    using output_ports=std::tuple<
            out_0_product,
            out_1_product,
            out_2_product,
            out_0_information,
            out_1_information,
            out_2_information
    >;
    using input_ports=std::tuple<in_0>;
};

template<class TIME>
class enzyme {
public:

    using rid=std::string;
    using cid=std::string;
    using sid=std::string;

    using rejected_type=map<pair<string, Way>, map<rid, Integer>>;

    using Product=typename enzyme_ports::product_type;
    using Information=typename enzyme_ports::information_type;

    using input_ports=typename enzyme_ports::input_ports;
    using output_ports=typename enzyme_ports::output_ports;

    using output_bags=typename make_message_bags<output_ports>::type;
    using input_bags=typename make_message_bags<input_ports>::type;

    /****** State and Properties *******/

    struct reaction_props_type {
        TIME rate;
        TIME reject_rate;
        double koff_STP;
        double koff_PTS;
        map<rid, MetaboliteAmounts> substrate_sctry; // the stoichiometry is separated by compartments
        map<rid, MetaboliteAmounts> products_sctry; // the stoichiometry is separated by compartments

        reaction_props_type() = default;

        reaction_props_type(const TIME& other_rate, const TIME& other_reject_rate, double other_koff_STP, double other_koff_PTS) {

            this->rate = other_rate;
            this->reject_rate = other_reject_rate;
            this->koff_STP = other_koff_STP;
            this->koff_PTS = other_koff_STP;

            this->substrate_sctry.clear();
            this->products_sctry.clear();
        }
    };

    /**
     *
     * @author Laouen Mayal Louan Belloli
     * @date 2 Jan 2019
     *
     * @struct enzyme::props_type enzyme.hpp
     *
     * @brief This struct stores all the variables related with the static properties 
     * of an enzyme atomic model.
     */
    struct props_type {
        std::string id;
        pmgbp::structs::space::EnzymeAddress location;
        TIME reject_rate; // Temporal hotFix to fast test with the same reject_rate for all the reactions
        TIME rate; // Temporal hotFix to fast test with the same rate for all the reactions
        std::map<rid, reaction_props_type> reactions;
        RoutingTable<string> routing_table;
    };

    struct reaction_state_type {
        map<cid, Integer> substrate_comps;
        map<cid, Integer> product_comps;

        reaction_state_type() {
            this->substrate_comps.clear();
            this->product_comps.clear();
        }
    };

    /**
     *
     * @author Laouen Mayal Louan Belloli
     * @date 2 Jan 2019
     *
     * @struct enzyme::state_type enzyme.hpp
     *
     * @brief This struct stores all the variables related with the state of a
     * enzyme atomic model.
     *
     */
    struct state_type {
        std::string id; // only te print then enzyme ID
        std::map<rid, reaction_state_type> reactions;
        TaskScheduler<TIME, output_bags> tasks;
    };

    state_type state;
    props_type props;

    /********** Constructors **************/

    enzyme() = default;

    /**
     * @brief Default constructor
     * @details Construct a new instance of a reaction model using the state passed
     * as parameter to initialize the new instance.
     *
     * @param state_other A reaction::state_type already initialized.
     */
    explicit enzyme(const props_type& props_other, const state_type& state_other) noexcept {
        this->props = props_other;
        this->state = state_other;
        this->logger.setModuleName("Enzyme_" + this->state.id);

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
    explicit enzyme(const char* xml_file, const char* id, pmgbp::structs::space::EnzymeAddress location) {
        // Set enzyme id
        this->props.id = id;
        this->props.location = location;
        this->state.id = id;

        // Set the module name as the enzyme and location
        std::ostringstream oss;
        oss << "Enzyme_" << this->props.id << ":" << this->props.location;
        this->logger.setModuleName(oss.str());

        // Initialize random generators
        this->initialize_random_engines();

        tinyxml2::XMLDocument doc;
        tinyxml2::XMLError opened = doc.LoadFile(xml_file);
        assert(opened == tinyxml2::XML_SUCCESS);

        // document root element
        tinyxml2::XMLElement* root = doc.RootElement();

        // Search the enzyme information
        tinyxml2::XMLElement* enzyme = root->FirstChildElement("spaces")
                ->FirstChildElement(this->props.location.compartment.c_str())
                ->FirstChildElement("enzymes")
                ->FirstChildElement(id);

        while (enzyme->FirstChildElement("address")->Attribute("esn") != this->props.location.reaction_set) {
            enzyme = enzyme->NextSiblingElement();
        }

        // Load reactions information
        tinyxml2::XMLElement* enzyme_reaction = enzyme->FirstChildElement("reactions")->FirstChildElement("reaction");
        while (enzyme_reaction != nullptr) {
            const char* reaction_id = enzyme_reaction->Attribute("id");
            tinyxml2::XMLElement* reaction = root->FirstChildElement("reactions")->FirstChildElement(reaction_id);
            this->load_reaction_props_and_state(reaction, reaction_id);

            enzyme_reaction = enzyme_reaction->NextSiblingElement();
        }
    }

    /************** PDEVS methods ********************/

    void internal_transition() {
        this->logger.info("Begin internal_transition");
        this->state.tasks.advance();
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin external_transition");

        std::ostringstream oss;
        oss << "External: ";
        pmgbp::tuple::print(oss, mbs);
        this->logger.debug(oss.str());

        this->state.tasks.update(e);

        // Inserting new accepted metabolites
        rejected_type rejected = {}; // first = STP, second = PTS
        this->bindMetabolites(mbs, rejected);

        // New task for the rejected metabolites to be send it.
        output_bags rejected_metabolites;
        sendBackRejected(rejected, rejected_metabolites);
        //TODO: reject_rate should be different for each reaction
        this->state.tasks.add(this->props.reject_rate, rejected_metabolites);

        // looking for new reactions
        output_bags products;
        this->lookForNewReactions(products);
        //TODO: rate should be different for each reaction
        this->state.tasks.add(this->props.rate, products);
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

        std::ostringstream oss;
        oss << "Output: ";
        pmgbp::tuple::print(oss, bags);
        this->logger.debug(oss.str());

        this->logger.info("End confluence_transition");
        return bags;
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");

        TIME result = this->state.tasks.time_advance();

        this->logger.info("End time_advance");
        return result;
    }

    /*************** print state *********************/

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename enzyme<TIME>::state_type& s) {
        os << "{";
        os << "\"model_class\":\"enzyme\",";
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

    RealRandom<double> real_random; // used for uniform random number generation
    Logger logger;

    /*********** Private attributes **********/

    /***************************************
    ********* helper functions *************
    ***************************************/

    void load_reaction_props_and_state(tinyxml2::XMLElement *reaction, rid reaction_id) {

        reaction_props_type new_reaction_props(
                TIME(reaction->FirstChildElement("rate")->GetText()),
                TIME(reaction->FirstChildElement("rejectRate")->GetText()),
                std::stod(reaction->FirstChildElement("koffSTP")->GetText()),
                std::stod(reaction->FirstChildElement("koffPTS")->GetText())
        );

        // TODO: Remove this temporal hotFix
        this->props.reject_rate = new_reaction_props.reject_rate;
        this->props.rate = new_reaction_props.rate;

        reaction_state_type new_reaction_state;

        // Read stoichiometry by compartments
        MetaboliteAmounts  substrate_sctry, products_sctry;
        cid compartment_id;
        sid specie_id;
        Integer specie_amount;

        tinyxml2::XMLElement* compartment_sctry = reaction->FirstChildElement("stoichiometryByCompartments")->FirstChildElement("compartmentStoichiometry");
        while (compartment_sctry != nullptr) {

            compartment_id = compartment_sctry->Attribute("cid");

            // Load substrate stoichiometry
            substrate_sctry.clear();

            tinyxml2::XMLElement* stoichiometry_specie = compartment_sctry->FirstChildElement("substrate");
            if (stoichiometry_specie != nullptr) {
                stoichiometry_specie = stoichiometry_specie->FirstChildElement();
            }

            while (stoichiometry_specie != nullptr) {
                specie_id = stoichiometry_specie->Attribute("id");
                specie_amount = Integer(std::stoi(stoichiometry_specie->Attribute("amount")));
                substrate_sctry.insert({specie_id, specie_amount});

                stoichiometry_specie = stoichiometry_specie->NextSiblingElement();
            }

            if (!substrate_sctry.empty()) {
                new_reaction_props.substrate_sctry.insert({compartment_id, substrate_sctry});
                new_reaction_state.substrate_comps.insert({compartment_id, 0});
            }

            // Load product stoichiometry
            products_sctry.clear();

            stoichiometry_specie = compartment_sctry->FirstChildElement("product");
            if (stoichiometry_specie != nullptr) {
                stoichiometry_specie = stoichiometry_specie->FirstChildElement();
            }

            while (stoichiometry_specie != nullptr) {
                specie_id = stoichiometry_specie->Attribute("id");
                specie_amount = Integer(std::stoi(stoichiometry_specie->Attribute("amount")));
                products_sctry.insert({specie_id, specie_amount});

                stoichiometry_specie = stoichiometry_specie->NextSiblingElement();
            }

            if (!products_sctry.empty()) {
                new_reaction_props.products_sctry.insert({compartment_id, products_sctry});
                new_reaction_state.product_comps.insert({compartment_id, 0});
            }

            compartment_sctry = compartment_sctry->NextSiblingElement();
        }

        this->props.reactions.insert({reaction_id, new_reaction_props});
        this->state.reactions.insert({reaction_id, new_reaction_state});

        // Add reaction metabolite addresses to the routing_table
        tinyxml2::XMLElement* entry = reaction->FirstChildElement("routingTable")->FirstChildElement();
        while (entry != nullptr) {
            if (this->props.routing_table.at(entry->Attribute("metaboliteId")) >= 0) {

                assert(this->props.routing_table.at(entry->Attribute("metaboliteId")) == std::stoi(entry->Attribute("port")));
            } else {

                this->props.routing_table.insert(
                        entry->Attribute("metaboliteId"),
                        std::stoi(entry->Attribute("port"))
                );
            }

            entry = entry->NextSiblingElement();
        }
    }

    void initialize_random_engines() {

        // real_random is initialized with a random generator engine
        random_device real_rd; // Random generator engine
        this->real_random.seed(real_rd());
    }

    void push_metabolite_to_correct_port(string metabolite_id, output_bags &bags, const Product &m) const {
        int port_number = this->props.routing_table.at(metabolite_id);
        switch (port_number) {
            case 0: std::get<0>(bags).messages.emplace_back(m); break;
            case 1: std::get<1>(bags).messages.emplace_back(m); break;
            case 2: std::get<2>(bags).messages.emplace_back(m); break;
            default:
                throw std::exception();
        }
    }

    void push_information_to_correct_port(string metabolite_id, output_bags &bags, const Information &m) const {
        int port_number = this->props.routing_table.at(metabolite_id);
        switch (port_number) {
            case 0: std::get<3>(bags).messages.emplace_back(m); break;
            case 1: std::get<4>(bags).messages.emplace_back(m); break;
            case 2: std::get<5>(bags).messages.emplace_back(m); break;
            default:
                throw std::exception();
        }
    }

    void bindMetabolites(typename make_message_bags<input_ports>::type mbs, rejected_type& rejected) {

        for (const auto &x : get_messages<typename enzyme_ports::in_0>(mbs)) {

            reaction_state_type& reaction_state = this->state.reactions.at(x.rid);
            const reaction_props_type& reaction_props = this->props.reactions.at(x.rid);

            if (x.reaction_direction == Way::STP) {

                // TODO: this step could be replaced by a single step using uniform distribution to calculate the rejected and accepted amount
                for (int i = 0; i < x.reaction_amount; ++i) {
                    if (acceptedMetabolites(reaction_props.koff_STP)) reaction_state.substrate_comps.at(x.from) += 1;
                    else increaseRejected(rejected, x.from, x.rid, Way::STP);
                }

            } else {

                // TODO: this step could be replaced by a single step using uniform distribution to calculate the rejected and accepted amount
                for (int i = 0; i < x.reaction_amount; ++i) {
                    if (acceptedMetabolites(reaction_props.koff_PTS)) reaction_state.product_comps.at(x.from) += 1;
                    else increaseRejected(rejected, x.from, x.rid, Way::PTS);
                }
            }
        }
    }

    bool acceptedMetabolites(double k) {
        return (this->real_random.drawNumber(0.0, 1.0) > k);
    }

    void increaseRejected(rejected_type& rejected, const string& compartment, const rid& reaction_id, Way direction) {

        // insert rejected information compartment and rejected reaction direction
        pair<string, Way> key = make_pair(compartment, direction);
        if(rejected.find(key) != rejected.end()) {
            if (rejected.at(key).find(reaction_id) != rejected.at(key).end()) {
                rejected.at(key).at(reaction_id) += 1;
            } else {
                rejected.at(key).insert({reaction_id, 1});
            }
        } else {
            map<rid, Integer> rejected_by_reactions;
            rejected_by_reactions.insert({reaction_id, 1});
            rejected.insert({key, rejected_by_reactions});
        }
    }

    void sendBackRejected(const rejected_type& rejected, output_bags &bags) const {

        Product message;
        for (const auto& it : rejected) {
            for (const auto& jt : it.second) {

                const reaction_props_type& reaction_props = this->props.reactions.at(jt.first);

                // Send the released enzymes
                Information informationMessage;
                informationMessage.enzyme_id = this->state.id;
                informationMessage.released_enzymes = jt.second;
                informationMessage.location = this->props.location;

                if (it.first.second == Way::STP) {
                    assert(reaction_props.substrate_sctry.find(it.first.first) != reaction_props.substrate_sctry.end());

                    // Send metabolites
                    for (const auto& metabolite : reaction_props.substrate_sctry.at(it.first.first)) {

                        message.clear();
                        message.metabolites.insert({metabolite.first, jt.second * metabolite.second});
                        this->push_metabolite_to_correct_port(metabolite.first, bags, message);
                    }

                    // Send the released enzymes
                    this->push_information_to_correct_port(reaction_props.substrate_sctry.at(it.first.first).cbegin()->first, bags, informationMessage);

                } else {
                    assert(reaction_props.products_sctry.find(it.first.first) != reaction_props.products_sctry.end());

                    // Send metabolites
                    for (const auto& metabolite : reaction_props.products_sctry.at(it.first.first)) {

                        message.clear();
                        message.metabolites.insert({metabolite.first, jt.second * metabolite.second});
                        this->push_metabolite_to_correct_port(metabolite.first, bags, message);
                    }

                    // Send the released enzymes
                    this->push_information_to_correct_port(reaction_props.products_sctry.at(it.first.first).cbegin()->first, bags, informationMessage);
                }

            }
        }
    }

    void lookForNewReactions(output_bags& bags) {

        for (const auto& reaction : this->state.reactions) {

            rid reaction_id = reaction.first;

            reaction_state_type& reaction_state = this->state.reactions[reaction_id];
            const reaction_props_type& reaction_props = this->props.reactions[reaction_id];

            Integer stp_ready = totalReadyFor(reaction_state.substrate_comps);
            Integer pts_ready = totalReadyFor(reaction_state.product_comps);

            if (stp_ready > 0) {

                this->removeMetabolites(reaction_state.substrate_comps, stp_ready);

                // Send the produced metabolites
                for (const auto &compartment_sctry : reaction_props.products_sctry) {

                    for (const auto &metabolite : compartment_sctry.second) {
                        Product metaboliteMessage;
                        metaboliteMessage.metabolites.insert({metabolite.first, stp_ready * metabolite.second});
                        this->push_metabolite_to_correct_port(metabolite.first, bags, metaboliteMessage);
                    }
                }

                // Send the released enzymes
                Information informationMessage;
                informationMessage.enzyme_id = this->state.id;
                informationMessage.released_enzymes = stp_ready;
                informationMessage.location = this->props.location;

                // Each compartment of the reactant stoichiometry will have relesed enzymes.
                // The reactant stoichiometry is used to route the released enzymes.
                for (const auto &compartment_sctry : reaction_props.substrate_sctry) {
                    this->push_information_to_correct_port(compartment_sctry.second.cbegin()->first, bags, informationMessage);
                }
            }

            if (pts_ready > 0) {

                this->removeMetabolites(reaction_state.product_comps, pts_ready);

                // Send the produced metabolites
                for (const auto &compartment_sctry : reaction_props.substrate_sctry) {

                    for (const auto &metabolite : compartment_sctry.second) {
                        Product message;
                        message.metabolites.insert({metabolite.first, pts_ready * metabolite.second});
                        this->push_metabolite_to_correct_port(metabolite.first, bags, message);
                    }
                }

                // Send the released enzymes
                Information informationMessage;
                informationMessage.enzyme_id = this->state.id;
                informationMessage.released_enzymes = pts_ready;
                informationMessage.location = this->props.location;

                // Each compartment of the reactant stoichiometry will have relesed enzymes.
                // The reactant stoichiometry is used to route the released enzymes.
                for (const auto &compartment_sctry : reaction_props.products_sctry) {
                    this->push_information_to_correct_port(compartment_sctry.second.cbegin()->first, bags, informationMessage);
                }
            }
        }
    }

    void removeMetabolites(map<string, Integer>& comp, Integer a) {

        for (auto &it : comp) {
            it.second -= a;
        }
    }

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

#endif //PMGBP_PDEVS_MODEL_ENZYME_HPP