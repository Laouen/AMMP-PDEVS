//
// Created by lao on 28/09/17.
//

#ifndef PMGBP_PDEVS_MODEL_ROUTER_HPP
#define PMGBP_PDEVS_MODEL_ROUTER_HPP

#include <map>

#include <cadmium/modeling/message_bag.hpp>
#include <TupleOperators.hpp>
#include <Logger.hpp>

#include <tinyxml2.h>

namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;

template <typename PORTS, class TIME>
class router {
public:
    using input_ports=typename PORTS::input_ports;
    using output_ports=typename PORTS::output_ports;

    using output_bags=typename make_message_bags<output_ports>::type;
    using input_bags=typename make_message_bags<input_ports>::type;

    struct state_type {
        string id;
        output_bags output;
        map<string, int> routing_table;
    };

    state_type state;

    router() noexcept = default;

    /**
     * @brief Parser constructor
     * @details Construct a new router atomic model instance by opening and parsing the xml
     * file in the path xml_file.
     *
     * @param xml_file path where the xml file containing all the parameters is located.
     * @param id model id.
     */
    explicit router(const char* xml_file, const char* id) {
        this->state.id = id;
        logger.setModuleName("Router" + this->state.id);

        tinyxml2::XMLDocument doc;
        tinyxml2::XMLError opened = doc.LoadFile(xml_file);
        assert(opened == tinyxml2::XML_SUCCESS);

        tinyxml2::XMLElement* root = doc.RootElement()
                ->FirstChildElement("routers")
                ->FirstChildElement(id);

        // Read routing table
        tinyxml2::XMLElement* routing_table = root->FirstChildElement("routingTable");
        tinyxml2::XMLElement* entry = routing_table->FirstChildElement();
        string key;
        int port_number;
        while (entry != nullptr) {
            key = entry->Attribute("metaboliteId");
            port_number = std::stoi(entry->Attribute("port"));
            this->state.routing_table.insert({key, port_number});

            entry = entry->NextSiblingElement();
        }
    }

    /********** P-DEVS functions **************/

    void internal_transition() {
        this->logger.info("Begin internal_transition");
        pmgbp::tuple::map<typename PORTS::output_type>(this->state.output, router::clear_bag);
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin external_transition");
        for (const auto &x : get_messages<typename PORTS::in>(mbs)) {
            this->push_to_correct_port(x);
        }
        this->logger.info("End external_transition");
    }

    void confluence_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin confluence_transition");
        internal_transition();
        external_transition(TIME::zero(), mbs);
        this->logger.info("End confluence_transition");
    }

    TIME time_advance() const {
        this->logger.info("Begin time_advance");
        TIME next_internal = TIME::infinity();
        if (!pmgbp::tuple::empty(this->state.output)) {
            next_internal = TIME::zero();
        }
        this->logger.info("End time_advance");
        return next_internal;
    }

    output_bags output() const {
        this->logger.info("Begin output");
        this->logger.info("End output");
        return this->state.output;
    }

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename router<PORTS,TIME>::state_type& s) {
        os << "Router_" << s.id << " ";
        const auto size = std::tuple_size<decltype(s.output)>::value;

        os << "Messages per port number: ";
        for (int port_number = 0; port_number < size; port_number++) {
            os << "(" << port_number << ", ";
            os << pmgbp::tuple::cget<typename PORTS::output_type>(s.output, port_number).size() << ")";
        }
    }

private:

    Logger logger;

    static void clear_bag(cadmium::bag<typename PORTS::output_type>& bag) {
        bag.clear();
    }

    void push_to_correct_port(const typename PORTS::output_type& message) {
        int port_num = this->state.routing_table.at(message.rid);
        pmgbp::tuple::get<typename PORTS::output_type>(this->state.output, port_num).emplace_back(message);
    }
};

}
}

#endif //PMGBP_PDEVS_MODEL_ROUTER_HPP
