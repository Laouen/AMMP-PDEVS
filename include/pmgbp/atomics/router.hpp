#ifndef PMGBP_PDEVS_MODEL_ROUTER_HPP
#define PMGBP_PDEVS_MODEL_ROUTER_HPP

#include <map>

#include <cadmium/modeling/message_bag.hpp>

#include <pmgbp/lib/TupleOperators.hpp>
#include <pmgbp/lib/Logger.hpp>

#include <tinyxml2.h>

namespace pmgbp {
namespace models {

using namespace std;
using namespace cadmium;

template <typename PORTS, class TIME>
class router_template {
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

    router_template() = default;

    /**
     * @brief Parser constructor
     * @details Construct a new router atomic model instance by opening and parsing the xml
     * file in the path xml_file.
     *
     * @param xml_file path where the xml file containing all the parameters is located.
     * @param id model id.
     */
    explicit router_template(const char* xml_file, const char* id) {
        this->state.id = id;
        logger.setModuleName("Router_" + this->state.id);
        logger.info("Loading from XML");

        tinyxml2::XMLDocument doc;
        tinyxml2::XMLError opened = doc.LoadFile(xml_file);
        assert(opened == tinyxml2::XML_SUCCESS);

        tinyxml2::XMLElement* root = doc.RootElement()
                ->FirstChildElement("routers")
                ->FirstChildElement(id);

        // Load routing table
        tinyxml2::XMLElement* routing_table = root->FirstChildElement("routingTable");
        for (tinyxml2::XMLElement* entry = routing_table->FirstChildElement(); entry != nullptr; entry = entry->NextSiblingElement()) {
            string enzyme_id = entry->Attribute("enzymeID");
            int port_number = std::stoi(entry->Attribute("port"));
            this->state.routing_table.insert({enzyme_id, port_number});
        }
    }

    /********** P-DEVS functions **************/

    void internal_transition() {
        this->logger.info("Begin internal_transition");
        pmgbp::tuple::map<typename PORTS::output_type>(this->state.output, router_template::clear_bag);
        this->logger.info("End internal_transition");
    }

    void external_transition(TIME e, input_bags mbs) {
        this->logger.info("Begin external_transition");
        for (const auto &x : get_messages<typename PORTS::in_0>(mbs)) {
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

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename router_template<PORTS,TIME>::state_type& s) {
        os << "{";
        os << "\"model_class\":\"router\",";
        os << "\"id\":\"" << s.id << "\",";
        os << "\"messages\": [";
        
        int messages_in_port;
        bool separate = false;
        const auto size = std::tuple_size<decltype(s.output)>::value;
        for (int port = 0; port < size; port++) {

            messages_in_port = pmgbp::tuple::cget<typename PORTS::output_type>(s.output, port).size();
            if (messages_in_port > 0) {
                if (separate) {
                    os << ",";
                }
                separate = true;
                
                os << "{";
                os << "\"port\":" << port << ",";
                os << "\"messages_amount\":" << messages_in_port;
                os << "}";
            }
        }

        os << "]";
        os << "}";

        return os;
    }

private:

    Logger logger;

    static void clear_bag(cadmium::bag<typename PORTS::output_type>& bag) {
        bag.clear();
    }

    void push_to_correct_port(const typename PORTS::output_type& message) {
        int port_num = this->state.routing_table.at(message.enzyme_id);
        pmgbp::tuple::get<typename PORTS::output_type>(this->state.output, port_num).emplace_back(message);
    }
};

struct router_ports {

    struct out_0: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_1: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_2: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_3: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_4: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_5: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_6: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_7: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_8: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_9: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_10: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_11: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_12: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_13: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_14: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_15: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_16: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_17: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_18: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_19: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_20: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_21: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_22: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_23: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_24: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_25: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_26: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_27: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_28: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_29: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_30: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_31: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_32: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_33: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_34: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_35: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_36: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_37: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_38: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_39: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_40: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_41: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_42: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_43: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_44: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_45: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_46: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_47: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_48: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_49: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_50: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_51: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_52: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_53: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_54: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_55: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_56: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_57: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_58: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_59: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_60: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_61: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_62: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_63: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_64: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_65: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_66: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_67: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_68: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_69: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_70: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_71: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_72: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_73: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_74: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_75: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_76: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_77: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_78: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_79: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_80: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_81: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_82: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_83: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_84: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_85: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_86: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_87: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_88: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_89: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_90: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_91: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_92: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_93: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_94: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_95: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_96: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_97: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_98: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_99: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_100: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_101: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_102: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_103: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_104: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_105: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_106: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_107: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_108: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_109: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_110: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_111: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_112: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_113: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_114: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_115: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_116: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_117: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_118: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_119: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_120: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_121: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_122: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_123: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_124: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_125: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_126: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_127: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_128: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_129: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_130: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_131: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_132: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_133: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_134: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_135: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_136: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_137: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_138: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_139: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_140: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_141: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_142: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_143: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_144: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_145: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_146: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_147: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_148: public cadmium::out_port<pmgbp::types::Reactant>{};
    struct out_149: public cadmium::out_port<pmgbp::types::Reactant>{};

    struct in_0: public cadmium::in_port<pmgbp::types::Reactant>{};

    using output_type=pmgbp::types::Reactant;
    using input_type=pmgbp::types::Reactant;

    using output_ports=std::tuple<
        out_0,
        out_1,
        out_2,
        out_3,
        out_4,
        out_5,
        out_6,
        out_7,
        out_8,
        out_9,
        out_10,
        out_11,
        out_12,
        out_13,
        out_14,
        out_15,
        out_16,
        out_17,
        out_18,
        out_19,
        out_20,
        out_21,
        out_22,
        out_23,
        out_24,
        out_25,
        out_26,
        out_27,
        out_28,
        out_29,
        out_30,
        out_31,
        out_32,
        out_33,
        out_34,
        out_35,
        out_36,
        out_37,
        out_38,
        out_39,
        out_40,
        out_41,
        out_42,
        out_43,
        out_44,
        out_45,
        out_46,
        out_47,
        out_48,
        out_49,
        out_50,
        out_51,
        out_52,
        out_53,
        out_54,
        out_55,
        out_56,
        out_57,
        out_58,
        out_59,
        out_60,
        out_61,
        out_62,
        out_63,
        out_64,
        out_65,
        out_66,
        out_67,
        out_68,
        out_69,
        out_70,
        out_71,
        out_72,
        out_73,
        out_74,
        out_75,
        out_76,
        out_77,
        out_78,
        out_79,
        out_80,
        out_81,
        out_82,
        out_83,
        out_84,
        out_85,
        out_86,
        out_87,
        out_88,
        out_89,
        out_90,
        out_91,
        out_92,
        out_93,
        out_94,
        out_95,
        out_96,
        out_97,
        out_98,
        out_99,
        out_100,
        out_101,
        out_102,
        out_103,
        out_104,
        out_105,
        out_106,
        out_107,
        out_108,
        out_109,
        out_110,
        out_111,
        out_112,
        out_113,
        out_114,
        out_115,
        out_116,
        out_117,
        out_118,
        out_119,
        out_120,
        out_121,
        out_122,
        out_123,
        out_124,
        out_125,
        out_126,
        out_127,
        out_128,
        out_129,
        out_130,
        out_131,
        out_132,
        out_133,
        out_134,
        out_135,
        out_136,
        out_137,
        out_138,
        out_139,
        out_140,
        out_141,
        out_142,
        out_143,
        out_144,
        out_145,
        out_146,
        out_147,
        out_148,
        out_149
    >;

    using input_ports=std::tuple<in_0>;
};

template<typename TIME>
class router : public router_template<pmgbp::models::router_ports, TIME> {
public:
    router() = default;
    explicit router(const char* xml_file, const char* id) : router_template<pmgbp::models::router_ports, TIME>(xml_file, id) {}
};

}
}

#endif //PMGBP_PDEVS_MODEL_ROUTER_HPP
