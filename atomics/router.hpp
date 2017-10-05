//
// Created by lao on 28/09/17.
//

#ifndef PMGBP_PDEVS_MODEL_ROUTER_HPP
#define PMGBP_PDEVS_MODEL_ROUTER_HPP

#include <map>

#include <cadmium/modeling/message_bag.hpp>
#include <TupleOperators.hpp>
#include <Logger.hpp>

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
        output_bags output;
        map<string, int> routing_table;
    };

    state_type state;

    router() noexcept = default;

    explicit router(const string& id) {
        this->logger.setModuleName("Router_" + id);
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

    friend std::ostringstream& operator<<(std::ostringstream& os, const typename router<PORTS,TIME>::state_type& i) {}

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
