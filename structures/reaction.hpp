//
// Created by lao on 28/09/17.
//

#ifndef PMGBP_PDEVS_REACTION_STRUCTURE_HPP
#define PMGBP_PDEVS_REACTION_STRUCTURE_HPP

#include <cadmium/modeling/message_bag.hpp>
#include <cadmium/modeling/ports.hpp>
#include "types.hpp"

namespace pmgbp {
namespace structs {
namespace reaction {

template<typename PRODUCT, typename REACTANT>
struct ports{

    struct out : public cadmium::out_port<PRODUCT> {};
    struct in : public cadmium::in_port<REACTANT> {};

    using output_type=PRODUCT;
    using input_type=REACTANT;

    using input_ports=std::tuple<typename in>;
    using output_ports=std::tuple<typename out>
};

/*******************************************/
/****************** Task *******************/
/*******************************************/

template<class OUTPUT_PORTS>
struct Task {

    typename cadmium::make_message_bags<OUTPUT_PORTS>::type message_bags;

    Task() = default;

    Task(const Task<OUTPUT_PORTS>& other) {
        this->message_bags = other.message_bags;
    }

    inline bool operator==(const Task<OUTPUT_PORTS>& o) const {
        return (kind == o.kind) && (message_bags == o.message_bags);
    }
};


/*******************************************/
/**************** End Task *****************/
/*******************************************/

}
}
}

#endif //PMGBP_PDEVS_REACTION_STRUCTURE_HPP
