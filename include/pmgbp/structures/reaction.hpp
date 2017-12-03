//
// Created by lao on 28/09/17.
//

#ifndef PMGBP_PDEVS_REACTION_STRUCTURE_HPP
#define PMGBP_PDEVS_REACTION_STRUCTURE_HPP

#include <cadmium/modeling/message_bag.hpp>
#include <cadmium/modeling/ports.hpp>

namespace pmgbp {
namespace structs {
namespace reaction {

template<typename PRODUCT, typename REACTANT>
struct ports{

    struct out : public cadmium::out_port<PRODUCT> {};
    struct in : public cadmium::in_port<REACTANT> {};

    using output_type=PRODUCT;
    using input_type=REACTANT;

    using input_ports=std::tuple<in>;
    using output_ports=std::tuple<out>;
};

}
}
}

#endif //PMGBP_PDEVS_REACTION_STRUCTURE_HPP
