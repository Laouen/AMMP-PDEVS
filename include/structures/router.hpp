//
// Created by lao on 04/10/17.
//

#ifndef PMGBP_PDEVS_ROUTER_STRUCTURES_HPP
#define PMGBP_PDEVS_ROUTER_STRUCTURES_HPP

#include <cadmium/modeling/ports.hpp>
#include <tuple>

namespace pmgbp {
    namespace structs {
        namespace router {

            template<typename REACTANT>
            struct ports{

                struct out : public cadmium::out_port<REACTANT> {};
                struct in : public cadmium::in_port<REACTANT> {};

                using output_type=REACTANT;
                using input_type=REACTANT;

                using input_ports=std::tuple<in>;
                using output_ports=std::tuple<out>;
            };
        }
    }
}

#endif //PMGBP_PDEVS_ROUTER_STRUCTURES_HPP
