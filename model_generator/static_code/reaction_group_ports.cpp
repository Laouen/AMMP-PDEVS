#include <cadmium/modeling/ports.hpp>
#include <pmgbp/structures/messages.hpp>

struct reaction_group_ports {
    struct out_0: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_1: public cadmium::out_port<pmgbp::types::Product>{};
    struct out_2: public cadmium::out_port<pmgbp::types::Product>{};

    struct in_0: public cadmium::in_port<pmgbp::types::Reactant>{};
};

cadmium::dynamic::modeling::Ports reaction_group_iports = { typeid(reaction_group_ports::in_0) };
cadmium::dynamic::modeling::Ports reaction_group_oports = {
        typeid(reaction_group_ports::out_0),
        typeid(reaction_group_ports::out_1),
        typeid(reaction_group_ports::out_2),
};