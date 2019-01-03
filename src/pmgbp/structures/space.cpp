#include <pmgbp/structures/space.hpp>

std::ostream& operator<<(std::ostream& os, const pmgbp::structs::space::Status& s) {

    switch(s) {
        case pmgbp::structs::space::Status::SENDING_BIOMASS:
            os << (std::string)("SENDING_BIOMASS");
            break;
        case pmgbp::structs::space::Status::SENDING_REACTIONS:
            os << (std::string)("SENDING_REACTIONS");
            break;
        case pmgbp::structs::space::Status::SELECTING_FOR_REACTION:
            os << (std::string)("SELECTING_FOR_REACTION");
            break;
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const pmgbp::structs::space::EnzymeAddress& s) {
    os << (std::string)("compartment: ") << s.compartment;
    os << (std::string)(" - reaction set: ") << s.reaction_set;

    return os;
}