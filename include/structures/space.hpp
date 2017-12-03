//
// Created by lao on 19/09/17.
//

#ifndef PMGBP_PDEVS_SPACE_STRUCTURES_HPP
#define PMGBP_PDEVS_SPACE_STRUCTURES_HPP

#include <vector>
#include <string>

#include <cadmium/modeling/message_bag.hpp>
#include <cadmium/modeling/ports.hpp>

#include <TupleOperators.hpp>

namespace pmgbp {
namespace structs {
namespace space {


enum class Status {
    SELECTING_FOR_REACTION = 2,
    SENDING_BIOMASS = 3,
    SENDING_REACTIONS = 4
};

struct ReactionAddress {
    std::string compartment;
    std::string reaction_set;

    ReactionAddress() {
        this->compartment = "";
        this->reaction_set = "";
    };

    ReactionAddress(const std::string& other_compartment, const std::string& other_reaction_set) {
        this->compartment = other_compartment;
        this->reaction_set = other_reaction_set;
    }

    ReactionAddress(const ReactionAddress& other) {
        this->compartment = other.compartment;
        this->reaction_set = other.reaction_set;
    }

    inline bool operator==(const ReactionAddress& o) const {
        return this->compartment == o.compartment && this->reaction_set == o.reaction_set;
    }

    inline bool operator<(const ReactionAddress& o) const {
        return (this->compartment < o.compartment) ||
                ((this->compartment == o.compartment) && (this->reaction_set < o.reaction_set));
    }

    bool empty() const {
        return this->compartment.empty() && this->reaction_set.empty();
    }

    void clear() {
        this->compartment = "";
        this->reaction_set = "";
    }
};


/**********************************************/
/***************** PORTS **********************/
/**********************************************/
template<class REACTANT, class PRODUCT>
struct ports {

    struct inner : public cadmium::out_port<REACTANT> {};
    struct in : public cadmium::in_port<PRODUCT> {};

    using output_type=REACTANT;
    using input_type=PRODUCT;

    using input_ports=std::tuple<in>;
    using output_ports=std::tuple<inner>;
};

/*******************************************/
/****************** Task *******************/
/*******************************************/

template<class OUTPUT_PORTS>
struct Task {
    Status kind;
    typename cadmium::make_message_bags<OUTPUT_PORTS>::type message_bags;

    Task() = default;

    Task(const Task<OUTPUT_PORTS>& other) {
        this->kind = other.kind;
        this->message_bags = other.message_bags;
    }

    explicit Task(Status other_kind) {
        kind = other_kind;
    }

    inline bool operator==(const Task<OUTPUT_PORTS>& o) const {

        bool result = (kind == o.kind);

        if ((kind == Status::SENDING_REACTIONS) || (kind == Status::SENDING_BIOMASS)) {
            result = result && pmgbp::tuple::equals(message_bags, o.message_bags);
        }

        return result;
    }
};

/*******************************************/
/**************** End Task *****************/
/*******************************************/

}
}
}

std::ostream& operator<<(std::ostream& os, const pmgbp::structs::space::Status& s);
std::ostream& operator<<(std::ostream& os, const pmgbp::structs::space::ReactionAddress& s);

#endif //PMGBP_PDEVS_SPACE_STRUCTURES_HPP
