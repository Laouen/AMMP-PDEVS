//
// Created by lao on 19/09/17.
//

#ifndef PMGBP_PDEVS_SPACE_STRUCTURES_HPP
#define PMGBP_PDEVS_SPACE_STRUCTURES_HPP

#include <vector>
#include <string>
#include <cadmium/modeling/message_bag.hpp>

enum class SpaceState {
    SELECTING_FOR_REACTION = 2,
    SENDING_BIOMASS = 3,
    SENDING_REACTIONS = 4
};

struct ReactionAddress {
    std::string compartment;
    std::string reaction_set;

    ReactionAddress() = default;

    ReactionAddress(const std::string& other_compartment, const std::string& other_reaction_set) {
        this->compartment = other_compartment;
        this->reaction_set = other_reaction_set;
    }

    inline bool operator==(const ReactionAddress& o) const {
        return this->compartment == o.compartment && this->reaction_set == o.reaction_set;
    }
};


/*******************************************/
/**************** STask_t ******************/
/*******************************************/

template<class OUT_PORTS>
struct SpaceTask {
    SpaceState kind;
    typename typename cadmium::make_message_bags<OUT_PORTS>::type message_bags;

    SpaceTask() = default;

    SpaceTask(const SpaceTask<OUT_PORTS>& other) {
        this->kind = other.kind;
        this->message_bags = other.message_bags;
    }

    SpaceTask(SpaceState other_kind) {
        kind = other_kind;
    }

    inline bool operator==(const SpaceTask<OUT_PORTS>& o) const {

        bool result = (kind == o.kind);

        if ((kind == SpaceState::SENDING_REACTIONS) || (kind == SpaceState::SENDING_BIOMASS)) {
            result = result && (message_bags == o.message_bags);
        }

        return result;
    }
};


/*******************************************/
/************** End STask_t ****************/
/*******************************************/

#endif //PMGBP_PDEVS_SPACE_STRUCTURES_HPP
