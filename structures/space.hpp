//
// Created by lao on 19/09/17.
//

#ifndef PMGBP_PDEVS_SPACE_STRUCTURES_HPP
#define PMGBP_PDEVS_SPACE_STRUCTURES_HPP

#include <vector>

enum class SpaceState {
    SELECTING_FOR_REACTION = 2,
    SENDING_BIOMASS = 3,
    SENDING_REACTIONS = 4
};

/*******************************************/
/**************** STask_t ******************/
/*******************************************/

template<class MSG>
struct SpaceTask {
    SpaceState    kind;
    std::vector<MSG>   msgs;

    SpaceTask() {}

    SpaceTask(const SpaceTask<MSG>& other) {
        kind = other.kind;
        msgs = other.msgs;
    }

    SpaceTask(SpaceState other_kind) {
        kind = other_kind;
    }

    inline bool operator==(const SpaceTask<MSG>& o)  const {

        bool result = (kind == o.kind);

        if ((kind == SpaceState::SENDING_REACTIONS) || (kind == SpaceState::SENDING_BIOMASS)) {
            result = result && (msgs == o.msgs);
        }

        return result;
    }
};


/*******************************************/
/************** End STask_t ****************/
/*******************************************/

#endif //PMGBP_PDEVS_SPACE_STRUCTURES_HPP
