//
// Created by lao on 21/09/17.
//

#ifndef PMGBP_PDEVS_MESSAGES_HPP
#define PMGBP_PDEVS_MESSAGES_HPP

#include <string>
#include "types.hpp"

namespace message {

    struct Products {
        std::string rid;
        std::string from;
        Way_t rection_direction;
        Integer amount;

        void clear() {
            this->rid = "";
            this->from = "";
            this->amount = 0;
        }
    };
}

#endif //PMGBP_PDEVS_MESSAGES_HPP
