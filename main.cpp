/**
 * Copyright (c) 2017, Laouen Mayal Louan Belloli
 * Carleton University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <chrono>
#include <algorithm>

#include <cadmium/modeling/coupled_model.hpp>
#include <cadmium/engine/pdevs_runner.hpp>

#include <NDTime.hpp>

#include "atomics/reaction.hpp"
#include "atomics/space.hpp"
#include "atomics/router.hpp"

#include "structures/types.hpp"
#include "structures/router.hpp"

using namespace std;
using namespace pmgbp;
using hclock=chrono::high_resolution_clock;

/****** Reaction atomic model definition *********/

using reaction_ports = pmgbp::structs::reaction::ports<types::Product, types::Reactant>;

template<typename TIME>
using reaction_test_uninitialized=models::reaction<reaction_ports, TIME>;

template<class TIME>
class reaction_test : public reaction_test_uninitialized<TIME> {
public:
    reaction_test(): reaction_test_uninitialized<TIME>("../parameters.xml", "reaction_test") {};
};

/*******************************************/

/****** Router atomic model definition *********/

using router_ports = pmgbp::structs::router::ports<types::Reactant>;

template<typename TIME>
using router_test_uninitialized=models::router<router_ports, TIME>;

template<class TIME>
class router_test : public router_test_uninitialized<TIME> {
public:
    router_test(): router_test_uninitialized<TIME>("../parameters.xml", "router_test") {};
};

/*******************************************/

/****** Reaction set coupled model definition *********/

struct reaction_set_out: public cadmium::out_port<pmgbp::types::Product>{};
struct reaction_set_in: public cadmium::in_port<pmgbp::types::Reactant>{};

using iports_reaction_set = std::tuple<reaction_set_in>;
using oports_reaction_set = std::tuple<reaction_set_out>;
using submodels_reaction_set = cadmium::modeling::models_tuple<router_test, reaction_test>;
using eics_reaction_set = std::tuple<cadmium::modeling::EIC<reaction_set_in, router_test, router_ports::in>>;
using eocs_reaction_set = std::tuple<cadmium::modeling::EOC<reaction_test, reaction_ports::out, reaction_set_out>>;
using ics_reaction_set = std::tuple<cadmium::modeling::IC<router_test, router_ports::out, reaction_test, reaction_ports::in>>;

template<typename TIME>
struct reaction_set: public cadmium::modeling::coupled_model<
        TIME,
        iports_reaction_set,
        oports_reaction_set,
        submodels_reaction_set,
        eics_reaction_set,
        eocs_reaction_set,
        ics_reaction_set
> {};

/*******************************************/

/****** Space atomic model definition *********/

using space_ports = pmgbp::structs::space::ports<types::Reactant, types::Product>;

template<typename TIME>
using space_test_uninitialized=models::space<space_ports, TIME>;

template<class TIME>
class space_test : public space_test_uninitialized<TIME> {
public:
    space_test(): space_test_uninitialized<TIME>("../parameters.xml", "space_test") {};
};

/*******************************************/

/****** TOP coupled model definition *********/

using iports_top= std::tuple<>;
using oports_top= std::tuple<>;
using submodels_top = cadmium::modeling::models_tuple<space_test, reaction_set>;
using eics_top = std::tuple<>;
using eocs_top = std::tuple<>;
using ics_top = std::tuple<cadmium::modeling::IC<space_test, space_ports::inner, reaction_set , reaction_set_in>>;

template<typename TIME>
using top_model=cadmium::modeling::coupled_model<TIME, iports_top, oports_top, submodels_top, eics_top, eocs_top, ics_top>;

/*******************************************/


/*************** Loggers *******************/

namespace {

    struct oss_sink_provider{
        static std::ostream& sink(){
            return std::cout;
        }
    };
}

using log_states=cadmium::logger::logger<cadmium::logger::logger_state, cadmium::logger::verbatim_formatter, oss_sink_provider>;
using log_msg=cadmium::logger::logger<cadmium::logger::logger_messages, cadmium::logger::verbatim_formatter, oss_sink_provider>;
using log_gt=cadmium::logger::logger<cadmium::logger::logger_global_time, cadmium::logger::verbatim_formatter, oss_sink_provider>;
using logger_top=cadmium::logger::multilogger<log_states, log_msg, log_gt>;

/*******************************************/


int main() {

    auto start = hclock::now(); //to measure simulation execution time

    cadmium::engine::runner<NDTime, top_model, logger_top> r{{0}};
    r.runUntil({3000});

    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
    cout << "Simulation took:" << elapsed << "sec" << endl;
    return 0;
}