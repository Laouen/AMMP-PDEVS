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

#include "atomic-models/reaction.hpp"
#include "libs/message_types.hpp"

using namespace std;

using hclock=chrono::high_resolution_clock;

/****** Reaction atomic model definition *********/

using out_p = reaction_defs<string>::out;

template<typename TIME>
using reaction_metabolite=reaction<msg_event::Metabolites,TIME>;

template<class TIME>
class reaction_test : public reaction_metabolite<TIME> {
public:
    reaction_test(): reaction_metabolite<TIME>(typename reaction_metabolite<TIME>::state_type()) {};
};

/*******************************************/

/****** Reaction coupled model definition *********/

using iports = std::tuple<>;
using oports = std::tuple<>;
using submodels=cadmium::modeling::models_tuple<reaction_test>;
using eics=std::tuple<>;
using eocs=std::tuple<>;
using ics=std::tuple<>;

template<typename TIME>
using top_model=cadmium::modeling::coupled_model<TIME, iports, oports, submodels, eics, eocs, ics>;

/*******************************************/

int main() {

    auto start = hclock::now(); //to measure simulation execution time

    cadmium::engine::runner<NDTime, top_model> r{0};
    r.runUntil(3000);

    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
    cout << "Simulation took:" << elapsed << "sec" << endl;
    return 0;
}