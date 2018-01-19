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

#include <cadmium/engine/pdevs_dynamic_runner.hpp>

#include <NDTime.hpp>
#include <include/model_json_exporter.hpp>

#include "top.hpp"


using namespace std;
using namespace pmgbp;
using hclock=chrono::high_resolution_clock;


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


int main(int argc, char ** argv) {

    #ifdef DIAGRAM

        if (argc < 2) {
            cout << "Usage: " + string(argv[0]) + " <diagram_output_file>" << endl;
            exit(0);
        }

        std::ofstream file;
        file.open(argv[1]);

        if (argc == 2) {
            export_model_to_json<NDTime, coupled_cell>(file);
        } else {
            export_model_to_json<NDTime, coupled_cell>(file, atoi(argv[2]));
        }
        file.close();
    
    #else
    
        auto start = hclock::now();
        cadmium::dynamic::engine::runner<NDTime, logger_top> r(coupled_cell, NDTime({0}));
        r.run_until({3000});
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
        cout << "Simulation took:" << elapsed << "sec" << endl;
    
    #endif

    return 0;
}