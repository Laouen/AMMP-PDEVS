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
#include <cadmium/logger/dynamic_common_loggers.hpp>

#include <NDTime.hpp>
#include <model_json_exporter.hpp>

#include <memore/logger.hpp>

#include "top.hpp"


using namespace std;
using namespace pmgbp;
using hclock=chrono::high_resolution_clock;


/*************** Loggers *******************/

using log_states=cadmium::logger::logger<cadmium::logger::logger_state, cadmium::dynamic::logger::formatter<NDTime>, cadmium::logger::cout_sink_provider>;
using log_msg=cadmium::logger::logger<cadmium::logger::logger_messages, cadmium::dynamic::logger::formatter<NDTime>, cadmium::logger::cout_sink_provider>;
using log_gt=cadmium::logger::logger<cadmium::logger::logger_global_time, cadmium::dynamic::logger::formatter<NDTime>, cadmium::logger::cout_sink_provider>;

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
            export_model_to_json<NDTime, generate_model()>(file);
        } else {
            export_model_to_json<NDTime, generate_model()>(file, atoi(argv[2]));
        }
        file.close();
    
    #else
        
        std::cout << "hclock" << std::endl;
        auto start = hclock::now();

        std::cout << "model" << std::endl;
        std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> top_model = generate_model();
        
        std::cout << "runner" << std::endl;
        cadmium::dynamic::engine::runner<NDTime, logger_top> r(top_model, NDTime({0}));
        
        std::cout << "run_until 3000" << std::endl;
        r.run_until({3000});
        std::cout << "finished" << std::endl;
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
        cout << "Simulation took:" << elapsed << "sec" << endl;
    
    #endif

    return 0;
}