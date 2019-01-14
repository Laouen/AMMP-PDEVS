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
#include <dynamic_json_exporter.hpp>

#include <memore/logger.hpp>

#include "top.hpp"


/*************** Loggers *******************/

#ifdef MEMORE

#include <memore/sink.hpp>
#include <mongocxx/instance.hpp>

/************** MeMoRe sink ******************/

// NOTE: it is necessary to create a unique instance of the mongocxx::instance in order to use the memore::recorder
mongocxx::instance instance{};

namespace {  

    memore::sink memore_sink("pmgbp", "pmgbp", "simulation_results");
    
    struct sink_provider {
        static memore::sink& sink() {
            return memore_sink;
        }
    };
}

#else

/************** cout sink ******************/
namespace {
    struct sink_provider {
        static std::ostream& sink() {
            return std::cout;
        }
    };
}
#endif

using log_states=cadmium::logger::logger<cadmium::logger::logger_state, memore::logger::formatter<NDTime>, sink_provider>;
using log_msg=cadmium::logger::logger<cadmium::logger::logger_messages, memore::logger::formatter<NDTime>, sink_provider>;
using log_gt=cadmium::logger::logger<cadmium::logger::logger_global_time, memore::logger::formatter<NDTime>, sink_provider>;

using logger_top=cadmium::logger::multilogger<log_states, log_msg, log_gt>;

/*******************************************/

using namespace std;
using namespace pmgbp;
using hclock=chrono::high_resolution_clock;

int main(int argc, char ** argv) {

    #ifdef DIAGRAM

        if (argc != 3) {
            std::cout << "Usage: " + std::string(argv[0]) + "<xml_parameters_path> <diagram_output_file>" << std::endl;
            exit(0);
        }

        std::string xml_parameters_path = std::string(argv[1]);
        const char * json_file = argv[2];

        std::ofstream file;
        file.open(json_file);
        std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> top_model = generate_model(xml_parameters_path);
        dynamic_export_model_to_json<NDTime>(file, top_model);
        file.close();

    #else

        if (argc != 3) {
            std::cout << "Usage: " + std::string(argv[0]) + " <xml_parameters_path> <simulation_db_identifier>" << std::endl;
            exit(0);
        }
        
        std::string xml_parameters_path = std::string(argv[1]);
        const char * simulation_db_identifier = argv[2];

        #ifdef MEMORE
        // New custom collection used so Django or other platform can set the desired collection name to retrieve results
        sink_provider::sink().new_collection(simulation_db_identifier);
        #endif
        
        // Initialize model
        auto start = hclock::now();
        
        std::cout << "generate_model" << std::endl;
        std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> top_model = generate_model(xml_parameters_path);
        std::cout << "create runner" << std::endl;
        cadmium::dynamic::engine::runner<NDTime, logger_top> r(top_model, NDTime({0}));        
        
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
        cout << "Model initialization took:" << elapsed << "sec" << endl;        

        // Run simulation
        start = hclock::now();

        std::cout << "run until 3000:00:00" << std::endl;
        r.run_until({3000});
        std::cout << "simulation finished" << std::endl;

        elapsed = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(hclock::now() - start).count();
        cout << "Simulation took:" << elapsed << "sec" << endl;
    
    #endif

    return 0;
}