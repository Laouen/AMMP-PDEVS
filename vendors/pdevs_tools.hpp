/**
 * Copyright (c) 2013-2015, Damian Vicino
 * Carleton University, Universite de Nice-Sophia Antipolis
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


#ifndef PDEVS_TOOLS_H
#define PDEVS_TOOLS_H
#include <boost/simulation/pdevs/coupled.hpp>
#include <boost/core/demangle.hpp>
#include <typeinfo>
#include <sstream>

namespace bsp = boost::simulation::pdevs;

namespace pdevs_tools {
/**
 * @brief The pdevs_coupling_diagram class provides the tools to generate a plant-uml description of the coupled models
 */
template<class TIME, class MSG>
class pdevs_coupling_diagram
{
    typename boost::simulation::pdevs::coupled<TIME, MSG>::description_type desc;
    int ticket;
public:
    explicit pdevs_coupling_diagram( boost::simulation::pdevs::coupled<TIME, MSG> pc) : desc(pc.get_description()), ticket(0){}

    std::string get_plant_uml() {
        ticket = 0; //reset for repeteability
//print preamble
        std::string ret =
"@startuml\n"
//top model
//"set namespaceSeparator ::\n"

//declaring ports on top model
"class top::in <<(I, green) >>\n"
"class top::out <<(O, orange) >>\n"
;
//check no input connected
assert(desc.external_input_coupling.size() == 0);
//defining internal models
for (auto m : desc.models){
    //detecting model kind
    std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(m);
    std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(m);
    if (m_atomic == nullptr){ //coupled model
        assert(m_coupled != nullptr);
        ret.append(process_puml_coupled(m_coupled));
    } else { //atomic model
        ret.append("class ");
        ret.append(model_name(m_atomic));
    }
    ret.append("\n");
}
//connect model to output
for (auto out : desc.external_output_coupling){
    std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(out);
    std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(out);
    ret.append("\"");
    if (m_atomic == nullptr){ //coupled model
        assert(m_coupled != nullptr);
        ret.append(clean_model_name(model_name(m_coupled)));
    } else { //atomic model
        ret.append(clean_model_name(model_name(m_atomic)));
    }
    ret.append("\"").append(" --> ").append("top::out\n");
}

//connect internal couplings on top level
//connect models internally
for (auto inter : desc.internal_coupling){
    //first item
    {
        std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(inter.first);
        std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(inter.first);
        ret.append("\"");
        if (m_atomic == nullptr){ //coupled model
            assert(m_coupled != nullptr);
            ret.append(clean_model_name(model_name(m_coupled)));
        } else { //atomic model
            ret.append(clean_model_name(model_name(m_atomic)));
        }
        ret.append("\"").append(" --|> ");
    }
    //second item
    {
        std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(inter.second);
        std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(inter.second);
        ret.append("\"");
        if (m_atomic == nullptr){ //coupled model
            assert(m_coupled != nullptr);
            ret.append(clean_model_name(model_name(m_coupled)));
        } else { //atomic model
            ret.append(clean_model_name(model_name(m_atomic)));
        }
        ret.append("\"");
    }
    ret.append("\n");
}

//print epilogue
ret.append(
"hide empty fields\n"
"hide empty methods\n"
"@enduml"
);

        return ret;
    }

private:
    std::string process_puml_coupled(std::shared_ptr<bsp::coupled<TIME, MSG>> m){
        //define ports
        std::string ret;
        ret.append("class ").append(model_name(m));
        ret.append("\n");
        //defining internal models and agregating them
        auto internal_desc = m->get_description();
        for (auto subm : internal_desc.models){
            //detecting model kind
            std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(subm);
            std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(subm);
            if (m_atomic == nullptr){ //coupled model
                assert(m_coupled != nullptr);
                ret.append(process_puml_coupled(m_coupled));
            } else { //atomic model
                ret.append("class ");
                ret.append(model_name(m_atomic));
            }
            ret.append("\n");
        }
        //process EIC, EOC, IC
        //connect model to input
        for (auto in : internal_desc.external_input_coupling){
            std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(in);
            std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(in);
            ret.append("\"");
            if (m_atomic == nullptr){ //coupled model
                assert(m_coupled != nullptr);
                ret.append(clean_model_name(model_name(m_coupled)));
            } else { //atomic model
                ret.append(clean_model_name(model_name(m_atomic)));
            }
            ret.append("\"").append(" <-- ").append(clean_model_name(model_name(m)));
            ret.append("\n");
        }

        //connect model to output
        for (auto out : internal_desc.external_output_coupling){
            std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(out);
            std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(out);
            ret.append("\"");
            if (m_atomic == nullptr){ //coupled model
                assert(m_coupled != nullptr);
                ret.append(clean_model_name(model_name(m_coupled)));
            } else { //atomic model
                ret.append(clean_model_name(model_name(m_atomic)));
            }
            ret.append("\"").append(" --> ").append(clean_model_name(model_name(m)));
            ret.append("\n");
        }

        //connect models internally
        for (auto inter : internal_desc.internal_coupling){
            //first item
            {
                std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(inter.first);
                std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(inter.first);
                ret.append("\"");
                if (m_atomic == nullptr){ //coupled model
                    assert(m_coupled != nullptr);
                    ret.append(clean_model_name(model_name(m_coupled)));
                } else { //atomic model
                    ret.append(clean_model_name(model_name(m_atomic)));
                }
                ret.append("\"").append(" --|> ");
            }
            //second item
            {
                std::shared_ptr<bsp::atomic<TIME, MSG>> m_atomic = std::dynamic_pointer_cast<bsp::atomic<TIME, MSG>>(inter.second);
                std::shared_ptr<bsp::coupled<TIME, MSG>> m_coupled = std::dynamic_pointer_cast<bsp::coupled<TIME, MSG>>(inter.second);
                ret.append("\"");
                if (m_atomic == nullptr){ //coupled model
                    assert(m_coupled != nullptr);
                    ret.append(clean_model_name(model_name(m_coupled)));
                } else { //atomic model
                    ret.append(clean_model_name(model_name(m_atomic)));
                }
                ret.append("\"");
            }
            ret.append("\n");
        }

        return ret;
    }

    std::string model_name(std::shared_ptr<bsp::atomic<TIME, MSG>> a){
        const void * address = static_cast<const void*>(a.get());
        std::stringstream ss;
        ss << "_" << address;
        std::string suffix = ss.str();

        std::string ret = boost::core::demangle(typeid(*a).name());
        std::string needle ("<");
        std::size_t found = ret.find(needle);
        if (found!=std::string::npos) //found a template parameter
            ret.insert(found, suffix);
        else //no template param
            ret.append(suffix);

        ret.append("<<(A, yellow) >>");
        return ret;
    }

    std::string model_name(std::shared_ptr<bsp::coupled<TIME, MSG>> c){
        const void * address = static_cast<const void*>(c.get());
        std::stringstream ss;
        ss << "_" << address;
        std::string suffix = ss.str();

        std::string ret = boost::core::demangle(typeid(*c).name());
        std::string needle ("<");
        std::size_t found = ret.find(needle);
        if (found!=std::string::npos) //found a template parameter
            ret.insert(found, suffix);
        else //no template param
            ret.append(suffix);
        ret.append("<<(c, pink) >>");
        return ret;
    }

    std::string clean_model_name(std::string mn){
        std::string ret(mn);
        std::string needle ("<");
        std::size_t found = ret.find(needle);
        if (found!=std::string::npos) //found a template parameter
            ret.erase(ret.begin()+found, ret.end());
        return ret;
    }
};

}


#endif // PDEVS_TOOLS_H
