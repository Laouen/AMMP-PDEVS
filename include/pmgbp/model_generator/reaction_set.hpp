#ifndef PMGBP_REACTION_SET_HPP
#define PMGBP_REACTION_SET_HPP

#include <string>
#include <vector>
#include <iostream>

#include <cadmium/modeling/dynamic_model_translator.hpp>
#include <cadmium/modeling/dynamic_coupled.hpp>
#include <cadmium/modeling/dynamic_model.hpp>

#include <pmgbp/model_generator/reaction_group.hpp>
#include <pmgbp/structures/types.hpp>
#include <pmgbp/atomics/router.hpp>
#include <pmgbp/atomics/reaction.hpp>

#include <NDTime.hpp>

std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> make_reaction_set(
	std::string cid,
    std::string rsn,
    std::vector<std::vector<std::string>> groups_reaction_ids,
    std::string parameters_xml) 
{

    std::string reaction_set_id = cid + '_' + rsn;
    std::string router_id = "router_" + reaction_set_id;
    std::string group_id;

    cadmium::dynamic::modeling::Models models;

    models.push_back(
        cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::router, NDTime, const char*, const char*>(
            router_id,
            parameters_xml.c_str(),
            reaction_set_id.c_str()
        )
    );

    cadmium::dynamic::modeling::EICs eics = {
        cadmium::dynamic::translate::make_EIC<pmgbp::models::reaction_ports::in_0, pmgbp::models::router_ports::in_0>(router_id)
    };

    cadmium::dynamic::modeling::EOCs eocs;
    cadmium::dynamic::modeling::ICs ics;

    for (int group_number = 0; group_number < groups_reaction_ids.size(); group_number++) {
        
        group_id = cid + '_' + rsn + '_' + std::to_string(group_number);
        models.push_back(make_reaction_group(group_id, groups_reaction_ids[group_number], parameters_xml));

        ics.push_back(make_router_reaction_ic(group_number, router_id, group_id));

        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_0, pmgbp::models::reaction_ports::out_0>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_1, pmgbp::models::reaction_ports::out_1>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_2, pmgbp::models::reaction_ports::out_2>(group_id));
    }

    cadmium::dynamic::modeling::Ports iports = { typeid(pmgbp::models::reaction_ports::in_0) };
    cadmium::dynamic::modeling::Ports oports = { 
        typeid(pmgbp::models::reaction_ports::out_0),
        typeid(pmgbp::models::reaction_ports::out_1),
        typeid(pmgbp::models::reaction_ports::out_2)
    };

    return std::make_shared<cadmium::dynamic::modeling::coupled<NDTime>>(
        reaction_set_id,
        models,
        iports,
        oports,
        eics,
        eocs,
        ics
    );
}

#endif //PMGBP_REACTION_SET_HPP
