#ifndef PMGBP_REACTION_SET_HPP
#define PMGBP_REACTION_SET_HPP

#include <string>
#include <vector>
#include <iostream>

#include <cadmium/modeling/dynamic_model_translator.hpp>
#include <cadmium/modeling/dynamic_coupled.hpp>
#include <cadmium/modeling/dynamic_model.hpp>

#include <pmgbp/model_generator/enzyme_group.hpp>
#include <pmgbp/structures/types.hpp>
#include <pmgbp/structures/space.hpp>
#include <pmgbp/atomics/router.hpp>
#include <pmgbp/atomics/enzyme.hpp>

#include <NDTime.hpp>

std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> make_enzyme_set(
	std::string cid,
    std::string esn,
    std::vector<std::vector<std::string>> groups_enzyme_ids,
    std::string parameters_xml) 
{

    std::string enzyme_set_id = cid + '_' + esn;
    std::string router_id = "router_" + enzyme_set_id;
    std::string group_id;

    pmgbp::structs::space::EnzymeAddress enzymes_location(cid, esn);

    cadmium::dynamic::modeling::Models models;

    models.push_back(
        cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::router, NDTime, const char*, const char*>(
            router_id,
            parameters_xml.c_str(),
            enzyme_set_id.c_str()
        )
    );

    cadmium::dynamic::modeling::EICs eics = {
        cadmium::dynamic::translate::make_EIC<pmgbp::models::enzyme_ports::in_0, pmgbp::models::router_ports::in_0>(router_id)
    };

    cadmium::dynamic::modeling::EOCs eocs;
    cadmium::dynamic::modeling::ICs ics;

    for (int group_number = 0; group_number < groups_enzyme_ids.size(); group_number++) {
        
        group_id = cid + '_' + esn + '_' + std::to_string(group_number);
        models.push_back(make_enzyme_group(group_id, enzymes_location, groups_enzyme_ids[group_number], parameters_xml));

        ics.push_back(make_router_enzyme_ic(group_number, router_id, group_id));

        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_0_product, pmgbp::models::enzyme_ports::out_0_product>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_1_product, pmgbp::models::enzyme_ports::out_1_product>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_2_product, pmgbp::models::enzyme_ports::out_2_product>(group_id));

        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_0_information, pmgbp::models::enzyme_ports::out_0_information>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_1_information, pmgbp::models::enzyme_ports::out_1_information>(group_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::enzyme_ports::out_2_information, pmgbp::models::enzyme_ports::out_2_information>(group_id));
    }

    cadmium::dynamic::modeling::Ports iports = { typeid(pmgbp::models::enzyme_ports::in_0) };
    cadmium::dynamic::modeling::Ports oports = { 
        typeid(pmgbp::models::enzyme_ports::out_0_product),
        typeid(pmgbp::models::enzyme_ports::out_1_product),
        typeid(pmgbp::models::enzyme_ports::out_2_product),
        typeid(pmgbp::models::enzyme_ports::out_0_information),
        typeid(pmgbp::models::enzyme_ports::out_1_information),
        typeid(pmgbp::models::enzyme_ports::out_2_information)
    };

    return std::make_shared<cadmium::dynamic::modeling::coupled<NDTime>>(
        enzyme_set_id,
        models,
        iports,
        oports,
        eics,
        eocs,
        ics
    );
}

#endif //PMGBP_REACTION_SET_HPP
