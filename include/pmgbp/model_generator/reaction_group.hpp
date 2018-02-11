#ifndef PMGBP_REACTION_GROUP_HPP
#define PMGBP_REACTION_GROUP_HPP

#include <string>
#include <vector>

#include <cadmium/modeling/dynamic_coupled.hpp>
#include <cadmium/modeling/dynamic_model.hpp>
#include <cadmium/modeling/dynamic_model_translator.hpp>

cadmium::dynamic::modeling::IC make_router_reaction_ic(int, std::string&, std::string&);

std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> make_reaction_group(std::string, std::vector<std::string>, std::string);

#endif //PMGBP_REACTION_GROUP_HPP
