#ifndef BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H
#define BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H

#include <string>
#include <iostrim>
#include <utility>
#include <map>

#include "parser/parser.hpp"
#include "data-structures/types.hpp"

// Hacer template
template<class TIME>
using vectorOfModels_t = vector< shared_ptr< model<TIME> > >;
template<class TIME>
using vectorOfModelPairs_t = vector< pair< shared_ptr< model<TIME> >, shared_ptr< model<TIME> > > >;
template<class TIME>
using modelsMap_t = map<string, shared_ptr< model<TIME> > >;
template<class TIME,class MSG>
using vectorOfCoupledModels_t = vector< shared_ptr< flattened_coupled<TIME, MSG> > >;
template<class TIME,class MSG>
using coupledModelsMap_t = map<string, shared_ptr< flattened_coupled<TIME, MSG> > >;

template<class TIME,class MSG>
class ModelEngine() {

private:
	bool model_ready;
	double cell_weight;
	string e, c, p, b;
	Parser_t doc;

	// addresses
	map<string, Address_t> enzyme_addresses;

	// models ptr
  	shared_ptr<flattened_coupled<TIME, MSG>> 		biomass_model
  	map< string, coupledModelsMap_t<TIME, MSG> > 	enzyme_models;
  	coupledModelsMap_t<TIME, MSG> 					compartment_models;
  	coupledModelsMap_t<TIME, MSG>  					enzyme_set_models;
  	vectorOfCoupledModels_t<TIME, MSG>				organelle_models;
  	shared_ptr<flattened_coupled<TIME, MSG>> 		extra_cellular_model;
  	shared_ptr<flattened_coupled<TIME, MSG>> 		periplasm_model;
  	shared_ptr<flattened_coupled<TIME, MSG>> 		cytoplasm_model;
  	shared_ptr<flattened_coupled<TIME, MSG>> 		cell_model;
  	shared_ptr< model<TIME> >						input;
  	shared_ptr< model<TIME> >						output_request;

public:
	explicit ModelEngine() : model_ready(false), cell_weight(0), e("e"), c("c"), p("p") {}
};

#endif // BOOST_SIMULATION_PDEVS_MODEL_ENGINE_H