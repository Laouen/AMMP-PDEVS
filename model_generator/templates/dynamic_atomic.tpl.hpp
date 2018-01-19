/****************************** atomic model {{model_name}} *****************************************/
using {{model_name}}_ports = pmgbp::structs::{{model_name}}::ports<pmgbp::types::{{output_type}},pmgbp::types::{{input_type}}>;

template<typename T>
using {{model_name}}_template = pmgbp::models::{{model_class}}<{{model_name}}_ports, T>;

std::shared_ptr<cadmium::dynamic::modeling::model> {{model_name}} = 
	cadmium::dynamic::translate::make_dynamic_atomic_model<{{model_name}}_template, {TIME}, {{ARGS}}>(
		"{{model_name}}",
		{{parameters}}
	);
/**************************************************************************************************/