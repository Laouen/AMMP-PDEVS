/****************************** atomic model {{model_name}} *****************************************/

std::shared_ptr<cadmium::dynamic::modeling::model> {{model_name}} =
	cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::{{model_class}}, {TIME}, {{ARGS}}>(
		"{{model_name}}",
		{{parameters}}
	);
/**************************************************************************************************/