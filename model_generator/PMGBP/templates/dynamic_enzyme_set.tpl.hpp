/*************************** coupled model {{cid}} {{esn}} *******************************************/

std::shared_ptr<cadmium::dynamic::modeling::coupled<{TIME}>> {{cid}}_{{esn}} = make_enzyme_set(
	"{{cid}}",
	"{{esn}}",
	{{enzyme_ids}},
	xml_parameter_path
);

/**************************************************************************************************/