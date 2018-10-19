/*************************** coupled model {{cid}} {{rsn}} *******************************************/

std::shared_ptr<cadmium::dynamic::modeling::coupled<{TIME}>> {{cid}}_{{rsn}} = make_reaction_set(
	"{{cid}}",
	"{{rsn}}",
	{{reaction_ids}},
	xml_parameter_path
);

/**************************************************************************************************/