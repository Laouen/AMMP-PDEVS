/*************************** coupled model {{group_id}} *******************************************/

std::shared_ptr<cadmium::dynamic::modeling::coupled<{TIME}>> {{group_id}} = make_reaction_group(
	"{{group_id}}",
	{{reaction_ids}},
	"{{parameters_xml}}"
);

/**************************************************************************************************/