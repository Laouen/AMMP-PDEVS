/*************************** coupled model {{model_name}} *******************************************/
struct {{model_name}}_ports {{{{
    {{ports}}
}}}};

cadmium::dynamic::modeling::Models sub_models_{{model_name}} = {{{{ {{sub_models}} }}}};
cadmium::dynamic::modeling::Ports iports_{{model_name}} = {{{{ {{iports}} }}}};
cadmium::dynamic::modeling::Ports oports_{{model_name}} = {{{{ {{oports}} }}}};
cadmium::dynamic::modeling::EICs eics_{{model_name}} = {{{{ {{eic}} }}}};
cadmium::dynamic::modeling::EOCs eocs_{{model_name}} = {{{{ {{eoc}} }}}};
cadmium::dynamic::modeling::ICs ics_{{model_name}} = {{{{ {{ic}} }}}};

std::shared_ptr<cadmium::dynamic::modeling::coupled<{TIME}>> {{model_name}} = std::make_shared<cadmium::dynamic::modeling::coupled<{TIME}>>(
    "{{model_name}}",
    sub_models_{{model_name}},
    iports_{{model_name}},
    oports_{{model_name}},
    eics_{{model_name}},
    eocs_{{model_name}},
    ics_{{model_name}}
);
/**************************************************************************************************/