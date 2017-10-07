/***************************** ports for model {model_name} ***************************************/
namespace pmgbp {{
namespace structs {{
namespace {model_name} {{

template<class OUTPUT_TYPE, class INPUT_TYPE>
struct ports {{

    {output_ports_definitions}

    {input_ports_definitions}

    using output_type=OUTPUT_TYPE;
    using input_type=INPUT_TYPE;

    using input_ports=std::tuple<{input_port_names}>;
    using output_ports=std::tuple<{output_port_names}>;
}};

}}
}}
}}
/**************************************************************************************************/