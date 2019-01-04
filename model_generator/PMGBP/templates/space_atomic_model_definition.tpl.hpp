/***************************** ports for model {model_name} ***************************************/

namespace pmgbp {{
namespace structs {{
namespace {model_name} {{

struct {model_name}_ports {{

    {input_ports_definitions}

    {output_ports_definitions}

    using reactant_type=pmgbp::types::{reactant_type};
    using product_type=pmgbp::types::{product_type};
    using information_type=pmgbp::types::{information_type};

    using input_ports=std::tuple<{input_port_names}>;
    using output_ports=std::tuple<{output_port_names}>;
}};


}}
}}
}}

template<typename TIME>
using {model_name}_definition = pmgbp::models::space<pmgbp::structs::{model_name}::{model_name}_ports, TIME>;

/**************************************************************************************************/