/****************************** atomic model {model_name} *****************************************/
using {model_name}_ports = pmgbp::structs::{model_name}::ports<{output_type},{input_type}>;

template<typename TIME>
using {model_name}_template = pmgbp::models::{model_class}<{model_name}_ports, TIME>;

template<typename TIME>
class {model_name} : public {model_name}_template<TIME> {{
public:
    {model_name}(): {model_name}_template<TIME>(
        {parameters}
    ) {{}}
}};
/**************************************************************************************************/