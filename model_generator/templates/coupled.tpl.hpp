/*************************** coupled model {model_name} *******************************************/
{ports}

using oports_{model_name} = std::tuple<{oports}>;
using iports_{model_name} = std::tuple<{iports}>;
using submodels_{model_name} = cadmium::modeling::models_tuple<{submodels}>;
using eics_{model_name} = std::tuple<{eic}>;
using eocs_{model_name} = std::tuple<{eoc}>;
using ics_{model_name} = std::tuple<{ic}>;

template<typename TIME>
struct reaction_set: public cadmium::modeling::coupled_model<
        TIME,
        oports_{model_name},
        iports_{model_name},
        submodels_{model_name},
        eics_{model_name},
        eocs_{model_name},
        ics_{model_name}
> {{}};
/**************************************************************************************************/