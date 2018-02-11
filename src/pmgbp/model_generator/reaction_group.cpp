#include <pmggbp/model_generator/reaction_group.hpp>

cadmium::dynamic::modeling::IC make_router_reaction_ic(int index, std::string& router_id, std::string& reaction_id) {
    switch(index) {
    case 0: return cadmium::dynamic::translate::make_IC<router_ports::out_0,reaction_ports::in_0>(router_id, reaction_id); break;
    case 1: return cadmium::dynamic::translate::make_IC<router_ports::out_1,reaction_ports::in_0>(router_id, reaction_id); break;
    case 2: return cadmium::dynamic::translate::make_IC<router_ports::out_2,reaction_ports::in_0>(router_id, reaction_id); break;
    case 3: return cadmium::dynamic::translate::make_IC<router_ports::out_3,reaction_ports::in_0>(router_id, reaction_id); break;
    case 4: return cadmium::dynamic::translate::make_IC<router_ports::out_4,reaction_ports::in_0>(router_id, reaction_id); break;
    case 5: return cadmium::dynamic::translate::make_IC<router_ports::out_5,reaction_ports::in_0>(router_id, reaction_id); break;
    case 6: return cadmium::dynamic::translate::make_IC<router_ports::out_6,reaction_ports::in_0>(router_id, reaction_id); break;
    case 7: return cadmium::dynamic::translate::make_IC<router_ports::out_7,reaction_ports::in_0>(router_id, reaction_id); break;
    case 8: return cadmium::dynamic::translate::make_IC<router_ports::out_8,reaction_ports::in_0>(router_id, reaction_id); break;
    case 9: return cadmium::dynamic::translate::make_IC<router_ports::out_9,reaction_ports::in_0>(router_id, reaction_id); break;
    case 10: return cadmium::dynamic::translate::make_IC<router_ports::out_10,reaction_ports::in_0>(router_id, reaction_id); break;
    case 11: return cadmium::dynamic::translate::make_IC<router_ports::out_11,reaction_ports::in_0>(router_id, reaction_id); break;
    case 12: return cadmium::dynamic::translate::make_IC<router_ports::out_12,reaction_ports::in_0>(router_id, reaction_id); break;
    case 13: return cadmium::dynamic::translate::make_IC<router_ports::out_13,reaction_ports::in_0>(router_id, reaction_id); break;
    case 14: return cadmium::dynamic::translate::make_IC<router_ports::out_14,reaction_ports::in_0>(router_id, reaction_id); break;
    case 15: return cadmium::dynamic::translate::make_IC<router_ports::out_15,reaction_ports::in_0>(router_id, reaction_id); break;
    case 16: return cadmium::dynamic::translate::make_IC<router_ports::out_16,reaction_ports::in_0>(router_id, reaction_id); break;
    case 17: return cadmium::dynamic::translate::make_IC<router_ports::out_17,reaction_ports::in_0>(router_id, reaction_id); break;
    case 18: return cadmium::dynamic::translate::make_IC<router_ports::out_18,reaction_ports::in_0>(router_id, reaction_id); break;
    case 19: return cadmium::dynamic::translate::make_IC<router_ports::out_19,reaction_ports::in_0>(router_id, reaction_id); break;
    case 20: return cadmium::dynamic::translate::make_IC<router_ports::out_20,reaction_ports::in_0>(router_id, reaction_id); break;
    case 21: return cadmium::dynamic::translate::make_IC<router_ports::out_21,reaction_ports::in_0>(router_id, reaction_id); break;
    case 22: return cadmium::dynamic::translate::make_IC<router_ports::out_22,reaction_ports::in_0>(router_id, reaction_id); break;
    case 23: return cadmium::dynamic::translate::make_IC<router_ports::out_23,reaction_ports::in_0>(router_id, reaction_id); break;
    case 24: return cadmium::dynamic::translate::make_IC<router_ports::out_24,reaction_ports::in_0>(router_id, reaction_id); break;
    case 25: return cadmium::dynamic::translate::make_IC<router_ports::out_25,reaction_ports::in_0>(router_id, reaction_id); break;
    case 26: return cadmium::dynamic::translate::make_IC<router_ports::out_26,reaction_ports::in_0>(router_id, reaction_id); break;
    case 27: return cadmium::dynamic::translate::make_IC<router_ports::out_27,reaction_ports::in_0>(router_id, reaction_id); break;
    case 28: return cadmium::dynamic::translate::make_IC<router_ports::out_28,reaction_ports::in_0>(router_id, reaction_id); break;
    case 29: return cadmium::dynamic::translate::make_IC<router_ports::out_29,reaction_ports::in_0>(router_id, reaction_id); break;
    case 30: return cadmium::dynamic::translate::make_IC<router_ports::out_30,reaction_ports::in_0>(router_id, reaction_id); break;
    case 31: return cadmium::dynamic::translate::make_IC<router_ports::out_31,reaction_ports::in_0>(router_id, reaction_id); break;
    case 32: return cadmium::dynamic::translate::make_IC<router_ports::out_32,reaction_ports::in_0>(router_id, reaction_id); break;
    case 33: return cadmium::dynamic::translate::make_IC<router_ports::out_33,reaction_ports::in_0>(router_id, reaction_id); break;
    case 34: return cadmium::dynamic::translate::make_IC<router_ports::out_34,reaction_ports::in_0>(router_id, reaction_id); break;
    case 35: return cadmium::dynamic::translate::make_IC<router_ports::out_35,reaction_ports::in_0>(router_id, reaction_id); break;
    case 36: return cadmium::dynamic::translate::make_IC<router_ports::out_36,reaction_ports::in_0>(router_id, reaction_id); break;
    case 37: return cadmium::dynamic::translate::make_IC<router_ports::out_37,reaction_ports::in_0>(router_id, reaction_id); break;
    case 38: return cadmium::dynamic::translate::make_IC<router_ports::out_38,reaction_ports::in_0>(router_id, reaction_id); break;
    case 39: return cadmium::dynamic::translate::make_IC<router_ports::out_39,reaction_ports::in_0>(router_id, reaction_id); break;
    case 40: return cadmium::dynamic::translate::make_IC<router_ports::out_40,reaction_ports::in_0>(router_id, reaction_id); break;
    case 41: return cadmium::dynamic::translate::make_IC<router_ports::out_41,reaction_ports::in_0>(router_id, reaction_id); break;
    case 42: return cadmium::dynamic::translate::make_IC<router_ports::out_42,reaction_ports::in_0>(router_id, reaction_id); break;
    case 43: return cadmium::dynamic::translate::make_IC<router_ports::out_43,reaction_ports::in_0>(router_id, reaction_id); break;
    case 44: return cadmium::dynamic::translate::make_IC<router_ports::out_44,reaction_ports::in_0>(router_id, reaction_id); break;
    case 45: return cadmium::dynamic::translate::make_IC<router_ports::out_45,reaction_ports::in_0>(router_id, reaction_id); break;
    case 46: return cadmium::dynamic::translate::make_IC<router_ports::out_46,reaction_ports::in_0>(router_id, reaction_id); break;
    case 47: return cadmium::dynamic::translate::make_IC<router_ports::out_47,reaction_ports::in_0>(router_id, reaction_id); break;
    case 48: return cadmium::dynamic::translate::make_IC<router_ports::out_48,reaction_ports::in_0>(router_id, reaction_id); break;
    case 49: return cadmium::dynamic::translate::make_IC<router_ports::out_49,reaction_ports::in_0>(router_id, reaction_id); break;
    case 50: return cadmium::dynamic::translate::make_IC<router_ports::out_50,reaction_ports::in_0>(router_id, reaction_id); break;
    case 51: return cadmium::dynamic::translate::make_IC<router_ports::out_51,reaction_ports::in_0>(router_id, reaction_id); break;
    case 52: return cadmium::dynamic::translate::make_IC<router_ports::out_52,reaction_ports::in_0>(router_id, reaction_id); break;
    case 53: return cadmium::dynamic::translate::make_IC<router_ports::out_53,reaction_ports::in_0>(router_id, reaction_id); break;
    case 54: return cadmium::dynamic::translate::make_IC<router_ports::out_54,reaction_ports::in_0>(router_id, reaction_id); break;
    case 55: return cadmium::dynamic::translate::make_IC<router_ports::out_55,reaction_ports::in_0>(router_id, reaction_id); break;
    case 56: return cadmium::dynamic::translate::make_IC<router_ports::out_56,reaction_ports::in_0>(router_id, reaction_id); break;
    case 57: return cadmium::dynamic::translate::make_IC<router_ports::out_57,reaction_ports::in_0>(router_id, reaction_id); break;
    case 58: return cadmium::dynamic::translate::make_IC<router_ports::out_58,reaction_ports::in_0>(router_id, reaction_id); break;
    case 59: return cadmium::dynamic::translate::make_IC<router_ports::out_59,reaction_ports::in_0>(router_id, reaction_id); break;
    case 60: return cadmium::dynamic::translate::make_IC<router_ports::out_60,reaction_ports::in_0>(router_id, reaction_id); break;
    case 61: return cadmium::dynamic::translate::make_IC<router_ports::out_61,reaction_ports::in_0>(router_id, reaction_id); break;
    case 62: return cadmium::dynamic::translate::make_IC<router_ports::out_62,reaction_ports::in_0>(router_id, reaction_id); break;
    case 63: return cadmium::dynamic::translate::make_IC<router_ports::out_63,reaction_ports::in_0>(router_id, reaction_id); break;
    case 64: return cadmium::dynamic::translate::make_IC<router_ports::out_64,reaction_ports::in_0>(router_id, reaction_id); break;
    case 65: return cadmium::dynamic::translate::make_IC<router_ports::out_65,reaction_ports::in_0>(router_id, reaction_id); break;
    case 66: return cadmium::dynamic::translate::make_IC<router_ports::out_66,reaction_ports::in_0>(router_id, reaction_id); break;
    case 67: return cadmium::dynamic::translate::make_IC<router_ports::out_67,reaction_ports::in_0>(router_id, reaction_id); break;
    case 68: return cadmium::dynamic::translate::make_IC<router_ports::out_68,reaction_ports::in_0>(router_id, reaction_id); break;
    case 69: return cadmium::dynamic::translate::make_IC<router_ports::out_69,reaction_ports::in_0>(router_id, reaction_id); break;
    case 70: return cadmium::dynamic::translate::make_IC<router_ports::out_70,reaction_ports::in_0>(router_id, reaction_id); break;
    case 71: return cadmium::dynamic::translate::make_IC<router_ports::out_71,reaction_ports::in_0>(router_id, reaction_id); break;
    case 72: return cadmium::dynamic::translate::make_IC<router_ports::out_72,reaction_ports::in_0>(router_id, reaction_id); break;
    case 73: return cadmium::dynamic::translate::make_IC<router_ports::out_73,reaction_ports::in_0>(router_id, reaction_id); break;
    case 74: return cadmium::dynamic::translate::make_IC<router_ports::out_74,reaction_ports::in_0>(router_id, reaction_id); break;
    case 75: return cadmium::dynamic::translate::make_IC<router_ports::out_75,reaction_ports::in_0>(router_id, reaction_id); break;
    case 76: return cadmium::dynamic::translate::make_IC<router_ports::out_76,reaction_ports::in_0>(router_id, reaction_id); break;
    case 77: return cadmium::dynamic::translate::make_IC<router_ports::out_77,reaction_ports::in_0>(router_id, reaction_id); break;
    case 78: return cadmium::dynamic::translate::make_IC<router_ports::out_78,reaction_ports::in_0>(router_id, reaction_id); break;
    case 79: return cadmium::dynamic::translate::make_IC<router_ports::out_79,reaction_ports::in_0>(router_id, reaction_id); break;
    case 80: return cadmium::dynamic::translate::make_IC<router_ports::out_80,reaction_ports::in_0>(router_id, reaction_id); break;
    case 81: return cadmium::dynamic::translate::make_IC<router_ports::out_81,reaction_ports::in_0>(router_id, reaction_id); break;
    case 82: return cadmium::dynamic::translate::make_IC<router_ports::out_82,reaction_ports::in_0>(router_id, reaction_id); break;
    case 83: return cadmium::dynamic::translate::make_IC<router_ports::out_83,reaction_ports::in_0>(router_id, reaction_id); break;
    case 84: return cadmium::dynamic::translate::make_IC<router_ports::out_84,reaction_ports::in_0>(router_id, reaction_id); break;
    case 85: return cadmium::dynamic::translate::make_IC<router_ports::out_85,reaction_ports::in_0>(router_id, reaction_id); break;
    case 86: return cadmium::dynamic::translate::make_IC<router_ports::out_86,reaction_ports::in_0>(router_id, reaction_id); break;
    case 87: return cadmium::dynamic::translate::make_IC<router_ports::out_87,reaction_ports::in_0>(router_id, reaction_id); break;
    case 88: return cadmium::dynamic::translate::make_IC<router_ports::out_88,reaction_ports::in_0>(router_id, reaction_id); break;
    case 89: return cadmium::dynamic::translate::make_IC<router_ports::out_89,reaction_ports::in_0>(router_id, reaction_id); break;
    case 90: return cadmium::dynamic::translate::make_IC<router_ports::out_90,reaction_ports::in_0>(router_id, reaction_id); break;
    case 91: return cadmium::dynamic::translate::make_IC<router_ports::out_91,reaction_ports::in_0>(router_id, reaction_id); break;
    case 92: return cadmium::dynamic::translate::make_IC<router_ports::out_92,reaction_ports::in_0>(router_id, reaction_id); break;
    case 93: return cadmium::dynamic::translate::make_IC<router_ports::out_93,reaction_ports::in_0>(router_id, reaction_id); break;
    case 94: return cadmium::dynamic::translate::make_IC<router_ports::out_94,reaction_ports::in_0>(router_id, reaction_id); break;
    case 95: return cadmium::dynamic::translate::make_IC<router_ports::out_95,reaction_ports::in_0>(router_id, reaction_id); break;
    case 96: return cadmium::dynamic::translate::make_IC<router_ports::out_96,reaction_ports::in_0>(router_id, reaction_id); break;
    case 97: return cadmium::dynamic::translate::make_IC<router_ports::out_97,reaction_ports::in_0>(router_id, reaction_id); break;
    case 98: return cadmium::dynamic::translate::make_IC<router_ports::out_98,reaction_ports::in_0>(router_id, reaction_id); break;
    case 99: return cadmium::dynamic::translate::make_IC<router_ports::out_99,reaction_ports::in_0>(router_id, reaction_id); break;
    case 100: return cadmium::dynamic::translate::make_IC<router_ports::out_100,reaction_ports::in_0>(router_id, reaction_id); break;
    case 101: return cadmium::dynamic::translate::make_IC<router_ports::out_101,reaction_ports::in_0>(router_id, reaction_id); break;
    case 102: return cadmium::dynamic::translate::make_IC<router_ports::out_102,reaction_ports::in_0>(router_id, reaction_id); break;
    case 103: return cadmium::dynamic::translate::make_IC<router_ports::out_103,reaction_ports::in_0>(router_id, reaction_id); break;
    case 104: return cadmium::dynamic::translate::make_IC<router_ports::out_104,reaction_ports::in_0>(router_id, reaction_id); break;
    case 105: return cadmium::dynamic::translate::make_IC<router_ports::out_105,reaction_ports::in_0>(router_id, reaction_id); break;
    case 106: return cadmium::dynamic::translate::make_IC<router_ports::out_106,reaction_ports::in_0>(router_id, reaction_id); break;
    case 107: return cadmium::dynamic::translate::make_IC<router_ports::out_107,reaction_ports::in_0>(router_id, reaction_id); break;
    case 108: return cadmium::dynamic::translate::make_IC<router_ports::out_108,reaction_ports::in_0>(router_id, reaction_id); break;
    case 109: return cadmium::dynamic::translate::make_IC<router_ports::out_109,reaction_ports::in_0>(router_id, reaction_id); break;
    case 110: return cadmium::dynamic::translate::make_IC<router_ports::out_110,reaction_ports::in_0>(router_id, reaction_id); break;
    case 111: return cadmium::dynamic::translate::make_IC<router_ports::out_111,reaction_ports::in_0>(router_id, reaction_id); break;
    case 112: return cadmium::dynamic::translate::make_IC<router_ports::out_112,reaction_ports::in_0>(router_id, reaction_id); break;
    case 113: return cadmium::dynamic::translate::make_IC<router_ports::out_113,reaction_ports::in_0>(router_id, reaction_id); break;
    case 114: return cadmium::dynamic::translate::make_IC<router_ports::out_114,reaction_ports::in_0>(router_id, reaction_id); break;
    case 115: return cadmium::dynamic::translate::make_IC<router_ports::out_115,reaction_ports::in_0>(router_id, reaction_id); break;
    case 116: return cadmium::dynamic::translate::make_IC<router_ports::out_116,reaction_ports::in_0>(router_id, reaction_id); break;
    case 117: return cadmium::dynamic::translate::make_IC<router_ports::out_117,reaction_ports::in_0>(router_id, reaction_id); break;
    case 118: return cadmium::dynamic::translate::make_IC<router_ports::out_118,reaction_ports::in_0>(router_id, reaction_id); break;
    case 119: return cadmium::dynamic::translate::make_IC<router_ports::out_119,reaction_ports::in_0>(router_id, reaction_id); break;
    case 120: return cadmium::dynamic::translate::make_IC<router_ports::out_120,reaction_ports::in_0>(router_id, reaction_id); break;
    case 121: return cadmium::dynamic::translate::make_IC<router_ports::out_121,reaction_ports::in_0>(router_id, reaction_id); break;
    case 122: return cadmium::dynamic::translate::make_IC<router_ports::out_122,reaction_ports::in_0>(router_id, reaction_id); break;
    case 123: return cadmium::dynamic::translate::make_IC<router_ports::out_123,reaction_ports::in_0>(router_id, reaction_id); break;
    case 124: return cadmium::dynamic::translate::make_IC<router_ports::out_124,reaction_ports::in_0>(router_id, reaction_id); break;
    case 125: return cadmium::dynamic::translate::make_IC<router_ports::out_125,reaction_ports::in_0>(router_id, reaction_id); break;
    case 126: return cadmium::dynamic::translate::make_IC<router_ports::out_126,reaction_ports::in_0>(router_id, reaction_id); break;
    case 127: return cadmium::dynamic::translate::make_IC<router_ports::out_127,reaction_ports::in_0>(router_id, reaction_id); break;
    case 128: return cadmium::dynamic::translate::make_IC<router_ports::out_128,reaction_ports::in_0>(router_id, reaction_id); break;
    case 129: return cadmium::dynamic::translate::make_IC<router_ports::out_129,reaction_ports::in_0>(router_id, reaction_id); break;
    case 130: return cadmium::dynamic::translate::make_IC<router_ports::out_130,reaction_ports::in_0>(router_id, reaction_id); break;
    case 131: return cadmium::dynamic::translate::make_IC<router_ports::out_131,reaction_ports::in_0>(router_id, reaction_id); break;
    case 132: return cadmium::dynamic::translate::make_IC<router_ports::out_132,reaction_ports::in_0>(router_id, reaction_id); break;
    case 133: return cadmium::dynamic::translate::make_IC<router_ports::out_133,reaction_ports::in_0>(router_id, reaction_id); break;
    case 134: return cadmium::dynamic::translate::make_IC<router_ports::out_134,reaction_ports::in_0>(router_id, reaction_id); break;
    case 135: return cadmium::dynamic::translate::make_IC<router_ports::out_135,reaction_ports::in_0>(router_id, reaction_id); break;
    case 136: return cadmium::dynamic::translate::make_IC<router_ports::out_136,reaction_ports::in_0>(router_id, reaction_id); break;
    case 137: return cadmium::dynamic::translate::make_IC<router_ports::out_137,reaction_ports::in_0>(router_id, reaction_id); break;
    case 138: return cadmium::dynamic::translate::make_IC<router_ports::out_138,reaction_ports::in_0>(router_id, reaction_id); break;
    case 139: return cadmium::dynamic::translate::make_IC<router_ports::out_139,reaction_ports::in_0>(router_id, reaction_id); break;
    case 140: return cadmium::dynamic::translate::make_IC<router_ports::out_140,reaction_ports::in_0>(router_id, reaction_id); break;
    case 141: return cadmium::dynamic::translate::make_IC<router_ports::out_141,reaction_ports::in_0>(router_id, reaction_id); break;
    case 142: return cadmium::dynamic::translate::make_IC<router_ports::out_142,reaction_ports::in_0>(router_id, reaction_id); break;
    case 143: return cadmium::dynamic::translate::make_IC<router_ports::out_143,reaction_ports::in_0>(router_id, reaction_id); break;
    case 144: return cadmium::dynamic::translate::make_IC<router_ports::out_144,reaction_ports::in_0>(router_id, reaction_id); break;
    case 145: return cadmium::dynamic::translate::make_IC<router_ports::out_145,reaction_ports::in_0>(router_id, reaction_id); break;
    case 146: return cadmium::dynamic::translate::make_IC<router_ports::out_146,reaction_ports::in_0>(router_id, reaction_id); break;
    case 147: return cadmium::dynamic::translate::make_IC<router_ports::out_147,reaction_ports::in_0>(router_id, reaction_id); break;
    case 148: return cadmium::dynamic::translate::make_IC<router_ports::out_148,reaction_ports::in_0>(router_id, reaction_id); break;
    case 149: return cadmium::dynamic::translate::make_IC<router_ports::out_149,reaction_ports::in_0>(router_id, reaction_id); break;
    default: assert(false && "invalid router out port number");
    }
}

std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> make_reaction_group(
    std::string group_id,
    std::vector<std::string> model_ids,
    std::string parameters_xml) {

    std::string router_id = "router_" + group_id;
    std::string reaction_id;

    cadmium::dynamic::modeling::Models sub_models = {
        cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::router, NDTime, std::string, std::string>(
            router_id,
            parameters_xml,
            router_id
        )
    };

    cadmium::dynamic::modeling::EICs eics = {
        cadmium::dynamic::translate::make_EIC<reaction_group_ports::in_0, router_ports::in_0>(router_id;)
    };

    cadmium::dynamic::modeling::EOCs eocs;
    cadmium::dynamic::modeling::ICs ics;

    int ic_index = 0;
    for (auto& model_id : model_ids) {
        reaction_id = "reaction_" + model_id; 
        
        sub_models.push_back(
            cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::reaction, NDTime, std::string, std::string>(
                reaction_id,
                parameters_xml,
                reaction_id
            )
        );

        ics.push_back(make_router_reaction_ic(ic_index, router_id, reaction_id));
        ic_index++;

        eocs.push_back(cadmium::dynamic::translate::make_EOC<reaction_ports::out_0, reaction_group_ports::out_0>(reaction_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<reaction_ports::out_1, reaction_group_ports::out_1>(reaction_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<reaction_ports::out_2, reaction_group_ports::out_2>(reaction_id));

    }

    return std::make_shared<cadmium::dynamic::modeling::coupled<NDTime>>(
        group_id,
        sub_models,
        reaction_group_iports,
        reaction_group_oports,
        eics,
        eocs,
        ics
    );
}