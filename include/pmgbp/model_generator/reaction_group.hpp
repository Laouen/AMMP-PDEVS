#ifndef PMGBP_REACTION_GROUP_HPP
#define PMGBP_REACTION_GROUP_HPP

#include <string>
#include <vector>

#include <cadmium/modeling/dynamic_model_translator.hpp>
#include <cadmium/modeling/dynamic_coupled.hpp>
#include <cadmium/modeling/dynamic_model.hpp>

#include <pmgbp/structures/types.hpp>
#include <pmgbp/atomics/router.hpp>
#include <pmgbp/atomics/reaction.hpp>

#include <NDTime.hpp>

cadmium::dynamic::modeling::IC make_router_reaction_ic(int index, std::string& router_id, std::string& reaction_id) {
    switch(index) {
    case 0: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_0,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 1: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_1,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 2: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_2,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 3: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_3,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 4: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_4,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 5: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_5,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 6: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_6,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 7: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_7,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 8: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_8,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 9: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_9,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 10: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_10,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 11: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_11,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 12: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_12,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 13: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_13,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 14: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_14,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 15: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_15,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 16: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_16,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 17: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_17,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 18: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_18,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 19: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_19,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 20: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_20,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 21: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_21,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 22: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_22,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 23: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_23,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 24: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_24,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 25: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_25,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 26: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_26,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 27: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_27,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 28: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_28,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 29: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_29,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 30: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_30,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 31: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_31,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 32: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_32,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 33: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_33,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 34: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_34,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 35: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_35,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 36: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_36,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 37: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_37,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 38: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_38,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 39: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_39,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 40: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_40,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 41: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_41,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 42: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_42,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 43: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_43,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 44: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_44,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 45: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_45,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 46: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_46,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 47: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_47,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 48: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_48,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 49: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_49,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 50: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_50,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 51: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_51,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 52: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_52,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 53: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_53,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 54: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_54,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 55: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_55,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 56: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_56,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 57: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_57,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 58: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_58,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 59: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_59,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 60: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_60,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 61: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_61,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 62: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_62,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 63: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_63,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 64: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_64,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 65: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_65,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 66: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_66,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 67: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_67,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 68: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_68,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 69: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_69,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 70: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_70,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 71: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_71,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 72: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_72,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 73: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_73,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 74: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_74,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 75: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_75,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 76: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_76,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 77: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_77,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 78: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_78,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 79: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_79,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 80: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_80,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 81: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_81,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 82: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_82,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 83: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_83,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 84: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_84,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 85: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_85,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 86: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_86,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 87: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_87,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 88: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_88,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 89: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_89,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 90: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_90,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 91: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_91,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 92: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_92,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 93: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_93,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 94: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_94,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 95: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_95,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 96: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_96,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 97: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_97,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 98: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_98,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 99: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_99,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 100: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_100,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 101: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_101,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 102: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_102,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 103: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_103,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 104: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_104,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 105: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_105,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 106: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_106,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 107: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_107,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 108: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_108,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 109: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_109,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 110: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_110,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 111: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_111,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 112: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_112,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 113: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_113,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 114: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_114,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 115: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_115,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 116: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_116,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 117: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_117,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 118: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_118,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 119: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_119,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 120: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_120,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 121: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_121,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 122: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_122,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 123: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_123,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 124: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_124,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 125: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_125,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 126: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_126,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 127: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_127,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 128: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_128,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 129: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_129,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 130: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_130,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 131: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_131,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 132: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_132,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 133: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_133,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 134: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_134,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 135: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_135,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 136: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_136,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 137: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_137,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 138: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_138,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 139: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_139,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 140: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_140,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 141: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_141,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 142: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_142,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 143: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_143,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 144: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_144,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 145: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_145,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 146: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_146,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 147: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_147,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 148: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_148,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    case 149: return cadmium::dynamic::translate::make_IC<pmgbp::models::router_ports::out_149,pmgbp::models::reaction_ports::in_0>(router_id, reaction_id); break;
    default: assert(false && "invalid router out port number");
    }
}

std::shared_ptr<cadmium::dynamic::modeling::coupled<NDTime>> make_reaction_group(
	std::string group_id,
    std::vector<std::string> reaction_ids,
    std::string parameters_xml)
{

    std::string router_id = "router_" + group_id;
    std::string reaction_id;

    cadmium::dynamic::modeling::Models models = {
        cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::router, NDTime, const char*, const char*>(
            router_id,
            parameters_xml.c_str(),
            group_id.c_str()
        )
    };

    cadmium::dynamic::modeling::EICs eics = {
        cadmium::dynamic::translate::make_EIC<pmgbp::models::reaction_ports::in_0, pmgbp::models::router_ports::in_0>(router_id)
    };

    cadmium::dynamic::modeling::EOCs eocs;
    cadmium::dynamic::modeling::ICs ics;

    for (int reaction_index = 0; reaction_index < reaction_ids.size(); reaction_index++) {
        
        reaction_id = reaction_ids[reaction_index]; 
        models.push_back(
            cadmium::dynamic::translate::make_dynamic_atomic_model<pmgbp::models::reaction, NDTime, const char*, const char*>(
                reaction_id,
                parameters_xml.c_str(),
                reaction_id.c_str()
            )
        );

        ics.push_back(make_router_reaction_ic(reaction_index, router_id, reaction_id));

        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_0, pmgbp::models::reaction_ports::out_0>(reaction_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_1, pmgbp::models::reaction_ports::out_1>(reaction_id));
        eocs.push_back(cadmium::dynamic::translate::make_EOC<pmgbp::models::reaction_ports::out_2, pmgbp::models::reaction_ports::out_2>(reaction_id));
    }

    cadmium::dynamic::modeling::Ports iports = { typeid(pmgbp::models::reaction_ports::in_0) };
    cadmium::dynamic::modeling::Ports oports = { 
        typeid(pmgbp::models::reaction_ports::out_0),
        typeid(pmgbp::models::reaction_ports::out_1),
        typeid(pmgbp::models::reaction_ports::out_2)
    };

    return std::make_shared<cadmium::dynamic::modeling::coupled<NDTime>>(
        group_id,
        models,
        iports,
        oports,
        eics,
        eocs,
        ics
    );
}

#endif //PMGBP_REACTION_GROUP_HPP
