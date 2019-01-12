#include <pmgbp/structures/types.hpp>

using namespace std;
using namespace pmgbp::types;

ostream& operator<<(ostream& os, const Address_t& to) {

    os << "[";
    auto i = to.cbegin();
    while(i != to.cend()){
        os << *i;
        ++i;
        if (i != to.cend()) os << ", ";
    }
    os << "]";
    return os;
}

ostream& operator<<(ostream& os, const vector<string>& m) {

    os << "[";
    auto i = m.cbegin();
    while(i != m.cend()){
        os << *i;
        ++i;
        if (i != m.cend()) os << ", ";
    }
    os << "]";
    return os;
}

ostream& operator<<(ostream& os, const MetaboliteAmounts& m) {

    os << "[";
    auto i = m.cbegin();
    while(i != m.cend()){
        os << i->second << "-" << i->first;
        ++i;
        if (i != m.cend()) os << ", ";
    }
    os << "]";
    return os;
}

ostream& operator<<(ostream& os, const BState_t& s) {

    switch(s) {
        case BState_t::IDLE:        os << "IDLE";       break;
        case BState_t::WAITING:     os << "WAITING";    break;
        case BState_t::NOT_ENOUGH:  os << "NOT_ENOUGH"; break;
        case BState_t::ENOUGH:      os << "ENOUGH";     break;
    }

    return os;
}

ostream& operator<<(ostream& os, const Way& s) {

    switch(s) {
        case Way::STP:
            os << "STP";
            break;
        case Way::PTS:
            os << "PTS";
            break;
    }

    return os;
}

ostream& operator<<(ostream& os, const ReactionInfo& r) {

    os << "id: " << r.id << endl;
    os << "substrates: " << r.substrate_sctry << endl;
    os << "products: " << r.products_sctry << endl;
    os << "KonSTP: " << r.konSTP << endl;
    os << "KonPTS" << r.konPTS << endl;
    os << "reversible: " << ((r.reversible) ? "true" : "false");
    return os;
}

ostream& operator<<(ostream& os, const Product& p) {
    os << "{";
    os << "\"Message_type\":\"product\",";
    os << "\"Metabolites\":[";

    bool separate = false;
    for (const auto& metabolite : p.metabolites) {
        if (separate) {
            os << ",";
        }

        os << "{";
        os << "\"Id\":\"" << metabolite.first << "\",";
        os << "\"Amount\":" << metabolite.second;
        os << "}";
    }

    os << "]";
    os << "}";
    return os;
}

ostream& operator<<(ostream& os, const Reactant& r) {
    os << "{";
    os << "\"Message_type\":\"reactant\",";
    os << "\"Reaction_amount\":" << r.reaction_amount;
    os << "\"From\":\"" << r.from << "\",";
    os << "\"For_reaction\":\"" << r.rid << "\",";
    os << "\"Way\":\"" << r.reaction_direction << "\",";
    os << "}";
    return os;
}

ostream& operator<<(ostream& os, const Information& i) {
    os << "{";
    os << "\"Message_type\":\"information\",";
    os << "\"Enzyme ID\":\"" << i.enzyme_id << "\",";
    os << "\"Released amount\":\"" << i.released_enzymes << "\",";
    os << "\"Location\":\"" << i.location << "\",";
    os << "}";
    return os;
}

ostream& operator<<(ostream& os, const Enzyme& e) {

    os << "id: " << e.id << endl;
    os << "location: " << e.location << endl;
    os << "amount: " << e.amount << endl;
    os << "handled reactions: " << endl << endl;

    for (const auto& reaction : e.handled_reactions) {
        os << "------------ Reaction Information -----------" << endl;
        os << reaction.second << endl;
        os << "---------------------------------------------" << endl;
    }

    return os;
}