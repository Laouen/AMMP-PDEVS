#ifndef PMGBP_TYPES_HPP
#define PMGBP_TYPES_HPP

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>

#include "space.hpp"

namespace pmgbp {
namespace types {

using namespace std;

template <class ENTRY>
struct RoutingTable {
    std::map<ENTRY, int> entries;

    int at(const ENTRY& entry) const {

        if (entries.find(entry) != entries.cend()) {
            return entries.at(entry);
        } else {
            return -1;
        }
    }

    void insert(const ENTRY& entry, int port) {
        entries.insert({entry, port});
    }
};

/******************************************/
/********** Enums and renames *************/
/******************************************/

enum class BState_t { ENOUGH = 0, NOT_ENOUGH = 1, IDLE = 2, WAITING = 3 };
enum class Way { STP, PTS };

using Integer = unsigned long long;
using MetaboliteAmounts = map<string, Integer>;

/******************************************/
/******** End enums and renames ***********/
/******************************************/

/******************************************/
/************** Constants *****************/
/******************************************/

const long double L = 6.0221413e+23L;
const long double MOL = 1e-6;

/******************************************/
/************ End Constants ***************/
/******************************************/

/*******************************************/
/***************** Message *****************/
/*******************************************/

using Address_t = list<string>;

struct Product {
    MetaboliteAmounts metabolites;

    bool operator==(const Product& other) const {
        return metabolites.size() == other.metabolites.size() &&
                std::equal(metabolites.begin(), metabolites.end(), other.metabolites.begin());
    }

    void clear() {
        metabolites.clear();
    }
};

struct Reactant {
    string rid;
    string enzyme_id;
    string from;
    Way reaction_direction;
    Integer reaction_amount;

    bool operator==(const Reactant& other) const {
        return rid == other.rid &&
                from == other.from &&
                reaction_amount == other.reaction_amount &&
                reaction_direction == other.reaction_direction;
    }

    void clear() {
        rid = "";
        from = "";
        reaction_amount = 0;
    }
};

struct Information {
    string enzyme_id;
    Integer released_enzymes;
    pmgbp::structs::space::EnzymeAddress location;

    bool operator==(const Information& other) const {
        return enzyme_id == other.enzyme_id && released_enzymes == other.released_enzymes && location == other.location;
    }
};

/*******************************************/
/************** End Messages ***************/
/*******************************************/


/*******************************************/
/************* Data info type **************/
/*******************************************/

struct ReactionInfo {

  string id;
  MetaboliteAmounts  substrate_sctry;
  MetaboliteAmounts  products_sctry;
  double konSTP = 1;
  double konPTS = 1;
  double koffPTS = 1;
  double koffSTP = 1;
  bool reversible = false;

  ReactionInfo() = default;

  ReactionInfo(
    string other_id,
    const MetaboliteAmounts& other_substrate_sctry,
    const MetaboliteAmounts& other_products_sctry,
    double other_konSTP,
    double other_konPTS,
    double other_koffSTP,
    double other_koffPTS,
    bool other_reversible
    ) : id(other_id),
        substrate_sctry(other_substrate_sctry),
        products_sctry(other_products_sctry),
        konSTP(other_konSTP),
        konPTS(other_konPTS),
        koffPTS(other_koffPTS),
        koffSTP(other_koffSTP),
        reversible(other_reversible) {}

  ReactionInfo(const ReactionInfo& other)
  : id(other.id),
    substrate_sctry(other.substrate_sctry),
    products_sctry(other.products_sctry),
    konSTP(other.konSTP),
    konPTS(other.konPTS),
    koffPTS(other.koffPTS),
    koffSTP(other.koffSTP),
    reversible(other.reversible) {}

  void clear() {
    id = "";
    substrate_sctry.clear();
    products_sctry.clear();
    konSTP = 0;
    konPTS = 0;
    reversible = false;
  }

  bool empty() {
    return (id == "") && substrate_sctry.empty() && products_sctry.empty();
  }
};


/*******************************************/
/*********** End Data info type ************/
/*******************************************/

/*******************************************/
/*********** Data enzyme type **************/
/*******************************************/

// TODO: Change name to EnzymeInformation
struct Enzyme {
    string id;
    pmgbp::structs::space::EnzymeAddress location;
    Integer amount;
    map<string, ReactionInfo> handled_reactions;

    Enzyme()
            : id(""), location(), amount(0), handled_reactions() {};

    Enzyme(
            const string& other_id,
            const pmgbp::structs::space::EnzymeAddress& other_location,
            const Integer& other_amount,
            const map<string, ReactionInfo>& other_handled_reactions
    )
            : id(other_id), location(other_location), amount(other_amount), handled_reactions(other_handled_reactions) {}

    Enzyme(const Enzyme& other)
            : id(other.id), location(other.location), amount(other.amount), handled_reactions(other.handled_reactions) {}

    void clear() {
        id = "";
        location.clear();
        amount = 0;
        handled_reactions.clear();
    }
};


/*******************************************/
/*********** Data enzyme type ************/
/*******************************************/

}
}

std::ostream& operator<<(std::ostream& os, const pmgbp::types::Address_t& to);
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& m);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::MetaboliteAmounts& m);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::BState_t& s);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::Way& s);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::ReactionInfo& r);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::Product& s);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::Reactant& s);
std::ostream& operator<<(std::ostream& os, const pmgbp::types::Enzyme& e);

#endif // PMGBP_TYPES_HPP