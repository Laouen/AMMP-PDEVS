#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>

#include "space.hpp"

using namespace std;

// TODO: correctly separate by namespaces like space::, reaction::, etc.
// TODO(new line): Models and model internal states also should be under the same namespaces

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
};

/******************************************/
/********** Enums and renames *************/
/******************************************/

enum class RState_t { REJECTING = 1, REACTING = 0 };
enum class BState_t { ENOUGH = 0, NOT_ENOUGH = 1, IDLE = 2, WAITING = 3 };
enum class Way_t { STP, PTS };

using Integer = unsigned long long;
using MetaboliteAmounts  = map<string, Integer>;

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
/**************** RTask_t ******************/
/*******************************************/

template<class TIME, class MSG>
struct RTask_t {
  RState_t    task_kind;
  TIME        time_left;
  Way_t       direction;
  Integer   amount;
  vector<MSG> toSend;

  RTask_t() {}

  RTask_t(const TIME& other_t, const Way_t& other_d, const Integer& other_a)
  : task_kind(RState_t::REACTING), time_left(other_t), direction(other_d), amount(other_a) {}

  RTask_t(const TIME& other_t, const vector<MSG>& other_ts)
  : task_kind(RState_t::REJECTING), time_left(other_t), toSend(other_ts) {};

  RTask_t(const RTask_t<TIME, MSG>& other) 
  : task_kind(other.task_kind), time_left(other.time_left), direction(other.direction), amount(other.amount), toSend(other.toSend) {}

  inline bool operator<(const RTask_t<TIME, MSG>& o) const {

    bool result;
    if (time_left != o.time_left) result = (time_left < o.time_left);
    else                          result = (task_kind < o.task_kind);

    return result;
  }

  inline bool operator==(const RTask_t<TIME, MSG>& o)  const { 

    if (task_kind == RState_t::REJECTING) {
      return (time_left == o.time_left) && (toSend == o.toSend);
    } else {
      return (time_left == o.time_left) && (direction == o.direction) && (amount == o.amount);
    }
  }

};

template<class TIME, class MSG>
using RTaskQueue_t = list< RTask_t<TIME, MSG> >;

/*******************************************/
/************** End RTask_t ****************/
/*******************************************/


/*******************************************/
/**************** Message_t ****************/
/*******************************************/

using Address_t = list<string>;

struct Message_t {
  
  Address_t to;
  
  // fields for reactions
  string from;
  Way_t react_direction;
  Integer react_amount;
  
  // fields for spaces
  MetaboliteAmounts metabolites;
  
  // field for show request
  bool show_request;
  
  // field for biomass reaction request
  bool biomass_request;

  Message_t(const Address_t& other_to, string other_from, const MetaboliteAmounts& other_m, Way_t other_rd, Integer other_ra, bool other_sr, bool other_br)
  : to(other_to), from(other_from), metabolites(other_m), react_direction(other_rd), react_amount(other_ra), show_request(other_sr), biomass_request(other_br) {}

  Message_t()
  :to(), from(""), metabolites(), react_direction(), react_amount(), show_request(false), biomass_request(false) {}

  Message_t(const Message_t& other)
  : to(other.to), from(other.from), metabolites(other.metabolites), react_direction(other.react_direction), react_amount(other.react_amount), show_request(other.show_request), biomass_request(other.biomass_request) {}

  Message_t(Message_t& other)
  : to(other.to), from(other.from), metabolites(other.metabolites), react_direction(other.react_direction), react_amount(other.react_amount), show_request(other.show_request), biomass_request(other.biomass_request) {}

  void clear() {
    to.clear();
    from = "";
    react_amount = 0;
    metabolites.clear();
    show_request = false;
    biomass_request = false;
  }

  bool empty() {
    return to.empty() && metabolites.empty() && (from == "") && (react_amount == 0);
  }
};

ostream& operator<<(ostream& os, const Message_t& msg);
ostream& operator<<(ostream& os, const Address_t& to);
ostream& operator<<(ostream& os, const vector<string>& m);
ostream& operator<<(ostream& os, const MetaboliteAmounts& m);
ostream& operator<<(ostream& os, const SpaceState& s);
ostream& operator<<(ostream& os, const BState_t& s);
ostream& operator<<(ostream& os, const Way_t& s);

/*******************************************/
/************** End Message_t **************/
/*******************************************/


/*******************************************/
/************* Data info type **************/
/*******************************************/

// TODO this type must be removed after the new implementation.
struct reaction_info_t {

  string id;
  ReactionAddress location;
  MetaboliteAmounts  substrate_sctry;
  MetaboliteAmounts  products_sctry;
  double konSTP;
  double konPTS;
  double koffPTS;
  double koffSTP;
  bool reversible;

  reaction_info_t()
  : id(), location(), konSTP(1), konPTS(1), reversible(false) {}

  reaction_info_t(
    string                  other_id,
    const Address_t&        other_location,
    const MetaboliteAmounts& other_substrate_sctry,
    const MetaboliteAmounts& other_products_sctry,
    double                  other_konSTP,
    double                  other_konPTS,
    double                  other_koffSTP,
    double                  other_koffPTS,
    bool                    other_reversible
    ) : location(other_location), substrate_sctry(other_substrate_sctry), products_sctry(other_products_sctry), konSTP(other_konSTP), konPTS(other_konPTS), koffPTS(other_koffPTS), koffSTP(other_koffSTP), reversible(other_reversible) {}

  reaction_info_t(const reaction_info_t& other)
  : id(other.id), location(other.location), substrate_sctry(other.substrate_sctry), products_sctry(other.products_sctry), konSTP(other.konSTP), konPTS(other.konPTS), koffPTS(other.koffPTS), koffSTP(other.koffSTP), reversible(other.reversible) {}

  void clear() {
    id = "";
    location.clear();
    substrate_sctry.clear();
    products_sctry.clear();
    konSTP = 0;
    konPTS = 0;
    reversible = false;
  }

  bool empty() {
    return (id == "") && location.empty() && substrate_sctry.empty() && products_sctry.empty();
  }
};

ostream& operator<<(ostream& os, const reaction_info_t& r);


/*******************************************/
/*********** End Data info type ************/
/*******************************************/

/*******************************************/
/*********** Data enzyme type **************/
/*******************************************/

struct Enzyme {
  string      id;
  Integer   amount;
  map<string, reaction_info_t> handled_reactions;

  Enzyme()
  : id(""), amount(0), handled_reactions() {};

  Enzyme(string other_id, const Integer& other_amount, const map<string, reaction_info_t>& other_handled_reactions)
  : id(other_id), amount(other_amount), handled_reactions(other_handled_reactions) {}

  Enzyme(const Enzyme& other)
  : id(other.id), amount(other.amount), handled_reactions(other.handled_reactions) {}

  void clear() {
    id = "";
    amount = 0;
    handled_reactions.clear();
  }
};


/*******************************************/
/*********** Data enzyme type ************/
/*******************************************/





/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H