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
  os << "Address: " << r.location << endl;
  os << "substrates: " << r.substrate_sctry << endl;
  os << "products: " << r.products_sctry << endl;
  os << "KonSTP: " << r.konSTP << endl;
  os << "KonPTS" << r.konPTS << endl;
  os << "reversible: " << ((r.reversible) ? "true" : "false");
  return os;
}

ostream& operator<<(ostream& os, const Product& p) {
  os << "{";
  os << "\"message_type\":\"product\",";
  os << "\"metabolites\":[";

  bool separate = false;
  for (const auto& metabolite : p.metabolites) {
    if (separate) {
      os << ",";
    }

    os << "{";
    os << "\"id\":\"" << metabolite.first << "\",";
    os << "\"amount\":" << metabolite.second;
    os << "}";
  }

  os << "]";
  os << "}";
  return os;
}

ostream& operator<<(ostream& os, const Reactant& r) {
  os << "{";
  os << "\"message_type\":\"reactant\",";
  os << "\"for_reaction\":\"" << r.rid << "\",";
  os << "\"from\":\"" << r.from << "\",";
  os << "\"Way\":\"" << r.reaction_direction << "\",";
  os << "\"reaction_amount\":" << r.reaction_amount;
  os << "}";
  return os;
}