#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <boost/simulation/pdevs/atomic.hpp>
#include <string>
#include <utility>
#include <map>
#include "../data-structures/reaction_input.hpp"


using namespace boost::simulation::pdevs;
using namespace std;

/*************************************
*********type definations*************
*************************************/

/******** Basic types ************/
typedef pair<string, double> molecule;

/******** vector types ************/
typedef vector<molecule> vectorOfMolecules;

template<class TIME, class MSG>
class controler : public atomic<TIME, MSG>
{
private:

  map< string, vector<string> > _reaction_reactants;
  vector<MSG>                   _filtered_input;
  bool                          _hasToSend;

public:

  explicit controler(const vector<reaction_input> reactions) noexcept :
  _hasToSend(false) {
    _filtered_input.clear();
    _reaction_reactants.clear();
    
    for (vector<reaction_input>::const_iterator it = reactions.cbegin(); it != reactions.cend(); ++it){
      _reaction_reactants[it->name] = {};

      for (vectorOfMolecules::const_iterator jt = it->reactants.cbegin(); jt != it->reactants.cend(); ++jt){
        _reaction_reactants[it->name].push_back(jt->first);
      }
      
      for (vectorOfMolecules::const_iterator jt = it->enzymes.cbegin(); jt != it->enzymes.cend(); ++jt){
        _reaction_reactants[it->name].push_back(jt->first);
      }
    }
  }

  void internal() noexcept { 
    _hasToSend = false;
    _filtered_input.clear();
  }

  TIME advance() const noexcept {
    
    return _hasToSend?TIME(0):atomic<TIME, MSG>::infinity;
  }

  vector<MSG> out() const noexcept {

    return _filtered_input; 
  }

  void external(const std::vector<MSG>& mb, const TIME& t) noexcept {
    
    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it){
      
      molecule curr_msg     = boost::any_cast<molecule>(*it); 
      int reactions_left    = reactionsWithThisReactant(curr_msg.first, _reaction_reactants);
      int reactant_amount   = curr_msg.second;

      if (reactions_left <=0) {

        string reaction_name              = "Output";
        int amount_to_send                = reactant_amount;
        molecule m_to_send                = make_pair(curr_msg.first, amount_to_send);
        pair<string, molecule> new_output = make_pair(reaction_name, m_to_send);

        _filtered_input.push_back(boost::any(new_output));
     
      } else {

        for (map< string, vector<string> >::iterator it = _reaction_reactants.begin(); it != _reaction_reactants.end(); ++it){
          if (hasReactant(it->second, curr_msg.first)) {

            string reaction_name              = it->first;
            int amount_to_send                = (reactions_left == 1) ? reactant_amount : (reactant_amount / reactions_left);
            molecule m_to_send                = make_pair(curr_msg.first, amount_to_send);
            pair<string, molecule> new_output = make_pair(reaction_name, m_to_send);

            _filtered_input.push_back(boost::any(new_output));

            reactant_amount -= amount_to_send;
            reactions_left  -= 1;
          }
        }    
      }
    }

    _hasToSend = (_filtered_input.size() > 0);
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
    internal();
  }

  int reactionsWithThisReactant(string r, map< string, vector<string> > reactions) {
    
    int result = 0;
    for (map< string, vector<string> >::iterator it = reactions.begin(); it != reactions.end(); ++it){
      if (hasReactant(it->second, r)) {
        result += 1;
      }
    }

    return result;
  }

  bool hasReactant(vector<string> rs, string r) {

    bool result = false;
    for (vector<string>::iterator it = rs.begin(); it != rs.end(); ++it) {
      if (*it == r) {
        result = true;
      }
    }

    return result;
  }

};

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
