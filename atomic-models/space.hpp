#ifndef BOOST_SIMULATION_PDEVS_SPACE_H
#define BOOST_SIMULATION_PDEVS_SPACE_H
#include <string>
#include <utility>
#include <map>
#include <limits>
#include <memory>
#include <math.h>
// specially to shuffle the current_enzyme vector.
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // reaction_info_t, SState_t, Integer_t
#include "../data-structures/randomNumbers.hpp" // IntegerRandom_t


using namespace boost::simulation::pdevs;
using namespace boost::simulation;
using namespace std;

template<class TIME, class MSG>
class space : public pdevs::atomic<TIME, MSG>
{
private:
  string                  _id;
  TIME                    _it;
  TIME                    _br;
  Address_t               _biomass_address;
  SetOfMolecules_t        _metabolites;
  map<string, enzyme_t>   _enzymes;
  double                  _volume;

  // task queue
  STaskQueue_t<TIME, MSG> _tasks;

  // used for uniform random numbers
  RealRandom_t<double>       _real_random;
  IntegerRandom_t<Integer_t> _integer_random;

public:

  // precondition: the setOfMetabolites_t given as parameter must contain all the metabolites in the space defined
  explicit space(
    const string                          other_id,
    const TIME                            other_it,
    const TIME                            other_br,
    const Address_t&                      other_biomass_address,
    const SetOfMolecules_t                other_metabolites,
    const map<string, reaction_info_t>&   other_enzymes,
    const double                          other_volume
    ) noexcept :
  _id(other_id),
  _it(other_it),
  _br(other_br),
  _biomass_address(other_biomass_address),
  _metabolites(other_metabolites),
  _enzymes(other_enzymes),
  _volume(other_volume) {

    // The random atributes must be initilized with a random generator
    random_device real_rd; // Random generator variable
    _real_random.seed(real_rd());
    random_device integer_rd;
    _integer_random.seed(integer_rd());

    // just to confirm, the space and metabolites start empty.
    _tasks.clear();
  }

  // TODO: check the correct ordering of the tasks, if there is a biomass request and a selec for reaction, the biomass request must win.
  void internal() noexcept {

    MSG cm;
    STask_t<TIME, MSG> sr, sb; // sr = selected_reactants, sb = selected_biomass
    bool srah = false; // this boolean says if a SELECTING_FOR_REACTION tasks has already happen or not.
    bool sbah = false; // this boolean says if a SELECTING_FOR_BIOMASS task has already happen or not.

    this->updateTaskTimeLefts(_tasks.front().time_left);

    // For all the tasks that are happening now. because The tasks time_lefts were updated, the current time is ZERO.
    for (typename STaskQueue_t<TIME, MSG>::iterator it = _tasks.begin(); !_tasks.empty() && (it->time_left == ZERO); it = _tasks.erase(it)) {
      
      if ((it->task_kind == SState_t::SELECTING_FOR_REACTION) && !srah) {

        srah = true; // no more than one selection in a given time T;

        // set a new task to send the selected metabolites.
        sr.task_kind  = SState_t::SENDING_REACTIONS;
        sr.time_left  = ZERO;

        this->selectMetalobitesToReact(sr.msgs);
        unifyMessages(sr.msgs);
      } else if ((it->task_kind == SState_t::SELECTING_FOR_BIOMAS) && !sbah) {
        // no more than one selection in a given time T;
        sbah = true;

        // look for metabolites to send
        cm.to = _biomass_address;
        addMultipleMetabolites(cm.metabolites, _metabolites);

        // set a new task for out() to send the selected metabolites.
        sb.time_left  = _br;
        sb.task_kind  = SState_t::SENDING_BIOMAS;
        sb.msgs.push_back(cm);
        
        // once the metabolite are all send to biomass, there is no more metabolites in the space.
        this->removeAllMetabolites(); // TODO this method must be implemented
      }
    }


    // inserting new tasks
    if (!sr.msgs.empty()) this->insertTask(sr);
    if (!sb.msgs.empty()) this->insertTask(sb);

    // setting new selection
    this->setNextSelection();
  }

  TIME advance() const noexcept {

    TIME result;
    if (!_tasks.empty()) result = _tasks.front().time_left;
    else                 result = pdevs::atomic<TIME, MSG>::infinity;

    return result;
  }

  vector<MSG> out() const noexcept {

    vector<MSG> result;
    MSG b_msg;
    TIME current_time  = _tasks.front().time_left;

    // for all the tasks that ocurr in the current time. These tasks are processed now.
    for (typename STaskQueue_t<TIME, MSG>::const_iterator it = _tasks.cbegin(); (it != _tasks.end()) && (it->time_left == current_time); ++it) {

      if ((it->task_kind == SState_t::SENDING_BIOMAS) || (it->task_kind == SState_t::SENDING_REACTIONS)) {

        result.insert(result.end(), it->msgs.cbegin(), it->msgs.cend()); //TODO: test this method of insert with constant iterators.

      } else if (it->task_kind == SState_t::SHOWING) {

        // Send all the current free metabolites in the space
        b_msg.to  = {"output", _id};
        addMultipleMetabolites(b_msg.metabolites, _metabolites);
        result.push_back(b_msg);
      }
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    STask_t<TIME, MSG> new_task;

    this->updateTaskTimeLefts(t);

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {


      if (it->show_request) {
        
        new_task.time_left = ZERO; // TODO: aks to Gabriel if it shouldn't be 0.
        new_task.task_kind = SState_t::SHOWING;
        this->insertTask(new_task);

      } else if (it->biomass_request) {

        new_task.time_left = _br; // NOTE: There was ZERO before, and ZERO it's an error (generate zero-time loops), in the other model there is the same error.
        new_task.task_kind = SState_t::SELECTING_FOR_BIOMAS;
        this->insertTask(new_task);
      
      } else {

        addMultipleMetabolites(_metabolites, it->metabolites);
      }
    }

    // if some metabolites have just arrived (the third part of the if has happen), a selection task must be programed.
    this->setNextSelection();
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
    internal(); 
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  /************** addMultipleMetabolites ****************************/

  // TODO: generate test of all the helper functions
  // This funtion takes all the metabolites from om with an amount grater than 0 and add them to m.
  void addMultipleMetabolites(SetOfMolecules_t& m, const SetOfMolecules_t& om) {
  
    for (SetOfMolecules_t::const_iterator it = om.cbegin(); it != om.cend(); ++it) {
      // if the metabolite isn't defined in the compartment is not from there and an error ocurre.
      m.at(it->first) += it->second;
    }
  }

  /************** SetNextSelection *********************************/

  // this function tells if there is or not metabolites in the space.
  bool thereIsMetabolites() const {

    bool result = false;
    for (SetOfMolecules_t::const_iterator it = _metabolites.cbegin(); it != _metabolites.cend(); ++it) {
      if (it->second > 0) {
        result = true;
        break;
      }
    }
    return result;
  }

  // this function looks if there is a selection task already programed.
  bool thereIsNextSelection() const {
    bool result = false;

    for (typename STaskQueue_t<TIME, MSG>::const_iterator it = _tasks.cbegin(); it != _tasks.cend(); ++it) {
      if ((it->task_kind == SState_t::SELECTING_FOR_REACTION) && (it->time_left <= _it)) {
        result = true;
        break;
      }
    }

    return result;
  }

  // this function look if there is metabolites to send and in this case, if the space have not alreafy programed a selection task to send metabolites, it will program one.
  void setNextSelection() {
    STask_t<TIME, MSG> new_selection;
    
    if ( this->thereIsMetabolites() && !this->thereIsNextSelection() ) {

      new_selection.time_left = _it;
      new_selection.task_kind = SState_t::SELECTING_FOR_REACTION;
      this->insertTask(new_selection);
    }
  }

  /***************** selectMetabolitesToReact ********************/

  bool thereAreEnaughFor(const SetOfMolecules_t& stcry) const {
    bool result = true;

    for (SetOfMolecules_t::const_iterator it = stcry.begin(); it != stcry.end(); ++it) {
      if((_metabolites.find(it->first) != _metabolites.end()) && (_metabolites.at(it->first) < it->second)) {
        result = false;
        break;
      }
    }

    return result;
  }

  double summAll(const map<string, double>& ons) const {
    double result = 0;

    for (map<string, double>::const_iterator i = ons.cbegin(); i != ons.cend(); ++i) {
      result += i->second;
    }

    return result;
  }

  void unfoldEnzymes(vector<string>& ce) const {

    for (map<string, enzyme_t>::const_iterator it = _enzymes.cbegin(); it != _enzymes.cend(); ++it) {
      ce.insert(ce.end(), it->second.amount, it->second.id);
    }
  }

  // TODO test this function specially
  void shuffleEnzymes(vector<string>& ce) const {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(ce.begin(), ce.end(), g);
  }

  // TODO test this function specially
  double bindingTreshold(const SetOfMolecules_t& sctry, double kon) const {

    // calculation of the consentrations [A][B][C]
    double consentration = 1.0;
    for (SetOfMolecules_t::const_iterator it = sctry.cbegin(); it != sctry.cend(); ++it) {
      if (_metabolites.find(it->first) != _metabolites.end()) {
        consentration *= _metabolites.at(it->first) / (L * _volume);  
      }
    }

    if (consentration = 0.0) 
      return 0.0;

    return exp(-(1.0 / (consentration*kon)));;
  }

  // TODO for testing
  void collectOns(map<string, reaction_info_t>& r, map<string, double>& s, map<string, double>& p) {

    for (map<string, reaction_info_t>::iterator it = r.begin(); it != r.end(); ++it) {
      
      // calculating the sons and pons
      if (this->thereAreEnaughFor(it->second.substrate_sctry)) s.insert({it->first, this->bindingTreshold(it->second.substrate_sctry, it->second.konSTP)});
      else s.insert({it->first, 0});
      if (it->second.reversible && this->thereAreEnaughFor(it->second.products_sctry)) p.insert({it->first, this->bindingTreshold(it->second.products_sctry, it->second.konPTS)});
      else p.insert({it->first, 0});
    }
  }
  // TODO for testing
  void normalize(map<string, double>& ons, double t) {
    for (map<string, double>::iterator i = ons.begin(); i != ons.end(); ++i) {
      i->second = i->second / t;
    }
  }


  void selectMetalobitesToReact(vector<MSG>& m) {
    MSG cm;
    double rv, total, partial;
    map<string, double> sons, pons;
    enzyme_t en;
    reaction_info_t re;
    vector<string> enzyme_IDs;

    this->unfoldEnzymes(enzyme_IDs); // all the enzyme are considered individualy not grouped by kind
    shuffleEnzymes(enzyme_IDs); // the enzymes are iterated randomly

    for (vector<string>::iterator it = enzyme_IDs.begin(); it != enzyme_IDs.end(); ++it) {
      en = _enzymes.at(*it);

      collectOns(en.reacts, sons, pons);


      // sons + pons can't be greater than 1. If that happen, they are normalized
      // if sons + pons is smaller than 1, there is a chanse that the enzyme does'nt react
      total = summAll(sons) + sumAll(pons);
      if (total > 1) {
        normalize(sons, total);
        normalize(pons, total);
      }
      
      // the interval [0,1] is  divided in pieces: 
      // [0,son1), [son1, son1+son2), ... , [son1 + ... + sonk, son1 + ... + sonk + pon1), ... ,[son1 + ... + sonk + pon1 + ... + ponk, 1)
      // depending on which of those sub-interval rv belongs, the enzyme triggers the correct reaction or do nothing.

      rv = _real_random.drawNumber(0.0, 1.0);
      partial = 0;
      
      for (map<string, double>::iterator i = sons.begin(); i != sons.end(); ++i) {
        
        partial += i->second;
        if (rv < partial) {
          // send message to trigger the reaction
          re = en.reacts.at(i->first);
          cm.clear();
          cm.to = re.location;
          cm.from = _id;
          cm.react_direction = Way_t::STP;
          cm.react_amount = 1;
          m.push_back(cm);
          break;
        } 
      }

      // update the metabolite amount in the space
      if (!re.empty()) {
        for (SetOfMolecules_t::iterator it = re.substrate_sctry.begin(); it != re.substrate_sctry.end(); ++it) {
          if (_metabolites.find(it->first) != _metabolites.end()) {
            _metabolites.at(it->first) -= it->second;
          }
        }
      } else { // if re is empty, then no one of the STP reactions have ben triggered and the search continue with the PTS reactions.
        for (map<string, double>::iterator i = pons.begin(); i != pons.end(); ++i) {
        
          partial += i->second;
          if (rv < partial) {
            // send message to trigger the reaction
            re = en.reacts.at(i->first);
            cm.clear();
            cm.to = re.location;
            cm.from = _id;
            cm.react_direction = Way_t::PTS;
            cm.react_amount = 1;
            m.push_back(cm);
            break;
          } 
        }
      
        // update the metabolite amount in the space
        if (!re.empty()) {
          for (SetOfMolecules_t::iterator it = re.products_sctry.begin(); it != re.products_sctry.end(); ++it) {
            if (_metabolites.find(it->first) != _metabolites.end()) {
              _metabolites.at(it->first) -= it->second;
            }
          }
        }
      }
    }
  }

  /*****************************************************************/

  void insertTask(const STask_t<TIME, MSG>& t) {

    typename STaskQueue_t<TIME, MSG>::iterator it = lower_bound(_tasks.begin(), _tasks.end(), t);
    _tasks.insert(it, t);
  }

  void updateTaskTimeLefts(TIME t){

    for (typename STaskQueue_t<TIME, MSG>::iterator it = _tasks.begin(); it != _tasks.end(); ++it) {
      it->time_left -= t;
    }
  }

  void unifyMessages(vector<MSG>& m) const {

    map<Address_t, MSG> unMsgs;

    for (typename vector<MSG>::iterator it = m.begin(); it != m.end(); ++it) {
      insertMessageUnifying(unMsgs, *it);
    }

    m.clear();

    for (typename map<Address_t, MSG>::iterator it = unMsgs.begin(); it != unMsgs.end(); ++it) {
      m.push_back(it->second);
    }
  }

  void insertMessageUnifying(map<Address_t, MSG>& ms, MSG& m) const {

    if (m.react_amount > 0) {
      if (ms.find(m.to) != ms.end()) {
        ms.at(m.to).react_amount += m.react_amount;
      } else {
        ms.insert({m.to, m}); // TODO: change all the initializer_list because they don't work on windows
      }
    }
  }
};
  

#endif // BOOST_SIMULATION_PDEVS_SPACE_H
