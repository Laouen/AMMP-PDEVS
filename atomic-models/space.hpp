#ifndef BOOST_SIMULATION_PDEVS_CONTROLER_H
#define BOOST_SIMULATION_PDEVS_CONTROLER_H
#include <string>
#include <utility>
#include <map>
#include <limits>
#include <memory>

#include <boost/simulation/pdevs/atomic.hpp> // boost simalator include

#include "../data-structures/types.hpp" // metabolite_info_t, enzyme_info_t, SState_t, Integer_t
#include "../data-structures/randomNumbers.hpp" // IntegerRandom_t


using namespace boost::simulation::pdevs;
using namespace boost::simulation;
using namespace std;

long double e = 2.71828182845904523536028747135266249775724709369995L;

template<class TIME, class MSG>
class space : public pdevs::atomic<TIME, MSG>
{
private:
  string                 _id;
  TIME                   _it;
  TIME                   _biomass_request_rate;
  Address_t              _biomass_address;
  SetOfMolecules_t       _metabolites;
  vector<enzyme_info_t>  _enzymes;
  double                 _volume;

  // task queue
  STaskQueue_t<TIME, MSG> _tasks;

  // used for uniform random numbers
  RealRandom_t<double>       _real_random;
  IntegerRandom_t<Integer_t> _integer_random;

  // constant variables
  TIME ZERO;

public:

  // This constructor start a new space without any metabolite, so, there isn't programed tasks when the space start.
  explicit space(
    const string                 other_id,
    const TIME                   other_it,
    const TIME                   other_biomass_request_rate,
    const Address_t&             other_biomass_address,
    const vector<enzyme_info_t>& other_enzymes,
    const double                 other_volume,
    ) noexcept :
  _id(other_id),
  _it(other_it),
  _biomass_request_rate(other_biomass_request_rate),
  _biomass_address(other_biomass_address),
  _enzymes(other_enzymes),
  _volume(other_volume),
  ZERO(0) {

    // The random atributes must be initilized with a random generator
    random_device real_rd; // Random generator variable
    _real_random.seed(real_rd());
    random_device integer_rd;
    _integer_random.seed(integer_rd());

    // just to confirm, the space and metabolites start empty.
    _tasks.clear();
    _metabolites.clear();
  }

  // This constructor start a new space with metabolites, the SetOfMolecules_t must have elements, is a precondition.
  // because of this, the space start with a selection task
  explicit space(
    const string                 other_id,
    const TIME                   other_it,
    const TIME                   other_biomass_request_rate,
    const Address_t&             other_biomass_address,
    const SetOfMolecules_t&      other_metabolites,
    const vector<enzyme_info_t>& other_enzymes,
    const double                 other_volume,
    ) noexcept :
  _id(other_id),
  _it(other_it),
  _biomass_request_rate(other_biomass_request_rate),
  _biomass_address(other_biomass_address),
  _metabolites(other_metabolites),
  _enzymes(other_enzymes),
  _volume(other_volume),
  ZERO(0) {
    assert((_metabolites.size() > 0) && "This constructor require a non empty map of metabolites.");

    // The random atributes must be initilized with a random generator
    random_device real_rd; // Random generator variable
    _real_random.seed(real_rd());
    random_device integer_rd;
    _integer_random.seed(integer_rd());

    // Becouse it start already with metabolites in the space, a selection task is programed.
    _tasks.clear();
    STask_t<TIME, MSG> new_task;
    new_task.time_left = _it;
    new_task.task_kind = SState_t::SELECTING_FOR_REACTION;
    this->insertTask(new_task);
  }

  // TODO: this for must explore the enzyme space taking carre of the amount of each kind and with the same probability. 
  // It also delete de enzyme to not be elegible again, and if some thing was send to this enzyme, it must be deleted from the space.
  // the function kon must be implemented following the formula in the thierry email.
  void internal() noexcept {

    STask_t<TIME, MSG> sr, sb; // sr = selected_reactants, sb = selected_biomass
    MSG cm;
    bool rs = false; // this boolean says if a SELECTIN_FOR_REACTION tasks has already happen or not.
    bool bs = false; // this boolean says if a SELECTIN_FOR_BIOMASS task has already happen or not.
    double son, pon, rv;

    this->updateTaskTimeLefts(_tasks.front().time_left);

    // For all the tasks that are happening now. because The tasks time_lefts were updated, the current time is zero.
    for (typename STaskQueue_t<TIME, MSG>::iterator it = _tasks.begin(); (it != _tasks.end()) && (it->time_left == ZERO); it = _tasks.erase(it)) {
      
      if ((it->task_kind == SState_t::SELECTING_FOR_REACTION) && !rs) {

        // set a new task to send the selected metabolites.
        sr.time_left  = ZERO;
        sr.task_kind  = SState_t::SENDING_REACTIONS;
        
        // look for metabolites to send, all the messages are allocated in coutput
        for (vector<enzyme_info_t>::iterator it = _enzymes.begin(); it != _enzymes.end(); ++it) {
          cm.clear();

          // calculating the son and pon
          if (this->thereAreEnaughFor(it->reactants_sctry)) son = this->kon(it->reactants_sctry, it->kon1);
          else son = 0;
          if (it->reversible && this->thereAreEnaughFor(it->products_sctry)) pon = this->kon(it->products_sctry, it->kon2);
          else pon = 0;

          // son + pon can't be greater than 1. If that happen, they are normalized.
          if (son + pon > 1) {
            son = son / (son + pon);
            pon = pon / (son + pon);
          }

          // the interval [0,1] is  divided in three pieces, [0 - son), [son - son+pon) and [son+pon - 1]
          // depending in which of these three sub-interval belongs rv, the enzyme me a product to subtract, subtract to produc or nothing.
          rv = _real_random.drawNumber(0.0, 1.0);
          if (rv < son) {
            // send message to enzyme for a subtract to product reaction
            cm.to = it->location;
            cm.react_direction = Way_t::STP;
            sr.msgs.push_back(cm);

            // update the taken metabolite from the space
            for (SetOfMolecules_t::iterator jt = it->reactants_sctry.begin(); jt != it->reactants_sctry.end(); ++jt) {
               _metabolites.at(jt->first) -= jt->second;
            }
          } else if (rv < kon) {
            // send message to enzyme for a subtract to product reaction
            cm.to = it->location;
            cm.react_direction = Way_t::PTS;
            sr.msgs.push_back(cm);

            // update the taken metabolite from the space
            for (SetOfMolecules_t::iterator jt = it->products_sctry.begin(); jt != it->products_sctry.end(); ++jt) {
               _metabolites.at(jt->first) -= jt->second;
            }
          }
        }
        
        // no more than one selection in a given time T;
        rs = true;
      } else if (it->task_kind == SState_t::SELECTING_FOR_BIOMAS && !bs) {

        // look for metabolites to send
        cm.to = _biomass_address;
        appendMetabolites(cm.metabolites, _metabolites);
        sb.msgs.push_back(cm);

        // set a new task for out() to send the selected metabolites.
        sb.time_left  = _biomass_request_rate;
        sb.task_kind  = SState_t::SENDING_BIOMAS;
        
        // once the metabolite are all send to biomass, there is no more metabolites in the space.
        _metabolites.clear();

        // no more than one selection in a given time T;
        bs = true;
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
        appendMetabolites(result.metabolites, _metabolites);
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

        new_task.time_left = _biomass_request_rate; // NOTE: There was ZERO before, and ZERO it's an error (generate zero-time loops), in the other model there is the same error.
        new_task.task_kind = SState_t::SELECTING_FOR_BIOMAS;
        this->insertTask(new_task);
      
      } else {

        appendMetabolites(_metabolites, it->metabolites);
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

  // TODO: generate test of all the helper functions
  // This funtion takes all the metabolites from om an amount grater than 0 and add them to m.
  void appendMetabolites(SetOfMolecules_t& m, const SetOfMolecules_t& om) {
  
    for (SetOfMolecules_t::const_iterator it = om.cbegin(); it != om.cend(); ++it) {
      addMetabolite(m, it->first, it->second);
    }
  }

  void addMetabolite(SetOfMolecules_t& m, string n, Integer_t a){

    if (a > 0) {
      if (m.find(n) != m.end()) {
        m.at(n) += a;
      } else {
        m.insert({n, a}); // TODO: change all the initializer_list because they don't work in windows
      }
    }
  }

  bool thereAreEnaughFor(const SetOfMolecules_t& stcry) const {
    bool result = true;

    for (SetOfMolecules_t::const_iterator it = stcry.begin(); it != stcry.end(); ++it) {
      if((_metabolites.find(it->first) == _metabolites.end()) || (_metabolites.at(it->first) < it->second)) {
        result = false;
        break;
      }
    }

    return result;
  }

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
  // this function look if there is metabolites to send and in this case, if the space have not alreafy programed a selection task to send metabolites, it will program one.
  void setNextSelection() {
    STask_t<TIME, MSG> new_selection;
    
    if ( this->thereIsMetabolites() && !this->thereIsNextSelection() ) {

      new_selection.time_left = _it;
      new_selection.task_kind = SState_t::SELECTING_FOR_REACTION;
      this->insertTask(new_selection);
    }
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

  void insertTask(const STask_t<TIME, MSG>& t) {

    typename STaskQueue_t<TIME, MSG>::iterator it = lower_bound(_tasks.begin(), _tasks.end(), t);
    _tasks.insert(it, t);
  }

  void updateTaskTimeLefts(TIME t){

    for (typename STaskQueue_t<TIME, MSG>::iterator it = _tasks.begin(); it != _tasks.end(); ++it) {
      it->time_left -= t;
    }
  }
};
  

#endif // BOOST_SIMULATION_PDEVS_CONTROLER_H
