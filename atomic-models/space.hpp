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

  void internal() noexcept {

    STask_t<TIME, MSG> sr, sb;
    MSG cm;
    vector<MSG> coutput;
    vector<Integer_t> distributed_reactants = {};
    bool reaction_selected      = false;
    bool biomass_selected       = false;
    // Updating time left

    this->updateTaskTimeLefts(_tasks.front().time_left);

    for (typename STaskQueue_t<TIME, MSG>::iterator it = _tasks.begin(); it->time_left == ZERO; it = _tasks.erase(it)) {
      if ((it->task_kind == SState_t::SELECTING_FOR_REACTION) && !reaction_selected) {

        // look for metabolites to send
        for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {

          if (this->weightedRandomBool(it->second.amount)){

            distributed_reactants.clear();
            distributed_reactants.resize(it->second.enzymes.size());
            randomDistribution(distributed_reactants, it->second.amount);
            cm.specie = it->first; 

            for (int i = 0; i < distributed_reactants.size(); ++i) {
              
              cm.amount          = distributed_reactants[i];
              cm.to              = it->second.enzymes[i];
              coutput.push_back(cm);
              it->second.amount  -= distributed_reactants[i];
            }
          }
        }

        // set a new task for out() to send the selected metabolites.
        sr.time_left  = ZERO;
        sr.task_kind  = SState_t::SENDING_REACTIONS;
        sr.to_send    = coutput;
        
        // no more than one selection in a given time T;
        reaction_selected = true;
      } else if (it->task_kind == SState_t::SELECTING_FOR_BIOMAS && !biomass_selected) {

        // look for metabolites to send
        for (map<string, metabolite_info_t>::iterator it = _metabolites.begin(); it != _metabolites.end(); ++it) {

          cm.specie  = it->first; 
          cm.amount  = it->second.amount;
          cm.to      = _biomass_address;
          coutput.push_back(cm);
          it->second.amount = 0;

        }

        // set a new task for out() to send the selected metabolites.
        sb.time_left  = _biomass_request_rate;
        sb.task_kind  = SState_t::SENDING_BIOMAS;
        sb.to_send    = coutput;
        
        // no more than one selection in a given time T;
        biomass_selected = true;
      }
    }


    // inserting new tasks
    if (!sr.to_send.empty()) this->insertTask(sr);
    if (!sb.to_send.empty()) this->insertTask(sb);

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
    MSG current_message;
    TIME current_time  = _tasks.front().time_left;

    for (typename STaskQueue_t<TIME, MSG>::const_iterator it = _tasks.cbegin(); it->time_left == current_time; ++it) {

      if ((it->task_kind == SState_t::SENDING_BIOMAS) || (it->task_kind == SState_t::SENDING_REACTIONS)) {

        for (typename vector<MSG>::const_iterator mt = it->to_send.cbegin(); mt != it->to_send.cend(); ++mt) {       
          result.push_back(*mt);
        }
      } else if (it->task_kind == SState_t::SHOWING) {

        // look for metabolites to send
        for (map<string, metabolite_info_t>::const_iterator it = _metabolites.cbegin(); it != _metabolites.cend(); ++it) {

          current_message.to  = {"output", _id}; 
          current_message.specie  = it->first; 
          current_message.amount  = it->second.amount;
          result.push_back(current_message);
        }
      }
    }

    return result;
  }

  void external(const vector<MSG>& mb, const TIME& t) noexcept {

    STask_t<TIME, MSG> new_task;

    this->updateTaskTimeLefts(t);

    for (typename vector<MSG>::const_iterator it = mb.cbegin(); it != mb.cend(); ++it) {

      if (it->show_request) {
        
        new_task.time_left = ZERO;
        new_task.task_kind = SState_t::SHOWING;
        this->insertTask(new_task);

      } else if (it->biomass_request) {

        new_task.time_left = _biomass_request_rate; // NOTE: There was ZERO before, and ZERO it's an error (generate zero-time loops), in the other model there is the same error.
        new_task.task_kind = SState_t::SELECTING_FOR_BIOMAS;
        this->insertTask(new_task);
      
      } else {

        this->addToMetabolites(it->specie, it->amount);
      }
    }

    // if some metabolites have just arrived, a selection task must be programed.
    this->setNextSelection();
  }

  virtual void confluence(const std::vector<MSG>& mb, const TIME& t) noexcept {

    external(mb, t);
    internal(); 
  }

  /***************************************
  ********* helper functions *************
  ***************************************/

  void addToMetabolites(string n, Integer_t a){

    if (a > 0) {
      if (_metabolites.find(n) != _metabolites.end()) {
        _metabolites.at(n) += a;
      } else {
        _metabolites.insert({n, a});
      }
    }
  }

  vector<Address_t> lookForEnzymes(string n) const {
    vector<Address_t> result;

    for (map<string, enzyme_info_t>::const_iterator it = _enzymes.cbegin(); it != _enzymes.cend(); ++it) {
      if (belong(n, it->second.reactants)) {
        result.push_back(it->second.location);
      }
    }

    return result;
  }

  bool belong(const string& n, const vector<string>& ls) const {
    bool result = false;

    for (vector<string>::const_iterator it = ls.cbegin(); it != ls.cend(); ++it) {
      if (n == *it) {
        result = true;
        break;
      }
    }
    return result;
  }

  bool weightedRandomBool(Integer_t a) {

    double proportion;

    if (a > 0)  proportion = _volume / (double)a;
    else        proportion = numeric_limits<double>::infinity();
    
    double threshold  = (double)1.0 / pow( (double)e, (double)_factor*proportion );


    return (_real_random.drawNumber(0.0, 1.0) < threshold);
  }

  void randomDistribution(vector<Integer_t>& ds, Integer_t a) {
    
    if (!ds.empty()) {

      for (int i = 0; i < ds.size(); ++i) {
        ds[i] = 0;
      }

      for (Integer_t i = 0; i < a; ++i) {
        ds[_integer_random.drawNumber(0, ds.size() - 1)] += 1;
      }
    }
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
