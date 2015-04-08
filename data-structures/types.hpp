#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H


#include <list>
#include <vector>
#include <map>

/******************************************/
/********** Type definations **************/
/******************************************/

using namespace std;


enum class RState { REJECTING, REACTING, SELECTING };
enum class SState { SENDING, SELECTING };
enum class Way { RTP, PTR };

using Integer = unsigned long long;
using SetOfMolecules  = map<string, Integer>;

template<class TIME>
struct Task {
  TIME                time_left;
  RState              task_kind;
  SetOfMolecules      rejected;  
  pair<Way, Integer>  reaction;

  Task() {}

  Task(const Task& other) {
    time_left = other.time_left;
    task_kind = other.task_kind;
    rejected  = other.rejected;
    reaction  = other.reaction;
  }

  inline bool operator<(const Task<TIME>& o) {

    bool result;
    if (time_left != o.time_left) {
      
      result = (time_left < o.time_left);
    } else {

      if ((task_kind == RState::SELECTING) && (o.task_kind != RState::SELECTING)) result = true;
      else if ((task_kind == RState::REJECTING) && (o.task_kind == RState::REACTING)) result = true;
      else result = false;
    }

    return result;
  }

  inline bool operator==(const Task<TIME>& o) {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if (task_kind == RState::REJECTING) {
      result = result && (rejected == o.rejected);
    }

    if (task_kind == RState::REACTING) {
      result = result && (reaction == o.reaction);
    }    

    return result;
  }

};

template<class TIME>
using TaskQueue = list< Task<TIME> >;


/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H