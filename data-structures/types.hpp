#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H

/******************************************/
/********** Type definations **************/
/******************************************/

using namespace std;

enum RState { REJECTING, REACTING, SELECTING, IDLE };
enum Way { RTP, PTR };

using SetOfMolecules  = map<string, int>;

template<class TIME>
struct Task {
  TIME            time_left;
  RState          task_kind;
  SetOfMolecules  rejected;  
  pair<Way, int>  reaction;

  inline bool operator<(const Task<TIME>& o) {
    if (time_left < o.time_left) return true;
    else if ((task_kind == SELECTING) && (o.task_kind != SELECTING)) return true;
    else if ((task_kind == REJECTING) && (o.task_kind != SELECTING) && (o.task_kind != REJECTING)) return true;
    return false;
  }
};

template<class TIME>
using TaskQueue = list< Task<TIME> >;

/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H