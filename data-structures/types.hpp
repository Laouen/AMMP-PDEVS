#ifndef BOOST_SIMULATION_TYPES_H
#define BOOST_SIMULATION_TYPES_H

/******************************************/
/********** Type definations **************/
/******************************************/

using namespace std;

enum RState { REJECTING, REACTING, SELECTING };
enum Way { RTP, PTR };

using SetOfMolecules  = map<string, int>;

template<class TIME>
struct Task {
  TIME            time_left;
  RState          task_kind;
  SetOfMolecules  rejected;  
  pair<Way, int>  reaction;

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

      if ((task_kind == SELECTING) && (o.task_kind != SELECTING)) result = true;
      else if ((task_kind == REJECTING) && (o.task_kind == REACTING)) result = true;
      else result = false;
    }

    return result;
  }

  inline bool operator==(const Task<TIME>& o) {

    bool result;

    result = (time_left == o.time_left);
    result = result && (task_kind == o.task_kind);

    if (task_kind == REJECTING) {
      result = result && (rejected == o.rejected);
    }

    if (task_kind == REACTING) {
      result = result && (reaction == o.reaction);
    }    

    return result;
  }

};

template<class TIME>
using TaskQueue = list< Task<TIME> >;

ostream& show(ostream& os, const SetOfMolecules& to) {

  os << "[";

  SetOfMolecules::const_iterator it = to.cbegin();
  while ( it != to.cend()) {
    os << "(" << it->first << "," << it->second << ")";
    ++it;
    if (it != to.cend()) os << ",";
  }
  os << "]";
  return os;
}

template<class TIME>
ostream& show(ostream& os, const Task<TIME>& to) {

  string kind, w;
  if (to.task_kind == SELECTING)        kind = "SELECTING";
  else if (to.task_kind == REJECTING)   kind = "REJECTING";
  else if (to.task_kind == REACTING)    kind = "REACTING";

  os << "Task Kind: " << kind << endl;
  os << "Time left: " << to.time_left << endl;

  if(to.task_kind == REJECTING) {
    show(os, to.rejected);
  } else if (to.task_kind == REACTING) {

    if (to.reaction.first == RTP)       w = "RTP";
    else if (to.reaction.first == PTR)  w = "PTR";

    os << "way: " << w << " amount: " << to.reaction.second;
  }

  return os;
}

template<class TIME>
ostream& show(ostream& os, const TaskQueue<TIME>& to) {

  os << "Current Tasks in the queue: ";
  for (typename TaskQueue<TIME>::const_iterator it = to.cbegin(); it != to.cend(); ++it) {
    os << endl << endl;
    show(cout, *it);
  }
  return os;
}

/******************************************/
/******** End type definations ************/
/******************************************/

#endif // BOOST_SIMULATION_TYPES_H