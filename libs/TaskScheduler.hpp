//
// Created by lao on 17/09/17.
//

#ifndef PMGBP_PDEVS_TASKSCHEDULER_HPP
#define PMGBP_PDEVS_TASKSCHEDULER_HPP

#include <list>
#include <algorithm>

template <class TIME, class ELEMENT>
struct Tasks {
    TIME time_left;
    std::list<ELEMENT> task_elements;

    Tasks(TIME time, ELEMENT task_element) {
        this->time_left = time;
        this->task_elements.push_back(task_element);
    }
};

template <class TIME, class ELEMENT>
class TaskScheduler {
public:
    using T = Tasks<TIME, ELEMENT>;

    void add(TIME time_left, ELEMENT element) {

        typename std::list<T>::iterator insert_it = this->tasks_queue.begin();
        while(insert_it != this->tasks_queue.end() && insert_it->time_left < time_left) {
            ++insert_it;
        }

        if(insert_it != this->tasks_queue.end() && insert_it->time_left == time_left) {
            insert_it->task_elements.push_back(element);
        } else {
            this->tasks_queue.insert(insert_it, T(time_left, element));
        }
    }

    void pop() {
        if(!this->tasks_queue.empty()) {
            this->tasks_queue.pop_front();
        }
    }

    std::list<ELEMENT> next() const {
        if(this->tasks_queue.empty()) {
            return std::list<ELEMENT>();
        }

        return this->tasks_queue.front().task_elements;
    }

    TIME time_advance() const {
        if(this->tasks_queue.empty()) {
            return TIME("inf");
        }
        return this->tasks_queue.front().time_left;
    }

    void update(TIME elapsed_time) {
        typename std::list<T>::iterator it;
        for(it = this->tasks_queue.begin(); it != this->tasks_queue.end(); ++it) {
            it->time_left -= elapsed_time;
        }
    }

    void advance() {
        TIME time_to_advance = this->time_advance();
        this->pop();
        this->update(time_to_advance);
    }

private:
    typename std::list<T> tasks_queue;
};

template<typename TIME, typename ELEMENT>
std::ostream& operator<<(std::ostream& os, const TaskScheduler<TIME, ELEMENT>& tasks) {
    os << "time left: " << tasks.time_advance();
    return os;
}

#endif //PMGBP_PDEVS_TASKSCHEDULER_HPP
