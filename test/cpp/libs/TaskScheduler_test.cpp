//
// Created by lao on 17/09/17.
//

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "../../../libs/TaskScheduler.hpp"
#include <NDTime.hpp>

BOOST_AUTO_TEST_SUITE( libs_task_cheduler )

    BOOST_AUTO_TEST_CASE( empty_scheduler_infinit_advance_time ) {

        TaskScheduler<NDTime, int> empty_scheduler;
        BOOST_CHECK_EQUAL(empty_scheduler.time_advance(), NDTime("inf"));
        BOOST_CHECK_EQUAL(empty_scheduler.time_advance(), std::string("inf"));

        empty_scheduler.add(NDTime(1), 2);
        empty_scheduler.pop();

        BOOST_CHECK_EQUAL(empty_scheduler.time_advance(), NDTime("inf"));
    }

    BOOST_AUTO_TEST_CASE( not_empty_inifinity_advance_time_return_slower_time_left ) {

        TaskScheduler<NDTime, int> scheduler;
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime("inf"));
        scheduler.add(NDTime(3), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(3));
        scheduler.add(NDTime(2), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(2));
        scheduler.add(NDTime(1), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(1));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(2));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(3));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime("inf"));
        scheduler.add(NDTime(1), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(1));
        scheduler.add(NDTime(1), 3);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(1));
        scheduler.add(NDTime(1), 4);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime(1));
    }

    BOOST_AUTO_TEST_CASE( next_return_correct_tasks ) {
        TaskScheduler<NDTime, int> scheduler;
        BOOST_CHECK_EQUAL(scheduler.next(), std::list<int>());
        scheduler.add(NDTime(3), 2);
        BOOST_CHECK_EQUAL(scheduler.next(), std::list<int>({2}));
        scheduler.add(NDTime(3), 2);
        BOOST_CHECK_EQUAL(scheduler.next(), std::list<int>({2,2}));
    }

BOOST_AUTO_TEST_SUITE_END()
