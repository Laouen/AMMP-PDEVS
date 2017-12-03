//
// Created by lao on 17/09/17.
//

#define BOOST_TEST_DYN_LINK
#include <NDTime.hpp>
#include <boost/test/unit_test.hpp>
#include "../../../libs/TaskScheduler.hpp"

BOOST_AUTO_TEST_SUITE( libs_task_cheduler )

    BOOST_AUTO_TEST_CASE( empty_scheduler_infinit_advance_time ) {

        TaskScheduler<NDTime, int> empty_scheduler;
        BOOST_CHECK_EQUAL(empty_scheduler.time_advance(), NDTime::infinity());

        empty_scheduler.add(NDTime({1}), 2);
        empty_scheduler.pop();

        BOOST_CHECK_EQUAL(empty_scheduler.time_advance(), NDTime::infinity());
    }

    BOOST_AUTO_TEST_CASE( not_empty_inifinity_advance_time_return_slower_time_left ) {

        TaskScheduler<NDTime, int> scheduler;
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime::infinity());
        scheduler.add(NDTime({3}), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({3}));
        scheduler.add(NDTime({2}), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({2}));
        scheduler.add(NDTime({1}), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({2}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({3}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime::infinity());
        scheduler.add(NDTime({1}), 2);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.add(NDTime({1}), 3);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.add(NDTime({1}), 4);
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
    }

    BOOST_AUTO_TEST_CASE( next_returns_correct_tasks ) {
        TaskScheduler<NDTime, int> scheduler;
        BOOST_CHECK(scheduler.next() == std::list<int>());
        scheduler.add(NDTime({3}), 2);
        BOOST_CHECK(scheduler.next() == std::list<int>({2}));
        scheduler.add(NDTime({3}), 2);
        BOOST_CHECK(scheduler.next() == std::list<int>({2,2}));
    }

    BOOST_AUTO_TEST_CASE( add_next_methods__elements_correctly_grouped_by_time_left ) {
        TaskScheduler<NDTime, int> scheduler;
        scheduler.add(NDTime({1}), 2);
        scheduler.add(NDTime({1}), 3);
        std::list<int> next = scheduler.next();
        BOOST_CHECK(std::find(next.begin(), next.end(), 2) != next.end());
        BOOST_CHECK(std::find(next.begin(), next.end(), 3) != next.end());
        BOOST_CHECK(next.size() == 2);

        scheduler.add(NDTime({2}), 4);
        scheduler.add(NDTime({2}), 5);
        scheduler.add(NDTime({1}), 1);

        next = scheduler.next();
        BOOST_CHECK(std::find(next.begin(), next.end(), 1) != next.end());
        BOOST_CHECK(std::find(next.begin(), next.end(), 2) != next.end());
        BOOST_CHECK(std::find(next.begin(), next.end(), 3) != next.end());
        BOOST_CHECK(next.size() == 3);

        scheduler.pop();

        next = scheduler.next();
        BOOST_CHECK(std::find(next.begin(), next.end(), 4) != next.end());
        BOOST_CHECK(std::find(next.begin(), next.end(), 5) != next.end());
        BOOST_CHECK(next.size() == 2);
    }

    BOOST_AUTO_TEST_CASE( update_pop__tests ) {
        TaskScheduler<NDTime, int> scheduler;

        scheduler.add(NDTime({1}), 2);
        scheduler.add(NDTime({2}), 2);
        scheduler.add(NDTime({3}), 2);

        scheduler.update(NDTime({0,30}));
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({0,30}));

        scheduler.update(NDTime({0,30}));
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({0}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({2}));
        scheduler.pop();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime::infinity());
    }

    BOOST_AUTO_TEST_CASE( advance__tests ) {
        TaskScheduler<NDTime, int> scheduler;

        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime::infinity());

        scheduler.add(NDTime({1}), 2);
        scheduler.add(NDTime({1}), 3);
        scheduler.add(NDTime({2}), 2);
        scheduler.add(NDTime({2}), 3);
        scheduler.add(NDTime({3}), 2);
        scheduler.add(NDTime({3}), 3);

        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.advance();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.advance();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime({1}));
        scheduler.advance();
        BOOST_CHECK_EQUAL(scheduler.time_advance(), NDTime::infinity());
    }

BOOST_AUTO_TEST_SUITE_END()
