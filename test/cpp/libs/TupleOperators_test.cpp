//
// Created by lao on 21/09/17.
//

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cadmium/modeling/ports.hpp>
#include <cadmium/modeling/message_bag.hpp>
#include "../../../libs/TupleOperators.hpp"

using namespace std;
using namespace cadmium;

void pow2_test(cadmium::bag<int>& a) {

    for (cadmium::bag<int>::iterator it = a.begin(); it != a.end(); it++) {
        *it = std::pow(*it, 2);
    }
}

struct ports {
    struct one : public out_port<int> {};
    struct two : public out_port<int> {};
    struct three : public out_port<int> {};
    struct four : public out_port<int> {};
};

using output_ports=tuple<typename ports::one, typename ports::two, typename ports::three, typename ports::four>;

BOOST_AUTO_TEST_SUITE( get_operator )

    BOOST_AUTO_TEST_CASE( get_message_bag_by_index ) {

        typename make_message_bags<output_ports>::type bags;

        BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, 1).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, 0).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, 2).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, 3).size(), 0);

        for (int i = 1; i <= 10; i++) {
            for (int j = 0; j < 4; j++) {
                pmgbp::get<int>(bags, j).emplace_back(1);
                BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, j).size(), i);
            }
        }
    }

    BOOST_AUTO_TEST_CASE( get_message_bag_wrong_index ) {

        typename make_message_bags<output_ports>::type bags;

        for (int i = 4; i <= 10; i++) {
            BOOST_CHECK_THROW(pmgbp::get<int>(bags, i), pmgbp::TupleOperatorException);
        }
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( empty_tuple_tests )

    BOOST_AUTO_TEST_CASE( all_cambination_of_non_empty_bags_return_false ) {

        typename make_message_bags<output_ports>::type bags;

        // All combinations of one non empty bag
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                pmgbp::get<int>(bags, k).clear();
            }
            pmgbp::get<int>(bags, i).emplace_back(i);
            BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), false);
        }

        // All combinations of two non empty bags
        for (int i = 0; i < 4; i++) {
            for (int j = i+1; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    pmgbp::get<int>(bags, k).clear();
                }
                pmgbp::get<int>(bags, i).emplace_back(i);
                pmgbp::get<int>(bags, j).emplace_back(j);
                BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), false);
            }
        }

        // All combination of three non empty bags
        for (int i = 0; i < 4; i++) {
            for (int j = i+1; j < 4; j++) {
                for (int x = j+1; x < 4; x++) {

                    for (int k = 0; k < 4; k++) {
                        pmgbp::get<int>(bags, k).clear();
                    }
                    pmgbp::get<int>(bags, i).emplace_back(i);
                    pmgbp::get<int>(bags, j).emplace_back(j);
                    pmgbp::get<int>(bags, x).emplace_back(x);
                    BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), false);
                }
            }
        }

        for (int k = 0; k < 4; k++) {
            pmgbp::get<int>(bags, k).clear();
            pmgbp::get<int>(bags, k).emplace_back(k);
        }

        BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), false);
    }


    BOOST_AUTO_TEST_CASE( empty_returns_true ) {
        typename make_message_bags<output_ports>::type bags;

        BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), true);

        for (int k = 0; k < 4; k++) {
            pmgbp::get<int>(bags, k).emplace_back(k);
        }

        BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), false);

        for (int k = 0; k < 4; k++) {
            pmgbp::get<int>(bags, k).clear();
        }

        BOOST_CHECK_EQUAL(pmgbp::empty<int>(bags), true);
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( merge_testing )

    BOOST_AUTO_TEST_CASE( merge_testing_pow_two_function ) {

        typename make_message_bags<output_ports>::type bags;

        for (int k = 0; k < 10; k++) {
            pmgbp::get<int>(bags, 0).emplace_back(k);
        }

        pmgbp::map<int>(bags, pow2_test);

        for(int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::get<int>(bags, 0)[k], std::pow(k, 2));
        }

    }

BOOST_AUTO_TEST_SUITE_END()