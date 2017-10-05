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

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, 0).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, 1).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, 2).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, 3).size(), 0);

        for (int i = 1; i <= 10; i++) {
            for (int j = 0; j < 4; j++) {
                pmgbp::tuple::get<int>(bags, j).emplace_back(1);
                BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, j).size(), i);
            }
        }
    }

    BOOST_AUTO_TEST_CASE( get_message_bag_wrong_index ) {

        typename make_message_bags<output_ports>::type bags;

        for (int i = 4; i <= 10; i++) {
            BOOST_CHECK_THROW(pmgbp::tuple::get<int>(bags, i), pmgbp::tuple::Exception);
        }
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( empty_tuple_tests )

    BOOST_AUTO_TEST_CASE( all_cambination_of_non_empty_bags_return_false ) {

        typename make_message_bags<output_ports>::type bags;

        // All combinations of one non empty bag
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                pmgbp::tuple::get<int>(bags, k).clear();
            }
            pmgbp::tuple::get<int>(bags, i).emplace_back(i);
            BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), false);
        }

        // All combinations of two non empty bags
        for (int i = 0; i < 4; i++) {
            for (int j = i+1; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    pmgbp::tuple::get<int>(bags, k).clear();
                }
                pmgbp::tuple::get<int>(bags, i).emplace_back(i);
                pmgbp::tuple::get<int>(bags, j).emplace_back(j);
                BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), false);
            }
        }

        // All combination of three non empty bags
        for (int i = 0; i < 4; i++) {
            for (int j = i+1; j < 4; j++) {
                for (int x = j+1; x < 4; x++) {

                    for (int k = 0; k < 4; k++) {
                        pmgbp::tuple::get<int>(bags, k).clear();
                    }
                    pmgbp::tuple::get<int>(bags, i).emplace_back(i);
                    pmgbp::tuple::get<int>(bags, j).emplace_back(j);
                    pmgbp::tuple::get<int>(bags, x).emplace_back(x);
                    BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), false);
                }
            }
        }

        for (int k = 0; k < 4; k++) {
            pmgbp::tuple::get<int>(bags, k).clear();
            pmgbp::tuple::get<int>(bags, k).emplace_back(k);
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), false);
    }


    BOOST_AUTO_TEST_CASE( empty_returns_true ) {
        typename make_message_bags<output_ports>::type bags;

        BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), true);

        for (int k = 0; k < 4; k++) {
            pmgbp::tuple::get<int>(bags, k).emplace_back(k);
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), false);

        for (int k = 0; k < 4; k++) {
            pmgbp::tuple::get<int>(bags, k).clear();
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::empty(bags), true);
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( map_testing )

    BOOST_AUTO_TEST_CASE( map_testing_pow_two_function ) {

        typename make_message_bags<output_ports>::type bags;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags, 0).emplace_back(k);
        }

        pmgbp::tuple::map(bags, pow2_test);

        for(int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags, 0)[k], std::pow(k, 2));
        }

    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( merge_testing )

    BOOST_AUTO_TEST_CASE( merge_testing_empty_empty ) {

        typename make_message_bags<output_ports>::type bags_left;
        typename make_message_bags<output_ports>::type bags_right;

        BOOST_CHECK(pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags_right));
        pmgbp::tuple::merge(bags_left, bags_right);
        BOOST_CHECK(pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags_right));
    }

    BOOST_AUTO_TEST_CASE( merge_testing_empty_not_empty ) {

        typename make_message_bags<output_ports>::type bags_left;
        typename make_message_bags<output_ports>::type bags_right;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags_right, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags_right, 2).emplace_back(k);
        }

        BOOST_CHECK(pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags_right));
        pmgbp::tuple::merge(bags_left, bags_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags_right));

        for (int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 2)[k], k);
        }

        BOOST_CHECK(pmgbp::tuple::get<int>(bags_left, 1).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags_left, 3).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags_right, 1).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags_right, 3).empty());

        typename make_message_bags<output_ports>::type bags2_left;
        typename make_message_bags<output_ports>::type bags2_right;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags2_right, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 2).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 3).emplace_back(k);
        }

        BOOST_CHECK(pmgbp::tuple::empty(bags2_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags2_right));
        pmgbp::tuple::merge(bags2_left, bags2_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags2_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags2_right));

        for (int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3)[k], k);

            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 1)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 3)[k], k);
        }

        typename make_message_bags<output_ports>::type bags3_left;
        typename make_message_bags<output_ports>::type bags3_right;

        pmgbp::tuple::get<int>(bags3_right, 0).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 2).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 3).emplace_back(1);

        BOOST_CHECK(pmgbp::tuple::empty(bags3_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_right));
        pmgbp::tuple::merge(bags3_left, bags3_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_right));

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 0)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 2)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 3)[0], 1);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 0)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 2)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 3)[0], 1);

        typename make_message_bags<output_ports>::type bags4_left;
        typename make_message_bags<output_ports>::type bags4_right;

        pmgbp::tuple::get<int>(bags4_right, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags4_right, 3).emplace_back(1);

        BOOST_CHECK(pmgbp::tuple::empty(bags4_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags4_right));
        pmgbp::tuple::merge(bags4_left, bags4_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags4_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags4_right));

        BOOST_CHECK(pmgbp::tuple::get<int>(bags4_right, 0).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags4_right, 2).empty());
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_right, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_right, 3)[0], 1);

        BOOST_CHECK(pmgbp::tuple::get<int>(bags4_left, 0).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags4_left, 2).empty());
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 3)[0], 1);
    }

    BOOST_AUTO_TEST_CASE( merge_testing_not_empty_empty ) {

        typename make_message_bags<output_ports>::type bags_left;
        typename make_message_bags<output_ports>::type bags_right;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags_left, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags_left, 2).emplace_back(k);
        }

        BOOST_CHECK(!pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags_right));
        pmgbp::tuple::merge(bags_left, bags_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags_right));

        for (int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 2)[k], k);
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 0).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 2).size(), 10);
        BOOST_CHECK(pmgbp::tuple::get<int>(bags_left, 1).empty());
        BOOST_CHECK(pmgbp::tuple::get<int>(bags_left, 3).empty());

        typename make_message_bags<output_ports>::type bags2_left;
        typename make_message_bags<output_ports>::type bags2_right;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags2_left, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 2).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 3).emplace_back(k);
        }

        BOOST_CHECK(!pmgbp::tuple::empty(bags2_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags2_right));
        pmgbp::tuple::merge(bags2_left, bags2_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags2_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags2_right));

        for (int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3)[k], k);
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3).size(), 10);

        typename make_message_bags<output_ports>::type bags3_left;
        typename make_message_bags<output_ports>::type bags3_right;

        pmgbp::tuple::get<int>(bags3_right, 0).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 2).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 3).emplace_back(1);

        BOOST_CHECK(pmgbp::tuple::empty(bags3_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_right));
        pmgbp::tuple::merge(bags3_left, bags3_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_left));
        BOOST_CHECK(!pmgbp::tuple::empty(bags3_right));

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 0)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 2)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 3)[0], 1);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 0).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 1).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 2).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 3).size(), 1);

        typename make_message_bags<output_ports>::type bags4_left;
        typename make_message_bags<output_ports>::type bags4_right;

        pmgbp::tuple::get<int>(bags4_left, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags4_left, 3).emplace_back(1);

        BOOST_CHECK(!pmgbp::tuple::empty(bags4_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags4_right));
        pmgbp::tuple::merge(bags4_left, bags4_right);
        BOOST_CHECK(!pmgbp::tuple::empty(bags4_left));
        BOOST_CHECK(pmgbp::tuple::empty(bags4_right));

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 3)[0], 1);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 0).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 1).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 2).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags4_left, 3).size(), 1);
    }

    BOOST_AUTO_TEST_CASE( merge_not_empty_not_empty ) {

        typename make_message_bags<output_ports>::type bags_left;
        typename make_message_bags<output_ports>::type bags_right;

        pmgbp::tuple::get<int>(bags_left, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags_left, 3).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right, 1).emplace_back(2);
        pmgbp::tuple::get<int>(bags_right, 3).emplace_back(2);

        pmgbp::tuple::merge(bags_left, bags_right);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 1)[1], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 3)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 3)[1], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 1)[0], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 3)[0], 2);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 1).size(), 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 3).size(), 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 0).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_left, 2).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 1).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 3).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 0).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags_right, 2).size(), 0);

        typename make_message_bags<output_ports>::type bags2_left;
        typename make_message_bags<output_ports>::type bags2_right;

        for (int k = 0; k < 10; k++) {
            pmgbp::tuple::get<int>(bags2_left, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 2).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_left, 3).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 0).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 2).emplace_back(k);
            pmgbp::tuple::get<int>(bags2_right, 3).emplace_back(k);
        }

        pmgbp::tuple::merge(bags2_left, bags2_right);

        for (int k = 0; k < 10; k++) {
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0)[k+10], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1)[k+10], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2)[k+10], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3)[k+10], k);

            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 0)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 1)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 2)[k], k);
            BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 3)[k], k);
        }

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 0).size(), 20);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 1).size(), 20);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 2).size(), 20);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_left, 3).size(), 20);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 0).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 1).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 2).size(), 10);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags2_right, 3).size(), 10);

        typename make_message_bags<output_ports>::type bags3_left;
        typename make_message_bags<output_ports>::type bags3_right;

        pmgbp::tuple::get<int>(bags3_left, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_left, 3).emplace_back(1);
        pmgbp::tuple::get<int>(bags3_right, 0).emplace_back(2);
        pmgbp::tuple::get<int>(bags3_right, 2).emplace_back(2);

        pmgbp::tuple::merge(bags3_left, bags3_right);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 0)[0], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 1)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 2)[0], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 3)[0], 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 0)[0], 2);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 2)[0], 2);

        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 0).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 1).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 2).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_left, 3).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 0).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 2).size(), 1);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 1).size(), 0);
        BOOST_CHECK_EQUAL(pmgbp::tuple::get<int>(bags3_right, 3).size(), 0);

    }



    BOOST_AUTO_TEST_CASE( equals_tests ) {
        typename make_message_bags<output_ports>::type bags_left;
        typename make_message_bags<output_ports>::type bags_right;

        pmgbp::tuple::get<int>(bags_left, 0).emplace_back(1);
        pmgbp::tuple::get<int>(bags_left, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags_left, 2).emplace_back(1);
        pmgbp::tuple::get<int>(bags_left, 3).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right, 0).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right, 2).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right, 3).emplace_back(1);

        BOOST_CHECK(pmgbp::tuple::equals(bags_left, bags_right));

        typename make_message_bags<output_ports>::type bags_left1;
        typename make_message_bags<output_ports>::type bags_right1;

        pmgbp::tuple::get<int>(bags_left1, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags_left1, 3).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right1, 1).emplace_back(1);
        pmgbp::tuple::get<int>(bags_right1, 3).emplace_back(1);

        BOOST_CHECK(pmgbp::tuple::equals(bags_left1, bags_right1));

        typename make_message_bags<output_ports>::type bags_left2;
        typename make_message_bags<output_ports>::type bags_right2;

        for (int k = 0; k < 10; k++) {

            pmgbp::tuple::get<int>(bags_left2, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags_left2, 3).emplace_back(k);
            pmgbp::tuple::get<int>(bags_right2, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags_right2, 3).emplace_back(k);
        }

        BOOST_CHECK(pmgbp::tuple::equals(bags_left2, bags_right2));

        typename make_message_bags<output_ports>::type bags_left3;
        typename make_message_bags<output_ports>::type bags_right3;

        for (int k = 0; k < 10; k++) {

            pmgbp::tuple::get<int>(bags_left3, 1).emplace_back(k);
            pmgbp::tuple::get<int>(bags_left3, 3).emplace_back(k);
            pmgbp::tuple::get<int>(bags_right3, 0).emplace_back(k);
        }

        BOOST_CHECK(!pmgbp::tuple::equals(bags_left3, bags_right3));
    }

BOOST_AUTO_TEST_SUITE_END()