//
// Created by lao on 21/09/17.
//

#ifndef PMGBP_PDEVS_TUPLE_OPERATORS_HPP
#define PMGBP_PDEVS_TUPLE_OPERATORS_HPP

#include <tuple>
#include <cadmium/modeling/message_bag.hpp>
#include <cassert>

namespace pmgbp {

    class TupleOperatorException : public std::exception {
    private:
        std::string message_;

    public:
        TupleOperatorException() { message_ = ""; };

        TupleOperatorException(const std::string& message) {
            message_ = message;
        };

        virtual const char* what() const throw() {
            return message_.c_str();
        }
    };

    template<int index, typename T, typename... Ts>
    struct get_tuple {
        cadmium::bag<T>& operator() (std::tuple<Ts...>& t, int position) {

            if (position == index) {
                return std::get<index>(t).messages;
            }

            if (position < index) {
                return get_tuple<index - 1, T, Ts...>{}(t, position);
            }
        }
    };

    template<typename T, typename... Ts>
    struct get_tuple<0, T, Ts...> {
        cadmium::bag<T>& operator() (std::tuple<Ts...>& t, int position) {
            if (position == 0) {
                return std::get<0>(t).messages;
            }
        }
    };

    template<typename T, typename... Ts>
    cadmium::bag<T>& get(std::tuple<Ts...>& t, int position) {

        const auto size = std::tuple_size<std::tuple<Ts...>>::value;
        if (position < size) {
            return get_tuple<size - 1, T, Ts...>{}(t, position);
        }

        throw TupleOperatorException("Index value out of range: "
                                     + std::to_string(position)
                                     + " > "
                                     + std::to_string(size - 1));
    }

    template<int index, typename T, typename... Ts>
    struct empty_tuple {
        bool operator() (std::tuple<Ts...>& t) {
            if (std::get<index>(t).messages.empty()) {
                return empty_tuple<index - 1, T, Ts...>{}(t);
            }
            return false;
        }
    };

    template<typename T, typename... Ts>
    struct empty_tuple<0, T, Ts...> {
        bool operator() (std::tuple<Ts...>& t) {
            return std::get<0>(t).messages.empty();
        }
    };

    template<typename T, typename... Ts>
    bool empty(std::tuple<Ts...>& t) {
        const auto size = std::tuple_size<std::tuple<Ts...>>::value;
        return empty_tuple<size - 1, T, Ts...>{}(t);
    }

    template<int index, typename T, typename... Ts>
    struct map_tuple {
        void operator() (std::tuple<Ts...>& t, void (*modifier)(cadmium::bag<T>&)) {
            modifier(std::get<index>(t).messages);
            map_tuple<index - 1, T, Ts...>{}(t, modifier);
        }
    };

    template<typename T, typename... Ts>
    struct map_tuple<0, T, Ts...> {
        void operator() (std::tuple<Ts...>& t, void (*modifier)(cadmium::bag<T>&)) {
            modifier(std::get<0>(t).messages);
        }
    };

    template<typename T, typename... Ts>
    void map(std::tuple<Ts...>& t, void (*modifier)(cadmium::bag<T>&)) {
        const auto size = std::tuple_size<std::tuple<Ts...>>::value;
        map_tuple<size - 1, T, Ts...>{}(t, modifier);
    }
}

#endif //PMGBP_PDEVS_TUPLE_OPERATORS_HPP
