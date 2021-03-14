//
// Created by lao on 21/09/17.
//

#ifndef PMGBP_PDEVS_TUPLE_OPERATORS_HPP
#define PMGBP_PDEVS_TUPLE_OPERATORS_HPP

#include <tuple>
#include <cadmium/modeling/message_bag.hpp>
#include <cassert>
#include <iostream>
#include <boost/core/demangle.hpp>

namespace pmgbp {
namespace tuple {

class Exception : public std::exception {
private:
    std::string message_;

public:
    Exception() { message_ = ""; };

    Exception(const std::string &message) {
        message_ = message;
    };

    virtual const char *what() const throw() {
        return message_.c_str();
    }
};


/*****************************************/
/*************** PRINT *******************/
/*****************************************/

template<int index, typename... Ts>
struct print_tuple {
    std::ostream& operator()(std::ostream& os, const std::tuple<Ts...> &t) {

        using port_type = typename std::tuple_element<index, std::tuple<Ts...>>::type::port;
        std::string port_name = boost::core::demangle(std::type_index(typeid(port_type)).name());
        port_name = port_name.substr(port_name.find_last_of(':') + 1);

        os << "port_name: " << port_name << " messages:[";
        for (const auto& m : std::get<index>(t).messages) {
            os << " " << m;
        }
        os << "] ";

        return print_tuple<index - 1, Ts...>{}(os, t);
    }
};

template<typename... Ts>
struct print_tuple<0, Ts...> {
    std::ostream& operator()(std::ostream& os, const std::tuple<Ts...> &t) {

        using port_type = typename std::tuple_element<0, std::tuple<Ts...>>::type::port;
        std::string port_name = boost::core::demangle(std::type_index(typeid(port_type)).name());
        port_name = port_name.substr(port_name.find_last_of(':') + 1);

        os << "port_name: " << port_name << " messages:[";
        for (const auto& m : std::get<0>(t).messages) {
            os << " " << m;
        }
        os << "] ";

        return os;
    }
};

template<typename... Ts>
std::ostream& print(std::ostream& os, const std::tuple<Ts...> &t) {
    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    return print_tuple<size - 1, Ts...>{}(os, t);
}

/*****************************************/
/**************** GET ********************/
/*****************************************/

template<int index, typename T, typename... Ts>
struct get_tuple {
    cadmium::bag<T>& operator()(std::tuple<Ts...> &t, int position) {

        if (position == index) {
            return std::get<index>(t).messages;
        }

        if (position < index) {
            return get_tuple<index - 1, T, Ts...>{}(t, position);
        }

        throw Exception("Index value out of range (position > index): "
                        + std::to_string(position)
                        + " > "
                        + std::to_string(index));
    }
};

template<typename T, typename... Ts>
struct get_tuple<0, T, Ts...> {
    cadmium::bag<T>& operator()(std::tuple<Ts...> &t, int position) {
        if (position == 0) {
            return std::get<0>(t).messages;
        }

        throw Exception("Index value out of range (position > index): " + std::to_string(position) + " > 0");
    }
};

template<typename T, typename... Ts>
cadmium::bag<T>& get(std::tuple<Ts...> &t, int position) {

    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    if (position < size) {
        return get_tuple<size - 1, T, Ts...>{}(t, position);
    }

    throw Exception("Index value out of range (position > bag_size: "
                                 + std::to_string(position)
                                 + " > "
                                 + std::to_string(size - 1));
}

template<int index, typename T, typename... Ts>
struct cget_tuple {
    const cadmium::bag<T>& operator()(const std::tuple<Ts...> &t, int position) {

        if (position == index) {
            return std::get<index>(t).messages;
        }

        if (position < index) {
            return cget_tuple<index - 1, T, Ts...>{}(t, position);
        }

        throw Exception("Index value out of range (position > index): "
                        + std::to_string(position)
                        + " > "
                        + std::to_string(index));
    }
};

template<typename T, typename... Ts>
struct cget_tuple<0, T, Ts...> {
    const cadmium::bag<T>& operator()(const std::tuple<Ts...> &t, int position) {
        if (position == 0) {
            return std::get<0>(t).messages;
        }

        throw Exception("Index value out of range (position > index): " + std::to_string(position) + " > 0");
    }
};

template<typename T, typename... Ts>
const cadmium::bag<T>& cget(const std::tuple<Ts...> &t, int position) {

    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    if (position < size) {
        return cget_tuple<size - 1, T, Ts...>{}(t, position);
    }

    throw Exception("Index value out of range: "
                    + std::to_string(position)
                    + " > "
                    + std::to_string(size - 1));
}
/*****************************************/
/**************** EMPTY ******************/
/*****************************************/


template<int index, typename... Ts>
struct empty_tuple {
    bool operator()(const std::tuple<Ts...> &t) {

        if (std::get<index>(t).messages.empty()) {
            return empty_tuple<index - 1, Ts...>{}(t);
        }

        return false;
    }
};

template<typename... Ts>
struct empty_tuple<0, Ts...> {
    bool operator()(const std::tuple<Ts...> &t) {
        return std::get<0>(t).messages.empty();
    }
};

template<typename... Ts>
bool empty(const std::tuple<Ts...> &t) {
    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    return empty_tuple<size - 1, Ts...>{}(t);
}


/*****************************************/
/**************** MAP ********************/
/*****************************************/

template<int index, typename T, typename... Ts>
struct map_tuple {
    void operator()(std::tuple<Ts...> &t, void (*modifier)(cadmium::bag<T> &)) {
        modifier(std::get<index>(t).messages);
        map_tuple<index - 1, T, Ts...>{}(t, modifier);
    }
};

template<typename T, typename... Ts>
struct map_tuple<0, T, Ts...> {
    void operator()(std::tuple<Ts...> &t, void (*modifier)(cadmium::bag<T> &)) {
        modifier(std::get<0>(t).messages);
    }
};

template<typename T, typename... Ts>
void map(std::tuple<Ts...> &t, void (*modifier)(cadmium::bag<T> &)) {
    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    map_tuple<size - 1, T, Ts...>{}(t, modifier);
}


/*****************************************/
/**************** MERGE ******************/
/*****************************************/


template <int index, typename... Ts>
struct merge_tuple {
    void operator()(std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
        std::get<index>(l).messages.insert(std::get<index>(l).messages.end(),
                                           std::get<index>(r).messages.begin(),
                                           std::get<index>(r).messages.end());

        merge_tuple<index - 1, Ts...>{}(l, r);
    }
};

template <typename... Ts>
struct merge_tuple<0, Ts...> {
    void operator()(std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
        std::get<0>(l).messages.insert(std::get<0>(l).messages.end(),
                                           std::get<0>(r).messages.begin(),
                                           std::get<0>(r).messages.end());
    }
};

template <typename... Ts>
void merge(std::tuple<Ts...> &l, const std::tuple<Ts...> &r) {
    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    merge_tuple<size - 1, Ts...>{}(l, r);
}

template <int index, typename... Ts>
struct equals_tuple {
    bool operator()(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
        bool current_equal = std::get<index>(l).messages == std::get<index>(r).messages;
        return current_equal && equals_tuple<index - 1, Ts...>{}(l, r);
    }
};

template <typename... Ts>
struct equals_tuple<0, Ts...> {
    bool operator()(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
        return std::get<0>(l).messages == std::get<0>(r).messages;
    }
};


template <typename... Ts>
bool equals(const std::tuple<Ts...>& l, const std::tuple<Ts...>& r) {
    const auto size = std::tuple_size<std::tuple<Ts...>>::value;
    return equals_tuple<size - 1, Ts...>{}(l, r);
}

}
}

#endif //PMGBP_PDEVS_TUPLE_OPERATORS_HPP
