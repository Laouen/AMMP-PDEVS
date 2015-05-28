#ifndef BRITIME_H
#define BRITIME_H

#include <iostream>
#include <string>
#include <assert.h>
#include <boost/rational.hpp>

using namespace std;

class BRITime {

public:
    boost::rational<int> _value;
    bool                 _is_inf;

    BRITime() : _is_inf(false) {}
    BRITime(int n) : _value(n), _is_inf(false) {}
    BRITime(int n, int d) : _value(n,d), _is_inf(false) {}

    BRITime& operator=(const BRITime& o) noexcept { 
        this->_value = o._value;
        this->_is_inf = o._is_inf;
	return *this;
    }

    BRITime& operator+=(const BRITime& o) noexcept {
        
        this->_value += o._value;
        if (o._is_inf) this->_is_inf = true;
        return *this;
    }

    BRITime& operator-=(const BRITime& o) noexcept { 
        assert(this->_is_inf || !o._is_inf);

        this->_value -= o._value;
        return *this;
    }

    string toString() {
        return to_string(_value.numerator()) + "/" + to_string(_value.denominator());  
    }

    static BRITime infinity() noexcept {
        BRITime f;
        f._is_inf=true;
        return f;
    }
};

inline BRITime operator+(BRITime lhs, const BRITime& rhs) noexcept {
    lhs += rhs;
    return lhs;
}

inline BRITime operator-(BRITime lhs, const BRITime& rhs) noexcept {
    lhs -= rhs;
    return lhs;
}

inline BRITime operator/(BRITime lhs, const BRITime& rhs) noexcept {
    assert(rhs._is_inf || (rhs._value != boost::rational<int>()));

    if (!lhs._is_inf && !rhs._is_inf) {
        lhs._value /= rhs._value;
    } else if (!lhs._is_inf) {
        lhs = BRITime(0,1);
    }
    return lhs;
}

inline bool operator==(const BRITime& lhs, const BRITime& rhs) noexcept {

    if (lhs._is_inf && rhs._is_inf) return true;
    else if (lhs._is_inf || rhs._is_inf) return false;
    return (lhs._value == rhs._value);
}

inline bool operator!=(const BRITime& lhs, const BRITime& rhs) noexcept {return !operator==(lhs,rhs);}

inline bool operator< (const BRITime& lhs, const BRITime& rhs) noexcept {
    if (lhs._is_inf) return false;
    else if (rhs._is_inf) return true;
    return (lhs._value < rhs._value);
}

inline bool operator> (const BRITime& lhs, const BRITime& rhs) noexcept {return  operator< (rhs,lhs);}

inline bool operator<=(const BRITime& lhs, const BRITime& rhs) noexcept {return !operator> (lhs,rhs);}

inline bool operator>=(const BRITime& lhs, const BRITime& rhs) noexcept {return !operator< (lhs,rhs);}

inline std::ostream& operator<<(std::ostream& os, const BRITime& t) noexcept {
    
    if (t._is_inf) os << "inf";
    else os << t._value; 
    return os;
}

inline std::istream& operator>> (std::istream& is, BRITime& rhs) noexcept {
    string a;
    int n,d;
    is >> a;
    if (a == "inf") rhs = BRITime::infinity();
    else {
        n = std::stoi(a.substr(0, a.find_last_of("/")));
        d = std::stoi(a.substr(a.find_last_of("/")+1));
        rhs = BRITime(n,d);
    }
    return is;
}


// Specialize numeric_limits
namespace std {
template<>
class numeric_limits<BRITime>{
public:
    static constexpr bool is_specialized = true;
    static BRITime min() noexcept { return BRITime{numeric_limits<int>::min(), 1}; }
    static BRITime max() noexcept { return BRITime{numeric_limits<int>::max(), 1}; }
    static BRITime lowest() noexcept { return BRITime{numeric_limits<int>::lowest(), 1}; }

    static constexpr int  digits = numeric_limits<int>::digits;
    static constexpr int  digits10 = numeric_limits<int>::digits10;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = true;
    static constexpr int radix = 2; // trash_value
    static BRITime epsilon() noexcept { return BRITime{1,1} - BRITime{numeric_limits<int>::max(), numeric_limits<int>::max() - 1}; }
    static BRITime round_error() noexcept { return BRITime(0); }

    static constexpr int  min_exponent = numeric_limits<int>::min(); // trash_value
    static constexpr int  min_exponent10 = min_exponent/radix; // trash_value
    static constexpr int  max_exponent = numeric_limits<int>::max(); // trash_value
    static constexpr int  max_exponent10 = max_exponent/radix; // trash_value

    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = denorm_indeterminate;
    static constexpr bool has_denorm_loss = false;
    static BRITime infinity() noexcept { return BRITime::infinity(); }
//    static constexpr BRITime quiet_NaN() noexcept { return T(); }
//    static constexpr BRITime signaling_NaN() noexcept { return T(); }
//    static constexpr BRITime denorm_min() noexcept { return T(); }

    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;

    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
//    static constexpr float_round_style round_style = round_toward_zero;

};
}

#endif // BRITIME_H
