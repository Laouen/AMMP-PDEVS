/**
 * Copyright (c) 2013, Carleton University, Universite de Nice-Sophia Antipolis
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FDTIME_H
#define FDTIME_H

#include <ostream>
#include <exception>
#include <boost/math/common_factor.hpp>

namespace cdpp {
/**
 * @brief The FDTime class represents a time as a q * m * 2^exp, where q is a quotient that needs agreement before operate
 *
 */
class FDTime
{
public:
    mutable int _q_num;
    mutable int _q_denom;
    mutable int _m;
    mutable int _exp;

public:
    bool is_inf;

    FDTime() noexcept : _q_num{1}, _q_denom{1}, _m{0}, _exp{0}, is_inf{false}{}
    FDTime(const FDTime& rhs) noexcept :
        _q_num{rhs._q_num}, _q_denom{rhs._q_denom}, _m{rhs._m}, _exp{rhs._exp}, is_inf{rhs.is_inf}{}
    explicit FDTime(const int& rhs) noexcept :
        _q_num{1}, _q_denom{1}, _m{rhs}, _exp{0}, is_inf{false}{}
    explicit FDTime(const int& q_n, const int& q_d, const int& m, const int& exp) noexcept :
        _q_num{q_n}, _q_denom{q_d}, _m{m}, _exp{exp}, is_inf{false}{
        if (_q_num == 0) throw std::exception();
    }

    FDTime& operator=(const FDTime& arg) noexcept {
        this->_q_num = arg._q_num;
        this->_q_denom = arg._q_denom;
        this->_m = arg._m;
        this->_exp = arg._exp;
        this->is_inf = arg.is_inf;
        return *this;
    }

    FDTime& operator+=(const FDTime& arg) noexcept {
        if (is_inf || arg.is_inf) { is_inf=true; }
        else {
            if (_q_num != arg._q_num || _q_denom != arg._q_denom ) {
                // adjust q and exp;
                // - find common denom
                int denom = boost::math::lcm(_q_denom, arg._q_denom);
                // - find common num
                int num = boost::math::gcd(_q_num * (denom/_q_denom), arg._q_num * (denom/arg._q_denom));
                // - fix m
                _m = _m * (denom / _q_denom) * (_q_num / num); //common num should never be 0.
                arg._m = arg._m * (denom / arg._q_denom) * (arg._q_num / num); //common num should never be 0.
                //simplify factor
                int simplify = boost::math::gcd(num, denom);
                // - fix the qs
                _q_num = arg._q_num  = num/simplify;
                _q_denom = arg._q_denom = denom/simplify;
            }
            if(_exp != arg._exp) {
                // - fix exps
                if (_exp < arg._exp) {
                    arg._m<<=(_exp - arg._exp);
                    arg._exp = _exp;
                } else {
                    _m<<=(arg._exp - _exp);
                    _exp = arg._exp;
                }
            }
            _m +=  arg._m;
        }
        return *this; // return the result by reference
    }

    //time subtraction does internal timerep subtract when not inf, else inf
    //we dont throw when infinity substraction happens.
    FDTime& operator-=(const FDTime& arg) noexcept {
        if (is_inf || arg.is_inf) { is_inf=true; }
        else {
            if (_q_num != arg._q_num || _q_denom != arg._q_denom ) {
                // adjust q and exp;
                // - find common denom
                int denom = boost::math::lcm(_q_denom, arg._q_denom);
                // - find common num
                int num = boost::math::gcd(_q_num * (denom/_q_denom), arg._q_num * (denom/arg._q_denom));
                // - fix m
                _m = _m * (denom / _q_denom) * (_q_num / num); //common num should never be 0.
                arg._m = arg._m * (denom / arg._q_denom) * (arg._q_num / num); //common num should never be 0.
                //simplify factor
                int simplify = boost::math::gcd(num, denom);
                // - fix the qs
                _q_num = arg._q_num  = num/simplify;
                _q_denom = arg._q_denom = denom/simplify;
            }
            if(_exp != arg._exp) {
                // - fix exps
                if (_exp < arg._exp) {
                    arg._m<<=(_exp - arg._exp);
                    arg._exp = _exp;
                } else {
                    _m<<=(arg._exp - _exp);
                    _exp = arg._exp;
                }
            }
            _m -=  arg._m;
        }
        return *this; // return the result by reference
    }

    static FDTime infinity() noexcept {
        FDTime f;
        f.is_inf=true;
        return f;
    }

};

inline FDTime operator+(FDTime lhs, const FDTime& rhs) noexcept  // first arg by value, second by const ref
{
    lhs += rhs; // reuse compound assignment
    return lhs; // return the result by value
}

inline FDTime operator-(FDTime lhs, const FDTime& rhs) noexcept  // first arg by value, second by const ref
{
    lhs -= rhs; // reuse compound assignment
    return lhs; // return the result by value
}

//infinity time is considered higher than active
//all operations return false when comparing passive times
inline bool operator==(const FDTime& lhs, const FDTime& rhs) noexcept
{
    if (lhs._q_num != rhs._q_num || lhs._q_denom != rhs._q_denom || lhs._exp != rhs._exp){
        if (lhs._q_num != rhs._q_num || lhs._q_denom != rhs._q_denom ) {
            // adjust q and exp;
            // - find common denom
            int denom = boost::math::lcm(lhs._q_denom, rhs._q_denom);
            // - find common num
            int num = boost::math::gcd(lhs._q_num * (denom/lhs._q_denom), rhs._q_num * (denom/rhs._q_denom));
            // - fix m
            lhs._m = lhs._m * (denom / lhs._q_denom) * (lhs._q_num / num); //common num should never be 0.
            rhs._m = rhs._m * (denom / rhs._q_denom) * (rhs._q_num / num); //common num should never be 0.
            //simplify factor
            int simplify = boost::math::gcd(num, denom);
            // - fix the qs
            lhs._q_num = rhs._q_num  = num/simplify;
            lhs._q_denom = rhs._q_denom = denom/simplify;
        }
        if(lhs._exp != rhs._exp) {
            // - fix exps
            if (lhs._exp < rhs._exp) {
                rhs._m<<=(lhs._exp - rhs._exp);
                rhs._exp = lhs._exp;
            } else {
                lhs._m<<=(rhs._exp - lhs._exp);
                lhs._exp = rhs._exp;
            }
        }
    }
    return (lhs.is_inf && rhs.is_inf) || (!lhs.is_inf && !rhs.is_inf && lhs._m == rhs._m) ;
}


inline bool operator!=(const FDTime& lhs, const FDTime& rhs) noexcept {return !operator==(lhs,rhs);}

inline bool operator< (const FDTime& lhs, const FDTime& rhs) noexcept
{
    if (lhs._q_num != rhs._q_num || lhs._q_denom != rhs._q_denom || lhs._exp != rhs._exp){
        throw 1;
    }
    return 	((! lhs.is_inf) && (!rhs.is_inf) && (lhs._m < rhs._m ))
            || (( rhs.is_inf) && (!lhs.is_inf));
}

inline bool operator> (const FDTime& lhs, const FDTime& rhs) noexcept {return  operator< (rhs,lhs);}

inline bool operator<=(const FDTime& lhs, const FDTime& rhs) noexcept {return !operator> (lhs,rhs);}

inline bool operator>=(const FDTime& lhs, const FDTime& rhs) noexcept {return !operator< (lhs,rhs);}

//inline bool is_infinity(const FDTime& rhs){ return rhs.is_inf ;}

inline std::ostream& operator<<(std::ostream& out, const FDTime& t) noexcept {
    if (t.is_inf) out << "Infinity";
    else out << t._q_num << "/" << t._q_denom << " x " << t._m << " x 2^" << t._exp;
    return out;
}

//overload >> does not support infinity, only Ts
inline std::istream& operator>> (std::istream& stream, FDTime& rhs) noexcept {
    int a;
    stream >> a;
    rhs= FDTime{a};
    return stream;
}


}

//overload of numeric limits
namespace std {
template<>
class numeric_limits<cdpp::FDTime>{
public:
    static constexpr bool is_specialized = true;
    static cdpp::FDTime min() noexcept { return cdpp::FDTime{numeric_limits<int>::min(), 1, numeric_limits<int>::max(), numeric_limits<int>::max()}; }
    static cdpp::FDTime max() noexcept { return cdpp::FDTime{numeric_limits<int>::max(), 1, numeric_limits<int>::max(), numeric_limits<int>::max()}; }
    static cdpp::FDTime lowest() noexcept { return cdpp::FDTime{numeric_limits<int>::max(), -1, numeric_limits<int>::max(), numeric_limits<int>::max()}; }

    static constexpr int  digits = numeric_limits<int>::digits;
    static constexpr int  digits10 = numeric_limits<int>::digits10;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = true;
    static constexpr int radix = 2;
    static cdpp::FDTime epsilon() noexcept { return cdpp::FDTime{1, numeric_limits<int>::max(), 1, -1 * numeric_limits<int>::max()}; }
    static cdpp::FDTime round_error() noexcept { return cdpp::FDTime(0); }

    static constexpr int  min_exponent = -1 * numeric_limits<int>::max();
    static constexpr int  min_exponent10 = min_exponent/radix;
    static constexpr int  max_exponent = numeric_limits<int>::max();
    static constexpr int  max_exponent10 = max_exponent/radix;

    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr bool has_denorm_loss = false;
    static cdpp::FDTime infinity() noexcept { return cdpp::FDTime::infinity(); }
//    static constexpr cdpp::FDTime quiet_NaN() noexcept { return T(); }
//    static constexpr cdpp::FDTime signaling_NaN() noexcept { return T(); }
//    static constexpr cdpp::FDTime denorm_min() noexcept { return T(); }

    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;

    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
//    static constexpr float_round_style round_style = round_toward_zero;

};
}



#endif // FDTIME_H
