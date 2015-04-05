/*

  RealLib, a library for efficient exact real computation
  Copyright (C) 2006 Branimir Lambov

  This library is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this library except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

/*

    MachineEstimate.h

    Machine precision estimates.
    Classes:

        MachineEstimate - implements interval arithmetic on machine
        precision floats. To be used for the first fast stage
        of the evaluation.
        The class is to be inlined for good performance.

 */

#ifndef FILE_MACHINE_ESTIMATE_H
#define FILE_MACHINE_ESTIMATE_H

#include <stdlib.h>
#include <limits.h>
#include <ostream>
#include <exception>
#include <cfloat>
#include <cmath>
#include <iomanip>
#ifdef _MSC_VER
#include <emmintrin.h>
#else
#include <xmmintrin.h>
#endif
#include "defs.h"
#include "RealEstimate.h"

namespace RealLib {

typedef __m128d MachineEstimateBaseType;
#define REALLIB_MACHINEESTIMATE_ALIGNMENT_REQUIRED

// class MachineEstimate's definitions start here
class MachineEstimate;

// operations
static inline    MachineEstimate operator - (const MachineEstimate &arg);
static inline    MachineEstimate recip(const MachineEstimate &arg);

static inline    MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline    MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline    MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline    MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs);

// fast multiplication
static inline    MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs);
// and division
static inline    MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs);

static inline    std::ostream& operator <<(std::ostream &os, const MachineEstimate &e);

class MachineEstimate {
private:
public:

    // we'll be using round-to-minus-infinity mode only,
    // emulating round-to-plus-infinity via -round(-value)
    // we'll keep the high bound negated for quicker operation

    // the data will be stored in one __m128d variable,
    // -high bound in first element (the one that you can apply _sd operations to)
    // low bound in second
    __m128d interval;        

    static __m128d signmask;        // mask of only 1 in the MSB of the first double
    static __m128d mdelta;            // minus smallest representable number in both sides
    static __m128d half;            // 0.5 in both 
    static __m128d mhalf;            // 0.5, -0.5
    static __m128d mquarter;        // 0.25, -0.25
    static __m128d zero;            // 0.0
    static __m128d mone;            // -1.0
    static __m128d one;            // 1.0, -1.0
    static __m128d sqrt_corr;        // 0.0, NaN
    static __m128d coeff_sin[8];    // coefficients for sine
    static __m128d coeff_cos[8];
    static __m128d coeff_log[20];
    static __m128d coeff_exp[12];
    static __m128d coeff_atan[20];
    static __m128i expmask;            // exponent mask, i.e. 0x7ff000...
    static __m128i expbias;            // exponent bias - 1, i.e. 0x3fe000...
    static __m128d explimit;        // exponent limit, highest exponent that we admit, 1020
    static __m128d sign;                // 0x800...
    static __m128d log2e;            // log_2 e = 1/ln2, 1.49..., ++
    static __m128d sqrt2;            // sqrt(2), 1.41..., ++
    static __m128d sqrtsqrt2;        // sqrt(sqrt(2)), 1.1892..., ++
    static __m128d truehigh;        // 0x00000000, 0xffffffff. used in log to get the right comparison result

    static __m128d pi, pi_over_2, pi_over_4, rpi4, pi2;    // pi constants
    //static __m128d onethird;
    //static __m128i epi32incorrect;
    static __m128d three, four;    // for multiplication, 3.0 and 4.0
    static __m128d ln2;                // 0.69..., +-
    static __m128d ln2c;                // as ln2, but ++ for faster multiplication
    static int SavedRoundingMode; 

    MachineEstimate(double l, double h) : interval(_mm_set_pd(l, -h)) {}
    MachineEstimate(__m128d src) : interval(src) {}

    std::ostream& PrintInterval(std::ostream &os) const {
        double d[2];
        _mm_storeu_pd(d, interval);
        return (os << "(" << d[1] << ", " << -d[0] << ")");
    }

    std::ostream& PrintIntervalWHex(std::ostream &os) const {
        double d[2];
        u64 l[2];
        _mm_storeu_pd(d, interval);
        _mm_storeu_pd((double*)l, GetInterval());
        return (os << "(" << d[1] << ", " << -d[0] << "), hex diff " << 
                std::hex << (l[0] - l[1]) << std::dec);
    }

    // gets the sum of high and low, negated in the first component and positive in the second
    // (i.e. a proper interval)
    __m128d Sum() const { return _mm_sub_pd(interval, _mm_shuffle_pd(interval, interval, 1)); }
    // gets the difference, negated in both components of the __m128d
    // (i.e. negate second component to make a proper interval)
    __m128d MinusDiff() const { return _mm_add_pd(interval, _mm_shuffle_pd(interval, interval, 1)); }
    __m128d GetInterval() const { return _mm_xor_pd(interval, signmask); }
    void GetInterval(double i[2]) const { _mm_storeu_pd(i, GetInterval()); }

public:
    static int BeginComputation();    // return the previous rounding mode
    static void FinishComputation(int rm = 0);

    class Computation {
        int RoundingMode;
        Computation() : RoundingMode(BeginComputation()) {}
        ~Computation() { FinishComputation(RoundingMode); }
    };

    bool IsValueValid() const {
#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
        return !(_MM_GET_EXCEPTION_STATE() & 
                (_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW));
#else
        return !!_finite(weak_AsDouble());
#endif
    }

    MachineEstimate(double v = 0.0) : interval(_mm_xor_pd(_mm_set1_pd(v), signmask)) {}
    MachineEstimate(const char *val) {
        double v = atof(val);
        __m128d z(_mm_set1_pd(v));    // load double
        z = _mm_xor_pd(z, signmask);    // negate first component
        interval = _mm_add_pd(z, mdelta);        // round down
    }

    // error functions
    MachineEstimate GetError() const { return _mm_mul_pd(mhalf, MinusDiff()); }
    MachineEstimate& SetError(const MachineEstimate &err) { 
        __m128d s = _mm_mul_pd(Sum(), half);
        // assuming a positive error!
        __m128d e = _mm_shuffle_pd(err.interval, err.interval, 0);    // negated
        interval = _mm_add_pd(s, e);
        return *this;
    }
    MachineEstimate& AddError(const MachineEstimate &err) {
        // assuming a positive error!
        __m128d e = _mm_shuffle_pd(err.interval, err.interval, 0);    // negated
        interval = _mm_add_pd(interval, e);
        return *this;
    }

    // a lower bound on the correct binary digits
    // uses the exponents of the value and error to calculate it quickly
    i32 GetRelativeError() const {
        int e;
        double d;
        if (_mm_comieq_sd(Sum(), zero)) return I32_MAX;
        _mm_store_sd(&d, _mm_div_sd(Sum(), MinusDiff()));
        if (frexp(d, &e) == 0) return I32_MIN;
        else return e;
    }

    // get a rough estimate of the precision
    // used to determine the length of the approximations to functions
    u32 GetPrecision() const
    { return 3; }
    MachineEstimate& SetPrecision(u32 prec)
    { return *this; }

    // comparisons
    // these come in two flavors, strong (true if real is in relation to rhs)
    bool IsNegative() const
    { return !!_mm_comigt_sd(interval, zero); }    // -high > 0
    bool IsPositive() const
    { return (- *this).IsNegative(); }
    bool IsNonZero() const
    { return IsPositive() || IsNegative(); }

    // equality test is undecidable (i.e. would yield false for any precision)
    // thus ==, <= and >= are not included
    // also !(x<y) does not mean y<=x
    bool operator < (const MachineEstimate &rhs) const
    { return (*this - rhs).IsNegative(); }
    bool operator > (const MachineEstimate &rhs) const
    { return (rhs - *this).IsNegative(); }
    bool operator != (const MachineEstimate &rhs) const
                    { return (*this - rhs).IsNonZero(); }

    // truncation
    // used to make sure only arguments within the domain of the function
    // are processed for the closed ends of the domain.
    //    To this end, truncates the approximation interval so that
    // the indicated real numbers are thrown out. If nothing remains,
    // raise a DomainException(origin).

    // warning: an error in the approximation of the bound will be added to the
    // error in the end result, i.e. if [0, 3] is truncated below [1, 0.5], the
    // result will be [0.5, 3.5]. To avoid problems, use exact bounds (e.g. double)!

    // #ifndef REALLIB_RELY_ON_SSE_EXCEPTIONS
    // truncation checks if the inputs are valid values
    // to avoid losing information about an error
    // otherwise we only do the check at the end


    // removes the part of the approximation interval that is negative
    MachineEstimate TruncateNegative(const char *origin = "Truncate") const
    { if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comigt_sd(interval, zero))    // -high > 0 means
        throw DomainException(origin);    // provably negative input
    else return _mm_max_pd(interval, MachineEstimate::sqrt_corr);        }

    // removes the part of the approximation that is below a certain lower bound
    MachineEstimate TruncateBelow(const MachineEstimate &l, const char *origin = "Truncate") const
    {    __m128d ll(_mm_shuffle_pd(l.interval, l.interval, 1));
    __m128d a(_mm_max_pd(interval, l.interval)); // adjust the lower side to the largest of this and limit
    __m128d lh(_mm_xor_pd(interval, signmask));
    __m128d b(_mm_move_sd(a, interval));             // and keep the high part from this
    if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comilt_sd(lh, ll))    // high < l.low
        throw DomainException(origin);
    else return b;    }
    //{ return (*this - l).TruncateNegative(origin) + l; }

    MachineEstimate TruncateBelow(double l, const char *origin = "Truncate") const
    {    __m128d ll(MachineEstimate(l).interval);
    __m128d a(_mm_max_pd(interval, ll)); // adjust the lower side to the largest of this and limit
    __m128d b(_mm_move_sd(a, interval));             // and keep the high part from this
    if (_mm_comigt_sd(interval, ll))    // high < l.low, lh and ll are negated here
        throw DomainException(origin);
    else return b;    }

    // removes the part of the approximation that is above a certain upper bound
    MachineEstimate TruncateAbove(const MachineEstimate &h, const char *origin = "Truncate") const
    {    __m128d ll(_mm_shuffle_pd(interval, interval, 1));
    __m128d a(_mm_max_sd(interval, h.interval)); // adjust the higher side to the smallest of this and limit
    // these are high sides, i.e. negated, hence max instead of min
    __m128d lh(_mm_xor_pd(h.interval, signmask));
    if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comilt_sd(lh, ll))    // low > h.high
        throw DomainException(origin);
    else return a;    }
    //{ return h - (h - *this).TruncateNegative(origin); }

    MachineEstimate TruncateAbove(double h, const char *origin = "Truncate") const
    {    __m128d ll(_mm_shuffle_pd(interval, interval, 1));
    __m128d lh(_mm_set1_pd(h));
    __m128d lm(_mm_xor_pd(lh, signmask));
    __m128d a(_mm_max_sd(interval, lm)); // adjust the higher side to the smallest of this and limit
    // these are high sides, i.e. negated, hence max instead of min
    if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comilt_sd(lh, ll))    // low > h.high
        throw DomainException(origin);
    else return a;    }

    // removes the part of the approximation outside the specified interval
    MachineEstimate TruncateTo(double l, double h, const char *origin = "Truncate") const
    { __m128d a = MachineEstimate(l, h).interval;
    __m128d b = _mm_max_pd(interval, a);        // smallest high and largest low
    __m128d c = _mm_shuffle_pd(b, b, 1);
    __m128d d = _mm_xor_pd(b, signmask);
    if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comigt_sd(c, d))        // low > high in b
        throw DomainException(origin);
    else return b;        }

    //{ MachineEstimate e(max(l, low), min(h, high));
    //if (e.high < e.low) throw DomainException(origin);
    //else return e;}

    MachineEstimate TruncateTo(const MachineEstimate &l, const MachineEstimate &h, const char *origin = "Truncate") const
    { __m128d a = _mm_shuffle_pd(h.interval, l.interval, 2);        // check this
    __m128d b = _mm_max_pd(interval, a);        // smallest high and largest low
    __m128d c = _mm_shuffle_pd(b, b, 1);
    __m128d d = _mm_xor_pd(b, signmask);
    if (REALLIB_MACHEST_INVALID_CHECK(*this) || _mm_comigt_sd(c, d))        // low > high in b
        throw DomainException(origin);
    else return b;        }

    //{ return ((h-l) - (*this - l).TruncateNegative(origin)).TruncateNegative(origin) + l; }

    // and weak (true if m_Value is in relation to rhs)
    // should only be used if the transformation being aplied
    // would not differentiate on the two cases, e.g. to choose
    // whether to evaluate sin(x) and sin(pi - x)

    bool weak_IsPositive() const
    { return MachineEstimate(Sum()).IsNegative(); }
    bool weak_IsNegative() const
    { return MachineEstimate(Sum()).IsPositive(); }
    // Estimate does not provide zero test, so we don't either MachineEstimate

    bool weak_lt(const MachineEstimate &rhs) const
    { return (*this - rhs).weak_IsNegative(); }
    bool weak_eq(const MachineEstimate &rhs) const
    { return !!_mm_comieq_sd((*this - rhs).interval, zero); }

    bool weak_gt(const MachineEstimate &rhs) const
    { return rhs.weak_lt(*this); }

    bool weak_le(const MachineEstimate &rhs) const
    { return !weak_gt(rhs); }
    bool weak_ne(const MachineEstimate &rhs) const
    { return !weak_eq(rhs); }
    bool weak_ge(const MachineEstimate &rhs) const
    { return !weak_lt(rhs); }

    // among the weak operations is also rounding
    // the returned MachineEstimate is assumed exact
    // only to be used on periodic functions!
    MachineEstimate weak_round() const
    { return floor(weak_AsDouble() + 0.5); }

    // weak normalize, i.e. return an exponent such that
    // a >> a.weak_normalize()
    // is in the range [0.5, 1).
    i32 weak_normalize() const {
        int e;
        frexp(weak_AsDouble(), &e);
        return e;
    }

    // weak conversion
    double weak_AsDouble() const
    { double d;
    _mm_store_sd(&d, Sum());
    return d * -0.5; }

    // output
    char *weak_AsDecimal(char *buffer, u32 buflen) const
    {
#ifdef _MSC_VER
        return _gcvt(weak_AsDouble(), buflen-7, buffer);
#else
        return gcvt(weak_AsDouble(), buflen-7, buffer);
#endif
    }

    MachineEstimate weak_Center()
    { return weak_AsDouble(); }

    // operations
    friend MachineEstimate operator - (const MachineEstimate &arg);
    friend MachineEstimate recip(const MachineEstimate &arg);

    friend MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs);
    friend MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs);
    friend MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs);
    friend MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs);

    // fast multiplication
    friend MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs);
    // and division
    friend MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs);

    // binary shift
    MachineEstimate operator << (i32 howmuch) const
    { double d = ldexp(1.0, howmuch);
    return MachineEstimate(_mm_mul_pd(interval, _mm_set1_pd(d))); }
    MachineEstimate operator >> (i32 howmuch) const
    { return *this << -howmuch; }

    MachineEstimate& operator += (const MachineEstimate &rhs)
                 { return *this = *this + rhs; }
    MachineEstimate& operator -= (const MachineEstimate &rhs)
                 { return *this = *this - rhs; }
    MachineEstimate& operator *= (const MachineEstimate &rhs)
                 { return *this = *this * rhs; }
    MachineEstimate& operator /= (const MachineEstimate &rhs)
                 { return *this = *this / rhs; }


    MachineEstimate& operator >>= (i32 rhs) 
                 { return *this = *this >> rhs; }
    MachineEstimate& operator <<= (i32 rhs) 
                 { return *this = *this << rhs; }
    MachineEstimate& operator *= (i32 rhs) 
                 { return *this = *this * rhs; }
    MachineEstimate& operator /= (i32 rhs) 
                 { return *this = *this / rhs; }

    // should probably be somewhere else
    // conversion to string
    // char *AsDecimal(char *buffer, u32 buflen);
    friend    
    std::ostream& operator <<(std::ostream &os, const MachineEstimate &e);

};

// operations
static inline
MachineEstimate operator - (const MachineEstimate &arg)
{    // if you simply reverse the bounds' order you get the negation
    return MachineEstimate(_mm_shuffle_pd(arg.interval, arg.interval, 1)); }

static inline
MachineEstimate recip(const MachineEstimate &arg)
{
    if (!arg.IsNonZero()) throw PrecisionException("recip");

    __m128d a = _mm_div_pd(MachineEstimate::mone, arg.interval);
    return _mm_shuffle_pd(a, a, 1);
}

static inline
MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ return MachineEstimate(_mm_add_pd(lhs.interval, rhs.interval)); }

static inline
MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ return lhs + (-rhs); }

static inline
__m128d myshuffle(const __m128d a, const __m128d b, int m)
{
    switch(m) {
    case 0:
        return _mm_shuffle_pd(a, b, 0);
    case 1:
        return _mm_shuffle_pd(a, b, 1);
    case 2:
        return _mm_shuffle_pd(a, b, 2);
    case 3:
    default:
        return _mm_shuffle_pd(a, b, 3);
    }
}

static inline 
__m128d ivchoice(__m128d a, const __m128d b)
{
    a = _mm_xor_pd(a, _mm_set_pd(-0.0, 0.0));
    a = myshuffle(a, a, _mm_movemask_pd(b));
    return a;
}

static inline
__m128d ivmul(__m128d a, const __m128d b)
{
    a = _mm_xor_pd(a, _mm_set_pd(-0.0, 0.0));
    a = myshuffle(a, a, _mm_movemask_pd(b));
    a = _mm_mul_pd(a, b);
    return a;
}

static inline
__m128d mulpn(const __m128d a, const __m128d b)
{
    return _mm_mul_pd(a, _mm_xor_pd(b, MachineEstimate::signmask));
}

static inline
MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ 
    // 4 muls, sub delta version
    /*__m128d a = _mm_mul_pd(lhs.interval, rhs.interval);
    __m128d b = _mm_shuffle_pd(rhs.interval, rhs.interval, 1);
    __m128d c = _mm_sub_pd(MachineEstimate::zero, b);
    __m128d d = _mm_mul_pd(lhs.interval, c);
    __m128d e = _mm_min_pd(a, d);
    __m128d f = _mm_max_pd(a, d);
    __m128d g = _mm_sub_pd(MachineEstimate::mdelta, f);
    __m128d h = _mm_shuffle_pd(g, e, 1);
    __m128d i = _mm_shuffle_pd(g, e, 2);
    __m128d j = _mm_min_pd(h, i);
    return MachineEstimate(j);*/
    // lhs = (-a, b) rhs = (-c, d)

    // 8 muls version
    /*__m128d a = _mm_shuffle_pd(lhs.interval, lhs.interval, 1);    // b, -a
    __m128d b = _mm_mul_pd(lhs.interval, rhs.interval);            // ac, bd
    __m128d c = _mm_sub_pd(MachineEstimate::zero, rhs.interval);    // c, -d
    __m128d d = _mm_mul_pd(a, rhs.interval);                        // -bc, -ad
    __m128d e = _mm_mul_pd(lhs.interval, c);                        // -ac, -bd
    __m128d f = _mm_mul_pd(a, c);                                    // bc, ad
    __m128d g = _mm_min_pd(d, e);                                    // min(-bc,-ac), min(-ad, -bd)
    __m128d h = _mm_min_pd(b, f);                                    // min(ac, bc), min(bd, ad)
    __m128d i = _mm_unpackhi_pd(g, h);                            // min(-ad, -bd), min(ac, bc)
    __m128d j = _mm_unpacklo_pd(g, h);                            // min(-bc,-ac), min(bd, ad)
    __m128d k = _mm_min_pd(i, j);                                    // min(-all), min(all)
    return MachineEstimate(k);*/

    // choice2 ideas for faster hardware
    /*__m128d a = _mm_shuffle_pd(lhs.interval, lhs.interval, 1);
    __m128d z = _mm_shuffle_pd(rhs.interval, rhs.interval, 1);
    __m128d b = choice2mul(z, lhs.interval);
    __m128d d = choice2mul(rhs.interval, a);
    __m128d g = _mm_min_pd(b, d);
    return g;*/

    /*    __m128d a = _mm_shuffle_pd(rhs.interval, rhs.interval, 1);
    __m128d z = _mm_xor_pd(lhs.interval, MachineEstimate::signmask);
    __m128d b = _mm_shuffle_pd(z, z, 1);
    __m128d c = _mm_cmple_pd(z, MachineEstimate::sign);
    __m128d d = _mm_xor_pd(a, MachineEstimate::sign);
    __m128d y = _mm_shuffle_pd(c, c, 1);
    __m128d e = _mm_and_pd(c, d);
    __m128d f = _mm_andnot_pd(c, rhs.interval);
    __m128d g = _mm_or_pd(e, f);
    __m128d h = _mm_and_pd(y, d);
    __m128d i = _mm_andnot_pd(y, rhs.interval);
    __m128d j = _mm_or_pd(h, i);
    __m128d k = _mm_mul_pd(g, z);
    __m128d l = _mm_mul_pd(j, b);
    __m128d m = _mm_min_pd(k, l);
    return m;*/

    // 4mul choices version
    __m128d a = _mm_xor_pd(lhs.interval, MachineEstimate::signmask);
    __m128d b = rhs.interval;                            // no op
    __m128d p = _mm_shuffle_pd(b, b, 1);
    __m128d c = _mm_cmplt_pd(a, MachineEstimate::sign);
    __m128d d = _mm_xor_pd(p, MachineEstimate::sign);    
    __m128d g = _mm_shuffle_pd(c, c, 1);
    __m128d e = _mm_and_pd(c, d);
    __m128d f = _mm_andnot_pd(c, b);
    __m128d k = _mm_andnot_pd(g, b);
    __m128d i = _mm_and_pd(g, d);
    __m128d h = _mm_or_pd(e, f);
    __m128d l = _mm_shuffle_pd(a, a, 1);
    __m128d m = _mm_or_pd(i, k);
    __m128d j = _mm_mul_pd(a, h);
    __m128d n = _mm_mul_pd(l, m);
    __m128d o = _mm_min_pd(j, n);
    return o;


    // 4mul choices alternative
    /*    __m128d a = _mm_xor_pd(lhs.interval, MachineEstimate::signmask);
    __m128d b = _mm_xor_pd(rhs.interval, MachineEstimate::signmask);                            // no op
    __m128d c = _mm_cmplt_pd(a, MachineEstimate::zero);
    __m128d d = _mm_shuffle_pd(b, b, 1);    // two ops
    __m128d e = _mm_and_pd(c, d);
    __m128d f = _mm_andnot_pd(c, b);
    __m128d g = c;//_mm_shuffle_pd(c, c, 1);
    __m128d h = _mm_or_pd(e, f);
    __m128d i = _mm_and_pd(g, b);
    __m128d j = _mm_mul_pd(lhs.interval, h);
    __m128d k = _mm_andnot_pd(g, d);
    __m128d l = _mm_xor_pd(lhs.interval, MachineEstimate::sign);//_mm_shuffle_pd(a, a, 1);
    __m128d m = _mm_or_pd(i, k);
    __m128d n = _mm_mul_pd(l, m);
    n = _mm_shuffle_pd(n, n, 1);
    __m128d o = _mm_min_pd(j, n);
    return o;*/
}

static inline
MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ 
    __m128d r = _mm_xor_pd(rhs.interval, MachineEstimate::signmask);
    __m128d p = _mm_shuffle_pd(r, r, 1);

    __m128d a = _mm_cmpgt_pd(r, MachineEstimate::zero);
    if (_mm_movemask_pd(a) == 1) 
        throw PrecisionException("div");
    __m128d b = _mm_shuffle_pd(lhs.interval, lhs.interval, 1);
    b = _mm_xor_pd(b, MachineEstimate::sign);
    __m128d c = _mm_and_pd(a, lhs.interval);
    __m128d d = _mm_andnot_pd(a, b);
    __m128d e = _mm_or_pd(c, d);
    __m128d h = _mm_div_pd(e, p);
    __m128d i = _mm_div_pd(e, r);
    __m128d j = _mm_min_pd(h, i);
    return j;

    /*__m128d x = rhs.interval;
    __m128d y = lhs.interval;
        if (_mm_movemask_pd(x)==3)
                throw PrecisionException("div");
        __m128d a = _mm_shuffle_pd(x, x, 1);
        __m128d b = _mm_shuffle_pd(y, y, 1);
        __m128d c = ivchoice(b, x);
        __m128d d = ivchoice(y, a);    // here c==-d
        __m128d e = _mm_div_pd(c, x);
        __m128d f = _mm_div_pd(d, a);
        __m128d g = _mm_min_pd(e, f);
        return g;*/

    /*__m128d a = _mm_xor_pd(lhs.interval, MachineEstimate::signmask);
    __m128d p = rhs.interval;                            // no op
    __m128d d = _mm_shuffle_pd(p, p, 1);
    __m128d c = _mm_cmplt_pd(a, MachineEstimate::zero);
    __m128d b = _mm_xor_pd(p, MachineEstimate::sign);    
    __m128d g = _mm_shuffle_pd(c, c, 1);
    __m128d e = _mm_and_pd(c, d);
    __m128d f = _mm_andnot_pd(c, b);
    __m128d k = _mm_andnot_pd(g, b);
    __m128d i = _mm_and_pd(g, d);
    __m128d h = _mm_or_pd(e, f);
    __m128d l = _mm_shuffle_pd(a, a, 1);
    __m128d m = _mm_or_pd(i, k);
    __m128d j = _mm_div_pd(h, l);
    __m128d n = _mm_div_pd(m, a);
    __m128d o = _mm_min_pd(j, n);
    return o;*/

    //    return lhs * recip(rhs);
    /*    if (!rhs.IsNonZero()) throw PrecisionException("recip");

    __m128d dv = _mm_div_pd(MachineEstimate::mone, rhs.interval);
    __m128d b = _mm_shuffle_pd(dv, dv, 1);
    __m128d a = _mm_mul_pd(lhs.interval, b);
    __m128d c = _mm_sub_pd(MachineEstimate::zero, dv);
    __m128d d = _mm_mul_pd(lhs.interval, c);
    __m128d e = _mm_min_pd(a, d);
    __m128d f = _mm_max_pd(a, d);
    __m128d g = _mm_sub_pd(MachineEstimate::mdelta, f);
    __m128d h = _mm_shuffle_pd(g, e, 1);
    __m128d i = _mm_shuffle_pd(g, e, 2);
    __m128d j = _mm_min_pd(h, i);
    return MachineEstimate(j);*/
}

// fast multiplication
// (using the fact that 32 bits can be represented exactly in double)
static inline
MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs)
{ 
    /*
    __m128d v;
    v = _mm_cvtsi32_sd(v, rhs);
    v = _mm_shuffle_pd(v, v, 0);
    __m128d a(_mm_mul_pd(lhs.interval, v));
    __m128d b(_mm_shuffle_pd(a, a, 1));
    return _mm_min_pd(a, b);
     */
    if (rhs >=0) return MachineEstimate(_mm_mul_pd(lhs.interval, _mm_set1_pd(double(rhs)))); 
    else return -lhs * -rhs;
}

// and division
static inline
MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs)
{
    /*
        __m128d v;
        v = _mm_cvtsi32_sd(v, rhs);
        v = _mm_shuffle_pd(v, v, 0);
        __m128d a(_mm_div_pd(lhs.interval, v));
        __m128d b((-lhs).interval);
        v = _mm_xor_pd(v, MachineEstimate::sign);
        b = _mm_div_pd(b, v);
        return _mm_min_pd(a, b);
     */
    if (rhs > 0) {
        __m128d v;
        v = _mm_cvtsi32_sd(v, rhs);
        v = _mm_shuffle_pd(v, v, 0);
        // v= _mm_set1_pd(double(rhs));
        return MachineEstimate(_mm_div_pd(lhs.interval, v)); 
    } else if (rhs < 0) {
        __m128d v;
        v = _mm_cvtsi32_sd(v, -rhs);
        v = _mm_div_pd(lhs.interval, v);
        v = _mm_shuffle_pd(v, v, 0);
        // v= _mm_set1_pd(double(rhs));
        return MachineEstimate(v); 
    } else throw DomainException("integer div");
}

// shorthands
static inline 
MachineEstimate operator * (i32 lhs, const MachineEstimate &rhs)
{ return rhs * lhs; }

static inline 
MachineEstimate operator / (i32 lhs, const MachineEstimate &rhs)
{ return recip(rhs) * lhs; }

// C++-style output
static inline
std::ostream& operator <<(std::ostream &os, const MachineEstimate &e)
{    return os.operator<<(e.weak_AsDouble()); }
//    return os << e.weak_AsDouble(); }

#define REALLIB_ME_PRECISE_SQRT
static inline MachineEstimate sqrt(const MachineEstimate &arg)
{
#ifdef REALLIB_ME_PRECISE_SQRT
    __m128d a = _mm_xor_pd(arg.interval, MachineEstimate::signmask);
    __m128d b = _mm_sqrt_pd(a);
    __m128d d = _mm_mul_pd(a, a);

    // according to our convention we do not care about the invalid part of the intervals
    // but a NaN in the high part would mean the value was provably negative
    if (_mm_comilt_sd(b, b)) // i.e. sqrt(high) == NaN
        throw DomainException("sqrt");    // provably negative input
    b = _mm_max_pd(b, MachineEstimate::sqrt_corr);        // to get rid of a possible NaN in the low part

    if (_mm_comineq_sd(d, a)) {
        __m128d c = _mm_sub_sd(MachineEstimate::mdelta, b);
        b = _mm_move_sd(b, c);
    } else b = _mm_xor_pd(b, MachineEstimate::signmask);
    return MachineEstimate(b);
#else
#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
    __m128d a = _mm_xor_pd(arg.interval, MachineEstimate::signmask);
    __m128d x = _mm_max_pd(a, MachineEstimate::sqrt_corr);
    __m128d b = _mm_sqrt_pd(x);
    __m128d c = _mm_sub_sd(MachineEstimate::mdelta, b);
    __m128d d = _mm_move_sd(b, c); 
    return MachineEstimate(d);
#else
    __m128d a = _mm_xor_pd(arg.interval, MachineEstimate::signmask);
    __m128d b = _mm_sqrt_pd(a);
    __m128d c = _mm_sub_sd(MachineEstimate::mdelta, b);

    // according to our convention we do not care about the invalid part of the intervals
    // but a NaN in the high part would mean the value was provably negative
    if (_mm_comilt_sd(b, b)) // i.e. sqrt(high) == NaN
        throw DomainException("sqrt");    // provably negative input
    b = _mm_max_pd(b, MachineEstimate::sqrt_corr);        // to get rid of a possible NaN in the low part

    __m128d d = _mm_move_sd(b, c); 
    return MachineEstimate(d);
#endif
#endif
}

static inline MachineEstimate abs(const MachineEstimate &arg)
{
    __m128d a = arg.interval;                    // a, b
    __m128d b = _mm_shuffle_pd(a, a, 1);
    __m128d c = _mm_min_sd(a, b);                // if b has a NaN this preserves it
    __m128d d = _mm_max_sd(b, a);                // if a has a NaN this preserves it
    __m128d e = _mm_max_sd(MachineEstimate::zero, d); // keep a possible NaN in a
    __m128d f = _mm_unpacklo_pd(c, e);        
    return MachineEstimate(f);
}

static inline MachineEstimate sq(const MachineEstimate &arg)
{
    MachineEstimate x(abs(arg));
    return MachineEstimate(_mm_mul_pd(x.interval, _mm_xor_pd(x.interval, MachineEstimate::signmask))); 
}

static inline MachineEstimate rsqrt(const MachineEstimate &arg)
{
    // could do better (Newton-Raphson?)
    return recip(sqrt(arg));
} 

MachineEstimate cos(const MachineEstimate &arg);
MachineEstimate sin(const MachineEstimate &arg);
MachineEstimate log(const MachineEstimate &arg);
MachineEstimate exp(const MachineEstimate &arg);

MachineEstimate asin(const MachineEstimate &arg);
MachineEstimate atan(const MachineEstimate &arg);
MachineEstimate acos(const MachineEstimate &arg);
MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x);

template <class TYPE>
TYPE pi(unsigned int prec);// {return TYPE();}

template <class TYPE>
TYPE ln2(unsigned int prec);// {return TYPE();}

template <>
MachineEstimate ln2(unsigned int prec);

template <>
MachineEstimate pi(unsigned int prec);

static inline MachineEstimate tan(const MachineEstimate &arg)
{
    return sin(arg)/cos(arg);
}

}    // namespace

#endif // FILE
