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
#include "defs.h"
#include "RealEstimate.h"

namespace RealLib {

static inline double max(double u, double v) { return u < v ? v : u; }
static inline double min(double u, double v) { return u > v ? v : u; }

typedef double MachineEstimateBaseType[2];

// class MachineEstimate's definitions start here
class MachineEstimate;

// operations
static inline MachineEstimate operator - (const MachineEstimate &arg);
static inline MachineEstimate recip(const MachineEstimate &arg);

static inline MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs);

// fast multiplication
static inline MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs);
// and division
static inline MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs);

std::ostream& operator <<(std::ostream &os, const MachineEstimate &e);

class MachineEstimate {
private:
public:
    double low, high;

    static double minusinf;
    static double plusinf;
    static double OnePEps;
    static double OneMEps;

    MachineEstimate(double l, double h) : low(l), high(h) {}
    double Sum() const { return high+low; }
    double Diff() const { return high-low; }

    static double RoundUp(double v)
    { return v > 0 ? v * OnePEps : v * OneMEps; }
    //{ return _nextafter(v, plusinf); }
    static double RoundDown(double v)
    { return v > 0 ? v * OneMEps : v * OnePEps; }
    //{ return _nextafter(v, minusinf); }

    static double RoundToZero(double v)
    { return v * OneMEps; }
    static double RoundFromZero(double v)
    { return v * OnePEps; }

    std::ostream& PrintInterval(std::ostream &os) const {
        return (os << "(" << low << ", " << high << ")");
    }

    std::ostream& PrintIntervalWHex(std::ostream &os) const {
        return PrintInterval(os) << ", hex diff " << std::hex
                << *((i64*)(&high)) - *((i64*)(&low)) << std::dec;
    }

    void GetInterval(double i[2]) const {
        i[1] = low;
        i[0] = high;
    }

    // positive multiples, positive result
    MachineEstimate MulPositive(double v) const {
        return MachineEstimate(RoundToZero(low * v), RoundFromZero(high * v));
    }

    MachineEstimate MulPositive(const MachineEstimate &v) const {
        return MachineEstimate(RoundToZero(low * v.low), RoundFromZero(high * v.high));
    }

    // only the right hand-side is known to be positive, the result may be negative
    MachineEstimate MulPositiveRHS(double v) const {
        return MachineEstimate(RoundDown(low * v), RoundUp(high * v));
    }

    MachineEstimate MulPositiveRHS(const MachineEstimate &v) const {
        return MachineEstimate(RoundDown(low * v.low), RoundUp(high * v.high));
    }

    // the result needs to be positive
    MachineEstimate AddPositive(const MachineEstimate &v) const {
        return MachineEstimate(RoundToZero(low + v.low), RoundFromZero(high + v.high));
    }

    // multiplication by double, no restrictions
    MachineEstimate MulDouble(double v) const {
        return v >= 0 ? MulPositiveRHS(v) :
                MachineEstimate(RoundDown(high * v), RoundUp(low * v));
    }

    // positive multiples and positive result of the addition
    template <class A2>
    MachineEstimate AddProductPositive(const MachineEstimate &m, A2 d) const {
        return AddPositive(m.MulPositive(d));
    }

    // positive multiples, positive result of the substraction
    template <class A2>
    MachineEstimate SubProductPositive(const MachineEstimate &m, A2 d) const {
        return AddPositive(-m.MulPositive(d));
    }

public:
    static void BeginComputation();
    static void FinishComputation();

    bool IsValueValid() const { return (_finite(low) && _finite(high)); }

    MachineEstimate(double v = 0.0) : low(v), high(v) {}
    MachineEstimate(const char *val) : low(RoundDown(atof(val))), high(RoundUp(atof(val))) {}

    // error functions
    MachineEstimate GetError() const { return RoundUp(high - low)/2; }
    MachineEstimate& SetError(const MachineEstimate &err) {
        double s = (high + low)*0.5;
        double e = max(fabs(err.low), fabs(err.high));
        low = RoundDown(s - e);
        high = RoundUp(s + e);
        return *this;
    }
    MachineEstimate& AddError(const MachineEstimate &err) {
        double e = max(fabs(err.low), fabs(err.high));
        low = RoundDown(low - e);
        high = RoundUp(high + e);
        return *this;
    }

    // a lower bound on the correct binary digits
    // uses the exponents of the value and error to calculate it quickly
    i32 GetRelativeError() const {
        int e;
        if (high==low) return I32_MAX;
        else if (frexp((high+low)/(high-low), &e) == 0) return I32_MIN;
        else return e;
    }

    /*
     MachineEstimate& AddRoundingError() 
     { m_Error = m_Error + RoundingError(m_Value, m_Value.AdditionRoundingError()); 
         return *this; }

        MachineEstimate TheRoundingError() const // rounding error is assumed to be no more than
                                                                                // one in the least significant bit of the mantissa
                                                                                // Note! Newton-Raphson reciprocal is incorrect in the
                                                                                // least significant word (handled by recip())
     { return RoundingError(m_Value); }
     */

    // get a rough estimate of the precision
    // used to determine the length of the approximations to functions
    u32 GetPrecision() const
    { return 3; }
    MachineEstimate& SetPrecision(u32 prec)
    { return *this; }

    // truncation
    // used to make sure only arguments within the domain of the function
    // are processed for the closed ends of the domain.
    //    To this end, truncates the approximation interval so that
    // the indicated real numbers are thrown out. If nothing remains,
    // raise a DomainException(origin).

    // warning: an error in the approximation of the bound will be added to the
    // error in the end result, i.e. if [0, 3] is truncated below [1, 0.5], the
    // result will be [0.5, 3.5]. To avoid problems, use exact bounds (e.g. double)!

    // removes the part of the approximation interval that is negative
    MachineEstimate TruncateNegative(const char *origin = "Truncate") const
    { if (high < 0) throw DomainException(origin);
    else return MachineEstimate(max(low, 0.0), high); }

    // removes the part of the approximation that is below a certain lower bound
    MachineEstimate TruncateBelow(double l, const char *origin = "Truncate") const
    { if (high < l) throw DomainException(origin);
    else return MachineEstimate(max(low, l), high); }

    MachineEstimate TruncateBelow(const MachineEstimate &l, const char *origin = "Truncate") const
    { return (*this - l).TruncateNegative(origin) + l; }

    // removes the part of the approximation that is above a certain upper bound
    MachineEstimate TruncateAbove(double h, const char *origin = "Truncate") const
    { if (low > h) throw DomainException(origin);
    else return MachineEstimate(low, min(high, h)); }

    MachineEstimate TruncateAbove(const MachineEstimate &h, const char *origin = "Truncate") const
    { return h - (h - *this).TruncateNegative(origin); }

    // removes the part of the approximation outside the specified interval
    MachineEstimate TruncateTo(double l, double h, const char *origin = "Truncate") const
    {
        MachineEstimate e(max(l, low), min(h, high));
        if (e.high < e.low) throw DomainException(origin);
        else return e;
    }

    MachineEstimate TruncateTo(const MachineEstimate &l, const MachineEstimate &h, const char *origin = "Truncate") const
    { return ((h-l) - (*this - l).TruncateNegative(origin)).TruncateNegative(origin) + l; }

    // comparisons
    // these come in two flavors, strong (true if real is in relation to rhs)
    bool IsPositive() const
    { return low > 0; }
    bool IsNegative() const
    { return high < 0; }
    bool IsNonZero() const
    { return IsPositive() || IsNegative(); }

    // equality test is undecidable (i.e. would yield false for any precision)
    // thus ==, <= and >= are not included
    // also !(x<y) does not mean y<=x
    bool operator < (const MachineEstimate &rhs) const
    { return high < rhs.low; }
    bool operator > (const MachineEstimate &rhs) const
    { return low > rhs.high; }
    bool operator != (const MachineEstimate &rhs) const
                        { return *this < rhs || *this > rhs; }

    // and weak (true if m_Value is in relation to rhs)
    // should only be used if the transformation being aplied
    // would not differentiate on the two cases, e.g. to choose
    // whether to evaluate sin(x) and sin(pi - x)

    bool weak_IsPositive() const
    { return high > -low; }
    bool weak_IsNegative() const
    { return low < -high; }
    bool weak_IsNonZero() const
    { return low == -high; }

    bool weak_lt(const MachineEstimate &rhs) const
    { return Sum() < rhs.Sum(); }
    bool weak_eq(const MachineEstimate &rhs) const
    { return Sum() == rhs.Sum(); }

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
    { return floor((Sum() + 1.0)*0.5); }

    // weak normalize, i.e. return an exponent such that 
    // a >> a.weak_normalize()
    // is in the range [0.5, 1).
    i32 weak_normalize() const {
        int e;
        frexp(Sum(), &e);
        return e-1; 
    }

    // weak conversion
    double weak_AsDouble() const
    { return Sum()*0.5; }
    // output
    char *weak_AsDecimal(char *buffer, u32 buflen) const
    {
#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable:4996)
        return _gcvt(weak_AsDouble(), buflen-7, buffer);
#pragma warning (pop)
#else
        return gcvt(weak_AsDouble(), buflen-7, buffer);
#endif
    }

    MachineEstimate weak_Center()
    { return weak_AsDouble(); }

    /*
   // exponent and mantissa operations
   // needed to get initial approximations via double
   double weak_MantissaAsDouble() const
      { return m_Value.MantissaAsDouble(); }

   i32 weak_Exponent() const
      { return m_Value.Exponent(); }
   MachineEstimate& AddToExponent(i32 exp)
      { m_Value.AddToExponent(exp); 
        return *this; }
   char *weak_MantissaAsDecimal(char *buf, u32 buflen) const
      { return m_Value.MantissaAsDecimal(buf, buflen); }
     */
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
    // warning: this returns HUGE_VAL instead of infinity!
    { return MachineEstimate(ldexp(low, howmuch), ldexp(high, howmuch)); }  
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
{ return MachineEstimate(-arg.high, -arg.low); }

static inline
MachineEstimate recip(const MachineEstimate &arg)
{
    if (!arg.IsNonZero()) throw PrecisionException("recip");
    double l = 1.0 / arg.high;
    double h = 1.0 / arg.low;
    return MachineEstimate(MachineEstimate::RoundDown(l), 
                           MachineEstimate::RoundUp(h));
}

static inline
MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ return MachineEstimate(MachineEstimate::RoundDown(lhs.low + rhs.low), 
                         MachineEstimate::RoundUp(lhs.high + rhs.high)); }

static inline
MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ return lhs + (-rhs); }

static inline
MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ double ll = lhs.low * rhs.low;
double lh = lhs.low * rhs.high;
double hl = lhs.high * rhs.low;
double hh = lhs.high * rhs.high;
return MachineEstimate(MachineEstimate::RoundDown(min(min(ll, lh), min(hl, hh))),
                       MachineEstimate::RoundUp(max(max(ll, lh), max(hl, hh))));
}

static inline
MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs)
{
    return lhs * recip(rhs);
}

// fast multiplication
static inline
MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs)
{ return MachineEstimate(MachineEstimate::RoundDown(lhs.low * rhs),
                         MachineEstimate::RoundUp(lhs.high * rhs)); }
// and division
static inline
MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs)
{ return MachineEstimate(MachineEstimate::RoundDown(lhs.low / rhs),
                         MachineEstimate::RoundUp(lhs.high / rhs)); }

// shorthands
static inline 
MachineEstimate operator * (i32 lhs, const MachineEstimate &rhs)
{ return rhs * lhs; }

static inline 
MachineEstimate operator / (i32 lhs, const MachineEstimate &rhs)
{ return recip(rhs) * lhs; }

static inline MachineEstimate sqrt(const MachineEstimate &arg)
{
    double l = std::sqrt(arg.low);
    if (!(l >= 0))  // not the same as (l<0): unordered comparison (i.e. l >= 0 is false if l is a NaN, but l<0 is also false)
        l = 0;
    return MachineEstimate(MachineEstimate::RoundDown(l), 
                           MachineEstimate::RoundUp(std::sqrt(arg.high)));
}

static inline MachineEstimate rsqrt(const MachineEstimate &arg)
{
    return recip(sqrt(arg));
}

static inline MachineEstimate abs(const MachineEstimate &arg)
{
    if (arg.low >= 0) return arg;
    else if (arg.high <= 0) return -arg;
    else return MachineEstimate(0, max(arg.high, -arg.low));
}

static inline MachineEstimate sq(const MachineEstimate &arg)
{
    MachineEstimate a(abs(arg));
    return MachineEstimate(MachineEstimate::RoundDown(a.low * a.low),
                           MachineEstimate::RoundUp(a.high * a.high));
}

#ifndef TRUST_STDLIB
MachineEstimate cos(const MachineEstimate &arg);
MachineEstimate sin(const MachineEstimate &arg);
MachineEstimate log(const MachineEstimate &arg);
MachineEstimate exp(const MachineEstimate &arg);

MachineEstimate asin(const MachineEstimate &arg);
MachineEstimate atan(const MachineEstimate &arg);
MachineEstimate acos(const MachineEstimate &arg);
MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x);

#else

static inline MachineEstimate sin(const MachineEstimate &arg)
{
    double c = arg.weak_AsDouble();
    MachineEstimate d = arg.GetError();
    c = std::sin(c);
    return MachineEstimate(c).AddError(d);
}

static inline MachineEstimate cos(const MachineEstimate &arg)
{
    double c = arg.weak_AsDouble();
    MachineEstimate d = arg.GetError();
    c = std::cos(c);
    return MachineEstimate(c).AddError(d);
}

static inline MachineEstimate log(const MachineEstimate &arg)
{
    return MachineEstimate(MachineEstimate::RoundDown(std::log(arg.low)),
                           MachineEstimate::RoundUp(std::log(arg.high)));
}

static inline MachineEstimate exp(const MachineEstimate &arg)
{
    return MachineEstimate(MachineEstimate::RoundDown(std::exp(arg.low)),
                           MachineEstimate::RoundUp(std::exp(arg.high)));
}

static inline MachineEstimate asin(const MachineEstimate &arg)
{
    MachineEstimate a(arg);
    if (a.low < -1.0) a.low = -1.0;
    if (a.high > 1.0) a.high = 1.0;
    return MachineEstimate(MachineEstimate::RoundDown(std::asin(a.low)),
                           MachineEstimate::RoundUp(std::asin(a.high)));
}
#endif

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

} // namespace

#endif // FILE
