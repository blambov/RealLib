/*

    RealEstimate.h

    Error math and value+error container.
    Classes:

        Estimate - combines a LongFloat with its ErrorEstimate.
        the error is absolute: |value - real| < error.
    ErrorEstimate - a simple evaluation of error as 32-bit 
        mantissa and 32-bit exponent. Operations give results
        that are always greater than or equal to the actual
        result. Note the exponent in a ErrorEstimate is in 
        bits, not in words as in a LongFloat, and the mantissa
        always has 1 in its most significant bit.

 */

#ifndef FILE_REAL_ESTIMATE_H
#define FILE_REAL_ESTIMATE_H

#include <stdlib.h>
#include <limits.h>
#include <exception>
#include "defs.h"
#include "LongFloat.h"
#include "ErrorEstimate.h"
#include "RealExceptions.h"

namespace RealLib {

// class Estimate's definitions start here
class Estimate;

// operations
Estimate operator - (const Estimate &arg);
Estimate recip(const Estimate &arg);

Estimate operator + (const Estimate &lhs, const Estimate &rhs);
Estimate operator - (const Estimate &lhs, const Estimate &rhs);
Estimate operator * (const Estimate &lhs, const Estimate &rhs);
Estimate operator / (const Estimate &lhs, const Estimate &rhs);

// fast multiplication
Estimate operator * (const Estimate &lhs, i32 rhs);
// and division
Estimate operator / (const Estimate &lhs, i32 rhs);

static inline
std::ostream& operator <<(std::ostream &os, const Estimate &e);

class Estimate {
private:
    LongFloat m_Value;
    ErrorEstimate m_Error;        // error, |m_Value - real| < m_Error

    Estimate(const LongFloat &val, const ErrorEstimate &err = ErrorEstimate());
    void CorrectZero();         // 0 should not be used in calculations. Substitute (0, e) with (e, 2e)
    // and (0, 0) with (2^-MAXINT, 2^-MAXINT * 2)

public:
    Estimate(double v = 0.0);
    //    Estimate(const Estimate &rhs)
    //        : m_Value(rhs.m_Value), m_Error(rhs.m_Error) {}
    Estimate(const char *val);

    // error functions
    Estimate GetError() const;
    Estimate& SetError(const Estimate &err);
    Estimate& AddError(const Estimate &err);


    // a lower bound on the correct binary digits
    // uses the exponents of the value and error to calculate it quickly
    i32 GetRelativeError() const;

    /*
     Estimate& AddRoundingError() 
     { m_Error = m_Error + RoundingError(m_Value, m_Value.AdditionRoundingError()); 
         return *this; }

    Estimate TheRoundingError() const // rounding error is assumed to be no more than
                                                                            // one in the least significant bit of the mantissa
                                                                            // Note! Newton-Raphson reciprocal is incorrect in the
                                                                            // least significant word (handled by recip())
     { return RoundingError(m_Value); }
     */

    // get a rough estimate of the precision
    // used to determine the length of the approximations to functions
    u32 GetPrecision() const
    { return m_Value.GetPrecision(); }
    Estimate& SetPrecision(u32 prec)
    { m_Value.SetPrecision(prec);
    return *this; }


    // truncation
    // used to make sure only arguments within the domain of the function
    // are processed for the closed ends of the domain. 
    //    To this end, truncates the approximation interval so that
    // the indicated real numbers are thrown out. If nothing remains,
    // raise a DomainException(origin). 

    // warning: an error in the approximation of the bound will be added to the
    // error in the end result, i.e. if (center 0, error 3) is truncated below (c 1, e 0.5), the
    // result will be (c 2, e 1.5) (i.e. the interval [0.5, 3.5]). 
    // To avoid problems, use double arguments

    // removes the part of the approximation interval that is negative
    Estimate TruncateNegative(const char *origin = "Truncate") const;

    // removes the part of the approximation that is below a certain lower bound
    Estimate TruncateBelow(const Estimate &l, const char *origin = "Truncate") const
    { return (*this - l).TruncateNegative(origin) + l; }

    // removes the part of the approximation that is above a certain upper bound
    Estimate TruncateAbove(const Estimate &h, const char *origin = "Truncate") const
    { return h - (h - *this).TruncateNegative(origin); }

    // removes the part of the approximation outside the specified interval
    Estimate TruncateTo(double l, double h, const char *origin = "Truncate") const
    { return (h - ((h-l) - (*this - l).TruncateNegative(origin)).TruncateNegative(origin)); }
    Estimate TruncateTo(const Estimate &l, const Estimate &h, const char *origin = "Truncate") const
    { return (h - ((h-l) - (*this - l).TruncateNegative(origin)).TruncateNegative(origin)); }

    // comparisons
    // these come in two flavors, strong (true if real is in relation to rhs)
    bool IsPositive() const;
    bool IsNegative() const;
    bool IsNonZero() const;

    // equality test is undecidable (i.e. would yield false for any precision)
    // thus ==, <= and >= are not included
    // also !(x<y) does not mean y<=x
    bool operator < (const Estimate &rhs) const
    { return (*this - rhs).IsNegative(); }
    bool operator > (const Estimate &rhs) const
    { return (*this - rhs).IsPositive(); }
    bool operator != (const Estimate &rhs) const
                    { return (*this - rhs).IsNonZero(); }

    // and weak (true if m_Value is in relation to rhs)
    // should only be used if the transformation being aplied
    // would not differentiate on the two cases, e.g. to choose
    // whether to evaluate sin(x) and sin(pi - x)

    bool weak_IsPositive() const;
    bool weak_IsNegative() const;
    bool weak_IsNonZero() const
    { return true; }
    // an Estimate cannot be weakly zero -- see the remark for CorrectZero()

    bool weak_lt(const Estimate &rhs) const;
    bool weak_eq(const Estimate &rhs) const;

    bool weak_gt(const Estimate &rhs) const
    { return rhs.weak_lt(*this); }

    bool weak_le(const Estimate &rhs) const
    { return !weak_gt(rhs); }
    bool weak_ne(const Estimate &rhs) const
    { return !weak_eq(rhs); }
    bool weak_ge(const Estimate &rhs) const
    { return !weak_lt(rhs); }

    // among the weak operations is also rounding
    // the returned Estimate is assumed exact
    // only to be used on periodic functions!
    Estimate weak_round() const;

    // weak normalize, i.e. return an exponent such that
    // a >> a.weak_normalize()
    // is in the range [0.5, 1).
    i32 weak_normalize() const
    { return m_Value.normalize(); }

    // weak conversion
    double weak_AsDouble() const
    { return m_Value.AsDouble(); }
    // output
    char *weak_AsDecimal(char *buffer, u32 buflen) const
    { return m_Value.AsDecimal(buffer, buflen); }

    Estimate weak_Center()
    { return Estimate(m_Value, ErrorEstimate()); }

    /*
     // exponent and mantissa operations
     // needed to get initial approximations via double
     double weak_MantissaAsDouble() const
            { return m_Value.MantissaAsDouble(); }

     i32 weak_Exponent() const
            { return m_Value.Exponent(); }
     Estimate& AddToExponent(i32 exp)
            { m_Value.AddToExponent(exp); 
                return *this; }
     char *weak_MantissaAsDecimal(char *buf, u32 buflen) const
            { return m_Value.MantissaAsDecimal(buf, buflen); }
     */
    // operations
    friend Estimate operator - (const Estimate &arg);
    friend Estimate recip(const Estimate &arg);

    friend Estimate operator + (const Estimate &lhs, const Estimate &rhs);
    friend Estimate operator - (const Estimate &lhs, const Estimate &rhs);
    friend Estimate operator * (const Estimate &lhs, const Estimate &rhs);
    friend Estimate operator / (const Estimate &lhs, const Estimate &rhs);

    // fast multiplication
    friend Estimate operator * (const Estimate &lhs, i32 rhs);
    // and division
    friend Estimate operator / (const Estimate &lhs, i32 rhs);

    // binary shift
    Estimate operator << (i32 howmuch) const;
    Estimate operator >> (i32 howmuch) const
    { return *this << -howmuch; }

    Estimate& operator += (const Estimate &rhs)
                 { return *this = *this + rhs; }
    Estimate& operator -= (const Estimate &rhs)
                 { return *this = *this - rhs; }
    Estimate& operator *= (const Estimate &rhs)
                 { return *this = *this * rhs; }
    Estimate& operator /= (const Estimate &rhs)
                 { return *this = *this / rhs; }


    Estimate& operator >>= (i32 rhs) 
                 { return *this = *this >> rhs; }
    Estimate& operator <<= (i32 rhs) 
                 { return *this = *this << rhs; }
    Estimate& operator *= (i32 rhs) 
                 { return *this = *this * rhs; }
    Estimate& operator /= (i32 rhs) 
                 { return *this = *this / rhs; }

    // should probably be somewhere else
    // conversion to string
    // char *AsDecimal(char *buffer, u32 buflen);
    friend    
    std::ostream& operator <<(std::ostream &os, const Estimate &e);

};

// shorthands
static inline 
Estimate operator * (i32 lhs, const Estimate &rhs)
{ return rhs * lhs; }

static inline 
Estimate operator / (i32 lhs, const Estimate &rhs)
{ return recip(rhs) * lhs; }

// C++-style output
static inline
std::ostream& operator <<(std::ostream &os, const Estimate &e)
{ return os << e.m_Value; }



}    // namespace

#endif
