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

#include "defs.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "ErrorEstimate.h"

namespace RealLib {

// ErrorEstimate implementation

ErrorEstimate::ErrorEstimate(const u32 man, const i32 exp)
: m_Man(man), m_Exp(exp)
{
    if (m_Man==0) m_Exp = minusinf;
    else while ((m_Man & 0x80000000) == 0) { m_Man <<= 1; m_Exp -= 1; }
}

// conversion
ErrorEstimate::ErrorEstimate(const LongFloat &src, RoundingMode round)
{
    switch (src.Kind()) {
    case LongFloat::Normal: {
        // convert mantissa to double to extract the 32 most significant bits
        exp_type exp = exp_type(src.Exponent()) * 32;
        double man = src.MantissaAsDouble();
        int mexp;
        man = ldexp(frexp(fabs(man), &mexp), 32) + round;
        // offset the exponent with the double's exponent
        exp += mexp;

        u32 m = u32(man);

        // we might have overflowed
        if (m == 0) {
            m = 1u<<31;
            ++exp;
        }

        m_Exp = i32_saturated(exp);
        m_Man = m;
        break; }
    case LongFloat::Zero:
        m_Man = 0;
        m_Exp = minusinf;
        break;
    default:
        m_Man = 0;
        m_Exp = plusinf;
        break;
    }
}

ErrorEstimate::ErrorEstimate(const double man)
{
    switch (_fpclass(man)) {
    case _FPCLASS_PD:
    case _FPCLASS_ND:
    case _FPCLASS_PN:
    case _FPCLASS_NN: {
        int mexp;
        // split mantissa and exponent
        m_Man = u32(ldexp(frexp(fabs(man), &mexp), 32) + 1);
        // correct for possible overflow
        if (m_Man == 0) {
            m_Man = 1u<<31;
            m_Exp = mexp + 1;
        } else m_Exp = mexp;
        break; }
    case _FPCLASS_PZ:
    case _FPCLASS_NZ:
        m_Man = 0;
        m_Exp = minusinf;
        break;
    default:    // nan, inf
        m_Man = 0;
        m_Exp = plusinf;
        break;
    }
}

// helpers for ErrorEstimate

// add mantissas of ErrorEstimates, round up
bool DoEEManAdd(u32 &man, u32 full, u32 part, u32 start)
{
    if (start >= 32) man = ++full;
    else {
        // see if we need to round up
        if (part & ((1<<start)-1)) 
            part = (part >> start) + 1;
        else part = (part >> start);

        man = full + part;
    }
    return (man < full);    // tricky check for carry
}

// sub mantissas of ErrorEstimates, round down
exp_type DoEEManSub(u32 &man, u32 full, u32 part, exp_type start, exp_type exp)
{
    // substract
    if (start >= 32) man = --full;
    else {
        if (start > 0)
            if (part & ((1<<start)-1)) 
                part = (part >> start) + 1;
            else part = (part >> start);

        man = full - part;
    }
    assert (man <= full);

    // normalize
    if (man == 0) return ErrorEstimate::minusinf;
    while ((man & (1u<<31))==0) {
        --exp;
        man <<= 1;
    }
    return exp;
}

ErrorEstimate ErrorEstimate::operator + (const ErrorEstimate &rhs) const
{
    u32 man;
    exp_type exp;
    bool carry;

    // handle special cases
    if (m_Exp == plusinf || rhs.m_Exp == minusinf) return *this;
    if (rhs.m_Exp == plusinf || m_Exp == minusinf) return rhs;

    // do addition
    if (m_Exp == rhs.m_Exp) {
        exp = m_Exp;
        man = m_Man + rhs.m_Man;
        carry = true;
    } else if (m_Exp > rhs.m_Exp) {
        exp = m_Exp;
        carry = DoEEManAdd(man, m_Man, rhs.m_Man, m_Exp - rhs.m_Exp);
    } else {
        exp = rhs.m_Exp;
        carry = DoEEManAdd(man, rhs.m_Man, m_Man, rhs.m_Exp - m_Exp);
    }

    // update if carry
    if (carry) {
        if (man & 1) ++man;
        man = (man >> 1) | (1u << 31);
        exp = exp + 1;
    }

    return ErrorEstimate(man, i32_saturated(exp));
}

ErrorEstimate ErrorEstimate::operator - (const ErrorEstimate &rhs) const
{
    u32 man;
    exp_type exp;

    // handle specials
    if (m_Exp == plusinf || rhs.m_Exp == minusinf) return *this;
    if (rhs.m_Exp == plusinf || m_Exp == minusinf) return rhs;

    // errors are always positive, a negative result would mean error
    assert(m_Exp >= rhs.m_Exp);
    exp = DoEEManSub(man, m_Man, rhs.m_Man, m_Exp - rhs.m_Exp, m_Exp);

    return ErrorEstimate(man, i32_saturated(exp));
}

ErrorEstimate ErrorEstimate::operator * (const ErrorEstimate &rhs) const
{
    exp_type e = exp_type(m_Exp) + exp_type(rhs.m_Exp) - 1;

    // handle overflow and special cases
    if (m_Exp == plusinf || rhs.m_Exp == plusinf || e >= plusinf) return ErrorEstimate(0, plusinf);
    if (m_Exp == minusinf || rhs.m_Exp == minusinf || e <= minusinf) return ErrorEstimate(0, minusinf);

    // multiply. the result will at least have 1 in 62nd position
    // at most 1 in 63rd
    u64 m = u64(m_Man) * u64(rhs.m_Man);
    // round up if necessary
    if (u32((m<<1) & 0xFFFFFFFF)) m = (m >> 31) + 1;
    else m = (m >> 31);

    // make room for the 63rd bit if it is not 0
    if (m>>32) {
        if (m & 1) ++m;
        m >>= 1;
        ++e;
    }

    return ErrorEstimate(u32(m), i32(e));
}

ErrorEstimate ErrorEstimate::operator << (i32 howmuch) const
{ 
    // simply add to exponent saturating
    exp_type e = exp_type(m_Exp) + exp_type(howmuch);

    //if (m_Exp == plusinf || e >= plusinf) return ErrorEstimate(0, plusinf);
    //if (m_Exp == minusinf || e <= minusinf) return ErrorEstimate(0, minusinf);

    return ErrorEstimate(m_Man, i32_saturated(e));
}


ErrorEstimate ErrorEstimate::recip() const
{
    if (m_Exp == plusinf) return ErrorEstimate(0, minusinf);
    if (m_Exp <= minusinf + 2) return ErrorEstimate(0, plusinf);

    // calculate
    i32 exp = -(m_Exp - 2);
    u32 man = u32(((u64(1) << 62) + m_Man - 1) / m_Man);

    // normalize
    if (!(man & (1u<<31))) {
        man = (man << 1) + 1;
        --exp;
    }

    return ErrorEstimate(man, exp);
}

ErrorEstimate& ErrorEstimate::operator ++ ()
{
    if (++m_Man == 0) {
        m_Man = 1u<<31;
        ++m_Exp;
    }
    return *this;
}

double ErrorEstimate::AsDouble() const
{
    // this can very easily over or underflow
    return ldexp(m_Exp == plusinf ? 1.0 : m_Man, m_Exp - 32);
}

LongFloat ErrorEstimate::AsLongFloat() const
{
    // m_Man is u32, so we can't use LongFloat's i32 constructor
    // go through double, but separate the exponent to a remainder
    // which can safely be applied to a double, and a quotient
    // which can very quickly be added to a LongFloat's exponent.
    i32 expq = m_Exp >> 5;
    i32 expr = m_Exp - (expq << 5);

    return LongFloat(ldexp(double(m_Man), expr)).AddToExponent(expq - 1);
}

bool ErrorEstimate::operator >= (const ErrorEstimate &rhs) const
{
    if (m_Exp > rhs.m_Exp) return true;
    else if (m_Exp == rhs.m_Exp && 
            (m_Exp == plusinf || m_Exp == minusinf || m_Man >= rhs.m_Man))
        return true;
    else return false;
}

bool ErrorEstimate::operator > (const ErrorEstimate &rhs) const
{
    if (m_Exp > rhs.m_Exp) return true;
    else if (m_Exp == rhs.m_Exp && 
            (m_Exp != plusinf && m_Exp != minusinf && m_Man > rhs.m_Man))
        return true;
    else return false;
}

ErrorEstimate max(const ErrorEstimate &a, const ErrorEstimate &b) 
{
    return a >= b ? a : b;
}

ErrorEstimate min(const ErrorEstimate &a, const ErrorEstimate &b) 
{
    return a >= b ? b : a;
}

ErrorEstimate RoundingError(const LongFloat &lf, i32 re) 
// rounding error is assumed to be no more than
// one in the least significant bit of the mantissa
// Note! Newton-Raphson reciprocal is incorrect in the
// least significant word (handled by recip())
{ 
    exp_type exp = exp_type(lf.Exponent()) * 32 + re + 1;
    if (exp <= ErrorEstimate::minusinf) return ErrorEstimate(0, ErrorEstimate::minusinf);
    if (exp >= ErrorEstimate::plusinf) return ErrorEstimate(0, ErrorEstimate::plusinf);
    else return ErrorEstimate(1u<<31, (lf.Exponent()) * 32 + re + 1); 
}

} // namespace
