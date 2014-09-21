/*

  RealLib, a library for efficient exact real computation
  Copyright (C) 2006  Branimir Lambov <barnie@brics.dk>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/ 

#include "defs.h"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "RealEstimate.h"
#include "DataManager.h"
#include "LongFloat.h"
#include "ErrorEstimate.h"

namespace RealLib {

RealLibException::RealLibException(const char *wht) throw()
{ 
     if (wht) { 
        strncpy(m_what, wht, 127);
        m_what[127] = 0;
     } else m_what[0] = 0;
} 


// Estimate implementation

Estimate::Estimate(const LongFloat &val, const ErrorEstimate &err)
: m_Value(val), m_Error(err)
{
    CorrectZero();
}

Estimate::Estimate(double v)
: m_Value(v), m_Error(ErrorEstimate())
{
    CorrectZero();
}

Estimate::Estimate(const char *val)
: m_Value(val), m_Error(RoundingError(m_Value, m_Value.DivisionRoundingError()))
{
    CorrectZero();
}

void Estimate::CorrectZero()
{
    if (m_Value.Kind() == LongFloat::Zero) 
        if (m_Error.m_Exp == ErrorEstimate::minusinf) 
        {
            m_Error.m_Exp = MINIMUM_EXPONENT+1;
            m_Error.m_Man = 1<<31;
            m_Value = m_Error.AsLongFloat() >> 1;
        }
        else
        {
            m_Value = m_Error.AsLongFloat();
            m_Error = m_Error << 1;
        }
}

Estimate Estimate::GetError() const
{
   return Estimate(m_Error.AsLongFloat());
}

Estimate& Estimate::SetError(const Estimate &err)
{
   ErrorEstimate ee(err.m_Value);
   m_Error = ee + err.m_Error;
   return *this;
}

Estimate& Estimate::AddError(const Estimate &err)
{
   m_Error = m_Error + ErrorEstimate(err.m_Value) + err.m_Error;
   return *this;
}

i32 Estimate::GetRelativeError() const
{
   exp_type ee = m_Error.m_Exp;
   exp_type ev = weak_normalize();

   return i32_saturated(ev - ee - 1);

}

Estimate Estimate::TruncateNegative(const char *origin) const
{
    if (IsNegative()) throw DomainException(origin);

    if (IsPositive()) return *this;

    // (the theory says
   // we can't always give the correct DomainException, so we shouldn't try)    

    //Estimate a = (*this + GetError()) >> 1;

    // a little magic here:
    // SetError sets an error of the higher bound of the interval, i.e. 
    // a + a.GetError()
    Estimate a(*this >> 1);
    // with this we make an Estimate with center (*this + *this.GetError())/2
    // and the same distance 
    return (a + a.GetError()).SetError(a);

   //return a.SetError(a);
}

// operations
Estimate operator - (const Estimate &arg)
{
    return Estimate(-arg.m_Value, arg.m_Error);
}

Estimate recip(const Estimate &arg)
{
   if (!arg.IsNonZero()) throw PrecisionException("recip");
    
    LongFloat r(arg.m_Value.recip());
    ErrorEstimate e(arg.m_Value, ErrorEstimate::Down);
    ErrorEstimate re(RoundingError(r, r.DivisionRoundingError()));
    return Estimate(r, (arg.m_Error / (e - arg.m_Error) / e) + re); 
                   // multiplication in denominator would have the wrong rounding mode
}

Estimate operator + (const Estimate &lhs, const Estimate &rhs)
{
   LongFloat s(lhs.m_Value + rhs.m_Value);
   return Estimate(s, lhs.m_Error + rhs.m_Error + RoundingError(s, s.AdditionRoundingError()));
}

Estimate operator - (const Estimate &lhs, const Estimate &rhs)
{
   LongFloat s(lhs.m_Value - rhs.m_Value);
   return Estimate(s, lhs.m_Error + rhs.m_Error + RoundingError(s, s.AdditionRoundingError()));
}

Estimate operator * (const Estimate &lhs, const Estimate &rhs)
{
    LongFloat r(lhs.m_Value * rhs.m_Value);
    ErrorEstimate e(lhs.m_Error * rhs.m_Error + lhs.m_Error * ErrorEstimate(rhs.m_Value) + rhs.m_Error * ErrorEstimate(lhs.m_Value));

    return Estimate(r, e + RoundingError(r, r.MultiplicationRoundingError()));
}

Estimate operator / (const Estimate &lhs, const Estimate &rhs)
{
   if (!rhs.IsNonZero()) throw PrecisionException("division");
                   // this also assures e - rhs.m_Error > 0
    
    LongFloat r(lhs.m_Value / rhs.m_Value);
    ErrorEstimate e(rhs.m_Value, ErrorEstimate::Down);
    ErrorEstimate n(ErrorEstimate(lhs.m_Value) * rhs.m_Error + ErrorEstimate(rhs.m_Value, ErrorEstimate::Up) * lhs.m_Error);
    return Estimate(r, n / (e - rhs.m_Error) / e + RoundingError(r, r.DivisionRoundingError())); 
                   // multiplication in denominator would have the wrong rounding mode
}

Estimate Estimate::operator << (i32 howmuch) const
{
   LongFloat v(m_Value << howmuch);
    return Estimate(v, (m_Error << howmuch) + RoundingError(v, v.AdditionRoundingError()));
}

Estimate operator * (const Estimate &lhs, i32 rhs)
{
    LongFloat r(lhs.m_Value * rhs);

    return Estimate(r, lhs.m_Error * ErrorEstimate(double(rhs)) + RoundingError(r, r.AdditionRoundingError()));
}

Estimate operator / (const Estimate &lhs, i32 rhs)
{ 
   if (rhs == 0) throw DomainException("division by int");
   LongFloat r(lhs.m_Value / rhs);
   
   return Estimate(r, lhs.m_Error / ErrorEstimate(double(rhs)) + RoundingError(r, r.AdditionRoundingError()));
}   

bool Estimate::IsPositive() const
{
   return !m_Value.IsNegative() && ErrorEstimate(m_Value) >= m_Error;
}

bool Estimate::IsNegative() const
{
   return m_Value.IsNegative() && ErrorEstimate(m_Value) >= m_Error;
}

bool Estimate::IsNonZero() const
{
   return ErrorEstimate(m_Value) >= m_Error;
}

bool Estimate::weak_IsNegative() const
{
   return m_Value.Kind() == LongFloat::Normal && m_Value.IsNegative();
}

bool Estimate::weak_IsPositive() const
{
   return m_Value.Kind() == LongFloat::Normal && !m_Value.IsNegative();
}

bool Estimate::weak_lt(const Estimate& rhs) const
{
   return m_Value < rhs.m_Value;
}

bool Estimate::weak_eq(const Estimate& rhs) const
{
   return m_Value == rhs.m_Value;
}

Estimate Estimate::weak_round() const
{
   return Estimate(m_Value.round());
}

/*
char* Estimate::AsDecimal(char *bufptr, u32 buflen)
{
   int prec = GetPrecision();
    // must at least accomodate exponent
    assert (buflen >= 10);
    char *buffer = bufptr;

    // the format chosen is "[-].<mantissa>e<+/-><exponent>"
    // where mantissa has a leading non-zero decimal

    if (m_Value.IsNegative()) {
        buffer++[0] = '-';
        --buflen;
    }

    // handle special values
    switch (m_Value.Kind()) {
    case LongFloat::Nan:
        strcpy(buffer, "NaN");
        return bufptr;
    case LongFloat::Infinity:
        strcpy(buffer, "Infinity");
        return bufptr;
    case LongFloat::Zero:
        strcpy(buffer, "Zero");
        return bufptr;
    }
    if (m_Error.gePow2(0)) {
        strcpy(buffer, "Zero");
        return bufptr;
    }

    // calculate exponent: the least power of 10 that is greater than
    // or equal to the value. Using decimal logarithm and its ceiling.
    Estimate pwr = ln(abs(*this)) / Constants.ln10.GetEstimate(prec);
    // the following operation is safe: the exponent is at most
    // 2^(32 * log(10, 2)) which fits in less than the 53 mantissa bits
    // of double
    int e = int(ceil(pwr.Value().AsDouble()));

    // GCC doesn't understand _itoa, use sprintf instead
    sprintf(buffer, "%+d", e);
    // now we know the length of the exponent. move it to the end of
    // the string
    int explen = strlen(buffer);
    strcpy(buffer + buflen - explen - 1, buffer);
    buffer[0] = '.';

    // divide the value by the exponent to form decimal mantissa
    Estimate div(pow(LongFloat(10), e));
    Estimate man = *this / div;

    // when the value is power of ten, we can get two possible
    // representations of the mantissa: 0.(9) and 1.(0). This
    // behavior is not an error as it is consistent with the
    // theory. The function used to convert the mantissa will
    // not display 1.(0) correctly. Thus we must handle the
    // case differently: just pretend it's 0.(9)
    if (man.Value().Exponent() > 0) {
        for (u32 i=1;i<buflen - explen - 2; ++i)
            buffer[i] = '9';
    } else 
        // convert it to decimal using LongFloat's function
        man.Value().MantissaAsDecimal(buffer+1, buflen - explen - 2);
    
    char *ptr = buffer + (buflen - explen - 2);
    *ptr = 'e';

    return bufptr;
}
*/

}   // namespace
