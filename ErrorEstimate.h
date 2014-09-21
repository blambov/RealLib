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

/*

  ErrorEstimate.h

  Error math.
  Classes:

  ErrorEstimate - a simple evaluation of error as 32-bit 
        mantissa and 32-bit exponent. Operations give results
        that are always greater than or equal to the actual
        result. Note the exponent in a ErrorEstimate is in 
        bits, not in words as in a LongFloat, and the mantissa
        always has 1 in its most significant bit.

*/

#ifndef FILE_ERROR_ESTIMATE_H
#define FILE_ERROR_ESTIMATE_H

#include <stdlib.h>
#include <limits.h>
#include <exception>
#include "defs.h"
#include "LongFloat.h"

#define I32_MIN INT_MIN     // to be used as -inf
#define I32_MAX INT_MAX     // to be used as +inf

namespace RealLib {

class ErrorEstimate {
public:
    u32 m_Man;      // >= 2^31
    i32 m_Exp;      // error = m_Man * 2 ^ (m_Exp - 32)

    enum Inf { minusinf = -I32_MAX, plusinf = I32_MAX } ;
    enum RoundingMode { Down = 0, Up = 1 };

    explicit ErrorEstimate(const u32 man = 0, const i32 exp = 0);
    // conversions (rounding up)
    ErrorEstimate(const double err);
    ErrorEstimate(const LongFloat &src, RoundingMode rnd = Up);

    ErrorEstimate RoundDownLongFloat(const LongFloat &src);

    // destructor, operator = not needed

    // operations
    // round-up addition
    ErrorEstimate operator + (const ErrorEstimate &rhs) const;
    // round-down substraction
    ErrorEstimate operator - (const ErrorEstimate &rhs) const;
    // round-up multiplication
    ErrorEstimate operator * (const ErrorEstimate &rhs) const;
    // round-up reciprocal
    ErrorEstimate recip() const;
    // round-up division
    ErrorEstimate operator / (const ErrorEstimate &rhs) const
       { return *this * rhs.recip(); }
    // round-up <<
    ErrorEstimate operator << (i32 howmuch) const;  // shift left in bits

   // comparisons
    bool operator >= (const ErrorEstimate &rhs) const;

    ErrorEstimate& operator ++ ();  // add the minimum

    // conversions
    double AsDouble() const;
    LongFloat AsLongFloat() const;
    
};

// two needed operations
ErrorEstimate max(const ErrorEstimate &a, const ErrorEstimate &b); 
ErrorEstimate min(const ErrorEstimate &a, const ErrorEstimate &b);

ErrorEstimate RoundingError(const LongFloat &lf, i32 re);

} // namespace

#endif
