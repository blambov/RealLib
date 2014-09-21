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

  LongFloat.h

  Arbitrary precision math.
  Classes:

    LongFloat - arbitrary precision floating point value.
        Consists of a 32-bit exponent, long mantissa, sign and
        special flag. Special values are Zero, Infinity and Nan.
        The exponent is in words, and the most significant word
        is never zero. The mantissa is little endian, i.e.
        m_Man[0] is the least significant word.
        Note: a LongFloat will never change its mantissa, 
        because it can be referenced more than once (for faster
        assignments and copying). All non-const operations
        create a new LongFloat and assign *this to it.
    Mantissa - wrapper around the memory block that stores the
        mantissa. Blocks for mantissas are pre-allocated to
        avoid frequent calls to malloc and free.


*/

#ifndef FILE_LONG_FLOAT_H
#define FILE_LONG_FLOAT_H

#include <iosfwd>

#include "defs.h"
#include "limits.h"

namespace RealLib {

// current working precision
extern i32 g_WorkingPrecision;

// initialize must be called before any work
bool InitializeLongFloat(i32 precision, u32 initialbuffersize = 100);

// clear the buffers and check for unreleased data
void FinalizeLongFloat();

#ifndef max
static inline i32 max(i32 a, i32 b) { return a > b ? a : b; }
static inline u32 max(u32 a, u32 b) { return a > b ? a : b; }
static inline i32 min(i32 a, i32 b) { return a < b ? a : b; }
static inline u32 min(u32 a, u32 b) { return a < b ? a : b; }
#endif

static inline i32 i32_saturated(i64 a)
{ a = a < I32_MAX ? a : I32_MAX;
return i32(a > -I32_MAX ? a : -I32_MAX); }

typedef u32 Alloc;

class Mantissa {
private:
    Alloc m_pAlloc;

public:
    // create a new mantissa, which is referenced exactly once
    Mantissa();
    // copy a mantissa
    Mantissa(Alloc pAlloc);
    Mantissa(const Mantissa &rhs);
    // dereference mantissa
    ~Mantissa();

    // destruct, then copy
    void operator = (const Mantissa &rhs);

    // non-const reference to a word for assignment
    u32& operator [] (i32 index);
    // const operation, extract value of index-th word
    u32 operator () (i32 index) const;

    // non const pointer get
    u32* getWritablePtr();
    // const pointer get
    const u32* getConstPtr() const;
};

// LongFloat
// the representation is little endian
// (i.e. -1 ^ m_Neg * sum(m_Man[i] * 2^ (32 * (m_Exp - g_WorkingPrecision + i))))

class LongFloat {
public:
    // special values
    typedef enum { Normal, Zero, Infinity, Nan } Special;
private:
    // mantissa
    Mantissa m_Man; 
    // exponent, power of 2^32
    i32 m_Exp;
    // sign, true if negative
    bool m_Neg;
    // != Normal only in special values
    Special m_Special;
    // needed precision. this is the point at which multiplication and division will stop
    u32 m_Prec;

    #ifndef NDEBUG
    // to facilitate debugging
    double doubleval;
    #endif

    // private constructor, directly setting members
    LongFloat(Special spec, bool neg, const Mantissa &man, exp_type exp, u32 Prec = g_WorkingPrecision);

public:
    // constructor of special values
    LongFloat(Special spec = Zero, bool neg = false);
    // copy constructor
    LongFloat(const LongFloat &src);
    // destructor
    ~LongFloat();

    // assignment
    LongFloat& operator = (const LongFloat &rhs);

    // conversions, exact
    LongFloat(const double src);
    LongFloat(const i32 man, const i32 exp);    // set to man * 2^(32*exp)
    // set to closest representable value
    LongFloat(const char *value);

    // return the nearest double value
    double AsDouble() const;
   // output
   char *AsDecimal(char *buffer, u32 buflen) const;

    // mantissa conversion. Don't care about any other attribute
    LongFloat MantissaAsLongFloat() const;
    double MantissaAsDouble() const;
    // does not add the leading '.', returns true if the rounding resulted in carry, 
    // i.e. the mantissa rounds to 1.(0)
    bool MantissaAsDecimal(char *buffer, u32 buflen) const; 

    LongFloat SignedMantissaAsLongFloat() const;

    // attributes retrieval
    Special Kind() const
    { return m_Special; }
    i32 Exponent() const
    { return m_Exp; }
    bool IsNegative() const
    { return m_Neg; }
    
    i32 GetPrecision() const
       { return m_Prec; }
    void SetPrecision(i32 prec)
       { m_Prec = prec < g_WorkingPrecision ? prec : g_WorkingPrecision; }

    // at LongFloat level we're directly operating Mantissa words
    // Mantissas can be shared between different LongFloats, thus
    // all operations produce new LongFloats rather than modify existing
    
    LongFloat operator - () const;

    LongFloat operator + (const LongFloat &rhs) const;
    LongFloat operator - (const LongFloat &rhs) const;
    LongFloat operator * (const LongFloat &rhs) const;
    LongFloat operator / (const LongFloat &rhs) const;

    LongFloat recip() const;
    
   // rounding error functions. needed to minimize the
   // dependancies of Estimate to the actual LongFloat
   // implementation
   // return the index of the first possibly incorrect bit 
   // (provided the inputs were correct)
    i32 AdditionRoundingError() const
    { return -g_WorkingPrecision*32-1; }
    i32 MultiplicationRoundingError() const
    { return -GetPrecision()*32-1; }
    i32 DivisionRoundingError() const;
                   // Newton-Raphson may be complete wrong in the lsw

    bool operator == (const LongFloat &rhs) const; 
    bool operator >= (const LongFloat &rhs) const; 

    LongFloat operator << (i32 howmuch) const;  // scale binary

    // rounding
    LongFloat RoundTowardZero() const;
    LongFloat round() const;    // to nearest integer
    
    // normalization (returns an exponent that would
   // set the mantissa in the range [0.5; 1) )
    i32 normalize() const;

    // fast stuff
    LongFloat operator * (i32 rhs) const;
    LongFloat operator / (i32 rhs) const;
    LongFloat addProduct(const LongFloat &a, const LongFloat &b) const;
    LongFloat addProduct(const LongFloat &a, i32 b) const;

    LongFloat& AddToExponent(i32 howmuch);

    // shorthands; all of these provide no performance benefit

    bool operator != (const LongFloat &rhs) const
    { return !(*this == rhs); }
    bool operator <= (const LongFloat &rhs) const
    { return rhs >= *this; }
    bool operator < (const LongFloat &rhs) const
    { return !(*this >= rhs); }
    bool operator > (const LongFloat &rhs) const
    { return !(rhs >= *this); }

    LongFloat operator >> (i32 howmuch) const
    { return *this << (-howmuch); }
    LongFloat& operator <<= (i32 howmuch)
    { return *this = *this << (howmuch); }
    LongFloat& operator >>= (i32 howmuch)
    { return *this = *this << (-howmuch); }

    LongFloat& operator += (const LongFloat &rhs) 
    { return *this = *this + rhs; }
    LongFloat& operator -= (const LongFloat &rhs)
    { return *this = *this - rhs; }
    LongFloat& operator *= (const LongFloat &rhs)
    { return *this = *this * rhs; }
    LongFloat& operator /= (const LongFloat &rhs)
    { return *this = *this / rhs; }

    LongFloat& operator *= (i32 rhs)
    { return *this = *this * rhs; }


};

// C++-style output
std::ostream& operator <<(std::ostream &os, const LongFloat &lf);

} // namespace

#endif
