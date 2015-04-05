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

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>

#include <iostream>
#include <iomanip>

#if !defined(__GNUC__) && !defined(_MSC_VER)
#include <alloca.h>
#else
#define ios_base ios
#endif

#include "LongFloat.h"
#include "DataManager.h"

#include "kernels.h"

namespace RealLib {

// the data manager that will be used in the library
DataManager *g_pDataMan = 0;

// current working precision
i32 g_WorkingPrecision = 0;
// there should be an alloc to be assigned to all special values
Alloc g_pZeroAlloc = 0;

// initialization of the LongFloat library
// create a data manager, convolution object and
// buffers, a zero alloc
bool InitializeLongFloat(i32 precision, u32 bufsize)
{
    assert(sizeof (i32) == 4);
    assert(sizeof (i64) == 8);

    assert(!g_pDataMan);    // make sure finalize has been called before
                            // reinitialization
    assert(precision >= 3);     // at least cover the precision of double

    g_pDataMan = new DataManager(precision, bufsize);
    if (!g_pDataMan) return false;
    if (!g_pDataMan->isValid()) {
        delete g_pDataMan;
        return false;
    }
    g_WorkingPrecision = precision;

    // make a zero alloc
    g_pZeroAlloc = g_pDataMan->newAlloc();

    // and fill it with zeroes
    u32 *ZA = (*g_pDataMan)[g_pZeroAlloc];
    for (i32 u=0; u<precision; ++u)
        ZA[u] = 0;

    InitializeKernels(precision);
    return true;
}

// delete all globals
void FinalizeLongFloat()
{
    assert(g_pDataMan && g_pDataMan->isValid());
    g_pDataMan->releaseAlloc(g_pZeroAlloc);

    delete g_pDataMan;
    g_pDataMan = 0;

    FinalizeKernels();
}

// Mantissa implementation

// a new mantissa is a new block
inline Mantissa::Mantissa()
: m_pAlloc(g_pDataMan->newAlloc())
{
}

// all copying is done by simply increasing the block's reference count
inline Mantissa::Mantissa(Alloc pAlloc)
: m_pAlloc(pAlloc)
{ 
    g_pDataMan->referenceAlloc(m_pAlloc); 
}

inline Mantissa::Mantissa(const Mantissa &rhs)
: m_pAlloc(rhs.m_pAlloc)
{ 
    g_pDataMan->referenceAlloc(m_pAlloc); 
}

// don't delete at destruct, just dereference
inline Mantissa::~Mantissa() 
{ 
    g_pDataMan->releaseAlloc(m_pAlloc);
}

// combine destruction and copying
inline void Mantissa::operator = (const Mantissa &rhs)
{
    g_pDataMan->releaseAlloc(m_pAlloc); 
    m_pAlloc = rhs.m_pAlloc;
    g_pDataMan->referenceAlloc(m_pAlloc); 
}

// The following two functions have identical implementations.
// The difference is one is const, the other one is for writing.
// Using the non-const operation on a multiply referenced
// mantissa is an error.
inline u32& Mantissa::operator [] (i32 index)
{ 
    assert(g_pDataMan->AllocCanBeChanged(m_pAlloc));    

    assert(index >= 0 && index < g_WorkingPrecision);
    return (*g_pDataMan)[m_pAlloc][index]; 
}

inline u32 Mantissa::operator () (i32 index) const
{ 
    assert(index >= 0 && index < g_WorkingPrecision);
    return (*g_pDataMan)[m_pAlloc][index]; 
}

// getWritablePtr and conversion to const u32* are implemented to speed up some operations
inline u32* Mantissa::getWritablePtr()
{
    assert(g_pDataMan->AllocCanBeChanged(m_pAlloc));    

    return (*g_pDataMan)[m_pAlloc];
}

inline const u32* Mantissa::getConstPtr() const
{
    return (*g_pDataMan)[m_pAlloc];
}

// LongFloat implementation

LongFloat sq(const LongFloat &arg)
{
    return arg * arg;
}

// pow
// use pow(a, 2*p) == pow(a*a, p)
LongFloat pow(LongFloat arg, i32 pwr)
{
    bool neg = false;
    if (pwr < 0) {
        pwr = -pwr;
        neg = true;
    }
    LongFloat acc(1);

    while (pwr) {
        if (pwr & 1) acc *= arg;
        pwr >>= 1;
        arg = sq(arg);
    }

    return neg ? acc.recip() : acc;
}



LongFloat::LongFloat(Special s, bool neg, const Mantissa &man, exp_type exp, u32 prec)
: m_Man(man), m_Exp(i32_saturated(exp)), m_Neg(neg), m_Special(s), m_Prec(prec)
{
#ifndef NDEBUG
    // debugging is easier with this
    doubleval = AsDouble();
#endif
}

// specials have an all-zero mantissa (g_pZeroAlloc)
LongFloat::LongFloat(Special s, bool neg)
: m_Man(g_pZeroAlloc), m_Exp(0), m_Neg(neg), m_Special(s), m_Prec(g_WorkingPrecision)
{
#ifndef NDEBUG
    doubleval = AsDouble();
#endif
}

// copy constructor
LongFloat::LongFloat(const LongFloat &src)
: m_Man(src.m_Man), m_Exp(src.m_Exp), m_Neg(src.m_Neg), m_Special(src.m_Special), m_Prec(src.m_Prec)
{
#ifndef NDEBUG
    doubleval = AsDouble();
#endif
}

// destructor. If not explicitly written, C++ will inline it in all files
// that use LongFloat. It is preferable to have ~LongFloat called rather
// than ~Mantissa.
LongFloat::~LongFloat()
{
}

// the previous remark also applies to assignment
LongFloat& LongFloat::operator = (const LongFloat &rhs)
{
    m_Man = rhs.m_Man;
    m_Neg = rhs.m_Neg;
    m_Exp = rhs.m_Exp;
    m_Special = rhs.m_Special;
    m_Prec = rhs.m_Prec;

#ifndef NDEBUG
    doubleval = rhs.doubleval;
#endif

    return *this;
}

// initialize from int, man * (2^32)^exp
LongFloat::LongFloat(const i32 man, const i32 exp)
: m_Exp(exp + 1), m_Neg(man < 0), m_Special(Normal), m_Prec(g_WorkingPrecision)
{
    if (man == 0) {
        m_Special = Zero;
        m_Man = Mantissa(g_pZeroAlloc);
        m_Exp = 0;
    } else {
        m_Man[g_WorkingPrecision - 1] = ::abs(man);
        for (i32 i=0; i<g_WorkingPrecision-1; ++i)
            m_Man[i] = 0;
    }

#ifndef NDEBUG
    doubleval = AsDouble();
#endif
}

// initialize from double
LongFloat::LongFloat(const double val)
: m_Neg(false), m_Prec(g_WorkingPrecision)
{
    double man = val;
    int exp;

#ifndef NDEBUG
    doubleval = val;
#endif

    switch (_fpclass(val)) {
    case _FPCLASS_NN:   // negative normal
    case _FPCLASS_ND:   // negative denormal
        m_Neg = true;
        man = -man;
    case _FPCLASS_PN:   // positive normal
    case _FPCLASS_PD: { // positive denormal
        // split mantissa and exponent
        man = frexp(man, &exp);
        // our exponents are based on 2^32, correct for this
        int exp32 = (exp+31) >> 5;
        exp = exp - (exp32 << 5);   // division by 32 can give negative
        // remainder
        assert(exp <= 0 && exp > -32);

        // shift one word and save the integer part
        man = ldexp(man, 32 + exp);
        m_Man[g_WorkingPrecision - 1] = u32(man);
        // shift the fractional part one word and save the integer part
        man = ldexp(man - floor(man), 32);
        m_Man[g_WorkingPrecision - 2] = u32(man);
        // shift the fractional part one word and save the integer part
        man = ldexp(man - floor(man), 32);
        m_Man[g_WorkingPrecision - 3] = u32(man);
        // by now fraction is zero
        assert(man == floor(man));

        for (i32 i=0; i<g_WorkingPrecision - 3; ++i)
            m_Man[i] = 0;

        m_Exp = exp32;
        m_Special = Normal;

        break; }
    // all others are special values
    case _FPCLASS_NZ:   // negative zero
        m_Neg = true;
    case _FPCLASS_PZ:   // positive zero
        m_Special = Zero;
        m_Man = Mantissa(g_pZeroAlloc);
        m_Exp = 0;
        return;
    case _FPCLASS_SNAN:   // not a number
    case _FPCLASS_QNAN:
        m_Special = Nan;
        m_Man = Mantissa(g_pZeroAlloc);
        m_Exp = 0;
        return;
    case _FPCLASS_NINF:   // negative infinity
        m_Neg = true;
    case _FPCLASS_PINF:   // positive infinity
        m_Special = Infinity;
        m_Man = Mantissa(g_pZeroAlloc);
        m_Exp = 0;
        return;
    }

}

// initialize from string
LongFloat::LongFloat(const char *val)
{
    assert(val);
    while (val[0] == ' ' || val[0] == '\t') ++val;

    bool neg = false;
    if (val[0] == '-' || val[0] == '+') neg = (val++)[0] == '-';

    LongFloat t(Zero, false);
    exp_type expd = 0;

    while (val[0] >= '0' && val[0] <= '9') {
        t = LongFloat(i32(val[0] - '0'), 0).addProduct(t, 10);
        ++val;
    }

    // on fractions we proceed as usual, only remember how many
    // fractional digits we've processed
    if (val[0] == '.') {
        ++val;
        while (val[0] >= '0' && val[0] <= '9') {
            t = LongFloat(i32(val[0] - '0'), 0).addProduct(t, 10);
            expd -= 1;
            ++val;
        }
    }

    // add the exponent to the one we've gathered until now
    if (val[0] == 'e' || val[0] == 'E') {
        ++val;

        bool eneg = false;
        if (val[0] == '+' || val[0] == '-') eneg = (val++)[0] == '-';

        exp_type expt = 0;
        while (val[0] >= '0' && val[0] <= '9') {
            expt = expt * 10 + (val[0] - '0');
            ++val;
        }
        expd += eneg ? -expt : expt;
    }

    // multiply by the power of ten that is the exponent
    if (expd) 
        t *= pow(LongFloat(10), i32_saturated(expd)); 

    *this = neg ? -t : t;
}

// convert to double
double LongFloat::AsDouble () const
{
    double v = 0.0;
    switch (m_Special) {
    case Zero:
        break;
    case Infinity:
        v = 1.0 / v;
        break;
    case Nan:
        v = (1.0 / v) * 0.0;
        break;
    default:
    case Normal: {
        u32 prec = g_WorkingPrecision;
        if (m_Exp > -1024/32)
            if (m_Exp < 1024/32)
                v = ldexp((double)m_Man(prec - 1), m_Exp * 32 - 32) +
                ldexp((double)m_Man(prec - 2), m_Exp * 32 - 64) +
                ldexp((double)m_Man(prec - 3), m_Exp * 32 - 96);
            else v = 1.0/v;
        else v = 0.0;
        break; }
    }

    if (m_Neg) v = -v;
    return v;
}

LongFloat LongFloat::MantissaAsLongFloat() const
{
    // we need to copy the special, otherwise an error
    // might occur for zeroes
    return LongFloat(m_Special, false, m_Man, 0, m_Prec);
}

LongFloat LongFloat::SignedMantissaAsLongFloat() const
{
    // we need to copy the special, otherwise an error
    // might occur for zeroes
    return LongFloat(m_Special, m_Neg, m_Man, 0, m_Prec);
}

double LongFloat::MantissaAsDouble() const
{
    u32 prec = g_WorkingPrecision;
    return ldexp((double)m_Man(prec - 1), - 32) +
            ldexp((double)m_Man(prec - 2), - 64) +
            ldexp((double)m_Man(prec - 3), - 96);
}

// MantissaAsDecimal: converts the fraction to a decimal string
bool LongFloat::MantissaAsDecimal(char *buffer, u32 buflen) const
{
    Mantissa man;
    u32 *manptr = man.getWritablePtr();
    const u32 *m_manptr = m_Man.getConstPtr();
    i32 prec = g_WorkingPrecision - 1;

    // multiply by ten and output the carry
    buffer[0] = char(ScaleMantissa(manptr, m_manptr, 10)) + '0';

    i32 i;
    for (i=1; i<i32(buflen)-1; ++i) {
        buffer[i] = char(ScaleMantissa(manptr, manptr, 10)) + '0';
    }
    buffer[i] = 0;

    if (manptr[prec] >= (1u<<31)) { // we should round up
        while (i--)
            if (buffer[i]++ == '9') buffer[i] = '0';
            else break;
    }

    return i==-1;
}

// negation
LongFloat LongFloat::operator - () const
{
    return LongFloat(m_Special, !m_Neg, m_Man, m_Exp, m_Prec);
}

// addition
LongFloat LongFloat::operator + (const LongFloat &rhs) const
{
    // is it really substraction? 
    if (m_Neg != rhs.m_Neg)
        return *this - (-rhs);

    // handle special values
    switch (m_Special) {
    case Zero:
        return rhs;
    case Nan:
        return *this;
    case Infinity:
        switch (rhs.m_Special) {
        case Nan:
            return rhs;
        case Zero:
        case Normal:
            return *this;
        case Infinity:
            return *this;
        }
        case Normal:
            switch (rhs.m_Special) {
            case Nan:
            case Infinity:
                return rhs;
            case Zero:
                return *this;
            case Normal: {
                // Normal + Normal: add mantissas and adjust for carry
                Mantissa man;
                u32 *manptr = man.getWritablePtr();
                bool carry;
                exp_type exp = max(m_Exp, rhs.m_Exp);

                if (m_Exp == rhs.m_Exp)
                    carry = AddMantissa(manptr, m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), 0);
                else if (m_Exp > rhs.m_Exp)
                    carry = AddMantissa(manptr, m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), i32_saturated(exp - rhs.m_Exp));
                else
                    carry = AddMantissa(manptr, rhs.m_Man.getConstPtr(), m_Man.getConstPtr(), i32_saturated(exp - m_Exp));

                if (carry)
                    exp += AdjustForCarry(manptr, carry);

                return LongFloat(Normal, m_Neg, man, exp, min(m_Prec, rhs.m_Prec));
            }
            }
    }
    return LongFloat(Nan, false);
}

// substraction
LongFloat LongFloat::operator - (const LongFloat &rhs) const
{
    // is it really addition?
    if (m_Neg != rhs.m_Neg)
        return *this + (-rhs);

    // handle special values
    switch (m_Special) {
    case Zero:
        return -rhs;
    case Nan:
        return *this;
    case Infinity:
        switch (rhs.m_Special) {
        case Nan:
            return rhs;
        case Zero:
        case Normal:
            return *this;
        case Infinity:
            return LongFloat(Nan, false);
        }
        case Normal:
            switch (rhs.m_Special) {
            case Nan:
                return rhs;
            case Infinity:
                return -rhs;
            case Zero:
                return *this;
            case Normal: {
                // Normal - Normal: sub mantissas, negate if necessary, normalize
                Mantissa man;
                u32 *manptr = man.getWritablePtr();
                bool carry;
                exp_type exp;
                bool neg; 

                if (m_Exp == rhs.m_Exp) {
                    carry = SubMantissa(manptr, m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), 0);
                    exp = m_Exp;
                    neg = m_Neg;
                } else if (m_Exp > rhs.m_Exp) {
                    carry = SubMantissa(manptr, m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), i32_saturated(exp_type(m_Exp) - rhs.m_Exp));
                    exp = m_Exp;
                    neg = m_Neg;
                } else {
                    carry = SubMantissa(manptr, rhs.m_Man.getConstPtr(), m_Man.getConstPtr(), i32_saturated(exp_type(rhs.m_Exp) - m_Exp));
                    exp = rhs.m_Exp;
                    neg = !rhs.m_Neg;
                }

                if (carry) {
                    NegMantissa(manptr);
                    neg = !neg;
                } 

                int corr = NormalizeMantissa(manptr);
                if (corr == g_WorkingPrecision) return LongFloat(Zero);
                else return LongFloat(Normal, neg, man, exp-corr, min(m_Prec, rhs.m_Prec));
            }
            }
    }
    return LongFloat(Nan, false);
}

// binary scale
LongFloat LongFloat::operator << (i32 howmuch) const
{
    // has no effect on special values
    if (m_Special != Normal) return *this;

    // split to exponent offset and BScale mantissa amount
    i32 exp;
    if (howmuch >= 0)
        exp = howmuch / 32;
    else 
        exp = (howmuch - 31) / 32;


    howmuch -= exp * 32;
    if (howmuch == 0) 
        return LongFloat(Normal, m_Neg, m_Man, m_Exp + exp, m_Prec);

    Mantissa man;
    u32 *manptr = man.getWritablePtr();

    u32 carry = BScaleMantissa(manptr, m_Man.getConstPtr(), howmuch);
    if (carry)
        exp += AdjustForCarry(manptr, carry);

    return LongFloat(Normal, m_Neg, man, exp_type(m_Exp) + exp, m_Prec);

}

// multiplication
LongFloat LongFloat::operator * (const LongFloat &rhs) const
        {
    bool neg = m_Neg ^ rhs.m_Neg;

    switch (m_Special) {
    case Nan:
        return *this;
    case Zero:
        switch (rhs.m_Special) {
        case Nan:
            return rhs;
        case Infinity:
            return LongFloat(Nan, false);
        case Normal:
        case Zero:
            return rhs.m_Neg ? *this : - *this;
        }
        case Infinity:
            switch (rhs.m_Special) {
            case Nan:
                return rhs;
            case Zero:
                return LongFloat(Nan, false);
            case Infinity:
            case Normal:
                return rhs.m_Neg ? *this : - *this;
            }
            case Normal:
                switch (rhs.m_Special) {
                case Nan:
                    return rhs;
                case Zero:
                case Infinity:
                    return m_Neg ? rhs : -rhs;
                case Normal: {
                    // Normal * Normal: add exponents, mul mantissas, adjust for carry
                    exp_type exp = exp_type(m_Exp) + rhs.m_Exp - 1;
                    Mantissa man;
                    u32* manptr = man.getWritablePtr();

                    u32 carry = MulMantissa(manptr, m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), g_WorkingPrecision - GetPrecision(), GetPrecision());

                    if (carry)
                        exp += AdjustForCarry(manptr, carry);

                    assert (man[g_WorkingPrecision - 1] != 0);
                    return LongFloat(Normal, neg, man, i32_saturated(exp), min(m_Prec, rhs.m_Prec));
                }
                }
    }
    return LongFloat(Nan);
        }

// multiplication by int
LongFloat LongFloat::operator * (i32 rhs) const
        {
    bool neg = m_Neg;

    if (m_Special != Normal) 
        if (rhs == 0 && m_Special == Infinity)
            return LongFloat(Nan, false);
        else return *this;

    if (rhs == 0)
        return LongFloat(Zero, neg);

    if (rhs < 0) {
        rhs = -rhs;
        neg = !neg;
    }

    Mantissa man;
    u32 *manptr = man.getWritablePtr();
    u32 carry = ScaleMantissa(manptr, m_Man.getConstPtr(), rhs);
    exp_type exp = m_Exp;

    if (carry)
        exp += AdjustForCarry(manptr, carry);

    return LongFloat(Normal, neg, man, exp, m_Prec);
        }

// multiplication by 1/int
LongFloat LongFloat::operator / (i32 rhs) const
        {
    bool neg = m_Neg;

    if (m_Special != Normal) 
        if (rhs == 0 && m_Special == Infinity)
            return LongFloat(Nan, false);
        else return *this;

    if (rhs == 0)
        return LongFloat(Zero, neg);

    if (rhs < 0) {
        rhs = -rhs;
        neg = !neg;
    }

    Mantissa man;
    u32 *manptr = man.getWritablePtr();
    exp_type exp = exp_type(m_Exp) + exp_type(InvScaleMantissa(manptr, m_Man.getConstPtr(), rhs));

    return LongFloat(Normal, neg, man, exp, m_Prec);
        }

// directly add to exponent
LongFloat& LongFloat::AddToExponent(i32 howmuch)
{
    m_Exp += howmuch;
    return *this;
}

// addProduct: could be done faster
LongFloat LongFloat::addProduct(const LongFloat &a, const LongFloat &b) const
{
    return *this + a*b;
}

LongFloat LongFloat::addProduct(const LongFloat &a, i32 b) const
{
    return *this + a*b;
}

// reciprocal
LongFloat LongFloat::recip() const
{
    switch (m_Special) {
    case Zero:
        return LongFloat(Infinity, m_Neg);
    case Infinity:
        return LongFloat(Zero, m_Neg);
    case Nan:
        return *this;
    default:
    case Normal:
        int il = GetPrecision();
        int is = g_WorkingPrecision - il;
        if (MultipliedByConvolution(il)) {
            // Newton-Raphson iterations
            // the algorithm:
            // separate mantissa and exponent
            // get the mantissa to double precision, calculate reciprocal
            // this gets us >52 correct binary digits
            // for each iteration error = error * error * |y|,
            // i.e. at least doubles the number of correct digits
            double init(1.0 / MantissaAsDouble());

            // shouldn't I be using a mantissa directly ??
            LongFloat two(2, 0);
            LongFloat my(Normal, true, m_Man, 0);
            LongFloat r(init);

            for (int i=52; i/32 < il; i *= 2) {
                r.SetPrecision((i+31)/16);    // use i*2 to calculate precision
                my.SetPrecision((i+31)/16);
                r = r * two.addProduct(r, my);
            }

            return LongFloat(Normal, m_Neg, r.m_Man, exp_type(r.m_Exp) - m_Exp, m_Prec);
        } else {
            // direct division with O(n*n) complexity
            Mantissa a, t1, t2;
            u32 *aptr = a.getWritablePtr();
            for (int i=is; i<il+is-1; ++i)
                aptr[i] = 0;
            aptr[il+is-1] = 1;
            // a = 1 * 2^-32

            int e = DivMantissa(aptr, aptr, m_Man.getConstPtr(), is, il, t1.getWritablePtr(), t2.getWritablePtr());
            // e is 0 only if m_Man is the same as a

            return LongFloat(Normal, m_Neg, a, 1 + e -exp_type(m_Exp), m_Prec); 
        }
    }
}

// division
LongFloat LongFloat::operator /(const LongFloat &rhs) const
{
    if (MultipliedByConvolution(GetPrecision()) || m_Special != Normal || rhs.m_Special != Normal) 
        return *this * rhs.recip();

    int il = GetPrecision();
    int is = g_WorkingPrecision - GetPrecision();

    Mantissa man, t1, t2;
    int e = DivMantissa(man.getWritablePtr(), m_Man.getConstPtr(), rhs.m_Man.getConstPtr(), is, il, t1.getWritablePtr(), t2.getWritablePtr());
    LongFloat lf(LongFloat(Normal, m_Neg != rhs.m_Neg, man, exp_type(m_Exp) - rhs.m_Exp + e, m_Prec));

    LongFloat l;
    //    assert((l=lf*rhs - *this).Kind() == Zero || (Exponent() - l.Exponent()) >= GetPrecision()-1);
    return lf;
}

i32 LongFloat::DivisionRoundingError() const
{ 
    return -(GetPrecision() - (MultipliedByConvolution(GetPrecision()) ? 1 : 1)) * 32 + 2; 
}


// comparison: checks difference negative
bool LongFloat::operator >= (const LongFloat &rhs) const
        {
    return !(*this - rhs).IsNegative();
        }

// 
bool LongFloat::operator == (const LongFloat &rhs) const
        {
    return (*this - rhs).Kind() == Zero;
        }

// RoundTowardZero: truncates the mantissa at the point where exponent is zero
LongFloat LongFloat::RoundTowardZero() const
{
    if (m_Special != Normal || m_Exp >= g_WorkingPrecision) return *this;
    if (m_Exp <= 0) return LongFloat(Zero, m_Neg);

    Mantissa man;
    int i;

    for (i=1; i <= m_Exp; ++i)
        man[g_WorkingPrecision - i] = m_Man(g_WorkingPrecision - i);
    while (i<=g_WorkingPrecision) man[g_WorkingPrecision - i++] = 0;

    return LongFloat(Normal, m_Neg, man, m_Exp, m_Prec);
}

// round: round to nearest integer
LongFloat LongFloat::round() const
{
    if (m_Special != Normal || m_Exp >= g_WorkingPrecision) return *this;
    if (m_Exp < 0) return LongFloat(Zero, m_Neg);
    if (m_Exp == 0) return 
            m_Man(g_WorkingPrecision - 1) >= (1u<<31) ?
                    LongFloat(m_Neg ? -1.0 : 1.0) : LongFloat(Zero, m_Neg);

    Mantissa man;
    u32 *manptr = man.getWritablePtr();
    int i;

    for (i=1; i <= m_Exp; ++i)
        manptr[g_WorkingPrecision - i] = m_Man(g_WorkingPrecision - i);
    int j = g_WorkingPrecision - i;
    while (i<=g_WorkingPrecision) manptr[g_WorkingPrecision - i++] = 0;        

    exp_type exp = m_Exp;
    if (j>=0 && m_Man(j) >= (1u<<31)) {
        while (++j < g_WorkingPrecision && !++manptr[j]) ;

        if (j == g_WorkingPrecision)
            exp += AdjustForCarry(manptr, 1);

    }

    return LongFloat(Normal, m_Neg, man, exp, m_Prec);
}

i32 LongFloat::normalize() const
{
    if (m_Special != Normal) return 0;

    exp_type e = Exponent() * 32;
    u32 m = m_Man(g_WorkingPrecision - 1);

    assert(m!=0);
    while (!(m & 0x80000000)) { 
        m <<= 1;
        e -= 1;
    }
    return i32_saturated(e);
}

char* LongFloat::AsDecimal(char *bufptr, u32 buflen) const
{
    // must at least accomodate exponent
    assert (buflen >= 10);
    char *buffer = bufptr;

    // the format chosen is "[-].<mantissa>e<+/-><exponent>"
    // where mantissa has a leading non-zero decimal

    LongFloat a(*this);
    if (IsNegative()) {
        buffer++[0] = '-';
        --buflen;
        a = -a;
    }

    // handle special values
    switch (Kind()) {
    case Nan:
        strcpy(buffer, "NaN");
        return bufptr;
    case Infinity:
        strcpy(buffer, "Infinity");
        return bufptr;
    case Zero:
        strcpy(buffer, "Zero");
        return bufptr;
    default:
        ;
    }

    // calculate exponent: the least power of 10 that is greater than
    // or equal to the value.
    int pwr = int(normalize() / LOG_2_10);
    // divide the value by the exponent to form decimal mantissa
    a = a / pow(LongFloat(10.0), pwr);
    // double arithmetic can be wrong...
    if (a > 1) { a /= 10; pwr++; }
    else
    {
        LongFloat b = a * 10;
        if (b < 1) { a = b; pwr--; }
    }

    // GCC doesn't understand _itoa, use sprintf instead
    sprintf(buffer, "%+d", pwr);
    // now we know the length of the exponent. move it to the end of
    // the string
    size_t explen = strlen(buffer);
    strcpy(buffer + buflen - explen - 1, buffer);
    buffer[0] = '.';

    // when the value is power of ten, we can get two possible
    // representations of the mantissa: 0.(9) and 1.(0). This
    // behavior is not an error as it is consistent with the
    // theory. The function used to convert the mantissa will
    // not display 1.(0) correctly. Thus we must handle the
    // case differently: just pretend it's 0.(9)
    if (a==1)
        memset(buffer+1, '9', buflen - explen - 2);
    else
        // convert the mantissa to decimal using LongFloat's function
        a.MantissaAsDecimal(buffer+1, buflen - explen - 2);

    char *ptr = buffer + (buflen - explen - 2);
    *ptr = 'e';

    return bufptr;
}

std::ostream& operator <<(std::ostream &os, const LongFloat &lf)
{
    using namespace std;

    int w = os.width(0);    // the width only applies to this operation
    int pr = os.precision();
    if (os.flags() & ios_base::scientific
            || os.flags() & ios_base::fixed) ++pr;
    int prall = pr+20;

    if ((os.flags() & ios_base::fixed) && lf.Kind() == LongFloat::Normal) {
        int z = int(ceil(lf.normalize() / LOG_2_10)) + 1;
        if (z > 0) prall += z;
    }

    bool mallocd = false;
    char *buffer;
    if (prall <= 4096) buffer = (char *)alloca(prall);
    else {
        buffer = (char *) malloc(prall);
        mallocd = true;
    }

    char *buf = buffer;

    // the format chosen is "[-].<mantissa>e<+/-><exponent>"
    // where mantissa has a leading non-zero decimal

    LongFloat a(lf);
    if (lf.IsNegative()) {
        if (os.flags() & ios_base::internal) {
            os << '-';
            if (w) --w;
        } else *(buf++) = '-';
        a = -a;
    } else if (os.flags() & ios_base::showpos) {
        if (os.flags() & ios_base::internal) {
            os << '+';
            if (w) --w;
        } else *(buf++) = '+';
        if (w) --w;
    }

    // handle special values
    switch (lf.Kind()) {
    case LongFloat::Nan:
        os << setw(w) << "NaN";
        if (mallocd) free(buffer);
        return os;
    case LongFloat::Infinity:
        os << setw(w) << "Infinity";
        if (mallocd) free(buffer);
        return os;
    case LongFloat::Zero:
        os << setw(w) << "Zero";
        if (mallocd) free(buffer);
        return os;
    default:
        ;
    }

    // calculate exponent: the least power of 10 that is greater than
    // or equal to the value.
    int pwr = int(a.normalize() / LOG_2_10);
    // divide the value by the exponent to form decimal mantissa
    a = a / pow(LongFloat(10.0), pwr);
    // normalize can result in the wrong decimal exponent 
    if (a >= 1)
    { a /= 10; pwr++; }
    else
    {
        LongFloat b = a * 10;
        if (b < 1) { a = b; pwr--; }
    }

    if (os.flags() & ios_base::fixed && pwr > 1) pr += pwr-1;
    if (!(os.flags() & ios_base::fixed || 
            os.flags() & ios_base::scientific) && pwr <= 0 && pwr > -4) pr += -pwr + 1;

    // convert the mantissa to decimal using LongFloat's function
    if (a.MantissaAsDecimal(buf, pr+1)) {
        buf[0] = '1';
        for (int i=1;i<pr;++i) buf[i] = '0';
        buf[pr] = 0;
        pwr++;
    }

    int pwrsv = pwr;

    if ((os.flags() & ios_base::scientific) ||
            (!(os.flags() & ios_base::fixed) && (pwr > pr || (pwr) <= -4))) {
        pwr = 1;
    } 

    if (pwr <= 0) {
        int pww = min(-(pwr-1), i32(pr));
        //memcpy(buf+pwr, buf, pr - pwr);
        for (int i=pr;i>pww;--i)
            buf[i] = buf[i-pww-1];
        memset(buf, '0', pww+1);
        buf[1] = '.';
    } else {
        //memcpy(buf+pwr+1, buf+pwr, pr - pwr);
        for (int i=pr;i>pwr;--i)
            buf[i] = buf[i-1];
        buf[pwr] = '.';
    }
    buf[pr + 1] = 0;

    if (!(os.flags() & ios_base::showpoint || 
            os.flags() & ios_base::fixed ||
            os.flags() & ios_base::scientific)) {
        while (buf[pr] == '0') --pr;
        if (buf[pr] == '.') --pr;
        buf[pr+1] = 0;
    }

    if (os.flags() & ios_base::scientific || pwrsv != pwr) {
        // GCC doesn't understand _itoa, use sprintf instead
        sprintf(buf+pr+1, "%c%+d", os.flags() & ios_base::uppercase ? 'E' : 'e', pwrsv - 1);
    } 
    os << setw(w) << buffer;

    if (mallocd) free(buffer);
    return os;
}

} // namespace
