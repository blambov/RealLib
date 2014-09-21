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
#ifdef _MSC_VER
#include <emmintrin.h>
#else
#include <xmmintrin.h>
#endif
#include "defs.h"
#include "RealEstimate.h"

namespace RealLib {

// class MachineEstimate's definitions start here
class MachineEstimate;

   // operations
static inline   MachineEstimate operator - (const MachineEstimate &arg);
static inline   MachineEstimate recip(const MachineEstimate &arg);
    
static inline   MachineEstimate operator + (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline   MachineEstimate operator - (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline   MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs);
static inline   MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs);

    // fast multiplication
static inline   MachineEstimate operator * (const MachineEstimate &lhs, i32 rhs);
    // and division
static inline   MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs);
    
static inline   std::ostream& operator <<(std::ostream &os, const MachineEstimate &e);
   
class MachineEstimate {
private:
public:

    // we'll be using round-to-minus-infinity mode only,
    // emulating round-to-plus-infinity via -round(-value)
    // we'll keep the high bounded negated for quicker operation

    // the data will be stored in one __m128d variable,
    // -high bound in first element (the one that you can apply _sd operations to)
    // low bound in second
    __m128d interval;       
    
    static __m128d signmask;        // mask of only 1 in the MSB of the first double
    static __m128d mdelta;          // minus smallest representable number in both sides
    static __m128d half;            // 0.5 in both 
    static __m128d mhalf;           // 0.5, -0.5
    static __m128d zero;            // 0.0
    static __m128d mone;            // -1.0
    static __m128d sqrt_corr;       // 0.0, NaN
    static __m128d coeff_sin[6];    // coefficients for sine
    static __m128d pi, rpiover6;
    static __m128d onethird;
    static int SavedRoundingMode;

    MachineEstimate(double l, double h) : interval(_mm_set_pd(l, -h)) {}
    MachineEstimate(__m128d src) : interval(src) {}

    std::ostream& PrintInterval(std::ostream &os) const {
        double d[2];
        _mm_storeu_pd(d, interval);
        return (os << "(" << d[1] << ", " << -d[0] << ")");
    }
    
    // gets the sum of high and low, negated in the first component and positive in the second
    // (i.e. a proper interval)
    __m128d Sum() const { return _mm_sub_pd(interval, _mm_shuffle_pd(interval, interval, 1)); }
    // gets the difference, negated in both components of the __m128d
    // (i.e. negate second component to make a proper interval)
    __m128d MinusDiff() const { return _mm_add_pd(interval, _mm_shuffle_pd(interval, interval, 1)); }
    __m128d GetInterval() const { return _mm_xor_pd(interval, signmask); }

public:
    static void BeginComputation();
    static void FinishComputation();

    MachineEstimate(double v = 0.0) : interval(_mm_xor_pd(_mm_set1_pd(v), signmask)) {}
    MachineEstimate(const char *val) {
        double v = atof(val);
        __m128d z(_mm_set1_pd(v));  // load double
        z = _mm_xor_pd(z, signmask);    // negate first component
        interval = _mm_add_pd(z, mdelta);       // round down
    }

    // error functions
    MachineEstimate GetError() const { return _mm_mul_pd(mhalf, MinusDiff()); }
    MachineEstimate& SetError(const MachineEstimate &err) { 
        __m128d s = _mm_mul_pd(Sum(), half);
        // assuming a positive error!
        __m128d e = _mm_shuffle_pd(err.interval, err.interval, 0);  // negated
        interval = _mm_add_pd(s, e);
        return *this;
    }
    MachineEstimate& AddError(const MachineEstimate &err) {
        // assuming a positive error!
        __m128d e = _mm_shuffle_pd(err.interval, err.interval, 0);  // negated
        interval = _mm_add_pd(interval, e);
        return *this;
    }

   // a lower bound on the correct binary digits
   // uses the exponents of the value and error to calculate it quickly
    i32 GetRelativeError() const {
        int e;
        double d;
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
   { return !!_mm_comigt_sd(interval, zero); }  // -high > 0 
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
{   // if you simply reverse the bounds' order you get the negation
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
MachineEstimate operator * (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ /*__m128d a = _mm_mul_pd(lhs.interval, rhs.interval);
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
  __m128d a = _mm_shuffle_pd(lhs.interval, lhs.interval, 1);    // b, -a
  __m128d b = _mm_mul_pd(lhs.interval, rhs.interval);           // ac, bd
  __m128d c = _mm_sub_pd(MachineEstimate::zero, rhs.interval);  // c, -d
  __m128d d = _mm_mul_pd(a, rhs.interval);                      // -bc, -ad
  __m128d e = _mm_mul_pd(lhs.interval, c);                      // -ac, -bd
  __m128d f = _mm_mul_pd(a, c);                                 // bc, ad
  __m128d g = _mm_min_pd(d, e);                                 // min(-bc,-ac), min(-ad, -bd)
  __m128d h = _mm_min_pd(b, f);                                 // min(ac, bc), min(bd, ad)
  __m128d i = _mm_unpackhi_pd(g, h);                            // min(-ad, -bd), min(ac, bc)
  __m128d j = _mm_unpacklo_pd(g, h);                            // min(-bc,-ac), min(bd, ad)
  __m128d k = _mm_min_pd(i, j);                                 // min(-all), min(all)
  return MachineEstimate(k);
}

static inline
MachineEstimate operator / (const MachineEstimate &lhs, const MachineEstimate &rhs)
{ 
    return lhs * recip(rhs);
/*  if (!rhs.IsNonZero()) throw PrecisionException("recip");

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
{ if (rhs >=0) return MachineEstimate(_mm_mul_pd(lhs.interval, _mm_set1_pd(double(rhs)))); 
  else return -lhs * -rhs;}

// and division
static inline
MachineEstimate operator / (const MachineEstimate &lhs, i32 rhs)
{   if (rhs==0) throw DomainException("integer div");
    else if (rhs > 0) return MachineEstimate(_mm_div_pd(lhs.interval, _mm_set1_pd(double(rhs)))); 
    else return -lhs / -rhs;}
    
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
{   return os.operator<<(e.weak_AsDouble()); }
//  return os << e.weak_AsDouble(); }


static inline MachineEstimate sqrt(const MachineEstimate &arg)
{
    __m128d a = _mm_xor_pd(arg.interval, MachineEstimate::signmask);
    __m128d b = _mm_sqrt_pd(a);
    __m128d c = _mm_sub_sd(MachineEstimate::mdelta, b);
    __m128d d = _mm_move_sd(b, c); 
    __m128d e = _mm_max_pd(d, MachineEstimate::sqrt_corr);  // to get rid of the NaN's
                        // according to our convention we do not care about the invalid part of the intervals
                        // but a NaN in the high part would mean the value was provably negative
    return MachineEstimate(e);
}

static inline MachineEstimate abs(const MachineEstimate &arg)
{
    __m128d a = arg.interval;                   // a, b
    __m128d b = _mm_shuffle_pd(a, a, 1);
    __m128d c = _mm_min_sd(a, b);
    __m128d d = _mm_max_sd(a, b);
    __m128d e = _mm_max_sd(d, MachineEstimate::zero);
    __m128d f = _mm_unpacklo_pd(c, e);
    return MachineEstimate(f);
}

static inline MachineEstimate sq(const MachineEstimate &arg)
{
    MachineEstimate x(abs(arg));
    return MachineEstimate(_mm_mul_pd(x.interval, _mm_xor_pd(x.interval, MachineEstimate::signmask))); 
}

// sine in its primary domain [0, pi/2]
MachineEstimate sinprimary(const MachineEstimate &arg);
// sine in the full domain
MachineEstimate sin(const MachineEstimate &arg);

}   // namespace

#endif // FILE
