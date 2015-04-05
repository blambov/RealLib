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

   RealFuncs.h

   This file defines the code that calculates functions.
   There are LongFloat and Estimate functions. LongFloat
   do not require a precision specification, while Estimate
   functions use the second argument to know when to stop
   iterating.
 */

#ifndef FILE_REAL_MATH_H
#define FILE_REAL_MATH_H

#include "MachineEstimate.h"
#include "RealEstimate.h"

namespace RealLib {

// whenever LongFloat is closed there have to be no Estimates around. This function
// is called by FinalizeRealLib to clean the Estimates for the constants.
void DestroyConstantEstimates();

// helper function
template <class T>
static inline T sq(const T &arg)
{ return arg * arg; }

// NewtonIteration: takes argument, previous estimate and previous precision
// outputs new estimate and updates prec (this is roughly bits, expected precision
// at next iteration)
typedef void (*NewtonIterator)(const Estimate &arg, Estimate &est, i32 &prec);

// PerformNewton: applies NewtonIterations to increasingly precise estimates
// prec is always one step ahead, i.e. initialPrec should be the precision you
// expect after one iteration. Using increasing precision, so that the complexity
// can be of the order of the complexity of the inverse function.
Estimate PerformNewton(Estimate arg, NewtonIterator iter, 
                       Estimate initialEstimate, i32 initialPrec);


typedef Estimate (*SeriesIterator)(const Estimate &arg, Estimate &workspace, i32 index);

// PowerSeries: sums the series 
void PowerSeriesDirect(const Estimate &arg, SeriesIterator iter,
                       Estimate &sum, Estimate &workspace, i32 indexstart, i32 indexend);

// absolute value
Estimate abs(const Estimate &arg);

// integer power by squaring and multiplication: O(M(n)*ln(pwr))
Estimate pow(Estimate arg, i32 pwr);

// reciprocal square root: evaluated using newton iterations; O(M(n))
Estimate rsqrt(const Estimate &arg);
// square root: sqrt(x) = x * rsqrt(x); O(ln(n))
Estimate sqrt(const Estimate &arg);

// exponent: using McLauren expansion with strength reduction; O(M(n)*sqrt(n))
Estimate exp_primary(const Estimate &arg); // only in primary range: [-1;1]
Estimate exp(const Estimate &arg);
// logarithm: newton iterations with exp_primary; same complexity as exp
Estimate log_primary(const Estimate &arg); // only in primary range: [1/e;e]
Estimate log(const Estimate &arg);

//template<class TYPE>
//TYPE ln2(unsigned int prec);

template <> Estimate ln2<Estimate>(unsigned int prec);

// pi: using Borwein iterations
Estimate rpi(unsigned int prec);

//template<class TYPE>
//TYPE pi(unsigned int prec);

template <> Estimate pi(unsigned int prec);
//Estimate pi(unsigned int prec);

// forward trigonometric: using McLaurin expansion for sine and identities; O(M(n) * sqrt(n))
Estimate sin_primary(const Estimate &arg); // only in primary range: [-pi/2;pi/2]
Estimate sin(const Estimate &arg);
static inline Estimate cos(const Estimate &arg) { return sin(pi<Estimate>(arg.GetPrecision())/2 - arg); }

template <class TYPE>
static inline TYPE cosfromsin(const TYPE &s) { return sqrt(1 - sq(s)); }

Estimate tan(const Estimate &arg);

// inverse trigonometric: newton iterations with sin_primary; same complexity as forward
Estimate asin_primary(const Estimate &arg); // only in primary range: [-sqrt(2);sqrt(2)]
Estimate asin(const Estimate &arg);

template <class Estimate>
static inline Estimate acos(const Estimate &arg) { return pi<Estimate>(arg.GetPrecision())/2 - asin(arg); }
template <class Estimate>
static inline Estimate atan(const Estimate &arg) { return asin(arg * rsqrt(1 + sq(arg))); }
template <class Estimate>
Estimate atan2(const Estimate &y, const Estimate &x);

/*




// 1 / sqrt(arg), using Newton-Raphson iteration
LongFloat rsqrt(const LongFloat &arg);
// arg ^ pow, integer pow
LongFloat pow(LongFloat arg, i32 pwr);

// shorthands (with guarded Zero/Infinity conditions)
static inline LongFloat sq(const LongFloat &arg)
{ return arg * arg; }

// arg * rsqrt(arg)
LongFloat sqrt(const LongFloat &arg);

// Estimate math

// needed constants
// listed in order of the calls in RealConstants constructor

// 1 / pi calculated using an iterative sequence by Borwein brothers
Estimate rpi(u32 prec);
// 1 / rpi
Estimate pi(u32 prec);
// ln2 calculated by an AGM method
Estimate ln2(u32 prec);
// ln using ln2 on exponent and AGM on mantissa
Estimate ln(const Estimate &arg, u32 prec);

// unary functions

Estimate abs(const Estimate &arg, u32 prec);

static inline Estimate sq(const Estimate &arg)
{ return arg * arg; }

// 1 / sqrt(arg), Newton-Raphson
Estimate rsqrt(const Estimate &arg, u32 prec);
// arg * rsqrt(arg)
Estimate sqrt(const Estimate &arg, u32 prec);

// exponent, Newton-Raphson using ln
Estimate exp(const Estimate &arg, u32 prec);

// cosine, Taylor series with appropriate strength reduction
Estimate cos(const Estimate &arg, u32 prec);
// arccosine, Newton-Raphson using cos
Estimate acos(const Estimate &arg, u32 prec);

// cos(pi/2 - arg)
Estimate sin(const Estimate &arg, u32 prec);
// sqrt(1 - sq(cos(arg))) / cos(arg)
Estimate tan(const Estimate &arg, u32 prec);
// pi/2 - acos(arg)
Estimate asin(const Estimate &arg, u32 prec);
// acos(rsqrt(1 + sq(arg)))
Estimate atan(const Estimate &arg, u32 prec);

// binary

// atan2 using AGM on complex
Estimate atan2(const Estimate &l, const Estimate &r, u32 prec);
 */

} // namespace

#endif
