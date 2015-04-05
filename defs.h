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

#ifndef FILE_DEFS_H
#define FILE_DEFS_H

#include <stdint.h>

namespace RealLib {

typedef int32_t i32;
typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;

#define I32_MAX INT_MAX
#define I32_MIN INT_MIN

//#define TRUST_STDLIB

#if defined(__GNUC__) || defined(__MWERKS__)
// doesn't define NDEBUG in non-debug builds
// comment this if you want a debug build
#define NDEBUG

#define NAMESPACE_STD std

#else  // MS Visual C++

#define NAMESPACE_STD 
#endif

#if defined(__MWERKS__)
typedef i32 exp_type;
#define MINIMUM_EXPONENT (-(1<<28))
#else
typedef i64 exp_type;
#define MINIMUM_EXPONENT (-I32_MAX)
#endif
} // namespace

#if defined(__GNUC__) || defined(__MWERKS__)
#include "GCChelper.h"
#endif

#if defined(__MWERKS__)
#define FUNCTION(x) &x
#else
#define FUNCTION(x) &(x)
#endif

// if this is set, multiplications of large number will
// be done via faster convolution
// recommended to keep this on
#define MULTIPLY_BY_CONVOLUTION

#ifdef MULTIPLY_BY_CONVOLUTION
// direct multiplication will be used for precisions below
// this threshold and convolution for larger precisions
#define CONVOLUTION_THRESHOLD 60
#else
#define CONVOLUTION_THRESHOLD INT_MAX
#endif

// without this the system will not limit its recursion depth
// may run slightly faster but will probably cause errors
// for longer computations on the class Real
// recommended to keep this on
#define LIMIT_RECURSION

#ifdef LIMIT_RECURSION
// the size of the chunk that is processed recursively
#define EvaluationDepth 500
#else
#define EvaluationDepth INT_MAX
#endif

#define LOG_2_10 3.32192809488736234787031942948939017586483139

#if !defined(NDEBUG) && defined(_MSC_VER)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#ifndef NDEBUG
// with this the system prints a message every time the
// precision is increased
#define REALLIB_SHOW_PRECISION_INCREASES
#endif

#if (_MSC_VER >= 1400)
#pragma warning (disable: 4996)
#endif

// the SSE2 version of MachineEstimate can rely on SSE2
// exception flags to recognize invalid operations
// and overflows
// with this turned off it will do more checks instead
// (slower performance)
#define REALLIB_RELY_ON_SSE_EXCEPTIONS

#ifndef NDEBUG
// relying on SSE exceptions does less checks
// and gives somewhat less precise information
// about the origin of the exception
// so use the other mode for debug builds

#undef REALLIB_RELY_ON_SSE_EXCEPTIONS
#endif

#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
#define REALLIB_MACHEST_INVALID_CHECK(x) false
#else
#define REALLIB_MACHEST_INVALID_CHECK(x) !(x).IsValueValid()
#endif

#endif // file
