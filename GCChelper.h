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

#ifndef FILE_GCCHELPER_H
#define FILE_GCCHELPER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace RealLib {

/* Return values for fpclass. */
#define _FPCLASS_SNAN   0x0001  /* Signaling "Not a Number" */
#define _FPCLASS_QNAN   0x0002  /* Quiet "Not a Number" */
#define _FPCLASS_NINF   0x0004  /* Negative Infinity */
#define _FPCLASS_NN 0x0008  /* Negative Normal */
#define _FPCLASS_ND 0x0010  /* Negative Denormal */
#define _FPCLASS_NZ 0x0020  /* Negative Zero */
#define _FPCLASS_PZ 0x0040  /* Positive Zero */
#define _FPCLASS_PD 0x0080  /* Positive Denormal */
#define _FPCLASS_PN 0x0100  /* Positive Normal */
#define _FPCLASS_PINF   0x0200  /* Positive Infinity */

static inline bool __isnan(const double x)
{
    return x != x;
}

static inline bool __isinf(const double x)
{
    return x == x + 1.0;
}

static inline int _fpclass(const double x)
{
    if (__isnan(x)) return _FPCLASS_QNAN;
    else if (__isinf(x)) return x > 0 ? _FPCLASS_PINF : _FPCLASS_NINF;
    else if (x == 0.0) return _FPCLASS_PZ;
    else return x > 0 ? _FPCLASS_PN : _FPCLASS_NN;
}

} // namespace

#endif
