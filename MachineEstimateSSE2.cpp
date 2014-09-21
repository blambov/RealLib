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

#include <cfloat>
#include <cassert>
#include "MachineEstimate.h"

namespace RealLib {

__m128d MachineEstimate::signmask = _mm_set_pd(0.0, -1.0 * 0.0);
__m128d MachineEstimate::mdelta = _mm_set1_pd(-DBL_MIN);
__m128d MachineEstimate::half = _mm_set1_pd(0.5);
__m128d MachineEstimate::mhalf = _mm_set_pd(-0.5, 0.5);
__m128d MachineEstimate::zero = _mm_set1_pd(0.0);
__m128d MachineEstimate::mone = _mm_set1_pd(-1.0);
__m128d MachineEstimate::sqrt_corr = _mm_set1_pd(0.0);
__m128d MachineEstimate::coeff_sin[6];
__m128d MachineEstimate::pi;
__m128d MachineEstimate::rpiover6;
__m128d MachineEstimate::onethird;
int MachineEstimate::SavedRoundingMode = _MM_ROUND_NEAREST;
static int initialized = 0;

// from Hart et. al., "Computer Approximations", Wiley
// table SIN 3043
#define SIN0a   +0.785398163397448307014
#define SIN1a   -0.80745512188280530192e-1
#define SIN2a   +0.2490394570188736117e-2
#define SIN3a   -0.36576204158455695e-4
#define SIN4a   +0.313361621661904e-6
#define SIN5a   -0.1757149292755e-8     // rounds tw +inf, all others tw -inf
#define SIN6a   +0.6877100349e-11
// table SIN 2922
#define SIN0    +0.52359877559829885532
#define SIN1    -0.2392459620393377657e-1
#define SIN2    +0.32795319441392666e-3
#define SIN3    -0.214071970654441e-5
#define SIN4    +0.815113605169e-8
#define SIN5    -0.2020852964e-10
#define SIN0s   "+0.52359877559829885532"
#define SIN1s   "-0.2392459620393377657e-1"
#define SIN2s   "+0.32795319441392666e-3"
#define SIN3s   "-0.214071970654441e-5"
#define SIN4s   "+0.815113605169e-8"
#define SIN5s   "-0.2020852964e-10"
// table SIN 2923
/*#define SIN0  +0.523598775598298 873071308
#define SIN1    -0.23924596203935045866796e-1
#define SIN2    +0.327953194428661969081e-3
#define SIN3    -0.2140719769181988118e-5
#define SIN4    +0.8151256504047484e-8
#define SIN5    -0.20315350937751e-10
#define SIN6    +0.35539710328e-13
*/

#define M_PI    3.14159265358979323846264338328
#ifndef _MSC_VER
#define _nextafter nextafter
#endif
void MachineEstimate::BeginComputation()
{
    SavedRoundingMode = _MM_GET_ROUNDING_MODE();
    _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);

    if (!initialized) {
        double z = 0.0;
        double minusinf = -1.0 / z;
        double plusinf = 1.0 / z;
        sqrt_corr = //_mm_mul_pd(
            _mm_set_pd(0.0, minusinf);//, zero);
                // note: we're nudging the first coefficient a bit more to accommodate the approximation error
        coeff_sin[0] = _mm_set_pd(_nextafter(SIN0, minusinf), - _nextafter(SIN0, plusinf));   
        coeff_sin[1] = _mm_set_pd(_nextafter(SIN1, minusinf), -SIN1);
        coeff_sin[2] = _mm_set_pd(SIN2, - _nextafter(SIN2, plusinf));
        coeff_sin[3] = _mm_set_pd(_nextafter(SIN3, minusinf), -SIN3);
        coeff_sin[4] = _mm_set_pd(_nextafter(SIN4, minusinf), -SIN4);
        coeff_sin[5] = _mm_set_pd(SIN5, - _nextafter(SIN5, plusinf));   
        pi = _mm_set_pd(M_PI, - _nextafter(M_PI, plusinf));
        rpiover6 = _mm_xor_pd(_mm_div_pd(_mm_set1_pd(6.0), pi), signmask);
        onethird = _mm_set_pd(1.0/3.0, _nextafter(1.0/3.0, plusinf));
        initialized = 1;
    }

}

void MachineEstimate::FinishComputation()
{
    assert(_MM_GET_ROUNDING_MODE()==_MM_ROUND_DOWN);
    _MM_SET_ROUNDING_MODE(SavedRoundingMode);
}

// sine in its primary range [-pi/6, pi/6] mapped in [1, 1]
static inline __m128d sinprimarymapped(__m128d x)
{
    __m128d x2 = _mm_mul_pd(x, x);  
    x = _mm_xor_pd(x, MachineEstimate::signmask);
    __m128d p10 = _mm_add_pd(MachineEstimate::coeff_sin[0], _mm_mul_pd(MachineEstimate::coeff_sin[1], x2));
    __m128d x4 = _mm_mul_pd(x2, x2);
    __m128d p32 = _mm_add_pd(MachineEstimate::coeff_sin[2], _mm_mul_pd(MachineEstimate::coeff_sin[3], x2));
    __m128d x8 = _mm_mul_pd(x4, x4);
    __m128d p54 = _mm_add_pd(MachineEstimate::coeff_sin[4], _mm_mul_pd(MachineEstimate::coeff_sin[5], x2));
    __m128d p3210 = _mm_add_pd(p10, _mm_mul_pd(p32, x4));
    __m128d s = _mm_add_pd(p3210, _mm_mul_pd(p54, x8));
    s = _mm_mul_pd(s, x);
    return s;
}

// sine in its primary range [0, pi/6] mapped in [0, 1]
static inline __m128d sinprimarymapped1(__m128d x)
{
    __m128d xp = _mm_xor_pd(x, MachineEstimate::signmask);
    __m128d x2 = _mm_mul_pd(x, xp); // upper bound remains negative, correct rounding
    x2 = _mm_xor_pd(x2, MachineEstimate::signmask);
    __m128d s = MachineEstimate::coeff_sin[5];
    s = _mm_add_pd(_mm_mul_pd(s, x2), MachineEstimate::coeff_sin[4]);
    s = _mm_add_pd(_mm_mul_pd(s, x2), MachineEstimate::coeff_sin[3]);
    s = _mm_add_pd(_mm_mul_pd(s, x2), MachineEstimate::coeff_sin[2]);
    s = _mm_add_pd(_mm_mul_pd(s, x2), MachineEstimate::coeff_sin[1]);
    s = _mm_add_pd(_mm_mul_pd(s, x2), MachineEstimate::coeff_sin[0]);       // the roundings here account for the
                                                                            // approximation error
    s = _mm_mul_pd(s, x);
    return s;
}

static inline __m128d sin_reduce_before(__m128d a)
{
    return _mm_mul_pd(a, MachineEstimate::onethird);
}

static inline __m128d sin_reduce_after(__m128d a)
{
    __m128d c = _mm_xor_pd(a, MachineEstimate::signmask);
    __m128d b = _mm_mul_pd(a, c);   // correct rounding
    b = _mm_shuffle_pd(b, b, 1); // negation effect
    b = _mm_mul_pd(b, _mm_set_pd(4.0, 4.0));    // no rounding here
    b = _mm_add_pd(b, _mm_set_pd(3.0, -3.0));   // remains positive (provided a is within [0, pi/6]), correct rounding
                        // careful! wrong signs in constants can make the rounding incorrect
    b = _mm_mul_pd(c, b);
    return b;       // correct rounding and correct signs for both */
}

// sine in its primary domain [0, pi/2] mapped in [0, 3]
MachineEstimate sinprimary(const MachineEstimate &arg)
{
    __m128d a = sin_reduce_before(arg.interval);
    a = sinprimarymapped1(a);
    return sin_reduce_after(a);
}

MachineEstimate sin(const MachineEstimate &arg)
{
    __m128d a = _mm_mul_pd(arg.interval, MachineEstimate::rpiover6);
    MachineEstimate ae(a);
    // make sure we only have the fractional part
    if (!(ae < 6.0 && ae > -6.0)) {
        // !!! what if we're on different sides of 1?
        double ac = ae.weak_AsDouble() / 12;
        ac = (floor(ac + 0.5)) * 12;
        ae -= MachineEstimate(ac);
    }

    // here we should be within [-6, 6] unless the error is far too big

    if (ae > 3.0) ae = 6.0 - ae;
    if (ae < -3.0) ae = -6.0 - ae;

    if (!(ae < 3.0 && ae > -3.0)) {
        // this is a bit of a problem... 
        // we handle it by just doing the strength reduction twice

        // !!! I think we need more work here

        ae = MachineEstimate(_mm_mul_pd(ae.interval, MachineEstimate::onethird));
        ae = sinprimary(ae);
        ae = ae * (MachineEstimate(3) - sq(ae)*MachineEstimate(4.0));
        return ae;
    }

    // !!! what happens for negative ae?
    return sinprimary(ae.interval);
}

}   // namespace

