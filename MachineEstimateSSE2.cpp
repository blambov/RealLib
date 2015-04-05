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
#include <cstring>
#include "MachineEstimate.h"

namespace RealLib {

__m128d MachineEstimate::signmask = _mm_set_pd(0.0, -1.0 * 0.0);
__m128d MachineEstimate::mdelta = _mm_set1_pd(-DBL_MIN);
__m128d MachineEstimate::half = _mm_set1_pd(0.5);
__m128d MachineEstimate::mhalf = _mm_set_pd(-0.5, 0.5);
__m128d MachineEstimate::mquarter = _mm_set_pd(-0.25, 0.25);
__m128d MachineEstimate::zero = _mm_setzero_pd();
__m128d MachineEstimate::mone = _mm_set1_pd(-1.0);
__m128d MachineEstimate::one = _mm_set_pd(1.0, -1.0);
__m128d MachineEstimate::sqrt_corr = _mm_set1_pd(0.0);
__m128d MachineEstimate::coeff_sin[8];
__m128d MachineEstimate::coeff_cos[8];
__m128d MachineEstimate::coeff_log[20];
__m128d MachineEstimate::coeff_exp[12];
__m128d MachineEstimate::coeff_atan[20];
__m128d MachineEstimate::sqrt2;
__m128d MachineEstimate::sqrtsqrt2;
__m128i MachineEstimate::expmask;
__m128i MachineEstimate::expbias;
__m128d MachineEstimate::explimit;
__m128d MachineEstimate::sign;
__m128d MachineEstimate::log2e;
__m128d MachineEstimate::pi;
__m128d MachineEstimate::ln2;
__m128d MachineEstimate::ln2c;
__m128d MachineEstimate::pi_over_2;
__m128d MachineEstimate::pi_over_4;
__m128d MachineEstimate::pi2;
//__m128d MachineEstimate::onethird;
__m128d MachineEstimate::rpi4;
//__m128i MachineEstimate::epi32incorrect = _mm_set1_epi32(0x80000000);
__m128d MachineEstimate::three = _mm_set1_pd(3.0);
__m128d MachineEstimate::four = _mm_set1_pd(4.0);
__m128d MachineEstimate::truehigh;
int MachineEstimate::SavedRoundingMode = _MM_ROUND_NEAREST;
static int initialized = 0;
static int nesting = 0;

static u32 ConstsLong[316] =
{
#define CONSTS_COS 0
        0x00000001, 0xbff00000,
        0xffffffff, 0x3fefffff,
        0x418cafdb, 0x40018bc4,
        0x418cafdb, 0x40018bc4,
        0x9e3185f6, 0xbfe9a7b2,
        0x9e3185f6, 0x3fe9a7b2,
        0x5d05165d, 0x3fbe0270,
        0x5d05165d, 0x3fbe0270,
        0xa8e30653, 0xbf82ce22,
        0xa8e30653, 0x3f82ce22,
        0x70426553, 0x3f3d5450,
        0x70426553, 0x3f3d5450,
        0x749f656f, 0xbeef2edc,
        0x749f656f, 0x3eef2edc,
        0x44c498c8, 0x3e979e4b,
        0x44c498c8, 0x3e979e4b,
#define CONSTS_RPI4 32
        // 1/(Pi*4), multiplication constant
        0x6dc9c883, 0x3fb45f30,
        0x6dc9c882, 0x3fb45f30,
#define CONSTS_PI    36
        // Pi, proper interval
        0x54442d19, 0xc00921fb,
        0x54442d18, 0x400921fb,
#define CONSTS_SIN 40
        0x382d7366, 0xc000c152,
        0x382d7365, 0x4000c152,
        0x791049dc, 0x3ff87fb0,
        0x791049dc, 0x3ff87fb0,
        0x3ea3fdb3, 0xbfd57e24,
        0x3ea3fdb3, 0x3fd57e24,
        0x23972846, 0x3fa1f529,
        0x23972846, 0x3fa1f529,
        0x62748c9e, 0xbf618133,
        0x62748c9e, 0x3f618133,
        0x4e962080, 0x3f165652,
        0x4e962080, 0x3f165652,
        0xe58a04cb, 0xbec4189c,
        0xe58a04cb, 0x3ec4189c,
        0x3772c742, 0x3e6a705b,
        0x3772c742, 0x3e6a705b,
        /*
0x382d7366, 0xc000c152,
0x382d7365, 0x4000c152,
0x791049da, 0x3ff87fb0,
0x791049da, 0x3ff87fb0,
0x3ea3fd41, 0xbfd57e24,
0x3ea3fd41, 0x3fd57e24,
0x2396ef85, 0x3fa1f529,
0x2396ef85, 0x3fa1f529,
0x6259c82a, 0xbf618133,
0x6259c82a, 0x3f618133,
0x42bd2aad, 0x3f165652,
0x42bd2aad, 0x3f165652,
0xe56c4e54, 0xbec4189a,
0xe56c4e54, 0x3ec4189a,
0xb4b3764d, 0x3e6a7057,
0xb4b3764d, 0x3e6a7057,*/
#define CONSTS_EXPMASK 72
        0x00000000, 0x7ff00000,
        0x00000000, 0x7ff00000,
#define CONSTS_EXPBIAS 76
        0x00000000, 0x3fe00000,
        0x00000000, 0x3fe00000,
#define CONSTS_LOG 80
#ifdef LESS_PRECISE
        0x68347bc0, 0x401528b6,
        0x68347cac, 0xc01528b6,    // -
        //0x68347c36, 0xc01528b6,    // +
        0x51b7b3a9, 0xc03e98c3, // +
        0x51b7b3a9, 0x403e98c3, // -
        0x32b71351, 0x4062d98b,
        0x32b71351, 0xc062d98b,
        0x1b5f9043, 0xc08317d4,
        0x1b5f9043, 0x408317d4,
        0xd568356e, 0x409e00f8,
        0xd568356e, 0xc09e00f8,
        0x50bce334, 0xc0b25ddd,
        0x50bce334, 0x40b25ddd,
        0xc8acce73, 0x40c1ac7a,
        0xc8acce73, 0xc0c1ac7a,
        0x465ebb3d, 0xc0cae1af,
        0x465ebb3d, 0x40cae1af,
        0xf709f6ba, 0x40d02ad4,
        0xf709f6ba, 0xc0d02ad4,
        0x2c97c937, 0xc0ce9c26,
        0x2c97c937, 0x40ce9c26,
        0x452abf6b, 0x40c68b7a,
        0x452abf6b, 0xc0c68b7a,
        0x1485d4ea, 0xc0b9506e,
        0x1485d4ea, 0x40b9506e,
        0x497359bb, 0x40a4efd7,
        0x497359bb, 0xc0a4efd7,
        0x40df444c, 0xc0880d63,
        0x40df444c, 0x40880d63,
        0x8424b720, 0x40612455,
        0x8424b720, 0xc0612455,
        0x64920854, 0xc026d616,
        0x64920854, 0x4026d616,
        0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,
#else
        0xb100afab, 0xbc2dcabb,
        0xed268e66, 0xbca3df51,
        0x00000004, 0x3ff00000,
        0x00000004, 0xbff00000,
        0xffffeeea, 0x3fdfffff,
        0xffffeeea, 0xbfdfffff,
        0x55550aa6, 0x3fd55555,
        0x55550aa6, 0xbfd55555,
        0x002a2505, 0x3fd00000,
        0x002a2505, 0xbfd00000,
        0x9a79b04d, 0x3fc99999,
        0x9a79b04d, 0xbfc99999,
        0x0771c502, 0x3fc55555,
        0x0771c502, 0xbfc55555,
        0xfc94cb71, 0x3fc24923,
        0xfc94cb71, 0xbfc24923,
        0x35ec7035, 0x3fc00022,
        0x35ec7035, 0xbfc00022,
        0x9d9f4587, 0x3fbc722e,
        0x9d9f4587, 0xbfbc722e,
        0xab47707f, 0x3fb98a39,
        0xab47707f, 0xbfb98a39,
        0xb4029d62, 0x3fb73291,
        0xb4029d62, 0xbfb73291,
        0x0deb07e2, 0x3fb7085a,
        0x0deb07e2, 0xbfb7085a,
        0xa63cf31c, 0x3fb582e2,
        0xa63cf31c, 0xbfb582e2,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
#define CONSTS_TRUEHIGH 152
        0xffffffff, 0xffffffff,
        0x00000000, 0x00000000,
#define CONSTS_SQRTSQRT2 156
        0xa31b716, 0x3ff306fe,
        0xa31b715, 0x3ff306fe,
#endif
#define CONSTS_LN2    160
        0xfefa39f0, 0xbfe62e42,
        0xfefa39ef, 0x3fe62e42,
#define CONSTS_LN2C    164
        0xfefa39f0, 0x3fe62e42,
        0xfefa39ef, 0x3fe62e42,
#define CONSTS_LOG2E 168
        0x652b82ff, 0x3ff71547,
        0x652b82fe, 0x3ff71547,
#define CONSTS_SIGN 172
        0x00000000, 0x80000000,
        0x00000000, 0x80000000,
#define CONSTS_EXPLIMIT 176
        0x0, 0x408ff000,
        0x0, 0x408ff000,
#define CONSTS_EXP 180
        0x00000001, 0xc0000000,
        0xffffffff, 0x3fffffff,
        0xfefa39f9, 0x3ff62e42,
        0xfefa39f9, 0x3ff62e42,
        0xff82bdb1, 0xbfdebfbd,
        0xff82bdb1, 0x3fdebfbd,
        0xd706fa97, 0x3fbc6b08,
        0xd706fa97, 0x3fbc6b08,
        0x6f5ef210, 0xbf93b2ab,
        0x6f5ef210, 0x3f93b2ab,
        0xf7c7e6fd, 0x3f65d87f,
        0xf7c7e6fd, 0x3f65d87f,
        0x5d0bd9c1, 0xbf34308f,
        0x5d0bd9c1, 0x3f34308f,
        0x722cb340, 0x3efffd04,
        0x722cb340, 0x3efffd04,
        0x43ec690a, 0xbec628a6,
        0x43ec690a, 0x3ec628a6,
        0xab63f5ed, 0x3e8b898b,
        0xab63f5ed, 0x3e8b898b,
        0xdf1599b6, 0xbe4c140c,
        0xdf1599b6, 0x3e4c140c,
        0xc4fc16f9, 0x3e15a8b6,
        0xc4fc16f9, 0x3e15a8b6,
#define CONSTS_ATAN 228
#ifdef LESS_PRECISE
        0x00000001, 0xbff00000,
        0xfffffd8f, 0x3fefffff,
        0x5550366f, 0x3fd55555,
        0x5550366f, 0xbfd55555,
        0x9607c4dd, 0xbfc99999,
        0x9607c4dd, 0x3fc99999,
        0x1414bba2, 0x3fc24924,
        0x1414bba2, 0xbfc24924,
        0xbb222785, 0xbfbc71b4,
        0xbb222785, 0x3fbc71b4,
        0xfb06875f, 0x3fb74500,
        0xfb06875f, 0xbfb74500,
        0xdf9fe793, 0xbfb3ab1c,
        0xdf9fe793, 0x3fb3ab1c,
        0x86d97f40, 0x3fb0f0e4,
        0x86d97f40, 0xbfb0f0e4,
        0xc7e1fb2e, 0xbfad27f5,
        0xc7e1fb2e, 0x3fad27f5,
        0x5754d059, 0x3fa831c1,
        0x5754d059, 0xbfa831c1,
        0x927fbb67, 0xbfa2578e,
        0x927fbb67, 0x3fa2578e,
        0x224d2216, 0x3f97c132,
        0x224d2216, 0xbf97c132,
        0x4b6883d1, 0xbf885ec9,
        0x4b6883d1, 0x3f885ec9,
        0xd9477c17, 0x3f721bb6,
        0xd9477c17, 0xbf721bb6,
        0xc66b658f, 0xbf510dc6,
        0xc66b658f, 0x3f510dc6,
        0x2d256709, 0x3f1e43d5,
        0x2d256709, 0xbf1e43d5,
        0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,
#else
        0x00000001, 0xbff00000,
        0xffffffff, 0x3fefffff,
        0x555553d2, 0x3fd55555,
        0x555553d2, 0xbfd55555,
        0x9998037a, 0xbfc99999,
        0x9998037a, 0x3fc99999,
        0x91f33a63, 0x3fc24924,
        0x91f33a63, 0xbfc24924,
        0x09057800, 0xbfbc71c7,
        0x09057800, 0x3fbc71c7,
        0x1aa24579, 0x3fb745d0,
        0x1aa24579, 0xbfb745d0,
        0xf9b84bf5, 0xbfb3b12a,
        0xf9b84bf5, 0x3fb3b12a,
        0x01a930e2, 0x3fb11089,
        0x01a930e2, 0xbfb11089,
        0x556e5d85, 0xbfae177e,
        0x556e5d85, 0x3fae177e,
        0xa80e2f1b, 0x3faad32f,
        0xa80e2f1b, 0xbfaad32f,
        0xa58c31d6, 0xbfa7ee71,
        0xa58c31d6, 0x3fa7ee71,
        0x8b0ccaa5, 0x3fa4f50b,
        0x8b0ccaa5, 0xbfa4f50b,
        0x6c6308fe, 0xbfa17309,
        0x6c6308fe, 0x3fa17309,
        0x2b4c52ee, 0x3f9a77d1,
        0x2b4c52ee, 0xbf9a77d1,
        0x7e19f3dd, 0xbf916913,
        0x7e19f3dd, 0x3f916913,
        0xfa32033c, 0x3f82da21,
        0xfa32033c, 0xbf82da21,
        0xd33c5aff, 0xbf6fb050,
        0xd33c5aff, 0x3f6fb050,
        0x6bed862f, 0x3f532726,
        0x6bed862f, 0xbf532726,
        0x510269d4, 0xbf2d637e,
        0x510269d4, 0x3f2d637e,
        0x64cd132e, 0x3ef5619e,
        0x64cd132e, 0xbef5619e,
#endif
#define CONSTS_SQRT2    308
        // sqrt(2), mult. constant
        0x667f3bcd, 0x3ff6a09e,
        0x667f3bcc, 0x3ff6a09e,
#define CONSTS_PI2    312
        // Pi, proper interval
        0x54442d19, 0x402921fb,
        0x54442d19, 0x402921fb


        // 16 atan coefficients
        // for angle up to pi/8
        // idea: use 1 step of cordic-like
        // reduction (multiply by 1+i)
        // if ratio is greater
        // than tan(pi/8)
        /*
0x5507f1ac, 0x3c813d00,
0x21a63091, 0xbc6d0a3e,
0x00000005, 0x3ff00000,
0x00000005, 0x3ff00000,
0xefbc06b0, 0x3d4531bd,
0xefbc06b0, 0x3d4531bd,
0x55578c7e, 0xbfd55555,
0x55578c7e, 0xbfd55555,
0xcede2c9d, 0xbde292bd,
0xcede2c9d, 0xbde292bd,
0xbd829637, 0x3fc99999,
0xbd829637, 0x3fc99999,
0xe43f0a84, 0xbea0eb43,
0xe43f0a84, 0xbea0eb43,
0x1a212e5d, 0xbfc248dd,
0x1a212e5d, 0xbfc248dd,
0x5fa3489f, 0xbf18c60e,
0x5fa3489f, 0xbf18c60e,
0x74d114e0, 0x3fbca1e3,
0x74d114e0, 0x3fbca1e3,
0xb4f3ebe3, 0xbf70dcbe,
0xb4f3ebe3, 0xbf70dcbe,
0x88cf479a, 0xbfb2f5b0,
0x88cf479a, 0xbfb2f5b0,
0x6ee43ee8, 0xbfa99dde,
0x6ee43ee8, 0xbfa99dde,
0x6ceef691, 0x3fc74859,
0x6ceef691, 0x3fc74859,
0x7d996ad2, 0xbfc29b9b,
0x7d996ad2, 0xbfc29b9b,
0x1bc49f36, 0x3fa569b3,
0x1bc49f36, 0x3fa569b3
         */
};

#define MY_PI    3.14159265358979323846264338328
#define MY_LN2 0.69314718055994530941723212145817656807550013436025

#ifndef _MSC_VER
#define _nextafter nextafter
#endif

// these two should be callable multiple times. 
int MachineEstimate::BeginComputation()
{
    if (!initialized) {
        double z = 0.0;
        double minusinf = -1.0 / z;
        double plusinf = 1.0 / z;
        sqrt_corr = //_mm_mul_pd(
                _mm_set_pd(0.0, minusinf);//, zero);
        // note: we're nudging the first coefficient a bit more to accommodate the approximation error
        initialized = 1;
        memcpy(coeff_cos, ConstsLong + CONSTS_COS, sizeof(__m128d)*8);
        memcpy(&rpi4, ConstsLong + CONSTS_RPI4, sizeof(__m128d)*1);
        memcpy(&pi, ConstsLong + CONSTS_PI, sizeof(__m128d)*1);
        memcpy(&pi2, ConstsLong + CONSTS_PI2, sizeof(__m128d)*1);
        pi_over_2 = _mm_mul_pd(pi, _mm_set_pd(0.5, 0.5));        // exact multiplication
        pi_over_4 = _mm_mul_pd(pi_over_2, _mm_set_pd(0.5, 0.5));        // exact multiplication
        memcpy(coeff_sin, ConstsLong + CONSTS_SIN, sizeof(__m128d)*8);
        memcpy(coeff_log, ConstsLong + CONSTS_LOG, sizeof(__m128d)*20);
        memcpy(coeff_exp, ConstsLong + CONSTS_EXP, sizeof(__m128d)*12);
        memcpy(coeff_atan, ConstsLong + CONSTS_ATAN, sizeof(__m128d)*20);
        memcpy(&expmask, ConstsLong + CONSTS_EXPMASK, sizeof(__m128d)*1);
        memcpy(&expbias, ConstsLong + CONSTS_EXPBIAS, sizeof(__m128d)*1);
        memcpy(&ln2, ConstsLong + CONSTS_LN2, sizeof(__m128d)*1);
        memcpy(&ln2c, ConstsLong + CONSTS_LN2C, sizeof(__m128d)*1);
        memcpy(&log2e, ConstsLong + CONSTS_LOG2E, sizeof(__m128d)*1);
        memcpy(&sign, ConstsLong + CONSTS_SIGN, sizeof(__m128d)*1);
        memcpy(&explimit, ConstsLong + CONSTS_EXPLIMIT, sizeof(__m128d)*1);
        memcpy(&sqrt2, ConstsLong + CONSTS_SQRT2, sizeof(__m128d)*1);
        memcpy(&sqrtsqrt2, ConstsLong + CONSTS_SQRTSQRT2, sizeof(__m128d)*1);
        memcpy(&truehigh, ConstsLong + CONSTS_TRUEHIGH, sizeof(__m128d)*1);
    }
#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
    if (nesting && 
            (_MM_GET_EXCEPTION_STATE() &
                    (_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW)))
        throw PrecisionException("SSE exception");
    _MM_SET_EXCEPTION_STATE(0);
#endif

    SavedRoundingMode = _MM_GET_ROUNDING_MODE();
    _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);

    ++nesting;
    return SavedRoundingMode;
}

void MachineEstimate::FinishComputation(int rm)
{
    if (!rm) rm = SavedRoundingMode;
    assert(_MM_GET_ROUNDING_MODE()==_MM_ROUND_DOWN);
    _MM_SET_ROUNDING_MODE(rm);
    --nesting;
}

static inline
__m128d MulPositiveConstant(const __m128d &lhs, const __m128d &rhs) 
// for this the constant rhs has to be prepared without a negation of the high bound
{     
    __m128d a = _mm_shuffle_pd(rhs, rhs, 1);    // d, c
    __m128d b = _mm_mul_pd(lhs, rhs);        // ac, -bd        
    // correct for positive lhs
    __m128d c = _mm_mul_pd(lhs, a);                    // bc, -ad        
    // correct for negative lhs
    __m128d d = _mm_min_pd(b, c);
    // also correct for lhs that includes zero
    return d;
}


MachineEstimate sin(const MachineEstimate &x)
{
    // idea: use the center and finally add the input error (because
    // sin is lipschitz 1).
    // first remove periodicity by dividing by pi*2,
    // taking the center and rounding (now x ranges [-pi,pi));
    // then divide by 3, use polynomial approximation
    // of sin(|x|)/|x| on up to pi/3 (8 coefficients suffice)
    // then multiply by x/3 to recover sign and get the full sin(x/3)
    // finally use sin (3*x) = sin(x) * (3 - 4*sin(x)*sin(x)) and
    // add the error of the argument

    // (weak_AsDouble uses a multiplication by 1/2 which we avoid
    //    by using 1/(pi*4) instead of 1/(pi*2); we also avoid the 
    //    division by 3 by incorporating it in the coefficients)
    __m128d z = MulPositiveConstant(x.interval, MachineEstimate::rpi4);
    __m128d v = _mm_shuffle_pd(z, z, 1);
    __m128d r = _mm_add_pd(z, v);

    // the rounding in the following could be causing problems
    z = _mm_sub_pd(z, v);
    // account for it either using doubled r (i.e. adjusted pi2) or through the
    // next line (slower, but more efficient error propagation)
    //r = _mm_add_pd(r, _mm_add_pd(z, _mm_shuffle_pd(z, z, 1)));

    // could the addition here cause trouble?
    int i = _mm_cvtsd_si32(_mm_add_sd(z, MachineEstimate::half));    // towards -inf, 
    // adjusted by 1/2 to get rounding to nearest
    v = _mm_cvtsi32_sd(v, i);    // exact result
#ifndef REALLIB_RELY_ON_SSE_EXCEPTIONS
    if (i == 0x80000000) // i.e. overflow occurred
        throw PrecisionException("sin");
#endif
    // make a positive z, abs(z-v)

    z = _mm_sub_sd(z, v);    // this is an exact operation
    bool b = (_mm_movemask_pd(z) & 1);
    z = _mm_or_pd(z, MachineEstimate::signmask);
    v = _mm_xor_pd(z, MachineEstimate::signmask);
    z = _mm_shuffle_pd(z, v, 0);
    __m128d w = _mm_shuffle_pd(v, v, 0);        // no sign here

    // z has (-pi, pi] mapped into (-0.5, 0.5]
    z = _mm_mul_pd(z, w);    // correctly rounded
    __m128d z2 = _mm_mul_pd(z, _mm_xor_pd(z, MachineEstimate::signmask)); // correctly rounded
    z = _mm_shuffle_pd(z, z, 1);                // to avoid negations of s0-s3
    __m128d s3 = _mm_mul_pd(MachineEstimate::coeff_sin[7], z);
    v = _mm_xor_pd(z2, MachineEstimate::signmask);
    __m128d s2 = _mm_mul_pd(MachineEstimate::coeff_sin[5], z);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_sin[6]);
    __m128d s1 = _mm_mul_pd(MachineEstimate::coeff_sin[3], z);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_sin[4]);
    __m128d s0 = _mm_mul_pd(MachineEstimate::coeff_sin[1], z);            // coeff[1] does not have a negated high
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_sin[2]);
    s3 = _mm_mul_pd(s3, v);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_sin[0]);                    // coeff[0] is a proper interval, i.e. high part negated
    z2 = _mm_mul_pd(z2, v);
    s2 = _mm_add_pd(s2, s3);
    z2 = _mm_xor_pd(z2, MachineEstimate::signmask);
    s1 = _mm_mul_pd(s1, v);
    s0 = _mm_add_pd(s0, s1);
    s2 = _mm_mul_pd(s2, z2);
    s0 = _mm_add_pd(s0, s2);
    s0 = _mm_mul_pd(s0, w);

    r = _mm_mul_pd(r, MachineEstimate::pi2);                    

    s2 = _mm_xor_pd(s0, MachineEstimate::signmask);
    s3 = _mm_mul_pd(s0, MachineEstimate::three);
    s1 = _mm_mul_pd(s2, MachineEstimate::four);    // rounding-free multiplication
    s3 = _mm_shuffle_pd(s3, s3, 1);
    s0 = _mm_mul_pd(s0, s2);
    s3 = _mm_add_pd(s3, r);                // add error from argument
    s0 = _mm_mul_pd(s0, s1);
    s0 = _mm_add_pd(s0, s3);            // it's important that this substraction
    // is the last operation
    // otherwise negative values could enter
    // the computation and require more
    // complicated multiplications

    return b ? -MachineEstimate(s0) : s0; 

}

MachineEstimate cos(const MachineEstimate &x)
{
    // very similar to sin, but the polynomial approximation
    // is not multiplied by x.
    __m128d z = _mm_mul_pd(abs(x).interval, MachineEstimate::rpi4);
    __m128d v = _mm_shuffle_pd(z, z, 1);
    __m128d r = _mm_add_pd(z, v);

    // the rounding in the following could be causing problems
    z = _mm_sub_pd(z, v);
    // account for it either using doubled r (i.e. adjusted pi2) or through the
    // next line (slower, but more efficient error propagation)
    //r = _mm_add_pd(r, _mm_add_pd(z, _mm_shuffle_pd(z, z, 1)));

    // could the addition here cause trouble?
    int i = _mm_cvtsd_si32(_mm_add_sd(z, MachineEstimate::half));    // towards -inf, 
    // adjusted by 1/2 to get rounding to nearest
    v = _mm_cvtsi32_sd(v, i);    // exact result
#ifndef REALLIB_RELY_ON_SSE_EXCEPTIONS
    if (i == 0x80000000) // i.e. overflow occurred
        throw PrecisionException("cos");
#endif
    // make a positive z, abs(z-v)

    z = _mm_sub_sd(z, v);            // *** rounding?
    z = _mm_or_pd(z, MachineEstimate::signmask);
    v = _mm_xor_pd(z, MachineEstimate::signmask);
    z = _mm_shuffle_pd(z, v, 0);

    // z has (-pi, pi] mapped into (-1, 1]
    z = _mm_mul_pd(z, _mm_xor_pd(z, MachineEstimate::signmask));    // correctly rounded
    __m128d z2 = _mm_mul_pd(z, _mm_xor_pd(z, MachineEstimate::signmask)); // correctly rounded
    z = _mm_shuffle_pd(z, z, 1);                // to avoid negations of s0-s3
    __m128d s3 = _mm_mul_pd(MachineEstimate::coeff_cos[7], z);
    v = _mm_xor_pd(z2, MachineEstimate::signmask);
    __m128d s2 = _mm_mul_pd(MachineEstimate::coeff_cos[5], z);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_cos[6]);
    __m128d s1 = _mm_mul_pd(MachineEstimate::coeff_cos[3], z);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_cos[4]);
    __m128d s0 = _mm_mul_pd(MachineEstimate::coeff_cos[1], z);            // coeff[1] does not have a negated high
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_cos[2]);
    s3 = _mm_mul_pd(s3, v);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_cos[0]);                    // coeff[0] is a proper interval, i.e. high part negated
    z2 = _mm_mul_pd(z2, v);
    s2 = _mm_add_pd(s2, s3);
    z2 = _mm_xor_pd(z2, MachineEstimate::signmask);
    s1 = _mm_mul_pd(s1, v);
    s0 = _mm_add_pd(s0, s1);
    s2 = _mm_mul_pd(s2, z2);
    s0 = _mm_add_pd(s0, s2);

    r = _mm_mul_pd(r, MachineEstimate::pi2);

    s2 = _mm_xor_pd(s0, MachineEstimate::signmask);
    s3 = _mm_mul_pd(s0, MachineEstimate::three);
    s1 = _mm_mul_pd(s2, MachineEstimate::four);    // rounding-free multiplication
    s3 = _mm_shuffle_pd(s3, s3, 1);
    s0 = _mm_mul_pd(s0, s2);
    s3 = _mm_add_pd(s3, r);                // add error from argument
    s0 = _mm_mul_pd(s0, s1);
    s0 = _mm_add_pd(s0, s3);            // it's important that this substraction
    // is the last operation
    // otherwise negative values could enter
    // the computation and require more
    // complicated multiplications

    return s0;

}

MachineEstimate log(const MachineEstimate &arg)
{
    // simultaneously on both components:
    // separate mantissa and exponent
    // use polynomial approximation on log2(mantissa)
    // (16 coefficients with adjusted signs so that 
    // the roundings are in the correct directions, not full precision)
    // add exponent (with removed bias)
    // finally multiply by ln 2

    // important: the computation of the polynomial is not a computation
    // on the interval; we compute the polynomial on the lower bound with
    // all roundings down to get a lower bound for the value on the interval
    // and similarly for the higher bound with roundings up.
    // this suffices, because the function is monotone

    // in other words if
    // p(x) = c_n x^n + ... + c_1 x + c_0
    // p(x)-e <= f(x) <= p(x)+e
    // monotone f(x)
    // and x_l <= x <= x_h
    // then
    // p(x_l)-e <= f(x_l) <= p(x_l)+e
    // p(x_h)-e <= f(x_h) <= p(x_h)+e
    // f(x_l) <= f(x) <= f(x_h)
    // and then p_d(x_l) <= p(x_l)-e <= f(x) <= p(x_h)+e <= p_u(x_h)
    // where p_d and p_u are, respectively, p computed via rounding down
    // with an adjusted c_0_d <= c_0-e, and via rounding up and
    // c_0_u >= c_0+e. All other coefficients are not changed (because
    // they are given exact)

    // we only need to ensure we're rounding the right way, i.e.
    // if we compute a negative product, rounding up requires that
    // both multiples were rounded down in absolute value
    // if we compute a positive product, both multiples must have
    // been rounded up in absolute value
    // this is a mess
    __m128d a(arg.interval);
    if (_mm_movemask_pd(a) != 1) {
        if (!(_mm_movemask_pd(a) & 1)) // -high > 0 means provably negative argument
            throw DomainException("log");
        else throw PrecisionException("log");
    }
    __m128i i;
    *((__m128d*)&i) = a;

#ifdef LESS_PRECISE
    i = _mm_and_si128(i, MachineEstimate::expmask);    // extract exponent
    __m128d z = _mm_andnot_pd(*((__m128d*)&MachineEstimate::expmask), a); // clear exponent
    i = _mm_sub_epi64(i, MachineEstimate::expbias); // adjust for bias
    z = _mm_or_pd(z, *((__m128d*)&MachineEstimate::expbias));        // set exponent to 0, i.e. value within 0.5 and 1
    i = _mm_srai_epi32(i, 20);                                // move to unit place
    __m128d v = _mm_xor_pd(z, MachineEstimate::signmask);
    i = _mm_shuffle_epi32(i, _MM_SHUFFLE(3, 1, 3, 1)); // move values to places 0 and 1
    z = _mm_mul_pd(v, z);
    __m128d e(_mm_cvtepi32_pd(i));                        // e now holds the exponent value for both sides
    __m128d s0, s1, s2, s3;

    e = _mm_xor_pd(e, MachineEstimate::signmask);

    __m128d z4(z);
    z = _mm_xor_pd(z, MachineEstimate::signmask);
    // do we need to shuffle z?
    //z = _mm_shuffle_pd(z, z, 1);

    s3 = _mm_mul_pd(v, MachineEstimate::coeff_log[15]);    // coeff[15] is -,+
    s2 = _mm_mul_pd(v, MachineEstimate::coeff_log[13]);
    z4 = _mm_mul_pd(z4, z);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_log[14]);                            // coeff[13] is -,+
    s1 = _mm_mul_pd(v, MachineEstimate::coeff_log[11]);
    __m128d z2 = _mm_xor_pd(z4, MachineEstimate::signmask);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[12]);
    s0 = _mm_mul_pd(v, MachineEstimate::coeff_log[9]);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[10]);
    s3 = _mm_mul_pd(s3, z);
    z4 = _mm_mul_pd(z4, z2);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[8]);
    z4 = _mm_xor_pd(z4, MachineEstimate::signmask);
    s1 = _mm_mul_pd(s1, z);
    s2 = _mm_add_pd(s2, s3);
    s3 = _mm_mul_pd(v, MachineEstimate::coeff_log[7]);    // +,-
    s0 = _mm_add_pd(s0, s1);
    s1 = _mm_mul_pd(v, MachineEstimate::coeff_log[3]);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_log[6]);    // +,-
    a = _mm_mul_pd(s2, z2);
    s2 = _mm_mul_pd(v, MachineEstimate::coeff_log[5]);
    a = _mm_add_pd(s0, a);
    s3 = _mm_mul_pd(s3, z);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[2]);
    s0 = _mm_mul_pd(v, MachineEstimate::coeff_log[1]);
    a = _mm_mul_pd(a, z4);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[4]);
    s1 = _mm_mul_pd(s1, z);
    s2 = _mm_add_pd(s2, s3);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[0]);
    s2 = _mm_mul_pd(s2, z2);
    s0 = _mm_add_pd(s0, s1);
    s0 = _mm_add_pd(s0, s2);
    a = _mm_add_pd(a, s0);
    a = _mm_add_pd(a, e);
    a = MulPositiveConstant(a, MachineEstimate::ln2c);
#else
    i = _mm_and_si128(i, MachineEstimate::expmask);    // extract exponent
    __m128d z = _mm_andnot_pd(*((__m128d*)&MachineEstimate::expmask), a); // clear exponent
    i = _mm_sub_epi64(i, MachineEstimate::expbias); // adjust for bias
    z = _mm_or_pd(z, *((__m128d*)&MachineEstimate::expbias));        // set exponent to 0, i.e. value within 0.5 and 1
    i = _mm_srai_epi32(i, 20);                                // move to unit place
    i = _mm_shuffle_epi32(i, _MM_SHUFFLE(3, 1, 3, 1)); // move values to places 0 and 1
    __m128d e(_mm_cvtepi32_pd(i));                        // e now holds the exponent value for both sides
    __m128d s0, s1, s2, s3;

    e = _mm_xor_pd(e, MachineEstimate::signmask);

    // two reductions using ln (x*y) = ln x + ln y
    // using q=2^(1/2^i) as multiples, we take xq if smaller
    // than one and just x if not
    // we reflect this by adding -1/2^i to the exponent part
    // for the values we've switched
    __m128d v = _mm_mul_pd(z, MachineEstimate::sqrt2);
    __m128d m = _mm_cmpgt_pd(v, MachineEstimate::one);            // one is an exact constant, sub works instead of shuffle+add
    m = _mm_xor_pd(m, MachineEstimate::truehigh);
    z = _mm_and_pd(m, z);
    v = _mm_andnot_pd(m, v);
    s0 = _mm_add_pd(e, MachineEstimate::mhalf);
    z = _mm_or_pd(z, v);
    e = _mm_and_pd(m, e);
    s0 = _mm_andnot_pd(m, s0);
    e = _mm_or_pd(e, s0);

    v = _mm_mul_pd(z, MachineEstimate::sqrtsqrt2);
    m = _mm_cmpgt_pd(v, MachineEstimate::one);            // one is an exact constant, sub works instead of shuffle+add
    m = _mm_xor_pd(m, MachineEstimate::truehigh);
    z = _mm_and_pd(m, z);
    v = _mm_andnot_pd(m, v);
    s0 = _mm_add_pd(e, MachineEstimate::mquarter);
    z = _mm_or_pd(z, v);
    e = _mm_and_pd(m, e);
    s0 = _mm_andnot_pd(m, s0);
    e = _mm_or_pd(e, s0);

    // the rounding up until here has been in the normal direction
    // we want a negative mantissa part in the end
    // after this substraction the high part is rounded towards -inf
    // and the low towards +inf. We keep this sort of rounding until the
    // end, so that multiplications with negative constants yield
    // correctly rounded products

    m = MachineEstimate::signmask;
    m = _mm_shuffle_pd(m, m, 1);
    // note this substraction has the wrong order: we
    // want to compute 1 - z. this order is to get the correct rounding
    // and the correct signs of the results
    z = _mm_sub_pd(z, MachineEstimate::one);            

    v = _mm_xor_pd(z, m);            // multiplication-friendly format (+,+)
    z = _mm_mul_pd(v, z);

    __m128d z4(z);
    z = _mm_xor_pd(z, m);

    s2 = _mm_mul_pd(v, MachineEstimate::coeff_log[13]);
    z4 = _mm_mul_pd(z4, z);
    s1 = _mm_mul_pd(v, MachineEstimate::coeff_log[11]);
    __m128d z2 = _mm_xor_pd(z4, m);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[12]);
    s0 = _mm_mul_pd(v, MachineEstimate::coeff_log[9]);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[10]);
    __m128d z8 = _mm_mul_pd(z4, z2);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[8]);
    z4 = _mm_xor_pd(z8, m);
    s1 = _mm_mul_pd(s1, z);
    s3 = _mm_mul_pd(v, MachineEstimate::coeff_log[7]);    // +,-
    s0 = _mm_add_pd(s0, s1);
    s1 = _mm_mul_pd(v, MachineEstimate::coeff_log[3]);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_log[6]);    // +,-
    a = _mm_mul_pd(s2, z2);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[2]);
    s2 = _mm_mul_pd(v, MachineEstimate::coeff_log[5]);
    a = _mm_add_pd(s0, a);
    s3 = _mm_mul_pd(s3, z);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[4]);
    s0 = _mm_mul_pd(v, MachineEstimate::coeff_log[1]);
    a = _mm_mul_pd(a, z4);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[0]);
    s1 = _mm_mul_pd(s1, z);
    s2 = _mm_add_pd(s2, s3);
    s2 = _mm_mul_pd(s2, z2);
    s0 = _mm_add_pd(s0, s2);
    e = MulPositiveConstant(e, MachineEstimate::ln2c);
    a = _mm_add_pd(a, s0);
    e = _mm_add_pd(e, s1);
    a = _mm_add_pd(e, a);
#endif    // less precise
    return a;
}

MachineEstimate exp(const MachineEstimate &arg)
{
    // simultaneously for both components:
    // multiply by log2(e)
    // separate integer and fractional part
    // use 12-coefficient polynomial approximation
    // for 2^(fractional part)*2
    // make 0.5*2^(integer part) by moving some bits to the right place
    // finally multiply the two and you have the result

    // *** sign job: need to do some thinking what happens if the input is negative
    // should be no problem    

    __m128d a = MulPositiveConstant(arg.interval, MachineEstimate::log2e);
    __m128d e = _mm_xor_pd(a, MachineEstimate::signmask);
    __m128i i = _mm_cvtpd_epi32(e);
    __m128d d = _mm_andnot_pd(MachineEstimate::sign, a);    // clear signs

    __m128d c = _mm_cvtepi32_pd(i);
    d = _mm_sub_pd(d, MachineEstimate::explimit);        // we don't care about the rounding here
    c = _mm_xor_pd(c, MachineEstimate::signmask);
    i = _mm_slli_epi32(i, 20);
    a = _mm_sub_pd(a, c);        // correct rounding
    // a has a positive interval-type value (i.e. (low, -high), but high need not be > low)
    if (_mm_movemask_pd(d) != 3) 
        throw PrecisionException("exp");

    e = _mm_mul_pd(a, MachineEstimate::coeff_exp[11]);    // +,+            // 11x
    __m128d f = _mm_xor_pd(a, MachineEstimate::signmask);
    d = _mm_mul_pd(a, MachineEstimate::coeff_exp[9]);                                 // 9x
    e = _mm_add_pd(e, MachineEstimate::coeff_exp[10]); // +,-            // +- 10 + 11x
    f = _mm_mul_pd(f, a);                                            // +- xx
    d = _mm_add_pd(d, MachineEstimate::coeff_exp[8]);                        // +- 8 + 9x
    __m128d g = _mm_xor_pd(f, MachineEstimate::signmask);    // ++ xx
    c = _mm_mul_pd(a, MachineEstimate::coeff_exp[7]);                        // 7x
    __m128d b = _mm_mul_pd(a, MachineEstimate::coeff_exp[5]);            // 5x
    e = _mm_mul_pd(e, g);                                                            // +- 10xx + 11xxx
    c = _mm_add_pd(c, MachineEstimate::coeff_exp[6]);                        // +- 6 + 7x
    d = _mm_add_pd(d, e);                                                            // +- 8 + 9x + 10xx + 11xxx
    e = _mm_mul_pd(a, MachineEstimate::coeff_exp[3]);                        // 3x
    f = _mm_mul_pd(f, g);                                            // +- xxxx
    b = _mm_add_pd(b, MachineEstimate::coeff_exp[4]);                        // 4 + 5x 
    e = _mm_add_pd(e, MachineEstimate::coeff_exp[2]);                        // 2 + 3x
    c = _mm_mul_pd(c, g);                                                            // 6xx + 7xxx
    __m128d h = _mm_xor_pd(f, MachineEstimate::signmask); // ++ xxxx
    e = _mm_mul_pd(e, g);                                                            // 2xx + 3xxx
    b = _mm_add_pd(b, c);                                                            // 4 + 5x + 6xx + 7xxx
    c = _mm_mul_pd(a, MachineEstimate::coeff_exp[1]);                        // 1x
    f = _mm_mul_pd(f, h);                                            // +- xxxxxxxx
    i = _mm_shuffle_epi32(i, _MM_SHUFFLE(1, 3, 0, 2));
    c = _mm_add_pd(c, MachineEstimate::coeff_exp[0]);                        // 0 + 1x
    b = _mm_mul_pd(b, h);                                                            // 4xxxx + 5xxxxx + 6xxxxxx + 7xxxxxxx
    i = _mm_add_epi32(i, MachineEstimate::expbias);
    // i currently holds 1.0 * 2^round(arg)
    c = _mm_add_pd(c, e);                                                            // 0 + 1x + 2xx + 3xxx
    f = _mm_xor_pd(f, MachineEstimate::signmask); // a^8    // ++ xxxxxxxx
    d = _mm_mul_pd(d, f);                                                            // 8xxxxxxxx + 9xxxxxxxxx + 10xxxxxxxxxx + 11xxxxxxxxxxx
    c = _mm_add_pd(c, b);                                                            // 0..7
    a = *((__m128d*)(&i));
    c = _mm_add_pd(c, d);            // final result 0..11

    a = _mm_mul_pd(a, c);
    return a;

}

// atan for 0.0 <= argument <= 1.0
static
MachineEstimate atanprimary(const MachineEstimate &x)
{
    // no problems with negative multiplications because the partial results
    // are always positive (the negative coefficients are always smaller than
    // the previous positive and are additionally multiplied by (0, 1])
    __m128d a = _mm_xor_pd(x.interval, MachineEstimate::signmask);
    __m128d b = _mm_mul_pd(a, x.interval);            // xx
    __m128d c = _mm_xor_pd(b, MachineEstimate::signmask); // xx ++
    __m128d s3 = _mm_mul_pd(c, MachineEstimate::coeff_atan[15]); // -+
    b = _mm_mul_pd(b, c);    // xxxx
    __m128d s2 = _mm_mul_pd(c, MachineEstimate::coeff_atan[13]);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_atan[14]);    // +-
    __m128d s1 = _mm_mul_pd(c, MachineEstimate::coeff_atan[11]);
    __m128d d = _mm_xor_pd(b, MachineEstimate::signmask);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_atan[12]);
    s3 = _mm_mul_pd(s3, d);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_atan[10]);
    __m128d s0 = _mm_mul_pd(c, MachineEstimate::coeff_atan[9]);
    s2 = _mm_add_pd(s2, s3);
    s3 = _mm_mul_pd(c, MachineEstimate::coeff_atan[7]);
    b = _mm_mul_pd(b, d);    // xxxxxxxx
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_atan[8]);
    s1 = _mm_mul_pd(s1, d);
    __m128d e = _mm_xor_pd(b, MachineEstimate::signmask);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_atan[6]);
    s2 = _mm_mul_pd(s2, e);
    s0 = _mm_add_pd(s0, s1);
    b = _mm_mul_pd(b, e);
    s1 = _mm_mul_pd(c, MachineEstimate::coeff_atan[3]);
    s0 = _mm_add_pd(s0, s2);
    __m128d f = _mm_xor_pd(b, MachineEstimate::signmask);
    s2 = _mm_mul_pd(c, MachineEstimate::coeff_atan[5]);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_atan[2]);
    s3 = _mm_mul_pd(s3, d);
    s0 = _mm_mul_pd(s0, f);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_atan[4]);
    s1 = _mm_mul_pd(s1, d);
    s2 = _mm_add_pd(s2, s3);
#ifdef LESS_PRECISE
    s3 = _mm_mul_pd(c, MachineEstimate::coeff_atan[1]);
    s0 = _mm_add_pd(s0, s1);
    s2 = _mm_mul_pd(s2, e);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_atan[0]);
    s0 = _mm_add_pd(s0, s2);
    s0 = _mm_add_pd(s0, s3);
#else
    s3 = _mm_mul_pd(c, MachineEstimate::coeff_atan[1]);
    s0 = _mm_add_pd(s0, s1);
    s1 = _mm_mul_pd(c, MachineEstimate::coeff_atan[19]);
    s2 = _mm_mul_pd(s2, e);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_atan[18]);
    b = _mm_mul_pd(b, f);
    s0 = _mm_add_pd(s0, s2);
    s2 = _mm_mul_pd(c, MachineEstimate::coeff_atan[17]);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_atan[0]);
    s1 = _mm_mul_pd(s1, d);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_atan[16]);
    b = _mm_xor_pd(b, MachineEstimate::signmask);
    s2 = _mm_add_pd(s2, s1);
    s0 = _mm_add_pd(s0, s3);
    s2 = _mm_mul_pd(s2, b);
    s0 = _mm_add_pd(s0, s2);
#endif

    s0 = _mm_mul_pd(s0, a);
    return s0;
}

MachineEstimate atan(const MachineEstimate &x)
{
    if (x.IsPositive()) {
        if (x > 1.0) return MachineEstimate::pi_over_2 - atanprimary(recip(x));
        else if (x < 1.0) return atanprimary(x);
        else // atan is Lipschitz 1. Just return pi/4 + twice the error from argument
            return _mm_add_pd(MachineEstimate::pi_over_4, x.MinusDiff());
    } else if (x.IsNegative()) {
        MachineEstimate y(-x);
        if (y > 1.0) return atanprimary(recip(y)) - MachineEstimate::pi_over_2;
        else if (y < 1.0) return -atanprimary(y);
        else // atan is Lipschitz 1. Just return pi/4 + twice the error from argument
            return -MachineEstimate(_mm_add_pd(MachineEstimate::pi_over_4, x.MinusDiff()));
    } else {
        return x.MinusDiff();
    }
}

MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x)
{
    if (x.IsPositive()) {
        if (y.IsPositive()) {
            MachineEstimate d(x-y);
            if (d.IsPositive()) return atanprimary(y/x);
            else if (d.IsNegative()) return MachineEstimate::pi_over_2 - atanprimary(x/y);
            else {
                MachineEstimate r(y/x);
                return _mm_add_pd(MachineEstimate::pi_over_4, r.MinusDiff());
            }
        } else if (y.IsNegative()) {
            MachineEstimate d(x+y);
            if (d.IsPositive()) return -atanprimary((-y)/x);
            else if (d.IsNegative()) return atanprimary(x/(-y)) - MachineEstimate::pi_over_2;
            else {
                MachineEstimate r(y/x);        // sign doesn't matter here
                return -MachineEstimate(_mm_add_pd(MachineEstimate::pi_over_4, r.MinusDiff()));
            }
        } else {
            MachineEstimate r(y/x);
            return r.MinusDiff();
        }
    } else 
        // no check needed because x close to zero would either not cause a problem
        // if y is big enough, or will cause a PrecisionException at the moment of
        // the division y/x or directly
        //if (x.IsNegative()) 
    {
        if (y.IsPositive()) {
            MachineEstimate d(-x-y);
            if (d.IsPositive()) return MachineEstimate::pi - atanprimary(y/(-x));
            else if (d.IsNegative()) return MachineEstimate::pi_over_2 + atanprimary((-x)/y);
            else {
                MachineEstimate r(y/x);
                return _mm_add_pd(
                        _mm_add_pd(MachineEstimate::pi, MachineEstimate::pi_over_4),
                        r.MinusDiff());
            }
        } else if (y.IsNegative()) {
            MachineEstimate d(y-x);
            if (d.IsPositive()) return atanprimary(y/x) - MachineEstimate::pi;
            else if (d.IsNegative()) return -(atanprimary(x/y) + MachineEstimate::pi_over_2);
            else {
                MachineEstimate r(y/x);        // sign doesn't matter here
                return -MachineEstimate(_mm_add_pd(
                        _mm_add_pd(MachineEstimate::pi, MachineEstimate::pi_over_4),
                        r.MinusDiff()));
            }
        } else {
            throw PrecisionException("atan2");
        }
    }
}

MachineEstimate asin(const MachineEstimate &x)
{
    return atan2(x, sqrt(1.0 - sq(x)));
}

MachineEstimate acos(const MachineEstimate &x)
{
    return atan2(sqrt(1.0 - sq(x)), x);
}

/*
MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x)
{
    // idea: work simultaneously on both parts
    // remove signs and set some flags to know how to correct
    // the result
    // check if input contains 0 for y with a non-positive x, raise
    // an exception if this is the case

    // see which input has the greater value, adjust additive and flag if needed
    // divide the smaller by the larger
    // (under consideration) if the value is greater than tan(pi/8), use complex
    //            multiplication by 1+i, examine sign of y and multiply by constant to
    //            recover ratio
    // use polynomial approximation with 20 coefficients (not full precision) or 
    // 16 (in the case under consideration).
    // finally use sign and additive to get result.

    __m128d a, b;
    __m128d sgn = _mm_xor_pd(x.interval, y.interval);                    // this is the sign for the operation on the computed value
                        // there can also be an additive and this only depends on the sign of x
                        // but its sign depends on the sign of y
    __m128d d = _mm_xor_pd(x.interval, MachineEstimate::signmask);
    __m128d e = _mm_cmplt_pd(d, MachineEstimate::zero);    
    __m128d f = _mm_and_pd(e, MachineEstimate::pi);            // +,-
    __m128d g = _mm_or_pd(f, MachineEstimate::zero);
    sgn = _mm_and_pd(sgn, MachineEstimate::sign);
    __m128d adt = _mm_xor_pd(g, sgn);
            // adt is additive, sgn is the sign

    __m128d c = _mm_andnot_pd(x.interval, MachineEstimate::sign);
    __m128d h = _mm_andnot_pd(y.interval, MachineEstimate::sign);
    a = _mm_cmpgt(c, h);
    // swap the values, change the sign and adjust the additive if x > y
    b = _mm_and_pd(h, a);
    d = _mm_andnot_pd(c, a);
    e = _mm_and_pd(c, a);
    f = _mm_andnot_pd(h, a);
    g = _mm_or_pd(b, d);            // greater value
    h = _mm_or_pd(e, f);            // smaller value
    f = _mm_xor_pd(MachineEstimate::pio2, sgn);            // +,-
    g = _mm_xor_pd(g, MachineEstimate::signmask);        // for rounding in the next operation
    sgn = _mm_xor_pd(sgn, a);
    h = _mm_div_pd(g, h);
    f = _mm_and_pd(f, a);
    sgn = _mm_and_pd(sgn, MachineEstimate::sign);
    adt = _mm_add_pd(adt, f);            // correct rounding

    // what we have here:
    // h: ratio within 0,1    (+,-)
    // sgn: sign to be xored
    // adt: constant to be added to the result (+,-)
    g = _mm_xor_pd(h, MachineEstimate::signmask);
    h = _mm_mul_pd(h, g);

    // compute poly on h

    // multiply by g

    // xor sign
    // add constant

    // check if high > low, ow. raise precisionexception

    // sign job: need to do some thinking what happens if the input is negative


}*/

/*
MachineEstimate logg(const MachineEstimate &arg)
{
    __m128d a(arg.interval);
    if (_mm_movemask_pd(a) != 1) {
        if (!(_mm_movemask_pd(a) & 1)) // -high > 0 means provably negative argument
            throw DomainException("log");
        else throw PrecisionException("log");
    }
    __m128i i;
 *((__m128d*)&i) = a;

    i = _mm_and_si128(i, MachineEstimate::expmask);    // extract exponent
    __m128d z = _mm_andnot_pd(*((__m128d*)&MachineEstimate::expmask), a); // clear exponent
    i = _mm_sub_epi64(i, MachineEstimate::expbias); // adjust for bias
    z = _mm_or_pd(z, *((__m128d*)&MachineEstimate::expbias));        // set exponent to 0, i.e. value within 0.5 and 1
    i = _mm_srai_epi32(i, 20);                                // move to unit place
    __m128d v = _mm_xor_pd(z, MachineEstimate::signmask);
    i = _mm_shuffle_epi32(i, _MM_SHUFFLE(3, 1, 3, 1)); // move values to places 0 and 1
    z = _mm_mul_pd(v, z);
    __m128d e(_mm_cvtepi32_pd(i));                        // e now holds the exponent value for both sides
    __m128d s0, s1, s2, s3;
    // think about the strategy for the computations
    // you have big numbers and they don't guarantee sign
    e = _mm_xor_pd(e, MachineEstimate::signmask);

    // separately compute positive and negative sum
    // and substract them at the end
    __m128d z4(z);
    z = _mm_xor_pd(z, MachineEstimate::signmask);
    s3 = _mm_mul_pd(z, MachineEstimate::coeff_log[15]);    // coeff[15] is -,+
    s2 = _mm_mul_pd(z, MachineEstimate::coeff_log[11]);
    z4 = _mm_mul_pd(z4, z);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_log[13]);                            // coeff[13] is -,+
    s1 = _mm_mul_pd(z, MachineEstimate::coeff_log[7]);
    __m128d z2 = _mm_xor_pd(z4, MachineEstimate::signmask);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[9]);
    s0 = _mm_mul_pd(z, MachineEstimate::coeff_log[3]);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[5]);
    s3 = _mm_mul_pd(s3, z2);
    z4 = _mm_mul_pd(z4, z2);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[1]);
    z4 = _mm_xor_pd(z4, MachineEstimate::signmask);
    s1 = _mm_mul_pd(s1, z2);
    s2 = _mm_add_pd(s2, s3);
    s3 = _mm_mul_pd(z, MachineEstimate::coeff_log[14]);    // +,-
    s0 = _mm_add_pd(s0, s1);
    s1 = _mm_mul_pd(z, MachineEstimate::coeff_log[6]);
    s3 = _mm_add_pd(s3, MachineEstimate::coeff_log[12]);    // +,-
    a = _mm_mul_pd(s2, z4);
    s2 = _mm_mul_pd(z, MachineEstimate::coeff_log[10]);
    a = _mm_add_pd(s0, a);
    s3 = _mm_mul_pd(s3, z2);
    s1 = _mm_add_pd(s1, MachineEstimate::coeff_log[4]);
    s0 = _mm_mul_pd(z, MachineEstimate::coeff_log[2]);
    a = _mm_mul_pd(a, v);
    s2 = _mm_add_pd(s2, MachineEstimate::coeff_log[8]);
    s1 = _mm_mul_pd(s1, z2);
    s2 = _mm_add_pd(s2, s3);
    s0 = _mm_add_pd(s0, MachineEstimate::coeff_log[0]);
    s2 = _mm_mul_pd(s2, z4);
    s0 = _mm_add_pd(s0, s1);
    s0 = _mm_add_pd(s0, s2);
    a = _mm_add_pd(a, s0);
    a = _mm_add_pd(a, e);
    return MulPositiveConstant(a, _mm_xor_pd(MachineEstimate::ln2, MachineEstimate::signmask));
}

 */

template <>
MachineEstimate ln2(unsigned int prec)
{
    return MachineEstimate(MachineEstimate::ln2);
}

template <>
MachineEstimate pi(unsigned int prec)
{
    return MachineEstimate(MachineEstimate::pi);
}



}    // namespace

