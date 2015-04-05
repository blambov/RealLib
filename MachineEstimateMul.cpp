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

#include <iomanip>
#include "MachineEstimate.h"

namespace RealLib {

double MachineEstimate::plusinf = 0.0;
double MachineEstimate::minusinf = 0.0;
double MachineEstimate::OnePEps = 1.0;
double MachineEstimate::OneMEps = 1.0;

static u32 ConstsLong[172] =
{
#define CONSTS_COS 0
        0x00000001, 0x3ff00000, // coeff[0] rounded up
        0xffffffff, 0x3fefffff, // coeff[0] rounded down
        0x418cafdb, 0x40018bc4, // coeff[1], exact
        0x9e3185f6, 0x3fe9a7b2, // coeff[2], exact
        0x5d05165d, 0x3fbe0270,
        0xa8e30653, 0x3f82ce22,
        0x70426553, 0x3f3d5450,
        0x749f656f, 0x3eef2edc,
        0x44c498c8, 0x3e979e4b,
#define CONSTS_RPI4 9
        // 1/(Pi*4), multiplication constant
        0x6dc9c883, 0x3fb45f30,
        0x6dc9c882, 0x3fb45f30,
#define CONSTS_PI    11
        // Pi, proper interval
        0x54442d19, 0x400921fb,
        0x54442d18, 0x400921fb,
#define CONSTS_SIN 13
        0x382d7366, 0x4000c152,
        0x382d7365, 0x4000c152,
        0x791049dc, 0x3ff87fb0,
        0x3ea3fdb3, 0x3fd57e24,
        0x23972846, 0x3fa1f529,
        0x62748c9e, 0x3f618133,
        0x4e962080, 0x3f165652,
        0xe58a04cb, 0x3ec4189c,
        0x3772c742, 0x3e6a705b,
#define CONSTS_EXPMASK 22
        0x00000000, 0x7ff00000,
#define CONSTS_EXPBIAS 23
        0x00000000, 0x3fe00000,
#define CONSTS_LOG 24
        0xed268e66, 0x3ca3df51,
        0xb100afab, 0xbc2dcabb,
        0x00000004, 0x3ff00000,
        0xffffeeea, 0x3fdfffff,
        0x55550aa6, 0x3fd55555,
        0x002a2505, 0x3fd00000,
        0x9a79b04d, 0x3fc99999,
        0x0771c502, 0x3fc55555,
        0xfc94cb71, 0x3fc24923,
        0x35ec7035, 0x3fc00022,
        0x9d9f4587, 0x3fbc722e,
        0xab47707f, 0x3fb98a39,
        0xb4029d62, 0x3fb73291,
        0x0deb07e2, 0x3fb7085a,
        0xa63cf31c, 0x3fb582e2,
#define CONSTS_SQRTSQRT2 39
        0xa31b716, 0x3ff306fe,
        0xa31b715, 0x3ff306fe,
#define CONSTS_LN2    41
        0xfefa39f0, 0x3fe62e42,
        0xfefa39ef, 0x3fe62e42,
#define CONSTS_LN2C    43
        0xfefa39f0, 0x3fe62e42,
        0xfefa39ef, 0x3fe62e42,
#define CONSTS_LOG2E 45
        0x652b82ff, 0x3ff71547,
        0x652b82fe, 0x3ff71547,
#define CONSTS_SIGN 47
        0x00000000, 0x80000000,
#define CONSTS_EXPLIMIT 48
        0x0, 0x408ff000,
#define CONSTS_EXP 49
        0x00000001, 0x40000000,
        0xffffffff, 0x3fffffff,
        0xfefa39f9, 0x3ff62e42,
        0xff82bdb1, 0x3fdebfbd,
        0xd706fa97, 0x3fbc6b08,
        0x6f5ef210, 0x3f93b2ab,
        0xf7c7e6fd, 0x3f65d87f,
        0x5d0bd9c1, 0x3f34308f,
        0x722cb340, 0x3efffd04,
        0x43ec690a, 0x3ec628a6,
        0xab63f5ed, 0x3e8b898b,
        0xdf1599b6, 0x3e4c140c,
        0xc4fc16f9, 0x3e15a8b6,
#define CONSTS_ATAN 62
        0x00000001, 0x3ff00000,
        0xffffffff, 0x3fefffff,
        0x555553d2, 0xbfd55555,
        0x9998037a, 0x3fc99999,
        0x91f33a63, 0xbfc24924,
        0x09057800, 0x3fbc71c7,
        0x1aa24579, 0xbfb745d0,
        0xf9b84bf5, 0x3fb3b12a,
        0x01a930e2, 0xbfb11089,
        0x556e5d85, 0x3fae177e,
        0xa80e2f1b, 0xbfaad32f,
        0xa58c31d6, 0x3fa7ee71,
        0x8b0ccaa5, 0xbfa4f50b,
        0x6c6308fe, 0x3fa17309,
        0x2b4c52ee, 0xbf9a77d1,
        0x7e19f3dd, 0x3f916913,
        0xfa32033c, 0xbf82da21,
        0xd33c5aff, 0x3f6fb050,
        0x6bed862f, 0xbf532726,
        0x510269d4, 0x3f2d637e,
        0x64cd132e, 0xbef5619e,
#define CONSTS_SQRT2    83
        // sqrt(2), mult. constant
        0x667f3bcd, 0x3ff6a09e,
        0x667f3bcc, 0x3ff6a09e,
#define CONSTS_PI2    85
        // Pi*2 rounded up
        0x54442d19, 0x402921fb
};

static const double *Consts = (double*)(ConstsLong);
static const double *CosConsts = Consts + CONSTS_COS;
static const double *SinConsts = Consts + CONSTS_SIN;
static const double *LogConsts = Consts + CONSTS_LOG;
static const double *ExpConsts = Consts + CONSTS_EXP;
static const double *AtanConsts = Consts + CONSTS_ATAN;
static const MachineEstimate me_rpi4(Consts[CONSTS_RPI4+1], Consts[CONSTS_RPI4]);
static const MachineEstimate me_pi_over_2(Consts[CONSTS_PI+1]*0.5, Consts[CONSTS_PI]*0.5);
static const MachineEstimate me_pi_over_4(Consts[CONSTS_PI+1]*0.25, Consts[CONSTS_PI]*0.25);
static const MachineEstimate me_pi(Consts[CONSTS_PI+1], Consts[CONSTS_PI]);
static const MachineEstimate me_ln2(Consts[CONSTS_LN2+1], Consts[CONSTS_LN2]);
static const double &me_upped_pi2 = Consts[CONSTS_PI2];
static const MachineEstimate me_sqrt_2(Consts[CONSTS_SQRT2+1], Consts[CONSTS_SQRT2]);
static const MachineEstimate me_sqrt_sqrt_2(Consts[CONSTS_SQRTSQRT2+1], Consts[CONSTS_SQRTSQRT2]);
static const MachineEstimate me_log2e(Consts[CONSTS_LOG2E+1], Consts[CONSTS_LOG2E]);

#ifdef _MSC_VER
#define nextafter _nextafter
#endif

// one positive and one negative multiple, positive result of the addition
// treating the elements separately
MachineEstimate AddProductPosNeg(const MachineEstimate add, const MachineEstimate &pos, double neg) 
{
    return MachineEstimate(
            MachineEstimate::RoundToZero(add.low + MachineEstimate::RoundFromZero(pos.low * neg)),
            MachineEstimate::RoundFromZero(add.high + MachineEstimate::RoundToZero(pos.high * neg)));
}


void MachineEstimate::BeginComputation()
{
    double zero = 0.0;
    plusinf = 1.0/zero;
    minusinf = -plusinf;
    OnePEps = nextafter(1.0, plusinf);
    OneMEps = nextafter(1.0, minusinf);
}

void MachineEstimate::FinishComputation()
{
}

std::ostream& operator <<(std::ostream &os, const MachineEstimate &e)
{    
    return os << e.weak_AsDouble();
    /*
    if (e.low < 0 && e.high > 0) return os << 0.0;
    else {
        double d(e.weak_AsDouble()); 
        if (d < 1) os << d;
        else {
            int op = os.precision();
            double v = d;
            int np = op;
            while (v >= 1 && np > 0) {
                --np;
                v /= 10;
            }
            os << std::setprecision(np) << d << std::setprecision(op);
        }
    }*/
}

#ifndef TRUST_STDLIB

MachineEstimate cos(const MachineEstimate &x)
{
    MachineEstimate z(x * me_rpi4);
    double d = z.Sum();    // possibly inexact, reflected in doubling r
    MachineEstimate r(z.Diff());    // exact
    z = MachineEstimate(fabs(d - floor(d+0.5)));    // exact operations
    // now z has (-pi,pi] mapped into (-0.5, 0.5]
    // the minimax we use does the cosine of (-pi/3, pi/3]
    // mapped into (-0.5, 0.5], squared
    z = z.MulPositive(z);
    MachineEstimate z2(z.MulPositive(z));
    MachineEstimate s0(CosConsts[1], CosConsts[0]);
    MachineEstimate s01(s0.SubProductPositive(z, CosConsts[2]));
    MachineEstimate s2(CosConsts[3]);
    MachineEstimate s23(s2.SubProductPositive(z, CosConsts[4]));
    MachineEstimate s4(CosConsts[5]);
    MachineEstimate s45(s4.SubProductPositive(z, CosConsts[6]));
    MachineEstimate s6(CosConsts[7]);
    MachineEstimate s67(s6.SubProductPositive(z, CosConsts[8]));
    MachineEstimate z4(z2.MulPositive(z2));
    s0 = s01.AddProductPositive(z2, s23);
    s4 = s45.AddProductPositive(z2, s67);
    s0 = s0.AddProductPositive(z4, s4);

    // now we use the formula cos(3x) = cos(x)*(4*cos(x)*cos(x) - 3)
    s4 = s0.MulPositive(s0);
    s2 = s0.MulPositive(3.0);
    s6 = s0.MulPositive(4.0);
    s0 = s4.MulPositive(s6) - s2;
    return s0.AddError((r.MulDouble(me_upped_pi2)));

}

MachineEstimate sin(const MachineEstimate &x)
{
    MachineEstimate z(x * me_rpi4);
    double d = z.Sum();    // possibly inexact, reflected in doubling error
    d = d - floor(d+0.5);    // exact operations
    MachineEstimate r(z.Diff());
    z = fabs(d);
    // now r has (-pi,pi] mapped into (-0.5, 0.5]
    // the minimax we use does the sine of (-pi/3, pi/3]
    // mapped into (-0.5, 0.5], squared
    z = z.MulPositive(z);
    MachineEstimate z2(z.MulPositive(z));
    MachineEstimate s0(SinConsts[1], SinConsts[0]);
    MachineEstimate s01(s0.SubProductPositive(z, SinConsts[2]));
    MachineEstimate s2(SinConsts[3]);
    MachineEstimate s23(s2.SubProductPositive(z, SinConsts[4]));
    MachineEstimate s4(SinConsts[5]);
    MachineEstimate s45(s4.SubProductPositive(z, SinConsts[6]));
    MachineEstimate s6(SinConsts[7]);
    MachineEstimate s67(s6.SubProductPositive(z, SinConsts[8]));
    MachineEstimate z4(z2.MulPositive(z2));
    s0 = s01.AddProductPositive(z2, s23);
    s4 = s45.AddProductPositive(z2, s67);
    s0 = s0.AddProductPositive(z4, s4);
    s0 = s0.MulDouble(d);

    // now we use the formula sin(3x) = sin(x)*(3 - 4*sin(x)*sin(x))
    s4 = s0.MulPositive(s0);
    s2 = s0.MulPositive(3.0);
    s6 = s0.MulPositive(4.0);
    s0 = s2.SubProductPositive(s4, s6);
    return s0.AddError((r.MulDouble(me_upped_pi2)));

}

static inline
void logreduction(const double additive, const double v, double &a, double &e)
{
    if (v < 1.0) {
        a = v;
        e -= additive;    // this is exact, because e only holds up to 12 bits until now
    }
}

static inline 
void logreduction(const MachineEstimate &factor,
        const double additive,
        MachineEstimate &arg,
        MachineEstimate &exp)
{
    MachineEstimate v(arg.MulPositive(factor));
    logreduction(additive, v.low, arg.low, exp.low);
    logreduction(additive, v.high, arg.high, exp.high);
}

MachineEstimate log(const MachineEstimate &arg)
{
    int el, eh;
    if (!(arg.IsPositive()))
        if (arg.IsNegative()) throw DomainException("log");
        else throw PrecisionException("log");

    MachineEstimate a(frexp(arg.low, &el), frexp(arg.high, &eh));    // these are both exact operations
    MachineEstimate e(el, eh);

    logreduction(me_sqrt_2, 0.5, a, e);
    logreduction(me_sqrt_sqrt_2, 0.25, a, e);

    // we're switching the two positions here so that we could use MulPositive on what are
    // actually negative values
    a = MachineEstimate(1.0 - a.high, 1.0 - a.low);    // this is exact as a is within (0.5, 1.0]

    MachineEstimate a2(a.MulPositive(a));

    MachineEstimate s2(LogConsts[13]);
    s2 = s2.AddProductPositive(a, LogConsts[14]);
    MachineEstimate s1(LogConsts[11]);
    s1 = s1.AddProductPositive(a, LogConsts[12]);
    MachineEstimate s0(LogConsts[9]);
    s0 = s0.AddProductPositive(a, LogConsts[10]);
    MachineEstimate a4(a2.MulPositive(a2));
    s0 = s0.AddProductPositive(s1, a2);
    MachineEstimate s3(LogConsts[7]);
    s0 = s0.AddProductPositive(s2, a4);
    s2 = MachineEstimate(LogConsts[5]);
    s3 = s3.AddProductPositive(a, LogConsts[8]);
    s1 = MachineEstimate(LogConsts[3]);
    s2 = s2.AddProductPositive(a, LogConsts[6]);
    s1 = s1.AddProductPositive(a, LogConsts[4]);
    s2 = s2.AddProductPositive(s3, a2);
    s3 = MachineEstimate(LogConsts[1], LogConsts[0]);
    s3 = s3.AddProductPositive(a, LogConsts[2]);
    MachineEstimate a8(a4.MulPositive(a4));
    s3 = s3.AddProductPositive(s1, a2);
    s3 = s3.AddProductPositive(s2, a4);
    s3 = s3.AddProductPositive(s0, a8);

    e = e.MulPositiveRHS(me_ln2);

    // this gets us just what we want 
    return e - s3;
}

MachineEstimate exp(const MachineEstimate &arg)
{
    MachineEstimate a(arg * me_log2e);

    if (a.high >= 1020.0 || a.low <= -1020.0) 
        throw PrecisionException("exp");

    MachineEstimate z(floor(a.low), floor(a.high));
    a = MachineEstimate(a.low - z.low, a.high - z.high);        // exact result

    MachineEstimate s1(ExpConsts[11]);
    MachineEstimate s0(ExpConsts[9]);
    s1 = s1.AddProductPositive(a, ExpConsts[12]);
    MachineEstimate a2(a.MulPositive(a));
    s0 = s0.AddProductPositive(a, ExpConsts[10]);
    MachineEstimate s3(ExpConsts[7]);
    s0 = s0.AddProductPositive(s1, a2);
    MachineEstimate s2(ExpConsts[5]);
    s3 = s3.AddProductPositive(a, ExpConsts[8]);
    MachineEstimate a4(a2.MulPositive(a2));
    s1 = ExpConsts[3];
    s2 = s2.AddProductPositive(a, ExpConsts[6]);
    s1 = s1.AddProductPositive(a, ExpConsts[4]);
    s2 = s2.AddProductPositive(s3, a2);
    s3 = MachineEstimate(ExpConsts[1], ExpConsts[0]);
    MachineEstimate a8(a4.MulPositive(a4));
    s3 = s3.AddProductPositive(a, ExpConsts[2]);
    s2 = s2.MulPositive(a4);
    s3 = s3.AddProductPositive(s1, a2);
    s2 = s2.AddProductPositive(s0, a8);
    s3 = s3 + s2;

    z = MachineEstimate(ldexp(s3.low, int(z.low)-1), ldexp(s3.high, int(z.high)-1));    // ldexp is exact
    return z;
}

// 0.0 <= arg <= 1.0
MachineEstimate atanprimary(const MachineEstimate &arg)
{
    MachineEstimate a(arg.MulPositive(arg));
    MachineEstimate s3(AtanConsts[15]);
    MachineEstimate a2(a.MulPositive(a));
    MachineEstimate s2(AtanConsts[13]);
    s3 = AddProductPosNeg(s3, a, AtanConsts[16]);
    MachineEstimate s1(AtanConsts[11]);
    s2 = AddProductPosNeg(s2, a, AtanConsts[14]);
    MachineEstimate s0(AtanConsts[9]);
    s1 = AddProductPosNeg(s1, a, AtanConsts[12]);
    s0 = AddProductPosNeg(s0, a, AtanConsts[10]);
    s2 = s2.AddProductPositive(a2, s3);
    MachineEstimate a4(a2.MulPositive(a2));
    s3 = AtanConsts[7];
    s0 = s0.AddProductPositive(a2, s1);
    s1 = AtanConsts[3];
    s3 = AddProductPosNeg(s3, a, AtanConsts[8]);
    MachineEstimate s(s0.AddProductPositive(a4, s2));
    MachineEstimate a8(a4.MulPositive(a4));
    s2 = AtanConsts[5];
    s1 = AddProductPosNeg(s1, a, AtanConsts[4]);
    s0 = MachineEstimate(AtanConsts[1], AtanConsts[0]);
    s2 = AddProductPosNeg(s2, a, AtanConsts[6]);
    s1 = s1.AddProductPositive(a4, s3);
    s0 = AddProductPosNeg(s0, a, AtanConsts[2]);
    s3 = AtanConsts[19];
    MachineEstimate a16(a8.MulPositive(a8));
    s0 = s0.AddProductPositive(a4, s2);
    s2 = AtanConsts[17];
    s3 = AddProductPosNeg(s3, a, AtanConsts[20]);
    s0 = s0.AddProductPositive(a2, s1);
    s2 = AddProductPosNeg(s2, a, AtanConsts[18]);
    s0 = s0.AddProductPositive(a8, s);
    s2 = s2.AddProductPositive(a2, s3);
    s0 = s0.AddProductPositive(a16, s2);
    return s0.MulPositive(arg);
}

MachineEstimate atan(const MachineEstimate &x)
{
    if (x.IsPositive()) {
        if (x > 1.0) return me_pi_over_2 - atanprimary(recip(x));
        else if (x < 1.0) return atanprimary(x);
        else // atan is Lipschitz 1. Just return pi/4 + twice the error from argument
            return MachineEstimate(me_pi_over_4).AddError(x.Diff());
    } else if (x.IsNegative()) {
        MachineEstimate y(-x);
        if (y > 1.0) return atanprimary(recip(y)) - me_pi_over_2;
        else if (y < 1.0) return -atanprimary(y);
        else // atan is Lipschitz 1. Just return pi/4 + twice the error from argument
            return -MachineEstimate(me_pi_over_4).AddError(x.Diff());
    } else {
        return MachineEstimate().SetError(x.Diff());
    }
}


MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x)
{
    if (!x.IsValueValid() || !y.IsValueValid()) 
        throw PrecisionException("atan2 argument");

    if (x.IsPositive()) {
        if (y.IsPositive()) {
            MachineEstimate d(x-y);
            if (d.IsPositive()) return atanprimary(y/x);
            else if (d.IsNegative()) return me_pi_over_2 - atanprimary(x/y);
            else {
                MachineEstimate r(y/x);
                return MachineEstimate(me_pi_over_4).AddError(r.Diff());
            }
        } else if (y.IsNegative()) {
            MachineEstimate d(x+y);
            if (d.IsPositive()) return -atanprimary((-y)/x);
            else if (d.IsNegative()) return atanprimary(x/(-y)) - me_pi_over_2;
            else {
                MachineEstimate r(y/x);        // sign doesn't matter here
                return -MachineEstimate(me_pi_over_4).AddError(r.Diff());
            }
        } else {
            MachineEstimate r(y/x);
            return MachineEstimate().SetError(r.Diff());
        }
    } else 
        // no check needed because x close to zero would either not cause a problem
        // if y is big enough, or will cause a PrecisionException at the moment of
        // the division y/x or directly
        //if (x.IsNegative()) 
    {
        if (y.IsPositive()) {
            MachineEstimate d(-x-y);
            if (d.IsPositive()) return me_pi - atanprimary(y/(-x));
            else if (d.IsNegative()) return me_pi_over_2 + atanprimary((-x)/y);
            else {
                MachineEstimate r(y/x);
                return (me_pi + me_pi_over_4).AddError(r.Diff());
            }
        } else if (y.IsNegative()) {
            MachineEstimate d(y-x);
            if (d.IsPositive()) return atanprimary(y/x) - me_pi;
            else if (d.IsNegative()) return -(atanprimary(x/y) + me_pi_over_2);
            else {
                MachineEstimate r(y/x);        // sign doesn't matter here
                return -(me_pi + me_pi_over_4).AddError(r.Diff());
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
#endif

template <>
MachineEstimate pi(unsigned int prec)
{
    return me_pi;
}

template <>
MachineEstimate ln2(unsigned int prec)
{
    return me_ln2;
}

}

