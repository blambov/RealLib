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
#include <stdio.h>

#include "MachineEstimate.h"
#include "RealFuncs.h"

namespace RealLib {

Estimate *g_pCachedEstimate_pi = NULL;
Estimate *g_pCachedEstimate_ln2 = NULL;

#define DELETE(x) if (x) { delete x; x = NULL; }
void DestroyConstantEstimates()
{
    DELETE(g_pCachedEstimate_pi);
    DELETE(g_pCachedEstimate_ln2);
}

static inline int recip(int x)
{ return 0*x; }

// pow
// use pow(a, 2*p) == pow(a*a, p)
template <class T>
T pow(T arg, i32 pwr)
{
    if (pwr < 0) {
        pwr = -pwr;
        arg = recip(arg);
    }
    T acc(1);

    while (pwr) {
        if (pwr & 1) acc *= arg;
        pwr >>= 1;
        arg = sq(arg);
    }

    return acc;
}

Estimate PerformNewton(Estimate arg, NewtonIterator iter, 
                       Estimate est, i32 prec)
{
    i32 targetprec = arg.GetPrecision() * 32 - 64;
    while (prec < targetprec) {
        arg.SetPrecision( (prec)/32 +2);
        est.SetPrecision( (prec)/32 +2);
        est.SetError(0.0);

        iter(arg, est, prec);
    }

    // this iteration should produce the actual error bound (caused by inexact
    // functions used in the iterations and the error in the input)
    arg.SetPrecision(targetprec / 32 + 2);
    est.SetPrecision(targetprec / 32 + 2);
    est.SetError(0.0);
    //   Estimate old = est;

    iter(arg, est, prec);

    return est;//.AddError(old - est);
}

void NewtonRSQRT(const Estimate &arg, Estimate &est, i32 &prec)
{
    est = est * (Estimate(1.5) - est * est * arg);
    prec = prec * 2 - 3;
}

Estimate rsqrt(const Estimate &a)
{
    Estimate arg(a);
    if (arg.IsNegative()) throw DomainException("rsqrt");

    if (!arg.IsPositive())
        throw PrecisionException("rsqrt");

    Estimate ei;
    {
        i32 exp = arg.weak_normalize();
        if (exp & 1) { exp--; }
        double d = 1.0 / ::sqrt((arg>>exp).weak_AsDouble());
        ei = Estimate(d) >> (exp / 2);
    }

    Estimate res(PerformNewton(arg/2, NewtonRSQRT, ei, 45 * 2 - 3));

    return res;

    // this should not be needed!
    Estimate err(arg.GetError() / res);

    return res.AddError(err / ((res - err) * res));
}

// separate function because of the handling of zeroes
Estimate sqrt(const Estimate &a)
{
    // ignore the negative part
    Estimate arg(a.TruncateNegative("sqrt"));

    Estimate ei;
    {
        i32 exp = arg.weak_normalize();
        if (exp & 1) { exp--; }
        double d = 1.0 / ::sqrt((arg>>exp).weak_AsDouble());
        ei = Estimate(d) >> (exp / 2);
    }

    Estimate res(arg * PerformNewton(arg/2, NewtonRSQRT, ei, 45 * 2 - 3));

    return res;

    // this should not be needed!
    if (arg.IsPositive())
        return res.AddError(arg.GetError() / res);
    else
        return res.SetError(res);
}

void PowerSeriesDirect(const Estimate &arg, SeriesIterator iter,
                       Estimate &sum, Estimate &workspace, i32 indexstart, i32 indexend)
{
    for (int i=indexstart; i<indexend; ++i)
        sum += iter(arg, workspace, i);

    sum.AddError(workspace);
}

Estimate abs(const Estimate &arg)
{
    if (arg.IsPositive()) return arg;
    if (arg.IsNegative()) return -arg;
    Estimate a;
    if (arg.weak_IsPositive()) a = arg;
    else a = -arg;
    a = a+a.GetError();
    return (a/2).SetError(a/2);
}

// exponent by McLauren expansion
Estimate IteratorEXP(const Estimate &arg, Estimate &workspace, i32 index)
{
    return workspace = workspace * arg / index;
}

// for the primary interval [-1;1]
Estimate exp_primary(const Estimate &arg)
{
    i32 i = i32(ceil(::pow(double(arg.GetPrecision()*32), 1./2.)));
    Estimate x(arg >> i);

    Estimate sum(x + 1.0);
    Estimate workspace(x);

    i32 indexend = i;

    PowerSeriesDirect(x, IteratorEXP, sum, workspace, 2, indexend);

    for (i32 k=0;k<i;++k)
        sum *= sum;
    return sum;
}

// logarithm by exponent and newton
void NewtonLN(const Estimate &arg, Estimate &est, i32 &prec)
{
    Estimate ex(exp_primary(est));
    est = est + (arg - ex) / ex;
    prec = prec * 2 - 2;
}

// primary interval [1/e; e]
Estimate log_primary(const Estimate &arg)
{
    double d(::log(arg.weak_AsDouble()));
    return PerformNewton(arg, NewtonLN, Estimate(d), 50 * 2 - 2);
}

template <>
Estimate ln2<Estimate>(unsigned int prec)
{
    if (g_pCachedEstimate_ln2 && g_pCachedEstimate_ln2->GetPrecision() >= prec)
        return *g_pCachedEstimate_ln2;
    prec = g_WorkingPrecision;

    Estimate v(2.0);

    v = log_primary(v);

    if (g_pCachedEstimate_ln2) *g_pCachedEstimate_ln2 = v;
    else g_pCachedEstimate_ln2 = new Estimate(v);

    return v;
}

Estimate exp(const Estimate &arg)
{
    Estimate l(ln2<Estimate>(arg.GetPrecision()));
    Estimate x(arg / l);
    Estimate de(arg - x*l);

    Estimate e(x.weak_round());
    // domainerror?

    x = (x-e) * l;

    Estimate y(exp_primary(x));
    return y << i32(e.weak_AsDouble());
}

Estimate log(const Estimate &arg)
{
    if (arg.IsNegative()) throw DomainException("log");
    if (!arg.IsPositive()) throw PrecisionException("log");

    Estimate l(ln2<Estimate>(arg.GetPrecision()));
    Estimate x(arg);
    i32 e = x.weak_normalize();

    return log_primary(x >> e) + l * e;
}

// rpi: 1/pi
// using Borwein iterations that quadruple the number of correct digits
// at each step
Estimate rpi(unsigned int prec)
{
    Estimate y(sqrt(Estimate(2)) - Estimate(1));
    Estimate alpha(Estimate(2) - y * Estimate(4));
    Estimate z(sqrt(sqrt(Estimate(1) - sq(sq(y)))));
    Estimate oa(0.0);
    Estimate one(1.0);

    y = (one - z) / (one + z);
    oa = alpha;
    alpha = sq(sq(y + one)) * alpha - 8 * y * (sq(y) + y + one);

    // OBS! This would not work if precision is greater that I32_MAX/128
    for (u32 i = 32; i < prec*128; i = i*4) {
        z = sqrt(sqrt(one - sq(sq(y))));
        y = (one - z) / (one + z);
        oa = alpha;
        alpha = sq(sq(y + one)) * alpha - y * i * (sq(y) + y + one);
    }

    return alpha.AddError(alpha - oa);
}

template <>
Estimate pi(unsigned int prec)
{ 
    if (g_pCachedEstimate_pi && g_pCachedEstimate_pi->GetPrecision() >= prec)
        return *g_pCachedEstimate_pi;
    prec = g_WorkingPrecision;

    Estimate v = recip(rpi(prec));

    if (g_pCachedEstimate_pi) *g_pCachedEstimate_pi = v;
    else g_pCachedEstimate_pi = new Estimate(v);

    return v;
}

// sine by McLauren expansion
Estimate IteratorSIN(const Estimate &arg, Estimate &workspace, i32 index)
{
    if (index < 30000)
        return workspace = workspace * arg / (2*index * (2*index+1));
    else
        return workspace = workspace * arg / (2*index) / (2*index+1);
}

// for the primary interval [-pi/2;pi/2]
Estimate sin_primary(const Estimate &arg)
{
    // here the strength reduction is done with base 3.
    // this means about 20.1897521 reductions per 32-bit word
    // additionally, 3^20 fits in u32, but not in i32, so use 3^19
    // for the division factor (1162261467)
    i32 i = i32(ceil(::pow(double(arg.GetPrecision())*20.2, 1./2.)))/2;
    i32 j = i % 19;
    Estimate x(arg / pow(3, j));
    for (;j<i;j+=19)
        x /= 1162261467;

    Estimate sum(x);
    Estimate workspace(x);

    i32 indexend = i*2;

    PowerSeriesDirect(-x*x, IteratorSIN, sum, workspace, 1, indexend);

    Estimate three(3);
    for (i32 k=0;k<i;++k)
        sum = sum * (three - 4*sum*sum);
    return sum;
}

Estimate sin(const Estimate &arg)
{
    Estimate pi2(pi<Estimate>(arg.GetPrecision())*2);
    Estimate x(arg / pi2);
    x -= x.weak_round();
    if (!(x<0.6125) || !(x>-0.6125))
        throw PrecisionException("sin"); // alternatively, we could just return [-1, 1]
    if (x.weak_gt(0.25)) x = 0.5 - x;
    else if (x.weak_lt(-0.25)) x = -0.5 - x;

    return sin_primary(x * pi2);
}

// logarithm by exponent and newton
void NewtonASIN(const Estimate &arg, Estimate &est, i32 &prec)
{
    Estimate ex(sin_primary(est));
    est = est + (arg - ex) * rsqrt(1 - sq(ex));
    prec = prec * 2 - 2;
}

Estimate asin_primary(const Estimate &arg)
{
    double d(::asin(arg.weak_AsDouble()));
    return PerformNewton(arg, NewtonASIN, Estimate(d), 50 * 2 - 2);
}

Estimate asin(const Estimate &arg)
{
    Estimate x(arg.TruncateTo(-1, 1, "asin"));
    /*bool sign = x.weak_IsNegative();

   if (sign) x = -x;
   Estimate x1(x - 1);
   if (x1.IsPositive()) throw DomainException("asin");*/
    /*
   if (!x1.IsNegative())
   {
      x = pi(arg.GetPrecision())/2 + x1.AddError(x1);
      return sign ? -x : x;
   }*/

    // we still have a problem with the rsqrt if arg is close to one.
    // in this case, use the cos-sin identities
    if (x.weak_gt(0.708))
        return pi<Estimate>(arg.GetPrecision())/2 - asin_primary(cosfromsin(x));
    if (x.weak_lt(-0.708))
        return pi<Estimate>(arg.GetPrecision())/-2 + asin_primary(cosfromsin(x));
    return asin_primary(x);
}

template <class Estimate>
Estimate atan2(const Estimate &y, const Estimate &x)
{
    // this is atan with result over the full range depending on the signs of both arguments
    // normal atan takes care of the sign of y (the sine value). (0, 0) is undefined, so
    // we can just return 0 with error pi -- so that some complex operations that use atan2
    // would work.
    //
    Estimate tpi(pi<Estimate>(x.GetPrecision()));

    if (x.IsPositive()) {
        return atan(y/x);
    } else if (x.IsNegative()) {
        if (y.IsPositive()) return tpi - atan(y/-x);
        else return -tpi - atan(y/-x);
    } else { // x cannot be distinguished from zero, but this does not stop us to give
        // estimates for the angle based on y and x's possible values
        if (y.IsPositive()) return (tpi/2).AddError(y/(2*x.GetError()));
        else if (y.IsNegative()) return (tpi/-2).AddError(y/(2*x.GetError()));
        else return Estimate().SetError(tpi);
    }
}

template
Estimate atan2(const Estimate &y, const Estimate &x);

template
MachineEstimate atan2(const MachineEstimate &y, const MachineEstimate &x);

Estimate tan(const Estimate &arg)
{ 
    Estimate pi2(pi<Estimate>(arg.GetPrecision())*2);
    Estimate x(arg / pi2);
    x -= x.weak_round();
    bool negc = false;
    if (x > 0.25) {
        x = 0.5 - x;
        negc = true;
    } else if (!(x < 0.25)) {
        throw PrecisionException("tan");
    } else if (x < -0.25) {
        x = -0.5 - x;
        negc = true;
    } else {
        if (!(x > -0.25)) throw PrecisionException("tan");
    }

    Estimate s(sin_primary(x * pi2));
    if (negc) return -s/cosfromsin(s);
    else return s/cosfromsin(s);

}

/*
// rsqrt
// use double of mantissa to get initial estimate.
// refine using Newton-Raphson iterations

#ifndef max
	static inline i32 max(i32 a, i32 b) { return a>b ? a : b; }
#endif

LongFloat rsqrt(const LongFloat &arg)
{
	switch (arg.Kind()) {
	case LongFloat::Zero:
		return LongFloat(LongFloat::Infinity, arg.IsNegative());
	case LongFloat::Infinity:
		return LongFloat(LongFloat::Zero, arg.IsNegative());
	case LongFloat::Nan:
		return arg;
	default:
	case LongFloat::Normal:
		double init(1.0 / ::sqrt(arg.MantissaAsDouble()));
		LongFloat y(arg.MantissaAsLongFloat());

		LongFloat three(3, 0);
		LongFloat half(0.5);
		LongFloat r(init);

		// log(err) = log(err) * 2 - 2 in theory. Make it 3 to accomodate rounding
		for (int i=50; i/32 < g_WorkingPrecision; i = i*2 - 3)
			r = half * r * (three - r * r * y);

		i32 exp = arg.Exponent();
		if (exp % 2) {
			++exp;
			r *= (1<<16);
		}

		return r.AddToExponent(-(exp/2));
	}


}

// sqrt
// to avoid NaN on sqrt(Zero/Infinity) handle them separately
LongFloat sqrt(const LongFloat &arg)
{
	switch (arg.Kind()) {
	case LongFloat::Zero:
	case LongFloat::Infinity:
	case LongFloat::Nan:
		return arg;
	default:
	case LongFloat::Normal:
		return arg * rsqrt(arg);
	}
}


// Estimate math
Estimate NonNegativeInterval(const Estimate& arg)
{
	if (!arg.Error().gePow2(0)) return arg;

	LongFloat a = (arg.Value() * (arg.Error().AsLongFloat() + 1)) >> 1;
	return Estimate(a, ErrorEstimate(a));
}

Estimate rsqrt(const Estimate &arg, u32 prec)
{
	prec;	// unused argument
	Estimate a(NonNegativeInterval(arg));
	if (a.Value().IsNegative()) return Estimate(LongFloat::Nan, 0);

	LongFloat r(rsqrt(a.Value()));

	if (a.Error().gePow2(-1)) 
		return Estimate(r, ErrorEstimate(0, ErrorEstimate::plusinf));

	return Estimate(r, a.Error() + a.Error() * (a.Error() << 1) + RoundingError());
}

// call only on arg > 0
Estimate sqrt(const Estimate &arg)
{
	return Estimate(sqrt(arg.Value()), arg.Error() + RoundingError()<<1);
}

Estimate sqrt(const Estimate &arg, u32 prec)
{
	Estimate a(NonNegativeInterval(arg));
	if (a.Value().IsNegative()) return Estimate(LongFloat::Nan, 0);

	return Estimate(sqrt(a.Value()), a.Error() + RoundingError()<<1);
}

// rpi: 1/pi
// using Borwein iterations that quadruple the number of correct digits
// at each step
Estimate rpi(u32 prec)
{
	Estimate y(sqrt(Estimate(2)) - Estimate(1));
	Estimate alpha(Estimate(2) - y * Estimate(4));
	Estimate z(sqrt(sqrt(Estimate(1) - sq(sq(y)))));
	Estimate pow2(8.0);
	Estimate oa(0.0);
	Estimate one(1.0);

	for (double i = 2; i < prec*32; i = i*4) {
		y = (one - z) / (one + z);
		if (y.Value().Kind() == LongFloat::Zero) return alpha;
		oa = alpha;
		alpha = sq(sq(y + one)) * alpha - pow2 * y * (sq(y) + y + one);
		pow2 = pow2 << 2;
		z = sqrt(sqrt(one - sq(sq(y))));
	}

	return Estimate(alpha.Value(), alpha.Error() + 
		ErrorEstimate((alpha - oa).Value()) * ErrorEstimate(alpha.Value()).recip());
}

// pi: 1/rpi
Estimate pi(u32 prec)
{
	return Constants.rpi.GetEstimate(prec).recip();
}

// perform AGM iteration and return 2 * AGM(a, b)
Estimate agm2(Estimate a, Estimate b, u32 prec)
{
	Estimate z((a + b) << -1);

	for (double i = 4; i<prec*32; i*=2) {	// two times per iteration
		b = sqrt(a * b);
		a = (z + b) << -1;
		b = sqrt(z * b);
		z = (a + b) << -1;

	}
	b = b - z;
	z = z << 1;
	if (b.Value().Kind() == LongFloat::Zero) return z;
	else return Estimate(z.Value(), z.Error() + 
		ErrorEstimate(b.Value()) * ErrorEstimate(z.Value()).recip());
}

// ln2 by AGM
Estimate ln2(u32 prec)
{
	i32 pow = prec * 16 + 10;
	Estimate z(Estimate(1) << (2 - pow));

	return (Constants.rpi.GetEstimate(prec) * (agm2(Estimate(1.0), z, prec) * pow)).recip();
}

// ln by AGM
Estimate ln(const Estimate &arg, u32 prec)
{
	// range control
	Estimate a(NonNegativeInterval(arg));
	if (a.Value().IsNegative()) return Estimate(LongFloat::Nan, 0);


	i32 pow = prec * 16 + 20;
	Estimate z((Estimate(4) / a.Value().MantissaAsLongFloat()) >> pow);
	i32 exp = a.Value().Exponent() * 32 - pow;

	z = (Constants.rpi.GetEstimate(prec) * (agm2(Estimate(1.0), z, prec))).recip();
	z = z + Constants.ln2.GetEstimate(prec) * exp;

	return Estimate(z.Value(), z.Error() + (a.Error() + a.Error() * (a.Error() << 1)) * ErrorEstimate(z.Value()).recip());
}

// exp
// use double math for initial estimate
// refine by Newton-Raphson and ln
Estimate exp(const Estimate &arg, u32 prec)
{
	// range control
	if (arg.Value().AsDouble() > ldexp(log(2), 32) ||
		arg.Value().AsDouble() < ldexp(-log(2), 32)) return Estimate(LongFloat::Nan, 0);

	Estimate a(abs(arg, prec));
	LongFloat ln2 = Constants.ln2.GetEstimate(prec).Value();
	LongFloat pwr = (a.Value() / ln2).RoundTowardZero();
	LongFloat flt = a.Value() - pwr * ln2;

	double dz = flt.AsDouble();
	LongFloat e(::exp(dz));
	LongFloat z(flt + LongFloat(1.0));

	for (u32 i=50; i<prec*32; i=i*2-3) {
		LongFloat l(ln(Estimate(e), i/16 + i/64 + 4).Value());
		e = e * (z - l);
	}
	Estimate r(Estimate(e) * (Estimate(z) - ln(Estimate(e), prec)));

	r = Estimate(r.Value(), r.Error() + ErrorEstimate(r.Value()-e) * ErrorEstimate(r.Value()).recip());

	r = r << i32(pwr.AsDouble());

	r = Estimate(r.Value(), r.Error() + ErrorEstimate(::exp((ErrorEstimate(a.Value()) * a.Error()).AsDouble()) - 1));

	return arg.Value().IsNegative() ? r.recip() : r;
}

// abs
Estimate abs(const Estimate &arg, u32 prec)
{
	return arg.Value().IsNegative() ? -arg : arg;
}

// complex AGM helpers

typedef RealLib::complex<Estimate> CmplxEst;

static Estimate abs(const CmplxEst& x)
{
	return sqrt(sq(x.re) + sq(x.im));
}

static
CmplxEst sqrt (const CmplxEst& x)
{
  Estimate r(abs(x));
  Estimate nr(LongFloat::Zero), ni(LongFloat::Zero);

  if (r.Value().Kind() != LongFloat::Normal)
    nr = ni = r;
  else if (!x.re.Value().IsNegative())
    {
      nr = sqrt ((r + x.re) << -1);
      ni = (x.im / nr) << -1;
    }
  else
    {
      ni = sqrt ((r - x.re) << -1);
      if (x.im.Value().IsNegative())
		ni = - ni;
      nr = (x.im / ni) << -1;
    }
  return CmplxEst (nr, ni); 
}

// perform comlpex AGM iteration and return 2 * AGM(a, b)
CmplxEst agm2complex(CmplxEst a, CmplxEst b, u32 prec)
{
	CmplxEst z((a + b) << -1);

	for (double i = 4; i<prec*32; i*=2) {
		b = sqrt(a * b);
		a = (z + b) << -1;
		b = sqrt(z * b);
		z = (a + b) << -1;

	}
	b = b - z;
	z = z << 1;

	ErrorEstimate err(ErrorEstimate(b.re.Value()) * ErrorEstimate(z.re.Value()).recip());
	ErrorEstimate eri(ErrorEstimate(b.im.Value()) * ErrorEstimate(z.im.Value()).recip());
	return CmplxEst(Estimate(z.re.Value(), z.re.Error() + err),
					Estimate(z.im.Value(), z.im.Error() + eri));
}

CmplxEst AdjustForAbsoluteError(const CmplxEst &arg, ErrorEstimate err)
{
	return CmplxEst(Estimate(arg.real().Value(), arg.real().Error() + err * ErrorEstimate(arg.real().Value()).recip()),
		            Estimate(arg.imag().Value(), arg.imag().Error() + err * ErrorEstimate(arg.imag().Value()).recip()));
}

CmplxEst lncomplex(const CmplxEst &arg, u32 prec)
{
	// no need for range control here
	ErrorEstimate err = arg.real().Error() * arg.real().Error() +
						arg.imag().Error() * arg.imag().Error();

	i32 exp = max(arg.real().Value().Exponent(),
				  arg.imag().Value().Exponent());

	i32 pow = prec * 16 + 20;
	CmplxEst z(arg.real().Value().SignedMantissaAsLongFloat().AddToExponent(arg.real().Value().Exponent() - exp),
			   arg.imag().Value().SignedMantissaAsLongFloat().AddToExponent(arg.imag().Value().Exponent() - exp));
	z = (CmplxEst(Estimate(4)) / z) << -pow;
	exp = exp * 32 - pow;

	z = CmplxEst(Constants.pi.GetEstimate(prec)) / (agm2complex(Estimate(1.0), z, prec));
	z = z + CmplxEst(Constants.ln2.GetEstimate(prec) * exp);

	return AdjustForAbsoluteError(z, err + err * (err << 1));
//	err = v.Error() + (err + err * (err << 1)) * ErrorEstimate(v.Value()).recip();
//	return CmplxEst(Estimate(z.real().Value(), err), Estimate(z.imag().Value(), err));
}

CmplxEst expcomplex(const CmplxEst& arg, u32 prec)
{
	// range control
	if (arg.real().Value().AsDouble() > ldexp(log(2), 32) ||
		arg.real().Value().AsDouble() < ldexp(-log(2), 32)) return Estimate(LongFloat::Nan, 0);
	bool sinneg = arg.real().Value().IsNegative();

	Estimate a(abs(arg.real(), prec));

	LongFloat ln2 = Constants.ln2.GetEstimate(prec).Value();
	LongFloat pwr = (a.Value() / ln2).RoundTowardZero();
	LongFloat flt = a.Value() - pwr * ln2;

	LongFloat pi2 = Constants.pi.GetEstimate(prec).Value() << 1;
	LongFloat alpha = arg.imag().Value() - pi2 * (a.Value() / pi2).round();
	if (alpha.IsNegative()) {
		sinneg = !sinneg;
		alpha = -alpha;
	}

	double dz = flt.AsDouble();
	double da = alpha.AsDouble();
	double de(::exp(dz));
	LongFloat s(::sin(da) * de);
	LongFloat c(::cos(da) * de);
	CmplxEst z(CmplxEst(Estimate(flt), Estimate(alpha)) + CmplxEst(Estimate(1.0)));
	CmplxEst e(c, s);

	for (u32 i=50; i<prec*32; i=i*2-3) {
		CmplxEst l(lncomplex(e, i/16 + i/64 + 4)); 
		e = e * (z - l);
		// clear error
		e = CmplxEst(e.real().Value(), e.imag().Value());
	}
	CmplxEst r(e * (z - lncomplex(e, prec)));
	ErrorEstimate err(r.real().Error() + ErrorEstimate(r.real().Value() - e.real().Value()) * ErrorEstimate(r.real().Value()).recip());
	ErrorEstimate eri(r.imag().Error() + ErrorEstimate(r.imag().Value() - e.imag().Value()) * ErrorEstimate(r.imag().Value()).recip());

	r = r << i32(pwr.AsDouble());

	ErrorEstimate abv(ErrorEstimate(arg.real().Value()) * arg.real().Error());
	ErrorEstimate abe(ErrorEstimate(arg.imag().Value()) * arg.imag().Error());
	abe = abv * abv + abe * abe;
	abe = ::exp(abe.AsDouble()) - 1.0;
	abv = ErrorEstimate(r.real().Value()) * ErrorEstimate(r.real().Value()) +
          ErrorEstimate(r.imag().Value()) * ErrorEstimate(r.imag().Value());
	abe = abe * abv;

	r = AdjustForAbsoluteError(r, abe);

	r = sinneg ? conj(r) : r;

	return arg.real().Value().IsNegative() ? CmplxEst(Estimate(1.0)) / r : r;
}

// atan2 by ln(complex(x, y)).im
Estimate atan2(const Estimate &y, const Estimate &x, u32 prec)
{
	CmplxEst z(x, y);
	z = lncomplex(z, prec);
	return z.imag();
}

// atan(arg) == atan2(arg, 1)
Estimate atan(const Estimate &arg, u32 prec)
{
	return atan2(arg, Estimate(1), prec);
}

// tan(arg) by exp complex
Estimate tan(const Estimate &arg, u32 prec)
{
	// put it in range:
	CmplxEst c(Estimate(), arg);
	c = expcomplex(c, prec);

	return c.imag() / c.real();
}

// cos(arg) = (1 - sq(tan(arg/2))) / (1 + sq(tan(arg/2)))
Estimate cos(const Estimate &arg, u32 prec)
{
	CmplxEst z(Estimate(), arg);
	z = expcomplex(z, prec);
	return z.real();
}

// sin(arg) = 2*tan(arg/2) / (1 + sq(tan(arg/2)))
Estimate sin(const Estimate &arg, u32 prec)
{
	CmplxEst z(Estimate(), arg);
	z = expcomplex(z, prec);
	return z.imag();
}

// asin(arg) = atan(arg / sqrt(1 - sq(arg)))
Estimate asin(const Estimate &arg, u32 prec)
{
	return atan(arg * rsqrt(Estimate(1) - sq(arg), prec), prec);
}

// acos(arg) = atan(sqrt(1 - sq(arg)) / arg)
Estimate acos(const Estimate &arg, u32 prec)
{
	return atan2(Estimate(1), arg * rsqrt(Estimate(1) - sq(arg), prec), prec);
}
 */
} // namespace
