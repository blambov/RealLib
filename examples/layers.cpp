#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include "RealEstimate.h"
#include "MachineEstimate.h"
#include "Real.h"
#include "RealFuncs.h"
#include "RealEncapsulation.h"

using namespace RealLib;
using namespace std;

int g_len;


template <class TYPE>
TYPE Harmonic(unsigned int prec, UserInt len)
{
    TYPE s(1.0);
    TYPE one(1.0);
    prec; // to avoid warning


    for (int i=1; i<=len; ++i) {
        s += one / i;
    }
    return s;
}

#define PREC 8
#define INITIALPREC 4

char *s_Harmonic = "s += one / i";
CreateIntRealFunction(Harmonic)

#define MakeFunction(name, operation, extrainit) \
template <class TYPE> \
TYPE name(unsigned int prec, UserInt len) \
{ \
    TYPE s(1.0); \
    TYPE one(1.0); \
	extrainit; \
    prec; \
 \
    for (int i=1; i<=len; ++i) { \
        operation; \
    } \
    return s; \
} \
	char * s_ ## name = #operation; \
CreateIntRealFunction(name)

namespace RealLib {
template<>
Real pi(unsigned int prec) {
	prec;
	return Pi;
}

template<>
double pi(unsigned int prec) {
	prec;
	return M_PI;
}
}

MakeFunction(Addition, s += one, )
MakeFunction(Multiplication, s *= v, TYPE v = one + TYPE(1.e-8))
MakeFunction(Reciprocal, s = one / s, s += one)
MakeFunction(Sine, s = sin(s), s = one)
MakeFunction(Logistic, s = coeff * s * (one - s), TYPE coeff(3.75); s = s / 2)
MakeFunction(Sqrt, s = sqrt(s), s = one + one)
MakeFunction(Abs, s = abs(s), s = one)
MakeFunction(Tan, s = tan(s), TYPE v = pi<TYPE>(prec)/4)
MakeFunction(Exp, s = exp(s - one), s = TYPE(0.0))
MakeFunction(ASine, s = v * asin(s), s = one / 2; TYPE v = 3 / pi<TYPE>(prec))
MakeFunction(ATan, s = v * atan(s), s = one; TYPE v = 4 / pi<TYPE>(prec))

int main()
{
    clock_t starttime, endtime;
    double t0, t1, t2, t3;
    char ch;

    {
        InitializeRealLib(INITIALPREC);
	int itc, c;

	cout << setprecision(70) << Harmonic(100000) << endl;
	//cin >> ws >> ch;
	
	ResetRealLib(INITIALPREC);

#define BENCH(function, endtype, text, itcount, res) { \
		itc = itcount ? itcount : 1; \
		ResetRealLib(INITIALPREC); \
		starttime = clock(); \
		try { \
			endtype rr;\
			c=0;\
			for (int i=0;i<itc;++i) \
			{ \
				rr = endtype(function); \
				if (rr > -2) ++c; \
			} \
			cout << setprecision(PREC); \
			cout << text << rr << endl; \
		} catch (PrecisionException e) { \
			cout << text << "PrecisionException in " << e.what() << endl; \
		} catch (DomainException e) { \
			cout << text << "DomainException in " << e.what() << endl; \
		} catch (...) { \
			cout << text << "unknown exception " << endl; \
		} \
		endtime = clock(); \
		res = double(endtime - starttime) / CLOCKS_PER_SEC / itc; }

#define GLUE(x, y) x ## y
		
//	BENCH(function<Real>(0), Real,         "Real:            ", ic3, t3);	\

#define BENCHALL(function, len, ic0, ic1, ic2, ic3) \
	cout << "Benchmarking " #function " (" << s_ ## function << "), " << len << " members" << endl; \
	BENCH(function<double>(0, len), double,     "double:          ", ic0, t0); \
	BENCH(function<MachineEstimate>(0, len), MachineEstimate,"MachineEstimate: ", ic1, t1); \
	BENCH(function(len), Real,						"Real, function:  ", ic3, t2);	\
	BENCH(function<Real>(0, len), Real,         "Real, direct:    ", ic2, t3); \
	\
        cout << setprecision(6); \
        cout << "time: double " << t0 << " mach " << t1 << " est " << t2 << " real " << t3 << endl \
			<< "mega ops per second: " << 1.e-6*g_len / t0 << ", " << 1.e-6*g_len / t1 \
			<< ", " << 1.e-6*g_len/t2 << ", " << 1.e-6*g_len/t3 << endl \
			<< "ratios 1:" << t1/t0 /*<< ":" << t2/t0 << ":" << t3/t0*/ \
			<< ", 1:" << t2/t1 /*<< ":" << t3/t1*/ << ", 1:" << t3/t2 << endl \
			<< endl; 
//	cin >> ch;

		BENCHALL(Addition, 300000, 1000, 1000, 2, 1000);
		BENCHALL(Multiplication, 300000, 1000, 300, 2, 300);
		BENCHALL(Reciprocal, 100000, 400, 200, 2, 200);
		BENCHALL(Sqrt, 60000, 4000, 600, 20, 600);
		BENCHALL(Abs, 60000, 4000, 2000, 20, 2000);
		BENCHALL(Sine, 30000, 400, 200, 1, 200);
		BENCHALL(Tan, 300, 4000, 2, 1, 2);
		BENCHALL(Exp, 30000, 400, 200, 1, 200);
		BENCHALL(ASine, 300, 400, 200, 1, 200);
		BENCHALL(ATan, 300, 400, 200, 1, 200);

		BENCHALL(Harmonic, 1000, 40000, 20000, 200, 20000);
		BENCHALL(Harmonic, 10000, 4000, 2000, 20, 2000);
		BENCHALL(Harmonic, 100000, 400, 200, 2, 200);
		BENCHALL(Harmonic, 1000000, 40, 200, 1, 200);

		BENCHALL(Logistic, 40, 100000, 10000, 1000, 500);
		BENCHALL(Logistic, 80, 100000, 1000, 1000, 500);
		BENCHALL(Logistic, 120, 100000, 1000, 1000, 500);
		BENCHALL(Logistic, 150, 100000, 1000, 1000, 1000);
		
		FinalizeRealLib();
    }

	 //cin >> ch;
    return 0;
}
