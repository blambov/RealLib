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

template <class TYPE>
TYPE* InitializeMulTest(int count)
{
    srand(1);
    TYPE *p = new TYPE[count];
    for (int i=0;i<count;++i) {
        int m(rand());
        int n(rand());
        int d(rand());
        switch (m%2) {
        case 0:
            p[i] = m*(double(n)/d)/65536;
            break;
        case 1:
            p[i] = -m*(double(n)/d)/65536;
            break;
        }
        //cout << i << ": " << p[i] << endl;
    }
    return p;
}

template <class TYPE>
TYPE MulWorst(unsigned int prec, UserInt len)
{
    TYPE m(1.0);
    int slen = 256;
    TYPE *p = InitializeMulTest<TYPE>(slen);
    for (int i=0;i<len;++i)
        for (int j=0;j<slen;j+=2)
            m += p[j]*p[j+1];
    delete p;
    return m;
}

const char *s_MulWorst = "multiply array";
CreateIntRealFunction(MulWorst)

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

#define PREC 6
#define INITIALPREC 4

const char *s_Harmonic = "s += one / i";
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
const char * s_ ## name = #operation; \
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
MakeFunction(Sine, s = sin(one += c), s = one; TYPE c(0.001))
MakeFunction(Cosine, s = cos(one += c), s = one; TYPE c(0.001))
MakeFunction(Logistic, s = coeff * s * (one - s), TYPE coeff(3.75); s = s / 2)
MakeFunction(Sqrt, s = sqrt(s), s = one + one)
MakeFunction(Abs, s = abs(s), s = one)
MakeFunction(Tan, s = tan(one += c), TYPE c(0.001)) //TYPE v = pi<TYPE>(prec)/4)
MakeFunction(Exp, s = exp(s - one), s = TYPE(0.0))
MakeFunction(ASine, s = asin(one -= c), TYPE c(0.000099))//s = v * asin(s), s = one / 2; TYPE v = 3 / pi<TYPE>(prec))
MakeFunction(ATan, s = atan(one += c), TYPE c(0.1))//v * atan(s), s = one; TYPE v = 4 / pi<TYPE>(prec))
MakeFunction(Log, s = log(one += c), TYPE c(0.1))

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
        //	BENCH(function<Real>(0, len), Real,         "Real, direct:    ", ic2, t3); \

#define BENCHALL(function, len, ic0, ic1, ic2, ic3) \
        cout << "Benchmarking " #function " (" << s_ ## function << "), " << len << " members" << endl; \
        BENCH(function<double>(0, len), double,     "double:          ", ic0, t0); \
        BENCH(function<MachineEstimate>(0, len), MachineEstimate,"MachineEstimate: ", ic1, t1); \
        BENCH(function(len), Real,						"Real, function:  ", ic3, t2);	\
        \
        cout << setprecision(6); \
        cout << "time: double " << t0 << " mach " << t1 << " est " << t2 << " real " << t3 << endl \
        << "mega ops per second: " << 1.e-6*len / t0 << ", " << 1.e-6*len / t1 \
        << ", " << 1.e-6*len/t2 << ", " << 1.e-6*len/t3 << endl \
        << "ratios 1:" << t1/t0 /*<< ":" << t2/t0 << ":" << t3/t0*/ \
        << ", 1:" << t2/t0 /*<< ":" << t3/t1*/ << ", 1:" << t3/t0 << endl \
                << endl;
        //	cin >> ch;

        BENCHALL(Addition, 300000, 1000, 1000, 2, 1000);
        BENCHALL(Multiplication, 300000, 1000, 300, 2, 300);
        BENCHALL(MulWorst, 30000, 100, 10, 1, 10);
        BENCHALL(Reciprocal, 100000, 400, 200, 2, 200);
        BENCHALL(Sqrt, 60000, 4000, 600, 20, 600);
        BENCHALL(Abs, 60000, 4000, 2000, 20, 2000);
        BENCHALL(Sine, 20000, 800, 400, 1, 200);
        BENCHALL(Cosine, 20000, 800, 400, 1, 200);
        BENCHALL(Tan, 20000, 800, 400, 1, 200);
        BENCHALL(Exp, 20000, 800, 400, 1, 200);
        BENCHALL(Log, 20000, 800, 400, 1, 200);
        BENCHALL(ASine, 20000, 800, 400, 1, 200);
        BENCHALL(ATan, 20000, 800, 400, 1, 200);

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
