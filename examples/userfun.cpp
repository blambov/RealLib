#include <iostream>
#include <iomanip>
#include <ctime>
#include "RealEstimate.h"
#include "Real.h"
#include "RealFuncs.h"

using namespace RealLib;
using namespace std;

Estimate etaylor(unsigned int prec)
{
   Estimate s(0.0);
   Estimate m(1.0);
   
	for (int i=1; i<=prec; ++i) {
		s += m;
		m /= i;
	}
   return s.AddError(m*2);
}

Estimate exptaylor(const Estimate &arg)
{
	unsigned int prec = arg.GetPrecision();
   Estimate s(0.0);
   Estimate m(1.0);
   
   if (abs(arg) > 1.0) throw DomainException("exptaylor");
   if (!(abs(arg) < 1.0)) throw PrecisionException("exptaylor");
   
	for (int i=1; i<=prec; ++i) {
		s += m;
		m = m*arg/i;
	}
   return s.AddError(m*2);
}

Real exptaylor(const Real &arg)
{
	return Real(exptaylor, arg);
}

int main()
{
	clock_t starttime = clock();
	InitializeRealLib(10);

	{
		Real one(Real(etaylor) / sq(exptaylor(Real(0.5))));

		try {
		cout << "etaylor / sq(exptaylor(0.5)):\t" << setprecision(50) << showpoint << one << endl;
		} catch (RealLibException e) {
		cout << "exception: " << e.what();
		}
	}

	unsigned pr = FinalizeRealLib();

	clock_t endtime = clock();
	cout << "prec: " << pr << " time elapsed: " << setprecision(6) <<
      double(endtime - starttime) / CLOCKS_PER_SEC << endl;

	cin.get();
	return 0;
}
