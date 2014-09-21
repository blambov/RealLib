#include <iostream>
#include <iomanip>
#include "Real.h"

using namespace RealLib;
using namespace std;

template <class TYPE>
TYPE myexp(const TYPE &a)
{
	unsigned int prec = a.GetPrecision();
   TYPE s(0.0);
   TYPE m(1.0);
   
	TYPE arg(a.TruncateTo(-1.0, 1.0, "myexp"));

	if (abs(arg).weak_le(0.75)) {
		TYPE err = (TYPE(1) >> (32 * prec)) / 3;
		for (unsigned i=1; abs(m) > err; ++i) {
			s += m;
			m = m*arg/i;
		}
	} else {
		if (prec < 6) prec = 6;
		unsigned int pc = prec * 23;
		for (unsigned i=1; i<=pc; ++i) {
			s += m;
			m = m*arg/i;
		}
	}
   return s.AddError(abs(m)*3); 
}
CreateUnaryRealFunction(myexp)

int main()
{
	InitializeRealLib();
	{
		Real a(myexp(Real(0.5)));
		Real b(exp(Real(1)));

		cout << fixed << setprecision(10);
		cout << "a(myexp(0.5)) =\t" << a << endl;
		cout << "a*a =\t\t" << a*a << endl;
		cout << "b(exp(1)) =\t" << b << endl;
		cout << "a*a/b =\t\t" << setprecision(300) << showpoint 
			<< a*a/b << endl;
	}
	cout << "precision used: " << FinalizeRealLib();

	return 0;
}
