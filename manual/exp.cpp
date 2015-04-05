#include <iostream>
#include <iomanip>
#include "Real.h"

using namespace RealLib;
using namespace std;

template <class TYPE>
TYPE myexp(const TYPE &arg)
{
    unsigned int prec = arg.GetPrecision();
    TYPE s(0.0);
    TYPE m(1.0);

    if (abs(arg) > 1.0) throw DomainException("myexp");
    if (!(abs(arg) < 1.0)) throw PrecisionException("myexp");

    for (unsigned i=1; i<=prec; ++i) {
        s += m;
        m = m*arg/i;
    }
    return s.AddError(m*3);
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
