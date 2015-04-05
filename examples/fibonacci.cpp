#include <iostream>
#include "Real.h"

using namespace RealLib;
using namespace std;

int main()
{
    InitializeRealLib();

    int N = 1000000000;

    Real sqf(sqrt(Real(5)));
    Real a((Real(1) + sqf)/2);
    Real b((sqf - Real(1))/2);

    a = exp(log(a)*N);
    b = exp(log(b)*N);

    if (N % 2 == 1) b = -b;
    a = (a-b)/sqf;

    cout << N << "th Fibonacci number: " << a;

    FinalizeRealLib();
    cin.get();
    return 0;
}
