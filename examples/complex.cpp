// real.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <ctime>
#include "Real.h"
#include "complex_fast.inl"

using namespace RealLib;
using namespace std;

template <class T>
ostream& operator<< (ostream &os, const complex<T>& c)
{
    return os << '(' << c.real() << ", " << c.imag() << ')';
}

int main()
{
    InitializeRealLib();

    {
        complex<Real> a(cos(Pi/(3*3)), sin(Pi/(3*3)));
        complex<Real> c(1);

        for (int i=0; i<54; ++i)
            c = c*a;

        cout << "c: " << c << endl;
    }

    int pr = FinalizeRealLib();

    cout << "prec: " << pr << endl;

    cin.get();
    return 0;
}

