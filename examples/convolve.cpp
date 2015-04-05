// real.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "Real.h"
#include "convolution.h"
#include "convolution.cpp"

// if you get an error saying Convolution<double> is defined twice, 
// exclude convolution.cpp (or .o) from the link

using namespace RealLib;
using namespace std;

#define size 256

int main()
{
    clock_t starttime, endtime;
    InitializeRealLib();

    starttime = clock();

    {
        int i;
        Real a[size], b[size];

        Convolution<Real> conv(size, Pi * 2);

        for (i=0;i<size/2;++i) {
            a[i] = i+1; b[i] = 5+i;
        }

        conv.Convolve(a, b);

        cout << setprecision(7);

        for (i=0;i<size-1;++i)
            cout << "result[" << setw(4)  << i << "]: " <<
            setw(8) << a[i] << endl;
    }

    int pr = FinalizeRealLib();
    endtime = clock();

    cout << "prec: " << pr << " time elapsed: " <<
            double(endtime - starttime) / CLOCKS_PER_SEC << endl;

    cin.get();
    return 0;
}

