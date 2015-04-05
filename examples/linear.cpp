#include <time.h>
#include <assert.h>
#include <valarray>
#include <iostream>
#include <iomanip>
#include "Real.h"
#include "RealEstimate.h"

using namespace RealLib;
using namespace std;

//#define Real double
//#define Estimate double
//#define recip(x) (1/(x))

#define EstVector valarray< TYPE >
#define EstMatrix valarray< EstVector >

template <class TYPE>
ostream& operator << (ostream& os, const valarray<TYPE> &v)
{
    int N = v.size();
    for (int i=0;i<N;++i)
        cout << v[i] << endl;
    return os;
}

template <class TYPE>
void doGaussian(EstMatrix &mat)
{
    size_t N = mat.size();
    size_t M = mat[0].size();
    int i,j;

    //cout << "pre" << endl << mat << endl;

    for (i=0;i<N;++i) {
        TYPE div(mat[i][i]);

        // no pivoting?

        mat[i] /= div;

        for (j=0;j<N;++j)
            if (j!=i) {
                TYPE m(mat[j][i]);
                mat[j] -= mat[i] * m;
            }
    }

    //   for (i=0;i<N;++i)
    //     mat[i] = mat[i] / mat[i][i];
}

Real anal(int i, int j, int n)
{
    Real I(i), J(j);
    Real v(I + J + 1);
    int k;

    for (k=n-j;k<=n+i;++k)
        v = v*k/(n+i-k+1);
    for (k=n-i;k<=n+j;++k)
        v = v*k/(n+j-k+1);
    for (k=i+1;k<=i+j;++k)
        v = v*k/(i+j-k+1);
    for (k=j+1;k<=i+j;++k)
        v = v*k/(i+j-k+1);
    return i+j & 1 ? -v : v;
}

int main()
{
    clock_t starttime = clock();
    InitializeRealLib();
    const int N = 100;
    //cout << fixed;

    try {
        int i,j;
        typedef valarray<Real> RealVec;
        typedef valarray<RealVec> RealMat;
        RealVec row(Real(), N+N);
        RealMat mat(row, N);

        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                mat[i][j] = recip(Real(i+j+1));
                mat[i][N+j] = 0.;
            }
            mat[i][N+i] = 1.;
        }
        for (i=0;i<N;++i)
            mat[i][i] += 1;

        doGaussian(mat);
        /*{
        ArrayInterface<Real, sizeof (Real)> ai(&mat[0][0], N*(N+N));
        Gaussian<Real, ArrayInterface<Real, sizeof (Real)> >(ai, N+N);
      } */ 


        for (i=0;i<N;++i)
            for (j=0;j<N;++j) {
                //i = j = N-1;

                cout << "inv[" << i << ", " << j << "]: " << mat[i][N+j];
                //cout << " analytical: " << anal(i, j, N);
                cout << endl;
            }

    }
    catch (RealLibException &e) {
        cout << "Exception caught: " << e;
    }

    unsigned pr = FinalizeRealLib();

    clock_t endtime = clock();
    printf("prec: %d time elapsed: %lf\n", pr, 
           double(endtime - starttime) / CLOCKS_PER_SEC);

    return 0;
}
