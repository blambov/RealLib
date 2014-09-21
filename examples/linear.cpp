#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <valarray>
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
void doGaussian(EstMatrix &mat)
{
   int N = mat.size();
   int M = mat[0].size();
   int i,j;

   for (i=0;i<N;++i) {
      j=i;
      TYPE div(mat[i][i]);

      // pivoting?
      /*
      while (div.IsNonZero()) 
         if (++j==N) throw PrecisionException("gaussian");
         else {
            div = mat[j][i];
         }
      if (i!=j) swap(mat[i], mat[j]);
      */

      div = recip(div);
      for (++j;j<N;++j)
         mat[j] = mat[j] - mat[i] * (mat[j][i] * div);
   }
}

template <class TYPE>
void solveGaussian(EstMatrix &mat)
{
   int N = mat.size();
   int M = mat[0].size();
   int i,j;
   assert(M == N+N);

   for (j=N-1; j>=0; --j) {
      for (i=N-1; i>j; --i)
         mat[j] -= mat[i] * mat[j][i];
      TYPE div(mat[j][j]);

      mat[j] = mat[j] / div;
   }
}

template <class TYPE, class ARRAY>
void Gaussian(ARRAY &arr, long M)
{
   int N = arr.size() / M;
   int i;
   assert(M*N == arr.size());

   // the STL in GCC does not allow changes in the length of a valarray on assignment,
   // therefore the rows of the matrix are initialized with zero vectors of length M
   EstVector vec(M);
   EstMatrix matrix(vec, N);
	i=0;
    for (int k=0;k<N;++k)
      for (int j=0;j<M;++j)
         matrix[k][j] = arr[i++];
/*  for (i=0;i<N;++i) {
      matrix[i] = EstVector(arr+i*M, M);
   }*/
   
   doGaussian(matrix);
   solveGaussian(matrix);

   i=0;
   for (int k=0;k<N;++k)
      for (int j=0;j<M;++j)
         arr[i++] = matrix[k][j];
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

void printmatrix(Real *mat, int N, int M)
{
   printf("matrix:");
   for (int i=0;i<N;++i) {
      printf("\nrow %d:", i);
      for (int j=0;j<M;++j)
         printf(" %7f", mat[i*M+j].AsDouble());
   }
   printf("\n");
}

CreateArrayAndIntRealFunction(Gaussian)

int main()
{
	clock_t starttime = clock();
	InitializeRealLib(12);
   const int N = 50;

   {
      int i,j;
      Real mat[N][N+N];
      Real src[N][N+N];
      for (i=0;i<N;++i) {
         for (j=0;j<N;++j) {
            mat[i][j] = src[i][j] = recip(Real(i+j+1));
            mat[i][N+j] = 0.;
         }
         mat[i][N+i] = 1.;
      }

		Gaussian((Real*)mat, N*(N+N), N+N);
      //Real(Gaussian, (Real*)mat, N*(N+N), (void*) (N+N));
      //Gaussian((Real*)mat, N*(N+N), (void*) (N+N));

      /*
      for (i=0;i<N;++i)
         for (j=0;j<N;++j)*/
      i = j = N-1;
            printf("inv[%d, %d]: %10.0lf analytical: %10.0lf\n", i, j, mat[i][N+j].AsDouble(), anal(i, j, N).AsDouble());
            //printf("inv[%d, %d]: %10.0lf analytical: %10.0lf\n", i, j, mat[i][N+j], anal(i, j, N));
      
   }

   unsigned pr = FinalizeRealLib();

	clock_t endtime = clock();
	printf("prec: %d time elapsed: %lf\n", pr, 
double(endtime - starttime) / CLOCKS_PER_SEC);

   getchar();
	return 0;
}
