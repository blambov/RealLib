#include <stdio.h>
#include <time.h>
#include "math.h"
#include "Real.h"

using namespace RealLib;

int main()
{
	clock_t starttime = clock();
	InitializeRealLib();
	{
		Real exr(exp(Real(1.0)));
		double ex(exp(1.0));
		double s(ex), m(1.0);

		for (int i=1; i<=1000; ++i) {
			s -= m;
			m /= i;
		}

		double diff(s);

		{	
			printf("exp(1.0):\t%20le\n", ex);
			printf("taylor: \t%20le\n", ex - s);
			printf("difference:\t%20le\n", diff);
			printf("diff by real:\t%20le\n", (ex - exr).AsDouble());
		}

	}

	unsigned pr = FinalizeRealLib();

	clock_t endtime = clock();
	printf("prec: %d time elapsed: %lf\n", pr, 
double(endtime - starttime) / CLOCKS_PER_SEC);

	getchar();
	return 0;
}
