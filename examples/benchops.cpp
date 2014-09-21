// real.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Real.h"

using namespace RealLib;

char *sizes[] = {"1e-100", "1e-1000", "1e-2500", "1e-10000", "1e-100000", "1e-1000000"}; 
double times[] = {200, 10, 1, 1, 0.1, 0.01};

int main()
{
	clock_t starttime, endtime;
	InitializeRealLib(11, 1000000);
	FILE *out = fopen("res.txt", "w+t");

	for (int i=0; i<4; ++i)
	{
		Real cmp(sizes[i]);
		Real half(0.5);
		int j;

		// warm up
		if (!(half+half - half*2 < cmp)) break;

		#define BENCHOP(x, tm) \
			starttime = clock(); \
			for (j=0; j<times[i]*tm; ++j) { \
				Real z(x); \
				if (!((z)-(z) < cmp)) break; \
			} \
			endtime = clock(); \
			fprintf(out, #x": %lf\n", double(endtime - starttime) / CLOCKS_PER_SEC/ ceil(times[i]*tm)); \
			printf(#x": %lf\n", double(endtime - starttime) / CLOCKS_PER_SEC / ceil(times[i]*tm));

		printf("prec: %s\n", sizes[i]);
		BENCHOP(half + half, 500);
		BENCHOP(half * half, 100);
		BENCHOP(recip(half), 10);
		BENCHOP(sqrt(half), 10);
		BENCHOP(Pi, 1.9);
		BENCHOP(Ln2, 1.8);
		BENCHOP(log(half), 0.75);
		BENCHOP(exp(half), 0.25);
		BENCHOP(sin(half), 0.4);
		BENCHOP(asin(half), 0.25);
	}

	fclose(out);
	unsigned pr = FinalizeRealLib();

	printf("prec: %d\n", pr);

	getchar();
	return 0;
}

