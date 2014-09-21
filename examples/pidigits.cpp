#include <stdio.h>
#include <time.h>
#include "Real.h"
#include "RealEstimate.h"

using namespace RealLib;

#define len 100000+5

int main()
{
	char* buf = new char[len];
	clock_t starttime = clock();
   FILE *f = fopen("pi.txt", "w");
	InitializeRealLib(5, 5);


	{
		sprintf(buf, "1e-%d", len);
		Real L(buf);
		printf("comparing with %s\n", buf);
		
		bool done = false;
		while (!done) {
			try {
				if ((Pi - Pi) < L) {
				   clock_t midtime = clock();
				   printf("computed in %lf seconds\n",
		            double(midtime - starttime) / CLOCKS_PER_SEC);
				   
				   starttime = midtime;
				}
				done = true;
			} catch (PrecisionException e) {
				printf("precision exception %s, ", e.what());
				int prec = FinalizeRealLib() * 2;
				printf("increasing precision to %d\n", prec);
				InitializeRealLib(prec, prec);
			}
		}

		Pi.AsDecimal(buf, len);
		
      // we're printing to file to avoid delays caused 
      // by displaying this quantity of information
		fputs(buf, f);
		buf[50] = 0;
		printf("pi is (to 50 digits): %s\n", buf);
      printf("look at pi.txt to see the full string\n");
	}

	unsigned pr = FinalizeRealLib();
	fclose(f);

	clock_t endtime = clock();
	printf("prec: %d printed in %lf seconds\n", pr, 
      double(endtime - starttime) / CLOCKS_PER_SEC);

	getchar();
	return 0;
}
