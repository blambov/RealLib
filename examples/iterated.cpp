#include <stdio.h>
#include <time.h>
#include "../Real.h"

using namespace RealLib;

int main()
{
    char buf[100];
    clock_t starttime = clock();
    InitializeRealLib();

    {
        Real x(0.5);
        Real c(3.75);
        int num = 1000;

        for (int i=1; i<=num; ++i)
            x = c*x*(1-x);

        printf("%dth: %s\n", num, x.AsDecimal(buf, 15));
    }

    unsigned pr = FinalizeRealLib();

    clock_t endtime = clock();
    printf("prec: %d time elapsed: %lf\n", pr,
           double(endtime - starttime) / CLOCKS_PER_SEC);

    getchar();
    return 0;
}
