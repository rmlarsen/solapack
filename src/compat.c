/*
 * Compatibility shims for legacy Fortran library functions that were
 * available in g77/libg2c but are not provided as linkable symbols
 * by modern gfortran.
 */

#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>

/* g77's rand() returned a double precision random number in [0,1). */
double rand_(void)
{
    return (double)random() / (double)RAND_MAX;
}

/* g77's dtime(tarray) returned elapsed CPU time since last call and
   filled tarray(1) with user time and tarray(2) with system time. */
float dtime_(float tarray[2])
{
    static struct tms last = {0, 0, 0, 0};
    struct tms now;
    float tick = 1.0f / (float)sysconf(_SC_CLK_TCK);

    times(&now);
    tarray[0] = (float)(now.tms_utime - last.tms_utime) * tick;
    tarray[1] = (float)(now.tms_stime - last.tms_stime) * tick;
    last = now;
    return tarray[0] + tarray[1];
}
