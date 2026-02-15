/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include "timer.h"
#include <stdio.h>
#include <sys/time.h>
#define N 100

static struct timeval first[N],second[N];

static struct timezone tzp[N];
static int nexttimer = 1;





/* Wall clock timer routines */
int GetNewTimer(void)
{
  int new;
  new = nexttimer;
  nexttimer = (nexttimer % N)+1;
  return new;
}

int StartNewTimer(void)
{
  int new;
  new = GetNewTimer();
  StartTimer(new);
  return new;
}

void StartTimer(int n)
{
  if (n<1 || n > N) 
    fprintf(stderr,"StartTimer: Timer number should be between 1 and %d\n",N);
  else {
    n = n-1; 
    gettimeofday (&first[n], &tzp[n]);
  }
}

float StopTimer(int n)
{
  if (n<1 || n > N) 
    fprintf(stderr,"StopTimer: Timer number should be between 1 and %d\n",N);
  else {
    n = n-1; 
    gettimeofday (&second[n], &tzp[n]);
    if (first[n].tv_usec > second[n].tv_usec) {  
      second[n].tv_usec += 1000000;
      second[n].tv_sec--;
    }
  }
  return (float) (second[n].tv_sec-first[n].tv_sec) +  
    (float) (second[n].tv_usec-first[n].tv_usec)/1000000.0;
}


/* FORTRAN Interface routines: */
void starttimer_(int *n)
{
  StartTimer(*n);
}

float stoptimer_(int *n)
{
  return StopTimer(*n);
}

void startnewtimer_(int *n)
{
  *n = StartNewTimer();
}

void getnewtimer_(int *n)
{
  *n = GetNewTimer();
}

