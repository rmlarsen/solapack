/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include "timer.h"
#include <stdio.h>
#include <sys/time.h>

#define MAX_TIMERS 100
#define USEC_PER_SEC 1000000

static struct timeval first[MAX_TIMERS], second[MAX_TIMERS];
static int nexttimer = 1;

/* Wall clock timer routines */
int GetNewTimer(void)
{
  int new;
  new = nexttimer;
  nexttimer = (nexttimer % MAX_TIMERS) + 1;
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
  if (n < 1 || n > MAX_TIMERS)
    fprintf(stderr, "StartTimer: Timer number should be between 1 and %d\n", MAX_TIMERS);
  else {
    n = n - 1;  /* Convert from 1-based (Fortran) to 0-based (C) index */
    gettimeofday(&first[n], NULL);
  }
}

float StopTimer(int n)
{
  if (n < 1 || n > MAX_TIMERS) {
    fprintf(stderr, "StopTimer: Timer number should be between 1 and %d\n", MAX_TIMERS);
    return 0.0f;
  }
  n = n - 1;  /* Convert from 1-based (Fortran) to 0-based (C) index */
  gettimeofday(&second[n], NULL);
  if (first[n].tv_usec > second[n].tv_usec) {
    second[n].tv_usec += USEC_PER_SEC;
    second[n].tv_sec--;
  }
  return (float)(second[n].tv_sec - first[n].tv_sec) +
    (float)(second[n].tv_usec - first[n].tv_usec) / (float)USEC_PER_SEC;
}


/* Fortran interface routines (1-based timer indices): */
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
