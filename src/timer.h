/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#ifndef SOLAPACK_TIMER_H
#define SOLAPACK_TIMER_H

void StartTimer(int n);
float StopTimer(int n);
int GetNewTimer(void);
int StartNewTimer(void);

/* Fortran interface */
void starttimer_(int *n);
float stoptimer_(int *n);
void startnewtimer_(int *n);
void getnewtimer_(int *n);

#endif
