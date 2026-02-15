/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#ifndef TIMER_H_DEF
#define TIMER_H_DEF

int GetNewTimer(void);
int StartNewTimer(void);
void StartTimer(int n);
float StopTimer(int n);

/* FORTRAN stubs */
void starttimer_(int *n);
float stoptimer_(int *n);
void startnewtimer_(int *n);
void getnewtimer_(int *n);

#endif
