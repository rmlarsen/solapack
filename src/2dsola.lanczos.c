/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include <stdio.h>
#include <stdlib.h>
#include "diskmalloc.h"
#include "maxsize.h"


#define SOLA2D sola2d_lanczos_

extern void SOLA2D(void *, void *, void *);

int main(void)
{
  void *U,*V,*target;
  int Vsize,disk=0;

  puts("Trying to allocate workspace with valloc.");
  U = (double *) safevalloc(sizeof(double)*MAXPTS*MAXITER);
  target = (double *) safevalloc(sizeof(double)*MAXPTS*MAXPTS);
  Vsize = sizeof(double)*MAXNLM*MAXITER+MAXPTS+MAXITER+MAXNLM; 
  /*   Vsize = sizeof(double)*MAXNLM*MAXPTS+MAXPTS+MAXITER+MAXNLM;  */
  printf("Vsize = %d\n",Vsize);
  V = (double *) valloc(Vsize);
  if (V == NULL)
  {
    puts("Malloc failed, attempting to allocate workspace using mmap.");
    puts("Allocating workspace using mmap.");
    V = (double *) diskmalloc(Vsize,"V.");
    if (V == NULL)
      exit(-1);
    disk = 1;
    /*    U = (double *) diskmalloc(sizeof(double)*MAXPTS*MAXITER,"U.");
    if (U == NULL)
    {
      diskfree(V);
      exit(-1);
    }
    target = (double *) diskmalloc(sizeof(double)*MAXPTS*MAXPTS,"TGT.");
    if (target == NULL)
    {
      diskfree(V);
      diskfree(U);
      exit(-1);
      } */
  }
  /*  else
  { 
    U = (double *) safevalloc(sizeof(double)*MAXPTS*MAXITER);
    target = (double *) safevalloc(sizeof(double)*MAXPTS*MAXPTS);
  }  */
  /*  scanf("%d",&disk);
  if (disk==1) {
    puts("Allocating workspace using mmap.");
    U = (double *) diskmalloc(sizeof(double)*MAXPTS*MAXITER,"U.");
    Vsize = sizeof(double)*MAXNLM*MAXITER+MAXPTS+MAXITER+MAXNLM;
    printf("Vsize = %d\n",Vsize);
    V = (double *) diskmalloc(Vsize,"V.");
    target = (double *) diskmalloc(sizeof(double)*MAXPTS*MAXPTS,"TGT.");
  }
  else {
  */

  SOLA2D(U,V,target);


  if (disk==1) {
    diskfree(U);
    diskfree(V);
    diskfree(target);
  }
  else {
    free(U);
    free(V);
    free(target);
  }
  return 0;
}


