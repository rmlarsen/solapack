/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include <stdio.h>
#include <stdlib.h>
#include "diskmalloc.h"
#include "maxsize.h"


#define SOLA2D sola2d_lanczos_

extern void SOLA2D(void *, void *, void *);

int main(void)
{
  void *U, *V, *target;
  size_t Vsize;
  int disk = 0;

  puts("Trying to allocate workspace with valloc.");
  U = safevalloc(sizeof(double) * MAXPTS * MAXITER);
  target = safevalloc(sizeof(double) * MAXPTS * MAXPTS);
  Vsize = sizeof(double) * MAXNLM * MAXITER + MAXPTS + MAXITER + MAXNLM;
  printf("Vsize = %zu\n", Vsize);
  V = valloc(Vsize);
  if (V == NULL)
  {
    puts("Malloc failed, attempting to allocate workspace using mmap.");
    V = diskmalloc(Vsize, "V.");
    if (V == NULL)
      exit(-1);
    disk = 1;
  }

  SOLA2D(U, V, target);

  if (disk) {
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
