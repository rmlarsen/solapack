/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#include <stdio.h>
#include <stdlib.h>
#include "diskmalloc.h"
#include "maxsize.h"


#define SOLA2D sola2d_lanczos_

extern void SOLA2D(void *, void *, void *);

int main(int argc, char *argv[])
{
  void *U, *V, *target;
  size_t Vsize;
  int V_is_mmaped = 0;

  puts("Trying to allocate workspace with valloc.");
  U = safevalloc(sizeof(double) * MAXPTS * MAXITER);
  target = safevalloc(sizeof(double) * MAXPTS * MAXPTS);
  Vsize = sizeof(double) * (size_t)MAXNLM * MAXITER + MAXPTS + MAXITER + MAXNLM;
  printf("Vsize = %zu\n", Vsize);
  V = valloc(Vsize);
  if (V == NULL)
  {
    puts("valloc failed, attempting to allocate workspace using mmap.");
    V = diskmalloc(Vsize, "V.");
    if (V == NULL)
      exit(-1);
    V_is_mmaped = 1;
  }

  SOLA2D(U, V, target);

  free(U);
  free(target);
  if (V_is_mmaped)
    diskfree(V);
  else
    free(V);
  return 0;
}
