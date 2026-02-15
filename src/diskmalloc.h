/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

#ifndef __DISKMALLOC_H_DEF
#define __DISKMALLOC_H_DEF


typedef struct memblock {
  void *addr;             /* Starting address of memory block */
  size_t size;            /* Length in bytes */
  int fildes;             /* Corresponding file */
  char *name;             /* Name of file */
  struct memblock *next;       /* Next block */
} memblock_t;


/* Prototypes */
void *safemalloc(size_t);
void *safevalloc(size_t);
void *diskmalloc(int , char *);
void diskfree(void *);


/* FORTRAN stubs */
void *diskmalloc_(size_t *);
void diskfree_(void **);
void diskfree1_(void *);
void *safemalloc_(size_t *);
void *safevalloc_(size_t *);



#endif
