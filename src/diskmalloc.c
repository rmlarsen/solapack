/*     Copyright Rasmus Munk Larsen, Stanford University, 2003 */

/* Fortran callable stub for malloc */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "diskmalloc.h"


static memblock_t *head = NULL;

void *safemalloc(size_t size)
{
  void *p;

  if ( (p = malloc(size)) == NULL) {
    fprintf(stderr, "safemalloc failed with size = %zu\n", size);
    perror("malloc");
    exit(1);
  }
  return p;
}

void *safevalloc(size_t size)
{
  void *p;

  if ( (p = valloc(size)) == NULL) {
    fprintf(stderr, "safevalloc failed with size = %zu\n", size);
    perror("valloc");
    exit(1);
  }
  return p;
}


void *diskmalloc(size_t size, const char *prefix)
{
  void *buf;
  memblock_t *bp;
  char *dir, *name;
  char defaultdir[] = "/tmp";
  char endchar = '\0';
  int fildes;

  if (size < 1)
    return NULL;

  if (prefix == NULL)
    prefix = "solapack.";

  if ( (dir = getenv("SCRATCH")) == NULL)
    dir = defaultdir;

  name = malloc(strlen(dir) + strlen(prefix) + 8);
  if (name == NULL) {
    perror("diskmalloc: malloc for name failed");
    return NULL;
  }
  sprintf(name, "%s/%sXXXXXX", dir, prefix);

  /* Create a temp file atomically with mkstemp */
  fildes = mkstemp(name);
  if (fildes == -1) {
    fprintf(stderr, "Could not create temp file %s: %s\n", name, strerror(errno));
    free(name);
    return NULL;
  }

  if ( lseek(fildes, size - 1, SEEK_SET) == -1) {
    fprintf(stderr, "Could not lseek in file %s: %s\n", name, strerror(errno));
    goto abort;
  }
  if ( write(fildes, &endchar, 1) == -1) {
    fprintf(stderr, "Could not allocate storage for file %s: %s\n",
            name, strerror(errno));
    goto abort;
  }

  /* Map in file */
  buf = mmap(NULL, size, (PROT_READ | PROT_WRITE),
                         MAP_SHARED, fildes, 0);
  if (buf == MAP_FAILED) {
    perror("mmap failed");
    goto abort;
  }

  /* Insert new block in info list */
  if (head == NULL) {
    head = (memblock_t *) malloc(sizeof(memblock_t));
    bp = head;
  }
  else {
    bp = head;
    while (bp->next != NULL) bp = bp->next;
    bp->next = (memblock_t *) malloc(sizeof(memblock_t));
    bp = bp->next;
  }
  if (bp == NULL) {
    perror("diskmalloc: malloc for memblock failed");
    munmap(buf, size);
    goto abort;
  }
  bp->addr = buf;
  bp->size = size;
  bp->fildes = fildes;
  bp->name = name;
  bp->next = NULL;

  return buf;

 abort:
  close(fildes);
  unlink(name);
  free(name);
  return NULL;
}


void diskfree(void *p)
{
  memblock_t *bp, *prev;

  /* Search through info list to find block with addr == p */
  prev = NULL;
  bp = head;
  while (bp != NULL) {
    if (bp->addr == p) {
      /* Unmap file and delete it from disk */
      munmap(bp->addr, bp->size);
      close(bp->fildes);
      unlink(bp->name);
      free(bp->name);
      if (bp == head)
	head = bp->next;
      else if (prev != NULL)
	prev->next = bp->next;
      free(bp);
      return;
    }
    prev = bp;
    bp = bp->next;
  }
  fprintf(stderr, "Warning: Trying to free block starting at %p that was never allocated.\n", p);
}


/* FORTRAN stubs */

void *diskmalloc_(size_t *size)
{
  return diskmalloc(*size, NULL);
}

void diskfree_(void **p)
{
  diskfree(*p);
}

void diskfree1_(void *p)
{
  diskfree(p);
}


void *safemalloc_(size_t *size)
{
  return safemalloc(*size);
}

void *safevalloc_(size_t *size)
{
  return safevalloc(*size);
}
