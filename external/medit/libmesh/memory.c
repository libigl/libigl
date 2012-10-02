/* file    : memory.c
 *   C code for memory debugging, to be used in lieu
 *   of standard memory functions
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "memory.h"


typedef struct memstack {
  size_t   size;
  void    *ptr;
  int      nxt;
  char     call[30];
} Memstack;

typedef Memstack   * pMemstack;

const int  MAXMEM = 300;
pMemstack  mstack;
int        stack,cur;


int M_memLeak() {
  int   i,c=0;

  for (i=1; i<=MAXMEM; i++)
    if (mstack[i].ptr)  c++;
  return(c);
}

/* print out allocated pointers */
void M_memDump() {
  size_t  size;
  int     i,c;
  static long mega = 1024 * 1024;
  static long kilo = 1024;

  fprintf(stdout,"\n  -- MEMORY USAGE\n");
  fprintf(stdout,"  Allocated pointers\n");
  size = 0;
  c    = 0;
  for (i=1; i<=MAXMEM; i++)
    if ( mstack[i].ptr ) {
      fprintf(stdout,"   %3d  %3d Pointer %10p  size ",++c,i,mstack[i].ptr);
      if (mstack[i].size > mega)
        fprintf(stdout,"   %10d Mbytes  ",(int)(mstack[i].size/mega));
      else if (mstack[i].size > kilo)
        fprintf(stdout,"   %10d Kbytes  ",(int)(mstack[i].size/kilo));
      else 
        fprintf(stdout,"   %10d  bytes  ",(int)(mstack[i].size));
      fprintf(stdout,"(%s)\n",mstack[i].call);
      size += mstack[i].size;
    }
  fprintf(stdout,"  Memory leaks    ");
  if ( size > mega )
    fprintf(stdout,"  %10d Mbytes  %d pointers\n",(int)(size/mega),c);
  else if ( size > kilo )
    fprintf(stdout,"  %10d Kbytes  %d pointers\n",(int)(size/kilo),c);
  else if ( size )
    fprintf(stdout,"  %10d bytes   %d pointers\n",(int)size,c);
}

/* Returns allocated memory space in bytes */
size_t M_memSize() {
  size_t size;
  int    i;

  size = 0;
  for (i=1; i<=MAXMEM; i++)
    if ( mstack[i].ptr )
      size += mstack[i].size;
  return size;
}

/* Allocates space for a block of at least size bytes,
   but does not initialize the space. */
void *M_malloc(size_t size,char *call) {
  int   i;

  /* check if first call */
  if ( !mstack ) {
    mstack = (Memstack *)calloc((1+MAXMEM),sizeof(Memstack));
    assert(mstack);
    for (i=1; i<MAXMEM; i++)
      mstack[i].nxt    = i+1;
    cur   = 1;
    stack = 0;
  }

  /* store pointer, size */
  if ( stack < MAXMEM ) {
    mstack[cur].ptr  = malloc(size);
    assert(mstack[cur].ptr);
    mstack[cur].size = size;
    /* i.e. mstack[cur].call = strdup(call) */
    /* mstack[cur].call = (char*)malloc((strlen(call)+1) * sizeof(char));
    assert(mstack[cur].call); */
    strncpy(mstack[cur].call,call,19);
    i = cur;
    cur = mstack[cur].nxt;
    ++stack;
#ifdef MEMDEBUG
    fprintf(stdout,"M_malloc: allocate %p of size %10d         (%g,%d)\n",
	    mstack[cur].ptr,size,stack,cur);
#endif
    return(mstack[i].ptr);
  }
  else {
    fprintf(stderr,"M_malloc: unable to store %10Zd bytes pointer. table full\n",
	    size);
    return(0);
  }
}


/* Allocates space for an array of nelem elements, each of size 
   elsize bytes, and initializes the space to zeros.  
   Actual amount of space allocated is >=  nelem * elsize bytes. */
void *M_calloc(size_t nelem, size_t elsize,char *call) {
  int    i;

  /* check if first call */
  if ( !mstack ) {
    mstack = (Memstack *)calloc((1+MAXMEM),sizeof(Memstack));
    assert(mstack);
    for (i=1; i<MAXMEM; i++)
      mstack[i].nxt    = i+1;
    cur   = 1;
    stack = 0;
  }

  /* store pointer, size */
  if ( stack < MAXMEM ) {
    mstack[cur].ptr  = calloc(nelem,elsize);
    if ( !mstack[cur].ptr )  return(0);

    /*assert(mstack[cur].ptr);*/
    mstack[cur].size = nelem * elsize;
    /* mstack[cur].call = (char*)malloc((strlen(call)+1) * sizeof(char));
    assert(mstack[cur].call); */
    strncpy(mstack[cur].call,call,19);
    i   = cur;
    cur = mstack[cur].nxt;
    ++stack;
#ifdef MEMDEBUG
    fprintf(stdout,"M_calloc: allocate %p of size %d         (%d,%d)\n",
	    mstack[cur].ptr,nelem*elsize,stack,cur);
#endif
    return(mstack[i].ptr);
  }
  else {
    fprintf(stderr,"M_calloc: unable to allocate %10Zd bytes. table full\n",
	    nelem*elsize);
    return(0);
  }
}

/* Changes the size of the block pointed to by ptr to size bytes 
   and returns a pointer to the (possibly moved) block. Existing 
   contents are unchanged up to the lesser of the new and old sizes. */
void *M_realloc(void *ptr, size_t size,char *call) {
  int    i;

  if ( !ptr )
    return 0;

  for (i=1; i<=MAXMEM; i++) {
    if (ptr == mstack[i].ptr) {
      /* free(mstack[i].call);
      mstack[cur].call = (char*)malloc((strlen(call)+1) * sizeof(char));
      assert(mstack[cur].call); */
      strncpy(mstack[i].call,call,19);
      mstack[i].ptr = realloc(mstack[i].ptr,size);
      if (size)
	assert(mstack[i].ptr);
      mstack[i].size = size;
#ifdef MEMDEBUG
      fprintf(stdout,"M_realloc: reallocate %p of size %d       (%d)\n",
	      mstack[i].ptr,mstack[i].size,size);
#endif
      return(mstack[i].ptr);
    }
  }
#ifdef MEMDEBUG
  fprintf(stderr,"M_realloc: pointer %p not found\n",ptr);
#endif
  return(0);
}

/* Deallocates the space pointed to by ptr (a pointer to a block 
   previously allocated by malloc() and makes the space available
   for further allocation.  If ptr is NULL, no action occurs. */
void M_free(void *ptr) {
  int   i;
  
  assert(ptr);
  for (i=1; i<=MAXMEM; i++) {
    if (mstack[i].ptr && ptr == mstack[i].ptr) {
      --stack;
      free(mstack[i].ptr);
      mstack[i].ptr  = 0;
      mstack[i].size = 0;
      mstack[i].nxt  = cur;
      mstack[i].call[0]  = '\0';
      cur = i;
#ifdef MEMDEBUG
      fprintf(stdout,"M_free: deallocate %p of size %d       (%d,%d)\n",
	      ptr,mstack[i].size,stack,cur);
#endif
      return;
    }
  }
#ifdef MEMDEBUG
  fprintf(stderr,"M_free: pointer %p not found\n",ptr);
#endif
}


/* dump memory requirements */
void primem(int np) {
  int memsize;

  memsize = M_memSize();
  if ( memsize ) {
    fprintf(stdout,"\n  -- MEMORY REQUIREMENTS\n");
    if (memsize > 1024*1024)
      fprintf(stdout,"  Total size :  %10Zd Mbytes",
	      (long int)(memsize/(1024.*1024.)));
    else if (memsize > 1024)
      fprintf(stdout,"  Total size :  %10Zd Kbytes",(long int)(memsize/1024.));
    else
      fprintf(stdout,"  Total size :  %10Zd bytes ",(long int)memsize);
    fprintf(stdout,"    (i.e. %d bytes/point)\n",memsize / np);
  }
}
