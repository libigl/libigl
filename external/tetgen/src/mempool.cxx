#include "../tetgen.h"
//// mempool_cxx //////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

// Initialize fast lookup tables for mesh maniplulation primitives.

int tetgenmesh::bondtbl[12][12] = {{0,},};
int tetgenmesh::enexttbl[12] = {0,};
int tetgenmesh::eprevtbl[12] = {0,};
int tetgenmesh::enextesymtbl[12] = {0,};
int tetgenmesh::eprevesymtbl[12] = {0,};
int tetgenmesh::eorgoppotbl[12] = {0,};
int tetgenmesh::edestoppotbl[12] = {0,};
int tetgenmesh::fsymtbl[12][12] = {{0,},};
int tetgenmesh::facepivot1[12] = {0,};
int tetgenmesh::facepivot2[12][12] = {{0,},};
int tetgenmesh::tsbondtbl[12][6] = {{0,},};
int tetgenmesh::stbondtbl[12][6] = {{0,},};
int tetgenmesh::tspivottbl[12][6] = {{0,},};
int tetgenmesh::stpivottbl[12][6] = {{0,},};

// Table 'esymtbl' takes an directed edge (version) as input, returns the
//   inversed edge (version) of it.

int tetgenmesh::esymtbl[12] = {9, 6, 11, 4, 3, 7, 1, 5, 10, 0, 8, 2};

// The following four tables give the 12 permutations of the set {0,1,2,3}.
//   An offset 4 is added to each element for a direct access of the points
//   in the tetrahedron data structure.

int tetgenmesh:: orgpivot[12] = {7, 7, 5, 5, 6, 4, 4, 6, 5, 6, 7, 4};
int tetgenmesh::destpivot[12] = {6, 4, 4, 6, 5, 6, 7, 4, 7, 7, 5, 5};
int tetgenmesh::apexpivot[12] = {5, 6, 7, 4, 7, 7, 5, 5, 6, 4, 4, 6};
int tetgenmesh::oppopivot[12] = {4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7};

// The twelve versions correspond to six undirected edges. The following two
//   tables map a version to an undirected edge and vice versa.

int tetgenmesh::ver2edge[12] = {0, 1, 2, 3, 3, 5, 1, 5, 4, 0, 4, 2};
int tetgenmesh::edge2ver[ 6] = {0, 1, 2, 3, 8, 5};

// Edge versions whose apex or opposite may be dummypoint.

int tetgenmesh::epivot[12] = {4, 5, 2, 11, 4, 5, 2, 11, 4, 5, 2, 11};


// Table 'snextpivot' takes an edge version as input, returns the next edge
//   version in the same edge ring.

int tetgenmesh::snextpivot[6] = {2, 5, 4, 1, 0, 3};

// The following three tables give the 6 permutations of the set {0,1,2}.
//   An offset 3 is added to each element for a direct access of the points
//   in the triangle data structure.

int tetgenmesh::sorgpivot [6] = {3, 4, 4, 5, 5, 3};
int tetgenmesh::sdestpivot[6] = {4, 3, 5, 4, 3, 5};
int tetgenmesh::sapexpivot[6] = {5, 5, 3, 3, 4, 4};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// inittable()    Initialize the look-up tables.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::inittables()
{
  int i, j;


  // i = t1.ver; j = t2.ver;
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      bondtbl[i][j] = (j & 3) + (((i & 12) + (j & 12)) % 12);
    }
  }


  // i = t1.ver; j = t2.ver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      fsymtbl[i][j] = (j + 12 - (i & 12)) % 12;
    }
  }


  for (i = 0; i < 12; i++) {
    facepivot1[i] = (esymtbl[i] & 3);
  }

  for (i = 0; i < 12; i++) {
    for (j = 0; j < 12; j++) {
      facepivot2[i][j] = fsymtbl[esymtbl[i]][j];
    }
  }

  for (i = 0; i < 12; i++) {
    enexttbl[i] = (i + 4) % 12;
    eprevtbl[i] = (i + 8) % 12;
  }

  for (i = 0; i < 12; i++) {
    enextesymtbl[i] = esymtbl[enexttbl[i]];
    eprevesymtbl[i] = esymtbl[eprevtbl[i]];
  }

  for (i = 0; i < 12; i++) {
    eorgoppotbl [i] = eprevtbl[esymtbl[enexttbl[i]]];
    edestoppotbl[i] = enexttbl[esymtbl[eprevtbl[i]]];
  }

  int soffset, toffset;

  // i = t.ver, j = s.shver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 6; j++) {
      if ((j & 1) == 0) {
        soffset = (6 - ((i & 12) >> 1)) % 6;
        toffset = (12 - ((j & 6) << 1)) % 12;
      } else {
        soffset = (i & 12) >> 1;
        toffset = (j & 6) << 1;
      }
      tsbondtbl[i][j] = (j & 1) + (((j & 6) + soffset) % 6);
      stbondtbl[i][j] = (i & 3) + (((i & 12) + toffset) % 12);
    }
  }


  // i = t.ver, j = s.shver
  for (i = 0; i < 12; i++) {
    for (j = 0; j < 6; j++) {
      if ((j & 1) == 0) {
        soffset = (i & 12) >> 1;
        toffset = (j & 6) << 1;
      } else {
        soffset = (6 - ((i & 12) >> 1)) % 6;
        toffset = (12 - ((j & 6) << 1)) % 12;
      }
      tspivottbl[i][j] = (j & 1) + (((j & 6) + soffset) % 6);
      stpivottbl[i][j] = (i & 3) + (((i & 12) + toffset) % 12);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// restart()    Deallocate all objects in this pool.                         //
//                                                                           //
// The pool returns to a fresh state, like after it was initialized, except  //
// that no memory is freed to the operating system.  Rather, the previously  //
// allocated blocks are ready to be used.                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::arraypool::restart()
{
  objects = 0l;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// poolinit()    Initialize an arraypool for allocation of objects.          //
//                                                                           //
// Before the pool may be used, it must be initialized by this procedure.    //
// After initialization, memory can be allocated and freed in this pool.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::arraypool::poolinit(int sizeofobject, int log2objperblk)
{
  // Each object must be at least one byte long.
  objectbytes = sizeofobject > 1 ? sizeofobject : 1;

  log2objectsperblock = log2objperblk;
  // Compute the number of objects in each block.
  objectsperblock = ((int) 1) << log2objectsperblock;
  objectsperblockmark = objectsperblock - 1;

  // No memory has been allocated.
  totalmemory = 0l;
  // The top array has not been allocated yet.
  toparray = (char **) NULL;
  toparraylen = 0;

  // Ready all indices to be allocated.
  restart();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// arraypool()    The constructor and destructor.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::arraypool::arraypool(int sizeofobject, int log2objperblk)
{
  poolinit(sizeofobject, log2objperblk);
}

tetgenmesh::arraypool::~arraypool()
{
  int i;

  // Has anything been allocated at all?
  if (toparray != (char **) NULL) {
    // Walk through the top array.
    for (i = 0; i < toparraylen; i++) {
      // Check every pointer; NULLs may be scattered randomly.
      if (toparray[i] != (char *) NULL) {
        // Free an allocated block.
        free((void *) toparray[i]);
      }
    }
    // Free the top array.
    free((void *) toparray);
  }

  // The top array is no longer allocated.
  toparray = (char **) NULL;
  toparraylen = 0;
  objects = 0;
  totalmemory = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getblock()    Return (and perhaps create) the block containing the object //
//               with a given index.                                         //
//                                                                           //
// This function takes care of allocating or resizing the top array if nece- //
// ssary, and of allocating the block if it hasn't yet been allocated.       //
//                                                                           //
// Return a pointer to the beginning of the block (NOT the object).          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

char* tetgenmesh::arraypool::getblock(int objectindex)
{
  char **newarray;
  char *block;
  int newsize;
  int topindex;
  int i;

  // Compute the index in the top array (upper bits).
  topindex = objectindex >> log2objectsperblock;
  // Does the top array need to be allocated or resized?
  if (toparray == (char **) NULL) {
    // Allocate the top array big enough to hold 'topindex', and NULL out
    //   its contents.
    newsize = topindex + 128;
    toparray = (char **) malloc((size_t) (newsize * sizeof(char *)));
    toparraylen = newsize;
    for (i = 0; i < newsize; i++) {
      toparray[i] = (char *) NULL;
    }
    // Account for the memory.
    totalmemory = newsize * (uintptr_t) sizeof(char *);
  } else if (topindex >= toparraylen) {
    // Resize the top array, making sure it holds 'topindex'.
    newsize = 3 * toparraylen;
    if (topindex >= newsize) {
      newsize = topindex + 128;
    }
    // Allocate the new array, copy the contents, NULL out the rest, and
    //   free the old array.
    newarray = (char **) malloc((size_t) (newsize * sizeof(char *)));
    for (i = 0; i < toparraylen; i++) {
      newarray[i] = toparray[i];
    }
    for (i = toparraylen; i < newsize; i++) {
      newarray[i] = (char *) NULL;
    }
    free(toparray);
    // Account for the memory.
    totalmemory += (newsize - toparraylen) * sizeof(char *);
    toparray = newarray;
    toparraylen = newsize;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  block = toparray[topindex];
  if (block == (char *) NULL) {
    // Allocate a block at this index.
    block = (char *) malloc((size_t) (objectsperblock * objectbytes));
    toparray[topindex] = block;
    // Account for the memory.
    totalmemory += objectsperblock * objectbytes;
  }

  // Return a pointer to the block.
  return block;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lookup()    Return the pointer to the object with a given index, or NULL  //
//             if the object's block doesn't exist yet.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::arraypool::lookup(int objectindex)
{
  char *block;
  int topindex;

  // Has the top array been allocated yet?
  if (toparray == (char **) NULL) {
    return (void *) NULL;
  }

  // Compute the index in the top array (upper bits).
  topindex = objectindex >> log2objectsperblock;
  // Does the top index fit in the top array?
  if (topindex >= toparraylen) {
    return (void *) NULL;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  block = toparray[topindex];
  if (block == (char *) NULL) {
    return (void *) NULL;
  }

  // Compute a pointer to the object with the given index.  Note that
  //   'objectsperblock' is a power of two, so the & operation is a bit mask
  //   that preserves the lower bits.
  return (void *)(block + (objectindex & (objectsperblock - 1)) * objectbytes);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// newindex()    Allocate space for a fresh object from the pool.            //
//                                                                           //
// 'newptr' returns a pointer to the new object (it must not be a NULL).     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::arraypool::newindex(void **newptr)
{
  // Allocate an object at index 'firstvirgin'.
  int newindex = objects;
  *newptr = (void *) (getblock(objects) +
    (objects & (objectsperblock - 1)) * objectbytes);
  objects++;

  return newindex;
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// memorypool()   The constructors of memorypool.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::memorypool::memorypool()
{
  firstblock = nowblock = (void **) NULL;
  nextitem = (void *) NULL;
  deaditemstack = (void *) NULL;
  pathblock = (void **) NULL;
  pathitem = (void *) NULL;
  alignbytes = 0;
  itembytes = itemwords = 0;
  itemsperblock = 0;
  items = maxitems = 0l;
  unallocateditems = 0;
  pathitemsleft = 0;
}

tetgenmesh::memorypool::memorypool(int bytecount, int itemcount, int wsize, 
                                   int alignment)
{
  poolinit(bytecount, itemcount, wsize, alignment);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ~memorypool()   Free to the operating system all memory taken by a pool.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::memorypool::~memorypool()
{
  while (firstblock != (void **) NULL) {
    nowblock = (void **) *(firstblock);
    free(firstblock);
    firstblock = nowblock;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// poolinit()    Initialize a pool of memory for allocation of items.        //
//                                                                           //
// A `pool' is created whose records have size at least `bytecount'.  Items  //
// will be allocated in `itemcount'-item blocks.  Each item is assumed to be //
// a collection of words, and either pointers or floating-point values are   //
// assumed to be the "primary" word type.  (The "primary" word type is used  //
// to determine alignment of items.)  If `alignment' isn't zero, all items   //
// will be `alignment'-byte aligned in memory.  `alignment' must be either a //
// multiple or a factor of the primary word size;  powers of two are safe.   //
// `alignment' is normally used to create a few unused bits at the bottom of //
// each item's pointer, in which information may be stored.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::poolinit(int bytecount,int itemcount,int wordsize,
                                      int alignment)
{
  // Find the proper alignment, which must be at least as large as:
  //   - The parameter `alignment'.
  //   - The primary word type, to avoid unaligned accesses.
  //   - sizeof(void *), so the stack of dead items can be maintained
  //       without unaligned accesses.
  if (alignment > wordsize) {
    alignbytes = alignment;
  } else {
    alignbytes = wordsize;
  }
  if ((int) sizeof(void *) > alignbytes) {
    alignbytes = (int) sizeof(void *);
  }
  itemwords = ((bytecount + alignbytes - 1) /  alignbytes)
            * (alignbytes / wordsize);
  itembytes = itemwords * wordsize;
  itemsperblock = itemcount;

  // Allocate a block of items.  Space for `itemsperblock' items and one
  //   pointer (to point to the next block) are allocated, as well as space
  //   to ensure alignment of the items. 
  firstblock = (void **) malloc(itemsperblock * itembytes + sizeof(void *)
                                + alignbytes); 
  if (firstblock == (void **) NULL) {
    terminatetetgen(NULL, 1);
  }
  // Set the next block pointer to NULL.
  *(firstblock) = (void *) NULL;
  restart();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// restart()   Deallocate all items in this pool.                            //
//                                                                           //
// The pool is returned to its starting state, except that no memory is      //
// freed to the operating system.  Rather, the previously allocated blocks   //
// are ready to be reused.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::restart()
{
  uintptr_t alignptr;

  items = 0;
  maxitems = 0;

  // Set the currently active block.
  nowblock = firstblock;
  // Find the first item in the pool.  Increment by the size of (void *).
  alignptr = (uintptr_t) (nowblock + 1);
  // Align the item on an `alignbytes'-byte boundary.
  nextitem = (void *)
    (alignptr + (uintptr_t) alignbytes -
     (alignptr % (uintptr_t) alignbytes));
  // There are lots of unallocated items left in this block.
  unallocateditems = itemsperblock;
  // The stack of deallocated items is empty.
  deaditemstack = (void *) NULL;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// alloc()   Allocate space for an item.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::memorypool::alloc()
{
  void *newitem;
  void **newblock;
  uintptr_t alignptr;

  // First check the linked list of dead items.  If the list is not 
  //   empty, allocate an item from the list rather than a fresh one.
  if (deaditemstack != (void *) NULL) {
    newitem = deaditemstack;                     // Take first item in list.
    deaditemstack = * (void **) deaditemstack;
  } else {
    // Check if there are any free items left in the current block.
    if (unallocateditems == 0) {
      // Check if another block must be allocated.
      if (*nowblock == (void *) NULL) {
        // Allocate a new block of items, pointed to by the previous block.
        newblock = (void **) malloc(itemsperblock * itembytes + sizeof(void *) 
                                    + alignbytes);
        if (newblock == (void **) NULL) {
          terminatetetgen(NULL, 1);
        }
        *nowblock = (void *) newblock;
        // The next block pointer is NULL.
        *newblock = (void *) NULL;
      }
      // Move to the new block.
      nowblock = (void **) *nowblock;
      // Find the first item in the block.
      //   Increment by the size of (void *).
      alignptr = (uintptr_t) (nowblock + 1);
      // Align the item on an `alignbytes'-byte boundary.
      nextitem = (void *)
        (alignptr + (uintptr_t) alignbytes -
         (alignptr % (uintptr_t) alignbytes));
      // There are lots of unallocated items left in this block.
      unallocateditems = itemsperblock;
    }
    // Allocate a new item.
    newitem = nextitem;
    // Advance `nextitem' pointer to next free item in block.
    nextitem = (void *) ((uintptr_t) nextitem + itembytes);
    unallocateditems--;
    maxitems++;
  }
  items++;
  return newitem;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// dealloc()   Deallocate space for an item.                                 //
//                                                                           //
// The deallocated space is stored in a queue for later reuse.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::dealloc(void *dyingitem)
{
  // Push freshly killed item onto stack.
  *((void **) dyingitem) = deaditemstack;
  deaditemstack = dyingitem;
  items--;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// traversalinit()   Prepare to traverse the entire list of items.           //
//                                                                           //
// This routine is used in conjunction with traverse().                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorypool::traversalinit()
{
  uintptr_t alignptr;

  // Begin the traversal in the first block.
  pathblock = firstblock;
  // Find the first item in the block.  Increment by the size of (void *).
  alignptr = (uintptr_t) (pathblock + 1);
  // Align with item on an `alignbytes'-byte boundary.
  pathitem = (void *)
    (alignptr + (uintptr_t) alignbytes -
     (alignptr % (uintptr_t) alignbytes));
  // Set the number of items left in the current block.
  pathitemsleft = itemsperblock;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// traverse()   Find the next item in the list.                              //
//                                                                           //
// This routine is used in conjunction with traversalinit().  Be forewarned  //
// that this routine successively returns all items in the list, including   //
// deallocated ones on the deaditemqueue. It's up to you to figure out which //
// ones are actually dead.  It can usually be done more space-efficiently by //
// a routine that knows something about the structure of the item.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void* tetgenmesh::memorypool::traverse()
{
  void *newitem;
  uintptr_t alignptr;

  // Stop upon exhausting the list of items.
  if (pathitem == nextitem) {
    return (void *) NULL;
  }
  // Check whether any untraversed items remain in the current block.
  if (pathitemsleft == 0) {
    // Find the next block.
    pathblock = (void **) *pathblock;
    // Find the first item in the block.  Increment by the size of (void *).
    alignptr = (uintptr_t) (pathblock + 1);
    // Align with item on an `alignbytes'-byte boundary.
    pathitem = (void *)
      (alignptr + (uintptr_t) alignbytes -
       (alignptr % (uintptr_t) alignbytes));
    // Set the number of items left in the current block.
    pathitemsleft = itemsperblock;
  }
  newitem = pathitem;
  // Find the next item in the block.
  pathitem = (void *) ((uintptr_t) pathitem + itembytes);
  pathitemsleft--;
  return newitem;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makeindex2pointmap()    Create a map from index to vertices.              //
//                                                                           //
// 'idx2verlist' returns the created map.  Traverse all vertices, a pointer  //
// to each vertex is set into the array.  The pointer to the first vertex is //
// saved in 'idx2verlist[in->firstnumber]'.                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makeindex2pointmap(point*& idx2verlist)
{
  point pointloop;
  int idx;

  if (b->verbose > 1) {
    printf("  Constructing mapping from indices to points.\n");
  }

  idx2verlist = new point[points->items + 1];

  points->traversalinit();
  pointloop = pointtraverse();
  idx =  in->firstnumber;
  while (pointloop != (point) NULL) {
    idx2verlist[idx++] = pointloop;
    pointloop = pointtraverse();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makesubfacemap()    Create a map from vertex to subfaces incident at it.  //
//                                                                           //
// The map is returned in two arrays 'idx2faclist' and 'facperverlist'.  All //
// subfaces incident at i-th vertex (i is counted from 0) are found in the   //
// array facperverlist[j], where idx2faclist[i] <= j < idx2faclist[i + 1].   //
// Each entry in facperverlist[j] is a subface whose origin is the vertex.   //
//                                                                           //
// NOTE: These two arrays will be created inside this routine, don't forget  //
// to free them after using.                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makepoint2submap(memorypool* pool, int*& idx2faclist,
                                  face*& facperverlist)
{
  face shloop;
  int i, j, k;

  if (b->verbose > 1) {
    printf("  Making a map from points to subfaces.\n");
  }

  // Initialize 'idx2faclist'.
  idx2faclist = new int[points->items + 1];
  for (i = 0; i < points->items + 1; i++) idx2faclist[i] = 0;

  // Loop all subfaces, counter the number of subfaces incident at a vertex.
  pool->traversalinit();
  shloop.sh = shellfacetraverse(pool);
  while (shloop.sh != (shellface *) NULL) {
    // Increment the number of incident subfaces for each vertex.
    j = pointmark((point) shloop.sh[3]) - in->firstnumber;
    idx2faclist[j]++;
    j = pointmark((point) shloop.sh[4]) - in->firstnumber;
    idx2faclist[j]++;
    // Skip the third corner if it is a segment.
    if (shloop.sh[5] != NULL) {
      j = pointmark((point) shloop.sh[5]) - in->firstnumber;
      idx2faclist[j]++;
    }
    shloop.sh = shellfacetraverse(pool);
  }

  // Calculate the total length of array 'facperverlist'.
  j = idx2faclist[0];
  idx2faclist[0] = 0;  // Array starts from 0 element.
  for (i = 0; i < points->items; i++) {
    k = idx2faclist[i + 1];
    idx2faclist[i + 1] = idx2faclist[i] + j;
    j = k;
  }

  // The total length is in the last unit of idx2faclist.
  facperverlist = new face[idx2faclist[i]];

  // Loop all subfaces again, remember the subfaces at each vertex.
  pool->traversalinit();
  shloop.sh = shellfacetraverse(pool);
  while (shloop.sh != (shellface *) NULL) {
    j = pointmark((point) shloop.sh[3]) - in->firstnumber;
    shloop.shver = 0; // save the origin.
    facperverlist[idx2faclist[j]] = shloop;
    idx2faclist[j]++;
    // Is it a subface or a subsegment?
    if (shloop.sh[5] != NULL) {
      j = pointmark((point) shloop.sh[4]) - in->firstnumber;
      shloop.shver = 2; // save the origin.
      facperverlist[idx2faclist[j]] = shloop;
      idx2faclist[j]++;
      j = pointmark((point) shloop.sh[5]) - in->firstnumber;
      shloop.shver = 4; // save the origin.
      facperverlist[idx2faclist[j]] = shloop;
      idx2faclist[j]++;
    } else {
      j = pointmark((point) shloop.sh[4]) - in->firstnumber;
      shloop.shver = 1; // save the origin.
      facperverlist[idx2faclist[j]] = shloop;
      idx2faclist[j]++;
    }
    shloop.sh = shellfacetraverse(pool);
  }

  // Contents in 'idx2faclist' are shifted, now shift them back.
  for (i = points->items - 1; i >= 0; i--) {
    idx2faclist[i + 1] = idx2faclist[i];
  }
  idx2faclist[0] = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedrondealloc()    Deallocate space for a tet., marking it dead.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::tetrahedrondealloc(tetrahedron *dyingtetrahedron)
{
  // Set tetrahedron's vertices to NULL. This makes it possible to detect
  //   dead tetrahedra when traversing the list of all tetrahedra.
  dyingtetrahedron[4] = (tetrahedron) NULL;

  // Dealloc the space to subfaces/subsegments.
  if (dyingtetrahedron[8] != NULL) {
    tet2segpool->dealloc((shellface *) dyingtetrahedron[8]);
  }
  if (dyingtetrahedron[9] != NULL) {
    tet2subpool->dealloc((shellface *) dyingtetrahedron[9]);
  }

  tetrahedrons->dealloc((void *) dyingtetrahedron);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedrontraverse()    Traverse the tetrahedra, skipping dead ones.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::tetrahedron* tetgenmesh::tetrahedrontraverse()
{
  tetrahedron *newtetrahedron;

  do {
    newtetrahedron = (tetrahedron *) tetrahedrons->traverse();
    if (newtetrahedron == (tetrahedron *) NULL) {
      return (tetrahedron *) NULL;
    }
  } while ((newtetrahedron[4] == (tetrahedron) NULL) ||
           ((point) newtetrahedron[7] == dummypoint));
  return newtetrahedron;
}

tetgenmesh::tetrahedron* tetgenmesh::alltetrahedrontraverse()
{
  tetrahedron *newtetrahedron;

  do {
    newtetrahedron = (tetrahedron *) tetrahedrons->traverse();
    if (newtetrahedron == (tetrahedron *) NULL) {
      return (tetrahedron *) NULL;
    }
  } while (newtetrahedron[4] == (tetrahedron) NULL); // Skip dead ones.
  return newtetrahedron;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// shellfacedealloc()    Deallocate space for a shellface, marking it dead.  //
//                       Used both for dealloc a subface and subsegment.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::shellfacedealloc(memorypool *pool, shellface *dyingsh)
{
  // Set shellface's vertices to NULL. This makes it possible to detect dead
  //   shellfaces when traversing the list of all shellfaces.
  dyingsh[3] = (shellface) NULL;
  pool->dealloc((void *) dyingsh);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// shellfacetraverse()    Traverse the subfaces, skipping dead ones. Used    //
//                        for both subfaces and subsegments pool traverse.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::shellface* tetgenmesh::shellfacetraverse(memorypool *pool)
{
  shellface *newshellface;

  do {
    newshellface = (shellface *) pool->traverse();
    if (newshellface == (shellface *) NULL) {
      return (shellface *) NULL;
    }
  } while (newshellface[3] == (shellface) NULL);          // Skip dead ones.
  return newshellface;
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// pointdealloc()    Deallocate space for a point, marking it dead.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::pointdealloc(point dyingpoint)
{
  // Mark the point as dead. This  makes it possible to detect dead points
  //   when traversing the list of all points.
  setpointtype(dyingpoint, DEADVERTEX);
  points->dealloc((void *) dyingpoint);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// pointtraverse()    Traverse the points, skipping dead ones.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

tetgenmesh::point tetgenmesh::pointtraverse()
{
  point newpoint;

  do {
    newpoint = (point) points->traverse();
    if (newpoint == (point) NULL) {
      return (point) NULL;
    }
  } while (pointtype(newpoint) == DEADVERTEX);            // Skip dead ones.
  return newpoint;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// maketetrahedron()    Create a new tetrahedron.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::maketetrahedron(triface *newtet)
{
  newtet->tet = (tetrahedron *) tetrahedrons->alloc();

  // Initialize the four adjoining tetrahedra to be "outer space".
  newtet->tet[0] = NULL;
  newtet->tet[1] = NULL;
  newtet->tet[2] = NULL;
  newtet->tet[3] = NULL;
  // Four NULL vertices.
  newtet->tet[4] = NULL;
  newtet->tet[5] = NULL;
  newtet->tet[6] = NULL;
  newtet->tet[7] = NULL;
  // No attached segments and subfaces yet.
  newtet->tet[8] = NULL; 
  newtet->tet[9] = NULL; 
  // Initialize the marker (clear all flags).
  setelemmarker(newtet->tet, 0);
  for (int i = 0; i < numelemattrib; i++) {
    setelemattribute(newtet->tet, i, 0.0);
  }
  if (b->varvolume) {
    setvolumebound(newtet->tet, -1.0);
  }

  // Initialize the version to be Zero.
  newtet->ver = 11;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makeshellface()    Create a new shellface with version zero. Used for     //
//                    both subfaces and subsegments.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makeshellface(memorypool *pool, face *newface)
{
  newface->sh = (shellface *) pool->alloc();

  // No adjointing subfaces.
  newface->sh[0] = NULL;
  newface->sh[1] = NULL;
  newface->sh[2] = NULL;
  // Three NULL vertices.
  newface->sh[3] = NULL;
  newface->sh[4] = NULL;
  newface->sh[5] = NULL;
  // No adjoining subsegments.
  newface->sh[6] = NULL;
  newface->sh[7] = NULL;
  newface->sh[8] = NULL;
  // No adjoining tetrahedra.
  newface->sh[9] = NULL;
  newface->sh[10] = NULL;
  if (checkconstraints) {
    // Initialize the maximum area bound.
    setareabound(*newface, 0.0);
  }
  // Clear the infection and marktest bits.
  ((int *) (newface->sh))[shmarkindex + 1] = 0;
  if (useinsertradius) {
    setfacetindex(*newface, 0);
  }
  // Set the boundary marker to zero.
  setshellmark(*newface, 0);

  newface->shver = 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makepoint()    Create a new point.                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makepoint(point* pnewpoint, enum verttype vtype)
{
  int i;

  *pnewpoint = (point) points->alloc();

  // Initialize the point attributes.
  for (i = 0; i < numpointattrib; i++) {
    (*pnewpoint)[3 + i] = 0.0;
  }
  // Initialize the metric tensor.
  for (i = 0; i < sizeoftensor; i++) {
    (*pnewpoint)[pointmtrindex + i] = 0.0;
  }
  setpoint2tet(*pnewpoint, NULL);
  setpoint2ppt(*pnewpoint, NULL);
  if (b->plc || b->refine) {
    // Initialize the point-to-simplex field.
    setpoint2sh(*pnewpoint, NULL);
    if (b->metric && (bgm != NULL)) {
      setpoint2bgmtet(*pnewpoint, NULL);
    }
  }
  // Initialize the point marker (starting from in->firstnumber).
  setpointmark(*pnewpoint, (int) (points->items) - (!in->firstnumber));
  // Clear all flags.
  ((int *) (*pnewpoint))[pointmarkindex + 1] = 0;
  // Initialize (set) the point type. 
  setpointtype(*pnewpoint, vtype);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initializepools()    Calculate the sizes of the point, tetrahedron, and   //
//                      subface. Initialize their memory pools.              //
//                                                                           //
// This routine also computes the indices 'pointmarkindex', 'point2simindex',//
// 'point2pbcptindex', 'elemattribindex', and 'volumeboundindex'.  They are  //
// used to find values within each point and tetrahedron, respectively.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::initializepools()
{
  int pointsize = 0, elesize = 0, shsize = 0;
  int i;

  if (b->verbose) {
    printf("  Initializing memorypools.\n");
    printf("  tetrahedron per block: %d.\n", b->tetrahedraperblock);
  }

  inittables();

  // There are three input point lists available, which are in, addin,
  //   and bgm->in. These point lists may have different number of 
  //   attributes. Decide the maximum number.
  numpointattrib = in->numberofpointattributes;
  if (bgm != NULL) {
    if (bgm->in->numberofpointattributes > numpointattrib) {
      numpointattrib = bgm->in->numberofpointattributes;
    }
  }
  if (addin != NULL) {
    if (addin->numberofpointattributes > numpointattrib) {
      numpointattrib = addin->numberofpointattributes;
    }
  }
  if (b->weighted || b->flipinsert) { // -w or -L.
    // The internal number of point attribute needs to be at least 1
    //   (for storing point weights).
    if (numpointattrib == 0) {    
      numpointattrib = 1;
    }
  }

  // Default varconstraint = 0;
  if (in->segmentconstraintlist || in->facetconstraintlist) {
    checkconstraints = 1;
  }
  if (b->plc || b->refine) {
    // Save the insertion radius for Steiner points if boundaries
    //   are allowed be split.
    if (!b->nobisect || checkconstraints) {
      useinsertradius = 1;
    }
  }

  // The index within each point at which its metric tensor is found. 
  // Each vertex has three coordinates.
  if (b->psc) {
    // '-s' option (PSC), the u,v coordinates are provided.
    pointmtrindex = 5 + numpointattrib;
    // The index within each point at which its u, v coordinates are found.
    // Comment: They are saved after the list of point attributes.
    pointparamindex = pointmtrindex - 2;
  } else {
    pointmtrindex = 3 + numpointattrib;
  }
  // For '-m' option. A tensor field is provided (*.mtr or *.b.mtr file).
  if (b->metric) {
    // Decide the size (1, 3, or 6) of the metric tensor.
    if (bgm != (tetgenmesh *) NULL) {
      // A background mesh is allocated. It may not exist though.
      sizeoftensor = (bgm->in != (tetgenio *) NULL) ? 
        bgm->in->numberofpointmtrs : in->numberofpointmtrs;
    } else {
      // No given background mesh - Itself is a background mesh.
      sizeoftensor = in->numberofpointmtrs;
    }
    // Make sure sizeoftensor is at least 1.
    sizeoftensor = (sizeoftensor > 0) ? sizeoftensor : 1;
  } else {
    // For '-q' option. Make sure to have space for saving a scalar value.
    sizeoftensor = b->quality ? 1 : 0;
  }
  if (useinsertradius) {
    // Increase a space (REAL) for saving point insertion radius, it is
    //   saved directly after the metric. 
    sizeoftensor++;
  }
  // The index within each point at which an element pointer is found, where
  //   the index is measured in pointers. Ensure the index is aligned to a
  //   sizeof(tetrahedron)-byte address.
  point2simindex = ((pointmtrindex + sizeoftensor) * sizeof(REAL)
                 + sizeof(tetrahedron) - 1) / sizeof(tetrahedron);
  if (b->plc || b->refine || b->voroout) {
    // Increase the point size by three pointers, which are:
    //   - a pointer to a tet, read by point2tet();
    //   - a pointer to a parent point, read by point2ppt()).
    //   - a pointer to a subface or segment, read by point2sh();
    if (b->metric && (bgm != (tetgenmesh *) NULL)) {
      // Increase one pointer into the background mesh, point2bgmtet().
      pointsize = (point2simindex + 4) * sizeof(tetrahedron);
    } else {
      pointsize = (point2simindex + 3) * sizeof(tetrahedron);
    }
  } else {
    // Increase the point size by two pointer, which are:
    //   - a pointer to a tet, read by point2tet();
    //   - a pointer to a parent point, read by point2ppt()). -- Used by btree.
    pointsize = (point2simindex + 2) * sizeof(tetrahedron);
  }
  // The index within each point at which the boundary marker is found,
  //   Ensure the point marker is aligned to a sizeof(int)-byte address.
  pointmarkindex = (pointsize + sizeof(int) - 1) / sizeof(int);
  // Now point size is the ints (indicated by pointmarkindex) plus:
  //   - an integer for boundary marker;
  //   - an integer for vertex type;
  //   - an integer for geometry tag (optional, -s option).
  pointsize = (pointmarkindex + 2 + (b->psc ? 1 : 0)) * sizeof(tetrahedron);

  // Initialize the pool of vertices.
  points = new memorypool(pointsize, b->vertexperblock, sizeof(REAL), 0);

  if (b->verbose) {
    printf("  Size of a point: %d bytes.\n", points->itembytes);
  }

  // Initialize the infinite vertex.
  dummypoint = (point) new char[pointsize];
  // Initialize all fields of this point.
  dummypoint[0] = 0.0;
  dummypoint[1] = 0.0;
  dummypoint[2] = 0.0;
  for (i = 0; i < numpointattrib; i++) {
    dummypoint[3 + i] = 0.0;
  }
  // Initialize the metric tensor.
  for (i = 0; i < sizeoftensor; i++) {
    dummypoint[pointmtrindex + i] = 0.0;
  }
  setpoint2tet(dummypoint, NULL);
  setpoint2ppt(dummypoint, NULL);
  if (b->plc || b->psc || b->refine) {
    // Initialize the point-to-simplex field.
    setpoint2sh(dummypoint, NULL);
    if (b->metric && (bgm != NULL)) {
      setpoint2bgmtet(dummypoint, NULL);
    }
  }
  // Initialize the point marker (starting from in->firstnumber).
  setpointmark(dummypoint, -1); // The unique marker for dummypoint.
  // Clear all flags.
  ((int *) (dummypoint))[pointmarkindex + 1] = 0;
  // Initialize (set) the point type. 
  setpointtype(dummypoint, UNUSEDVERTEX); // Does not matter.

  // The number of bytes occupied by a tetrahedron is varying by the user-
  //   specified options. The contents of the first 12 pointers are listed
  //   in the following table:
  //     [0]  |__ neighbor at f0 __|
  //     [1]  |__ neighbor at f1 __|
  //     [2]  |__ neighbor at f2 __|
  //     [3]  |__ neighbor at f3 __|
  //     [4]  |_____ vertex p0 ____|
  //     [5]  |_____ vertex p1 ____|
  //     [6]  |_____ vertex p2 ____|
  //     [7]  |_____ vertex p3 ____|
  //     [8]  |__ segments array __| (used by -p)
  //     [9]  |__ subfaces array __| (used by -p)
  //    [10]  |_____ reserved _____|
  //    [11]  |___ elem marker ____| (used as an integer)

  elesize = 12 * sizeof(tetrahedron); 

  // The index to find the element markers. An integer containing varies
  //   flags and element counter. 
  assert(sizeof(int) <= sizeof(tetrahedron));
  assert((sizeof(tetrahedron) % sizeof(int)) == 0);
  elemmarkerindex = (elesize - sizeof(tetrahedron)) / sizeof(int);

  // The actual number of element attributes. Note that if the
  //   `b->regionattrib' flag is set, an additional attribute will be added.
  numelemattrib = in->numberoftetrahedronattributes + (b->regionattrib > 0);

  // The index within each element at which its attributes are found, where
  //   the index is measured in REALs. 
  elemattribindex = (elesize + sizeof(REAL) - 1) / sizeof(REAL);
  // The index within each element at which the maximum volume bound is
  //   found, where the index is measured in REALs.
  volumeboundindex = elemattribindex + numelemattrib;
  // If element attributes or an constraint are needed, increase the number
  //   of bytes occupied by an element.
  if (b->varvolume) {
    elesize = (volumeboundindex + 1) * sizeof(REAL);
  } else if (numelemattrib > 0) {
    elesize = volumeboundindex * sizeof(REAL);
  }


  // Having determined the memory size of an element, initialize the pool.
  tetrahedrons = new memorypool(elesize, b->tetrahedraperblock, sizeof(void *),
                                16);

  if (b->verbose) {
    printf("  Size of a tetrahedron: %d (%d) bytes.\n", elesize,
           tetrahedrons->itembytes);
  }

  if (b->plc || b->refine) { // if (b->useshelles) {
    // The number of bytes occupied by a subface.  The list of pointers
    //   stored in a subface are: three to other subfaces, three to corners,
    //   three to subsegments, two to tetrahedra.
    shsize = 11 * sizeof(shellface);
    // The index within each subface at which the maximum area bound is
    //   found, where the index is measured in REALs.
    areaboundindex = (shsize + sizeof(REAL) - 1) / sizeof(REAL);
    // If -q switch is in use, increase the number of bytes occupied by
    //   a subface for saving maximum area bound.
    if (checkconstraints) { 
      shsize = (areaboundindex + 1) * sizeof(REAL);
    } else {
      shsize = areaboundindex * sizeof(REAL);
    }
    // The index within subface at which the facet marker is found. Ensure
    //   the marker is aligned to a sizeof(int)-byte address.
    shmarkindex = (shsize + sizeof(int) - 1) / sizeof(int);
    // Increase the number of bytes by two or three integers, one for facet
    //   marker, one for shellface type, and optionally one for pbc group.
    shsize = (shmarkindex + 2) * sizeof(shellface);
    if (useinsertradius) {
      // Increase the number of byte by one integer for storing facet index.
      //    set/read by setfacetindex() and getfacetindex.
      shsize = (shmarkindex + 3) * sizeof(shellface);
    }

    // Initialize the pool of subfaces. Each subface record is eight-byte
    //   aligned so it has room to store an edge version (from 0 to 5) in
    //   the least three bits.
    subfaces = new memorypool(shsize, b->shellfaceperblock, sizeof(void *), 8);

    if (b->verbose) {
      printf("  Size of a shellface: %d (%d) bytes.\n", shsize,
             subfaces->itembytes);
    }

    // Initialize the pool of subsegments. The subsegment's record is same
    //   with subface.
    subsegs = new memorypool(shsize, b->shellfaceperblock, sizeof(void *), 8);

    // Initialize the pool for tet-subseg connections.
    tet2segpool = new memorypool(6 * sizeof(shellface), b->shellfaceperblock, 
                                 sizeof(void *), 0);
    // Initialize the pool for tet-subface connections.
    tet2subpool = new memorypool(4 * sizeof(shellface), b->shellfaceperblock, 
                                 sizeof(void *), 0);

    // Initialize arraypools for segment & facet recovery.
    subsegstack = new arraypool(sizeof(face), 10);
    subfacstack = new arraypool(sizeof(face), 10);
    subvertstack = new arraypool(sizeof(point), 8);

    // Initialize arraypools for surface point insertion/deletion.
    caveshlist = new arraypool(sizeof(face), 8);
    caveshbdlist = new arraypool(sizeof(face), 8);
    cavesegshlist = new arraypool(sizeof(face), 4);

    cavetetshlist = new arraypool(sizeof(face), 8);
    cavetetseglist = new arraypool(sizeof(face), 8);
    caveencshlist = new arraypool(sizeof(face), 8);
    caveencseglist = new arraypool(sizeof(face), 8);
  }

  // Initialize the pools for flips.
  flippool = new memorypool(sizeof(badface), 1024, sizeof(void *), 0);
  unflipqueue = new arraypool(sizeof(badface), 10);

  // Initialize the arraypools for point insertion.
  cavetetlist = new arraypool(sizeof(triface), 10);
  cavebdrylist = new arraypool(sizeof(triface), 10);
  caveoldtetlist = new arraypool(sizeof(triface), 10);
  cavetetvertlist = new arraypool(sizeof(point), 10);
}

////                                                                       ////
////                                                                       ////
//// mempool_cxx //////////////////////////////////////////////////////////////

