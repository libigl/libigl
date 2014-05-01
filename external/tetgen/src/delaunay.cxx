#include "../tetgen.h"
//// delaunay_cxx /////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// transfernodes()    Read the vertices from the input (tetgenio).           //
//                                                                           //
// Transferring all points from input ('in->pointlist') to TetGen's 'points'.//
// All points are indexed (the first point index is 'in->firstnumber'). Each //
// point's type is initialized as UNUSEDVERTEX. The bounding box (xmax, xmin,//
// ...) and the diameter (longest) of the point set are calculated.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::transfernodes()
{
  point pointloop;
  REAL x, y, z, w;
  int coordindex;
  int attribindex;
  int mtrindex;
  int i, j;

  if (b->psc) {
    assert(in->pointparamlist != NULL);
  }

  // Read the points.
  coordindex = 0;
  attribindex = 0;
  mtrindex = 0;
  for (i = 0; i < in->numberofpoints; i++) {
    makepoint(&pointloop, UNUSEDVERTEX);
    // Read the point coordinates.
    x = pointloop[0] = in->pointlist[coordindex++];
    y = pointloop[1] = in->pointlist[coordindex++];
    z = pointloop[2] = in->pointlist[coordindex++];
    // Read the point attributes. (Including point weights.)
    for (j = 0; j < in->numberofpointattributes; j++) {
      pointloop[3 + j] = in->pointattributelist[attribindex++];
    }
    // Read the point metric tensor.
    for (j = 0; j < in->numberofpointmtrs; j++) {
      pointloop[pointmtrindex + j] = in->pointmtrlist[mtrindex++];
    }
    if (b->weighted) { // -w option
      if (in->numberofpointattributes > 0) {
        // The first point attribute is its weight.
        //w = in->pointattributelist[in->numberofpointattributes * i];
        w = pointloop[3];
      } else {
        // No given weight available. Default choose the maximum
        //   absolute value among its coordinates.        
        w = fabs(x);
        if (w < fabs(y)) w = fabs(y);
        if (w < fabs(z)) w = fabs(z);
      }
      if (b->weighted_param == 0) {
        pointloop[3] = x * x + y * y + z * z - w; // Weighted DT.
      } else { // -w1 option
        pointloop[3] = w;  // Regular tetrahedralization.
      }
    }
    // Determine the smallest and largest x, y and z coordinates.
    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmax = z;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
      zmin = (z < zmin) ? z : zmin;
      zmax = (z > zmax) ? z : zmax;
    }
    if (b->psc) {
      // Read the geometry parameters.
      setpointgeomuv(pointloop, 0, in->pointparamlist[i].uv[0]);
      setpointgeomuv(pointloop, 1, in->pointparamlist[i].uv[1]);
      setpointgeomtag(pointloop, in->pointparamlist[i].tag);
      if (in->pointparamlist[i].type == 0) {
        setpointtype(pointloop, RIDGEVERTEX);
      } else if (in->pointparamlist[i].type == 1) {
        setpointtype(pointloop, FREESEGVERTEX);
      } else if (in->pointparamlist[i].type == 2) {
        setpointtype(pointloop, FREEFACETVERTEX);
      } else if (in->pointparamlist[i].type == 3) {
        setpointtype(pointloop, FREEVOLVERTEX);
      }
    }
  }

  // 'longest' is the largest possible edge length formed by input vertices.
  x = xmax - xmin;
  y = ymax - ymin;
  z = zmax - zmin;
  longest = sqrt(x * x + y * y + z * z);
  if (longest == 0.0) {
    printf("Error:  The point set is trivial.\n");
    terminatetetgen(this, 3);
  }

  // Two identical points are distinguished by 'lengthlimit'.
  if (b->minedgelength == 0.0) {
    b->minedgelength = longest * b->epsilon;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// hilbert_init()    Initialize the Gray code permutation table.             //
//                                                                           //
// The table 'transgc' has 8 x 3 x 8 entries. It contains all possible Gray  //
// code sequences traveled by the 1st order Hilbert curve in 3 dimensions.   //
// The first column is the Gray code of the entry point of the curve, and    //
// the second column is the direction (0, 1, or 2, 0 means the x-axis) where //
// the exit point of curve lies.                                             //
//                                                                           //
// The table 'tsb1mod3' contains the numbers of trailing set '1' bits of the //
// indices from 0 to 7, modulo by '3'. The code for generating this table is //
// from: http://graphics.stanford.edu/~seander/bithacks.html.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::hilbert_init(int n)
{
  int gc[8], N, mask, travel_bit;
  int e, d, f, k, g;
  int v, c;
  int i;

  N = (n == 2) ? 4 : 8;
  mask = (n == 2) ? 3 : 7;

  // Generate the Gray code sequence.
  for (i = 0; i < N; i++) {
    gc[i] = i ^ (i >> 1);
  }

  for (e = 0; e < N; e++) {
    for (d = 0; d < n; d++) {
      // Calculate the end point (f).
      f = e ^ (1 << d);  // Toggle the d-th bit of 'e'.
      // travel_bit = 2**p, the bit we want to travel. 
      travel_bit = e ^ f;
      for (i = 0; i < N; i++) {
        // // Rotate gc[i] left by (p + 1) % n bits.
        k = gc[i] * (travel_bit * 2);
        g = ((k | (k / N)) & mask);
        // Calculate the permuted Gray code by xor with the start point (e).
        transgc[e][d][i] = (g ^ e);
      }
      assert(transgc[e][d][0] == e);
      assert(transgc[e][d][N - 1] == f);
    } // d
  } // e

  // Count the consecutive '1' bits (trailing) on the right.
  tsb1mod3[0] = 0;
  for (i = 1; i < N; i++) {
    v = ~i; // Count the 0s.
    v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
    tsb1mod3[i] = c % n;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// hilbert_sort3()    Sort points using the 3d Hilbert curve.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::hilbert_split(point* vertexarray,int arraysize,int gc0,int gc1,
                              REAL bxmin, REAL bxmax, REAL bymin, REAL bymax, 
                              REAL bzmin, REAL bzmax)
{
  point swapvert;
  int axis, d;
  REAL split;
  int i, j;


  // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which 
  //   correspoding to x-, or y- or z-axis.
  axis = (gc0 ^ gc1) >> 1; 

  // Calulate the split position along the axis.
  if (axis == 0) {
    split = 0.5 * (bxmin + bxmax);
  } else if (axis == 1) {
    split = 0.5 * (bymin + bymax);
  } else { // == 2
    split = 0.5 * (bzmin + bzmax);
  }

  // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction
  //   of the axis is to the positive of the axis, otherwise, it is -1.
  d = ((gc0 & (1<<axis)) == 0) ? 1 : -1;


  // Partition the vertices into left- and right-arrays such that left points
  //   have Hilbert indices lower than the right points.
  i = 0;
  j = arraysize - 1;

  // Partition the vertices into left- and right-arrays.
  if (d > 0) {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i][axis] >= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j][axis] < split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  } else {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i][axis] <= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j][axis] > split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }

  return i;
}

void tetgenmesh::hilbert_sort3(point* vertexarray, int arraysize, int e, int d, 
                               REAL bxmin, REAL bxmax, REAL bymin, REAL bymax, 
                               REAL bzmin, REAL bzmax, int depth)
{
  REAL x1, x2, y1, y2, z1, z2;
  int p[9], w, e_w, d_w, k, ei, di;
  int n = 3, mask = 7;

  p[0] = 0;
  p[8] = arraysize;

  // Sort the points according to the 1st order Hilbert curve in 3d.
  p[4] = hilbert_split(vertexarray, p[8], transgc[e][d][3], transgc[e][d][4], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[2] = hilbert_split(vertexarray, p[4], transgc[e][d][1], transgc[e][d][2], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[1] = hilbert_split(vertexarray, p[2], transgc[e][d][0], transgc[e][d][1], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax);
  p[3] = hilbert_split(&(vertexarray[p[2]]), p[4] - p[2], 
                       transgc[e][d][2], transgc[e][d][3], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[2];
  p[6] = hilbert_split(&(vertexarray[p[4]]), p[8] - p[4], 
                       transgc[e][d][5], transgc[e][d][6], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[5] = hilbert_split(&(vertexarray[p[4]]), p[6] - p[4], 
                       transgc[e][d][4], transgc[e][d][5], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
  p[7] = hilbert_split(&(vertexarray[p[6]]), p[8] - p[6], 
                       transgc[e][d][6], transgc[e][d][7], 
                       bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[6];

  if (b->hilbert_order > 0) {
    // A maximum order is prescribed. 
    if ((depth + 1) == b->hilbert_order) {
      // The maximum prescribed order is reached.
      return;
    }
  }

  // Recursively sort the points in sub-boxes.
  for (w = 0; w < 8; w++) {
    // w is the local Hilbert index (NOT Gray code).
    // Sort into the sub-box either there are more than 2 points in it, or
    //   the prescribed order of the curve is not reached yet.
    //if ((p[w+1] - p[w] > b->hilbert_limit) || (b->hilbert_order > 0)) {
    if ((p[w+1] - p[w]) > b->hilbert_limit) {
      // Calculcate the start point (ei) of the curve in this sub-box.
      //   update e = e ^ (e(w) left_rotate (d+1)).
      if (w == 0) {
        e_w = 0;
      } else {
        //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
        k = 2 * ((w - 1) / 2); 
        e_w = k ^ (k >> 1); // = gc(k).
      }
      k = e_w;
      e_w = ((k << (d+1)) & mask) | ((k >> (n-d-1)) & mask);
      ei = e ^ e_w;
      // Calulcate the direction (di) of the curve in this sub-box.
      //   update d = (d + d(w) + 1) % n
      if (w == 0) {
        d_w = 0;
      } else {
        d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
      }
      di = (d + d_w + 1) % n;
      // Calculate the bounding box of the sub-box.
      if (transgc[e][d][w] & 1) { // x-axis
        x1 = 0.5 * (bxmin + bxmax);
        x2 = bxmax;
      } else {
        x1 = bxmin;
        x2 = 0.5 * (bxmin + bxmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (bymin + bymax);
        y2 = bymax;
      } else {
        y1 = bymin;
        y2 = 0.5 * (bymin + bymax);
      }
      if (transgc[e][d][w] & 4) { // z-axis
        z1 = 0.5 * (bzmin + bzmax);
        z2 = bzmax;
      } else {
        z1 = bzmin;
        z2 = 0.5 * (bzmin + bzmax);
      }
      hilbert_sort3(&(vertexarray[p[w]]), p[w+1] - p[w], ei, di, 
                    x1, x2, y1, y2, z1, z2, depth+1);
    } // if (p[w+1] - p[w] > 1)
  } // w
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// brio_multiscale_sort()    Sort the points using BRIO and Hilbert curve.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::brio_multiscale_sort(point* vertexarray, int arraysize, 
                                      int threshold, REAL ratio, int *depth)
{
  int middle;

  middle = 0;
  if (arraysize >= threshold) {
    (*depth)++;
    middle = arraysize * ratio;
    brio_multiscale_sort(vertexarray, middle, threshold, ratio, depth);
  }
  // Sort the right-array (rnd-th round) using the Hilbert curve.
  hilbert_sort3(&(vertexarray[middle]), arraysize - middle, 0, 0, // e, d
                xmin, xmax, ymin, ymax, zmin, zmax, 0); // depth.
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// randomnation()    Generate a random number between 0 and 'choices' - 1.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

unsigned long tetgenmesh::randomnation(unsigned int choices)
{
  unsigned long newrandom;

  if (choices >= 714025l) {
    newrandom = (randomseed * 1366l + 150889l) % 714025l;
    randomseed = (newrandom * 1366l + 150889l) % 714025l;
    newrandom = newrandom * (choices / 714025l) + randomseed;
    if (newrandom >= choices) {
      return newrandom - choices;
    } else {
      return newrandom;
    }
  } else {
    randomseed = (randomseed * 1366l + 150889l) % 714025l;
    return randomseed % choices;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// randomsample()    Randomly sample the tetrahedra for point loation.       //
//                                                                           //
// Searching begins from one of handles:  the input 'searchtet', a recently  //
// encountered tetrahedron 'recenttet',  or from one chosen from a random    //
// sample.  The choice is made by determining which one's origin is closest  //
// to the point we are searching for.                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::randomsample(point searchpt,triface *searchtet)
{
  tetrahedron *firsttet, *tetptr;
  point torg;
  void **sampleblock;
  uintptr_t alignptr;
  long sampleblocks, samplesperblock, samplenum;
  long tetblocks, i, j;
  REAL searchdist, dist;

  if (b->verbose > 2) {
    printf("      Random sampling tetrahedra for searching point %d.\n",
           pointmark(searchpt));
  }

  if (!nonconvex) {
    if (searchtet->tet == NULL) {
      // A null tet. Choose the recenttet as the starting tet.
      *searchtet = recenttet;
      // Recenttet should not be dead.
      assert(recenttet.tet[4] != NULL);
    }

    // 'searchtet' should be a valid tetrahedron. Choose the base face
    //   whose vertices must not be 'dummypoint'.
    searchtet->ver = 3;
    // Record the distance from its origin to the searching point.
    torg = org(*searchtet);
    searchdist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
                 (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
                 (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);

    // If a recently encountered tetrahedron has been recorded and has not
    //   been deallocated, test it as a good starting point.
    if (recenttet.tet != searchtet->tet) {
      recenttet.ver = 3;
      torg = org(recenttet);
      dist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
             (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
             (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);
      if (dist < searchdist) {
        *searchtet = recenttet;
        searchdist = dist;
      }
    }
  } else {
    // The mesh is non-convex. Do not use 'recenttet'.
    assert(samples >= 1l); // Make sure at least 1 sample.
    searchdist = longest;
  }

  // Select "good" candidate using k random samples, taking the closest one.
  //   The number of random samples taken is proportional to the fourth root
  //   of the number of tetrahedra in the mesh. 
  while (samples * samples * samples * samples < tetrahedrons->items) {
    samples++;
  }
  // Find how much blocks in current tet pool.
  tetblocks = (tetrahedrons->maxitems + b->tetrahedraperblock - 1) 
            / b->tetrahedraperblock;
  // Find the average samples per block. Each block at least have 1 sample.
  samplesperblock = 1 + (samples / tetblocks);
  sampleblocks = samples / samplesperblock;
  sampleblock = tetrahedrons->firstblock;
  for (i = 0; i < sampleblocks; i++) {
    alignptr = (uintptr_t) (sampleblock + 1);
    firsttet = (tetrahedron *)
               (alignptr + (uintptr_t) tetrahedrons->alignbytes
               - (alignptr % (uintptr_t) tetrahedrons->alignbytes));
    for (j = 0; j < samplesperblock; j++) {
      if (i == tetblocks - 1) {
        // This is the last block.
        samplenum = randomnation((int)
                      (tetrahedrons->maxitems - (i * b->tetrahedraperblock)));
      } else {
        samplenum = randomnation(b->tetrahedraperblock);
      }
      tetptr = (tetrahedron *)
               (firsttet + (samplenum * tetrahedrons->itemwords));
      torg = (point) tetptr[4];
      if (torg != (point) NULL) {
        dist = (searchpt[0] - torg[0]) * (searchpt[0] - torg[0]) +
               (searchpt[1] - torg[1]) * (searchpt[1] - torg[1]) +
               (searchpt[2] - torg[2]) * (searchpt[2] - torg[2]);
        if (dist < searchdist) {
          searchtet->tet = tetptr;
          searchtet->ver = 11; // torg = org(t);
          searchdist = dist;
        }
      } else {
        // A dead tet. Re-sample it.
        if (i != tetblocks - 1) j--;
      }
    }
    sampleblock = (void **) *sampleblock;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// locate()    Find a tetrahedron containing a given point.                  //
//                                                                           //
// Begins its search from 'searchtet', assume there is a line segment L from //
// a vertex of 'searchtet' to the query point 'searchpt', and simply walk    //
// towards 'searchpt' by traversing all faces intersected by L.              //
//                                                                           //
// On completion, 'searchtet' is a tetrahedron that contains 'searchpt'. The //
// returned value indicates one of the following cases:                      //
//   - ONVERTEX, the search point lies on the origin of 'searchtet'.         //
//   - ONEDGE, the search point lies on an edge of 'searchtet'.              //
//   - ONFACE, the search point lies on a face of 'searchtet'.               //
//   - INTET, the search point lies in the interior of 'searchtet'.          //
//   - OUTSIDE, the search point lies outside the mesh. 'searchtet' is a     //
//              hull face which is visible by the search point.              //
//                                                                           //
// WARNING: This routine is designed for convex triangulations, and will not //
// generally work after the holes and concavities have been carved.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::locateresult tetgenmesh::locate(point searchpt, 
                                                 triface* searchtet)
{
  point torg, tdest, tapex, toppo;
  enum {ORGMOVE, DESTMOVE, APEXMOVE} nextmove;
  REAL ori, oriorg, oridest, oriapex;
  enum locateresult loc = OUTSIDE;
  int t1ver;
  int s;

  if (searchtet->tet == NULL) {
    // A null tet. Choose the recenttet as the starting tet.
    searchtet->tet = recenttet.tet;
  }

  // Check if we are in the outside of the convex hull.
  if (ishulltet(*searchtet)) {
    // Get its adjacent tet (inside the hull).
    searchtet->ver = 3;
    fsymself(*searchtet);
  }

  // Let searchtet be the face such that 'searchpt' lies above to it.
  for (searchtet->ver = 0; searchtet->ver < 4; searchtet->ver++) {
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);
    ori = orient3d(torg, tdest, tapex, searchpt); 
    if (ori < 0.0) break;
  }
  assert(searchtet->ver != 4);

  // Walk through tetrahedra to locate the point.
  while (true) {

    toppo = oppo(*searchtet);
    
    // Check if the vertex is we seek.
    if (toppo == searchpt) {
      // Adjust the origin of searchtet to be searchpt.
      esymself(*searchtet);
      eprevself(*searchtet);
      loc = ONVERTEX; // return ONVERTEX;
      break;
    }

    // We enter from one of serarchtet's faces, which face do we exit?
    oriorg = orient3d(tdest, tapex, toppo, searchpt); 
    oridest = orient3d(tapex, torg, toppo, searchpt);
    oriapex = orient3d(torg, tdest, toppo, searchpt);

    // Now decide which face to move. It is possible there are more than one
    //   faces are viable moves. If so, randomly choose one.
    if (oriorg < 0) {
      if (oridest < 0) {
        if (oriapex < 0) {
          // All three faces are possible.
          s = randomnation(3); // 's' is in {0,1,2}.
          if (s == 0) {
            nextmove = ORGMOVE;
          } else if (s == 1) {
            nextmove = DESTMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Two faces, opposite to origin and destination, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = ORGMOVE;
          } else {
            nextmove = DESTMOVE;
          }
        }
      } else {
        if (oriapex < 0) {
          // Two faces, opposite to origin and apex, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = ORGMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Only the face opposite to origin is viable.
          nextmove = ORGMOVE;
        }
      }
    } else {
      if (oridest < 0) {
        if (oriapex < 0) {
          // Two faces, opposite to destination and apex, are viable.
          //s = randomnation(2); // 's' is in {0,1}.
          if (randomnation(2)) {
            nextmove = DESTMOVE;
          } else {
            nextmove = APEXMOVE;
          }
        } else {
          // Only the face opposite to destination is viable.
          nextmove = DESTMOVE;
        }
      } else {
        if (oriapex < 0) {
          // Only the face opposite to apex is viable.
          nextmove = APEXMOVE;
        } else {
          // The point we seek must be on the boundary of or inside this
          //   tetrahedron. Check for boundary cases.
          if (oriorg == 0) {
            // Go to the face opposite to origin.
            enextesymself(*searchtet);
            if (oridest == 0) {
              eprevself(*searchtet); // edge oppo->apex
              if (oriapex == 0) {
                // oppo is duplicated with p.
                loc = ONVERTEX; // return ONVERTEX;
                break;
              }
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            if (oriapex == 0) {
              enextself(*searchtet); // edge dest->oppo
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oridest == 0) {
            // Go to the face opposite to destination.
            eprevesymself(*searchtet);
            if (oriapex == 0) {
              eprevself(*searchtet); // edge oppo->org
              loc = ONEDGE; // return ONEDGE;
              break;
            }
            loc = ONFACE; // return ONFACE;
            break;
          }
          if (oriapex == 0) {
            // Go to the face opposite to apex
            esymself(*searchtet);
            loc = ONFACE; // return ONFACE;
            break;
          }
          loc = INTETRAHEDRON; // return INTETRAHEDRON;
          break;
        }
      }
    }
    
    // Move to the selected face.
    if (nextmove == ORGMOVE) {
      enextesymself(*searchtet);
    } else if (nextmove == DESTMOVE) {
      eprevesymself(*searchtet);
    } else {
      esymself(*searchtet);
    }
    // Move to the adjacent tetrahedron (maybe a hull tetrahedron).
    fsymself(*searchtet);
    if (oppo(*searchtet) == dummypoint) {
      loc = OUTSIDE; // return OUTSIDE;
      break;
    }

    // Retreat the three vertices of the base face.
    torg = org(*searchtet);
    tdest = dest(*searchtet);
    tapex = apex(*searchtet);

  } // while (true)

  return loc;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flippush()    Push a face (possibly will be flipped) into flipstack.      //
//                                                                           //
// The face is marked. The flag is used to check the validity of the face on //
// its popup.  Some other flips may change it already.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flippush(badface*& fstack, triface* flipface)
{
  if (!facemarked(*flipface)) {
    badface *newflipface = (badface *) flippool->alloc();
    newflipface->tt = *flipface;
    markface(newflipface->tt);
    // Push this face into stack.
    newflipface->nextitem = fstack;
    fstack = newflipface;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// incrementalflip()    Incrementally flipping to construct DT.              //
//                                                                           //
// Faces need to be checked for flipping are already queued in 'flipstack'.  //
// Return the total number of performed flips.                               //
//                                                                           //
// Comment:  This routine should be only used in the incremental Delaunay    //
// construction.  In other cases, lawsonflip3d() should be used.             // 
//                                                                           //
// If the new point lies outside of the convex hull ('hullflag' is set). The //
// incremental flip algorithm still works as usual.  However, we must ensure //
// that every flip (2-to-3 or 3-to-2) does not create a duplicated (existing)//
// edge or face. Otherwise, the underlying space of the triangulation becomes//
// non-manifold and it is not possible to flip further.                      //
// Thanks to Joerg Rambau and Frank Lutz for helping in this issue.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::incrementalflip(point newpt, int hullflag, flipconstraints *fc)
{
  badface *popface;
  triface fliptets[5], *parytet;
  point *pts, *parypt, pe;
  REAL sign, ori;
  int flipcount = 0;
  int t1ver;
  int i;

  if (b->verbose > 2) {
    printf("      Lawson flip (%ld faces).\n", flippool->items);
  }

  if (hullflag) {
    // 'newpt' lies in the outside of the convex hull. 
    // Mark all hull vertices which are connecting to it.
    popface = flipstack;
    while (popface != NULL) {
      pts = (point *) popface->tt.tet;
      for (i = 4; i < 8; i++) {
        if ((pts[i] != newpt) && (pts[i] != dummypoint)) {
          if (!pinfected(pts[i])) {
            pinfect(pts[i]);
            cavetetvertlist->newindex((void **) &parypt);
            *parypt = pts[i];
          }
        } 
      }
      popface = popface->nextitem;
    }
  }

  // Loop until the queue is empty.
  while (flipstack != NULL) {

    // Pop a face from the stack.
    popface = flipstack;
    fliptets[0] = popface->tt;
    flipstack = flipstack->nextitem; // The next top item in stack.
    flippool->dealloc((void *) popface);

    // Skip it if it is a dead tet (destroyed by previous flips).
    if (isdeadtet(fliptets[0])) continue;
    // Skip it if it is not the same tet as we saved.
    if (!facemarked(fliptets[0])) continue;

    unmarkface(fliptets[0]);    

    if ((point) fliptets[0].tet[7] == dummypoint) {
      // It must be a hull edge.
      fliptets[0].ver = epivot[fliptets[0].ver];
      // A hull edge. The current convex hull may be enlarged.
      fsym(fliptets[0], fliptets[1]);
      pts = (point *) fliptets[1].tet;
      ori = orient3d(pts[4], pts[5], pts[6], newpt);
      if (ori < 0) {
        // Visible. The convex hull will be enlarged.
        // Decide which flip (2-to-3, 3-to-2, or 4-to-1) to use.
        // Check if the tet [a,c,e,d] or [c,b,e,d] exists.
        enext(fliptets[1], fliptets[2]); 
        eprev(fliptets[1], fliptets[3]); 
        fnextself(fliptets[2]); // [a,c,e,*]
        fnextself(fliptets[3]); // [c,b,e,*]
        if (oppo(fliptets[2]) == newpt) {
          if (oppo(fliptets[3]) == newpt) {
            // Both tets exist! A 4-to-1 flip is found.
            terminatetetgen(this, 2); // Report a bug.
          } else {
            esym(fliptets[2], fliptets[0]);
            fnext(fliptets[0], fliptets[1]); 
            fnext(fliptets[1], fliptets[2]); 
            // Perform a 3-to-2 flip. Replace edge [c,a] by face [d,e,b].
            // This corresponds to my standard labels, where edge [e,d] is
            //   repalced by face [a,b,c], and a is the new vertex. 
            //   [0] [c,a,d,e] (d = newpt)
            //   [1] [c,a,e,b] (c = dummypoint)
            //   [2] [c,a,b,d]
            flip32(fliptets, 1, fc);
          }
        } else {
          if (oppo(fliptets[3]) == newpt) {
            fnext(fliptets[3], fliptets[0]);
            fnext(fliptets[0], fliptets[1]); 
            fnext(fliptets[1], fliptets[2]); 
            // Perform a 3-to-2 flip. Replace edge [c,b] by face [d,a,e].
            //   [0] [c,b,d,a] (d = newpt)
            //   [1] [c,b,a,e] (c = dummypoint)
            //   [2] [c,b,e,d]
            flip32(fliptets, 1, fc);
          } else {
            if (hullflag) {
              // Reject this flip if pe is already marked.
              pe = oppo(fliptets[1]);
              if (!pinfected(pe)) {
                pinfect(pe);
                cavetetvertlist->newindex((void **) &parypt);
                *parypt = pe;
                // Perform a 2-to-3 flip.
                flip23(fliptets, 1, fc);
              } else {
                // Reject this flip.
                flipcount--;
              }
            } else {
              // Perform a 2-to-3 flip. Replace face [a,b,c] by edge [e,d].
              //   [0] [a,b,c,d], d = newpt.
              //   [1] [b,a,c,e], c = dummypoint.
              flip23(fliptets, 1, fc);
            }
          }
        }
        flipcount++;
      } 
      continue;
    } // if (dummypoint)

    fsym(fliptets[0], fliptets[1]);
    if ((point) fliptets[1].tet[7] == dummypoint) {
      // A hull face is locally Delaunay.
      continue;
    }
    // Check if the adjacent tet has already been tested.
    if (marktested(fliptets[1])) {
      // It has been tested and it is Delaunay.
      continue;
    }

    // Test whether the face is locally Delaunay or not.
    pts = (point *) fliptets[1].tet; 
    if (b->weighted) {
      sign = orient4d_s(pts[4], pts[5], pts[6], pts[7], newpt,
                        pts[4][3], pts[5][3], pts[6][3], pts[7][3],
                        newpt[3]);
    } else {
      sign = insphere_s(pts[4], pts[5], pts[6], pts[7], newpt);
    }


    if (sign < 0) {
      point pd = newpt;
      point pe = oppo(fliptets[1]);
      // Check the convexity of its three edges. Stop checking either a
      //   locally non-convex edge (ori < 0) or a flat edge (ori = 0) is
      //   encountered, and 'fliptet' represents that edge.
      for (i = 0; i < 3; i++) {
        ori = orient3d(org(fliptets[0]), dest(fliptets[0]), pd, pe);
        if (ori <= 0) break;
        enextself(fliptets[0]);
      }
      if (ori > 0) {
        // A 2-to-3 flip is found.
        //   [0] [a,b,c,d], 
        //   [1] [b,a,c,e]. no dummypoint.
        flip23(fliptets, 0, fc);
        flipcount++;
      } else { // ori <= 0
        // The edge ('fliptets[0]' = [a',b',c',d]) is non-convex or flat,
        //   where the edge [a',b'] is one of [a,b], [b,c], and [c,a].
        // Check if there are three or four tets sharing at this edge.        
        esymself(fliptets[0]); // [b,a,d,c]
        for (i = 0; i < 3; i++) {
          fnext(fliptets[i], fliptets[i+1]);
        }
        if (fliptets[3].tet == fliptets[0].tet) {
          // A 3-to-2 flip is found. (No hull tet.)
          flip32(fliptets, 0, fc); 
          flipcount++;
        } else {
          // There are more than 3 tets at this edge.
          fnext(fliptets[3], fliptets[4]);
          if (fliptets[4].tet == fliptets[0].tet) {
            if (ori == 0) {
              // A 4-to-4 flip is found. (Two hull tets may be involved.)
              // Current tets in 'fliptets':
              //   [0] [b,a,d,c] (d may be newpt)
              //   [1] [b,a,c,e]
              //   [2] [b,a,e,f] (f may be dummypoint)
              //   [3] [b,a,f,d]
              esymself(fliptets[0]); // [a,b,c,d] 
              // A 2-to-3 flip replaces face [a,b,c] by edge [e,d].
              //   This creates a degenerate tet [e,d,a,b] (tmpfliptets[0]).
              //   It will be removed by the followed 3-to-2 flip.
              flip23(fliptets, 0, fc); // No hull tet.
              fnext(fliptets[3], fliptets[1]);
              fnext(fliptets[1], fliptets[2]);
              // Current tets in 'fliptets':
              //   [0] [...]
              //   [1] [b,a,d,e] (degenerated, d may be new point).
              //   [2] [b,a,e,f] (f may be dummypoint)
              //   [3] [b,a,f,d]
              // A 3-to-2 flip replaces edge [b,a] by face [d,e,f].
              //   Hull tets may be involved (f may be dummypoint).
              flip32(&(fliptets[1]), (apex(fliptets[3]) == dummypoint), fc);
              flipcount++;
            }
          }
        }
      } // ori
    } else {
      // The adjacent tet is Delaunay. Mark it to avoid testing it again.
      marktest(fliptets[1]);
      // Save it for unmarking it later.
      cavebdrylist->newindex((void **) &parytet);
      *parytet = fliptets[1];
    }

  } // while (flipstack)

  // Unmark saved tetrahedra.
  for (i = 0; i < cavebdrylist->objects; i++) {
    parytet = (triface *) fastlookup(cavebdrylist, i);
    unmarktest(*parytet);
  }
  cavebdrylist->restart();

  if (hullflag) {
    // Unmark infected vertices.
    for (i = 0; i < cavetetvertlist->objects; i++) {
      parypt = (point *) fastlookup(cavetetvertlist, i);
      puninfect(*parypt);
    }
    cavetetvertlist->restart();
  }


  return flipcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// initialdelaunay()    Create an initial Delaunay tetrahedralization.       //
//                                                                           //
// The tetrahedralization contains only one tetrahedron abcd, and four hull  //
// tetrahedra. The points pa, pb, pc, and pd must be linearly independent.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::initialdelaunay(point pa, point pb, point pc, point pd)
{
  triface firsttet, tetopa, tetopb, tetopc, tetopd;
  triface worktet, worktet1;

  if (b->verbose > 2) {
    printf("      Create init tet (%d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd));
  }

  // Create the first tetrahedron.
  maketetrahedron(&firsttet);
  setvertices(firsttet, pa, pb, pc, pd);
  // Create four hull tetrahedra.
  maketetrahedron(&tetopa);
  setvertices(tetopa, pb, pc, pd, dummypoint);
  maketetrahedron(&tetopb);
  setvertices(tetopb, pc, pa, pd, dummypoint);
  maketetrahedron(&tetopc);
  setvertices(tetopc, pa, pb, pd, dummypoint);
  maketetrahedron(&tetopd);
  setvertices(tetopd, pb, pa, pc, dummypoint);
  hullsize += 4;

  // Connect hull tetrahedra to firsttet (at four faces of firsttet).
  bond(firsttet, tetopd);
  esym(firsttet, worktet);
  bond(worktet, tetopc); // ab
  enextesym(firsttet, worktet);
  bond(worktet, tetopa); // bc 
  eprevesym(firsttet, worktet);
  bond(worktet, tetopb); // ca

  // Connect hull tetrahedra together (at six edges of firsttet).
  esym(tetopc, worktet); 
  esym(tetopd, worktet1);
  bond(worktet, worktet1); // ab
  esym(tetopa, worktet);
  eprevesym(tetopd, worktet1);
  bond(worktet, worktet1); // bc
  esym(tetopb, worktet);
  enextesym(tetopd, worktet1);
  bond(worktet, worktet1); // ca
  eprevesym(tetopc, worktet);
  enextesym(tetopb, worktet1);
  bond(worktet, worktet1); // da
  eprevesym(tetopa, worktet);
  enextesym(tetopc, worktet1);
  bond(worktet, worktet1); // db
  eprevesym(tetopb, worktet);
  enextesym(tetopa, worktet1);
  bond(worktet, worktet1); // dc

  // Set the vertex type.
  if (pointtype(pa) == UNUSEDVERTEX) {
    setpointtype(pa, VOLVERTEX);
  }
  if (pointtype(pb) == UNUSEDVERTEX) {
    setpointtype(pb, VOLVERTEX);
  }
  if (pointtype(pc) == UNUSEDVERTEX) {
    setpointtype(pc, VOLVERTEX);
  }
  if (pointtype(pd) == UNUSEDVERTEX) {
    setpointtype(pd, VOLVERTEX);
  }

  setpoint2tet(pa, encode(firsttet));
  setpoint2tet(pb, encode(firsttet));
  setpoint2tet(pc, encode(firsttet));
  setpoint2tet(pd, encode(firsttet));

  // Remember the first tetrahedron.
  recenttet = firsttet;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// incrementaldelaunay()    Create a Delaunay tetrahedralization by          //
//                          the incremental approach.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


void tetgenmesh::incrementaldelaunay(clock_t& tv)
{
  triface searchtet;
  point *permutarray, swapvertex;
  REAL v1[3], v2[3], n[3];
  REAL bboxsize, bboxsize2, bboxsize3, ori;
  int randindex; 
  int ngroup = 0;
  int i, j;

  if (!b->quiet) {
    printf("Delaunizing vertices...\n");
  }

  // Form a random permuation (uniformly at random) of the set of vertices.
  permutarray = new point[in->numberofpoints];
  points->traversalinit();

  if (b->no_sort) {
    if (b->verbose) {
      printf("  Using the input order.\n"); 
    }
    for (i = 0; i < in->numberofpoints; i++) {
      permutarray[i] = (point) points->traverse();
    }
  } else {
    if (b->verbose) {
      printf("  Permuting vertices.\n"); 
    }
    srand(in->numberofpoints);
    for (i = 0; i < in->numberofpoints; i++) {
      randindex = rand() % (i + 1); // randomnation(i + 1);
      permutarray[i] = permutarray[randindex];
      permutarray[randindex] = (point) points->traverse();
    }
    if (b->brio_hilbert) { // -b option
      if (b->verbose) {
        printf("  Sorting vertices.\n"); 
      }
      hilbert_init(in->mesh_dim);
      brio_multiscale_sort(permutarray, in->numberofpoints, b->brio_threshold, 
                           b->brio_ratio, &ngroup);
    }
  }

  tv = clock(); // Remember the time for sorting points.

  // Calculate the diagonal size of its bounding box.
  bboxsize = sqrt(norm2(xmax - xmin, ymax - ymin, zmax - zmin));
  bboxsize2 = bboxsize * bboxsize;
  bboxsize3 = bboxsize2 * bboxsize;

  // Make sure the second vertex is not identical with the first one.
  i = 1;
  while ((distance(permutarray[0],permutarray[i])/bboxsize)<b->epsilon) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) identical (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(this, 10);
    }
  }
  if (i > 1) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Make sure the third vertex is not collinear with the first two.
  // Acknowledgement:  Thanks Jan Pomplun for his correction by using 
  //   epsilon^2 and epsilon^3 (instead of epsilon). 2013-08-15.
  i = 2;
  for (j = 0; j < 3; j++) {
    v1[j] = permutarray[1][j] - permutarray[0][j];
    v2[j] = permutarray[i][j] - permutarray[0][j];
  }
  cross(v1, v2, n);
  while ((sqrt(norm2(n[0], n[1], n[2])) / bboxsize2) < 
         (b->epsilon * b->epsilon)) {
    i++;
    if (i == in->numberofpoints - 1) {
      printf("Exception:  All vertices are (nearly) collinear (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(this, 10);
    }
    for (j = 0; j < 3; j++) {
      v2[j] = permutarray[i][j] - permutarray[0][j];
    }
    cross(v1, v2, n);
  }
  if (i > 2) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[2];
    permutarray[2] = swapvertex;
  }

  // Make sure the fourth vertex is not coplanar with the first three.
  i = 3;
  ori = orient3dfast(permutarray[0], permutarray[1], permutarray[2], 
                     permutarray[i]);
  while ((fabs(ori) / bboxsize3) < (b->epsilon * b->epsilon * b->epsilon)) {
    i++;
    if (i == in->numberofpoints) {
      printf("Exception:  All vertices are coplanar (Tol = %g).\n",
             b->epsilon);
      terminatetetgen(this, 10);
    }
    ori = orient3dfast(permutarray[0], permutarray[1], permutarray[2], 
                       permutarray[i]);
  }
  if (i > 3) {
    // Swap to move the non-identical vertex from index i to index 1.
    swapvertex = permutarray[i];
    permutarray[i] = permutarray[3];
    permutarray[3] = swapvertex;
  }

  // Orient the first four vertices in permutarray so that they follow the
  //   right-hand rule.
  if (ori > 0.0) {
    // Swap the first two vertices.
    swapvertex = permutarray[0];
    permutarray[0] = permutarray[1];
    permutarray[1] = swapvertex;
  }

  // Create the initial Delaunay tetrahedralization.
  initialdelaunay(permutarray[0], permutarray[1], permutarray[2],
                  permutarray[3]);

  if (b->verbose) {
    printf("  Incrementally inserting vertices.\n");
  }
  insertvertexflags ivf;
  flipconstraints fc;

  // Choose algorithm: Bowyer-Watson (default) or Incremental Flip
  if (b->incrflip) {
    ivf.bowywat = 0;
    ivf.lawson = 1;
    fc.enqflag = 1;
  } else {
    ivf.bowywat = 1;
    ivf.lawson = 0;
  }


  for (i = 4; i < in->numberofpoints; i++) {
    if (pointtype(permutarray[i]) == UNUSEDVERTEX) {
      setpointtype(permutarray[i], VOLVERTEX);
    }
    if (b->brio_hilbert || b->no_sort) { // -b or -b/1
      // Start the last updated tet.
      searchtet.tet = recenttet.tet;
    } else { // -b0
      // Randomly choose the starting tet for point location.
      searchtet.tet = NULL;
    }
    ivf.iloc = (int) OUTSIDE;
    // Insert the vertex.
    if (insertpoint(permutarray[i], &searchtet, NULL, NULL, &ivf)) {
      if (flipstack != NULL) {
        // Perform flip to recover Delaunayness.
        incrementalflip(permutarray[i], (ivf.iloc == (int) OUTSIDE), &fc);
      }
    } else {
      if (ivf.iloc == (int) ONVERTEX) {
        // The point already exists. Mark it and do nothing on it.
        swapvertex = org(searchtet);
        assert(swapvertex != permutarray[i]); // SELF_CHECK
        if (b->object != tetgenbehavior::STL) {
          if (!b->quiet) {
            printf("Warning:  Point #%d is coincident with #%d. Ignored!\n",
                   pointmark(permutarray[i]), pointmark(swapvertex));
          }
        }
        setpoint2ppt(permutarray[i], swapvertex);
        setpointtype(permutarray[i], DUPLICATEDVERTEX);
        dupverts++;
      } else if (ivf.iloc == (int) NEARVERTEX) {
        swapvertex = point2ppt(permutarray[i]);
        if (!b->quiet) {
          printf("Warning:  Point %d is replaced by point %d.\n",
                 pointmark(permutarray[i]), pointmark(swapvertex));
          printf("  Avoid creating a very short edge (len = %g) (< %g).\n",
                 permutarray[i][3], b->minedgelength);
          printf("  You may try a smaller tolerance (-T) (current is %g)\n", 
                 b->epsilon);
          printf("  or use the option -M0/1 to avoid such replacement.\n");
        }
        // Remember it is a duplicated point.
        setpointtype(permutarray[i], DUPLICATEDVERTEX);
        // Count the number of duplicated points.
        dupverts++;
      }
    }
  }



  delete [] permutarray;
}

////                                                                       ////
////                                                                       ////
//// delaunay_cxx /////////////////////////////////////////////////////////////

