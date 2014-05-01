#include "../tetgen.h"
//// refine_cxx ///////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// makefacetverticesmap()    Create a map from facet to its vertices.        //
//                                                                           //
// All facets will be indexed (starting from 0).  The map is saved in two    //
// global arrays: 'idx2facetlist' and 'facetverticeslist'.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::makefacetverticesmap()
{
  arraypool *facetvertexlist, *vertlist, **paryvertlist;
  face subloop, neighsh, *parysh, *parysh1;
  point pa, *ppt, *parypt;
  verttype vt;
  int facetindex, totalvertices;
  int i, j, k;

  if (b->verbose) {
    printf("  Creating the facet vertices map.\n");
  }

  facetvertexlist = new arraypool(sizeof(arraypool *), 10);
  facetindex = totalvertices = 0;

  subfaces->traversalinit();
  subloop.sh = shellfacetraverse(subfaces);
  while (subloop.sh != NULL) {
    if (!sinfected(subloop)) {
      // A new facet. Create its vertices list.
      vertlist = new arraypool(sizeof(point *), 8);
      ppt = (point *) &(subloop.sh[3]);
      for (k = 0; k < 3; k++) {
        vt = pointtype(ppt[k]);
        if ((vt != FREESEGVERTEX) && (vt != FREEFACETVERTEX)) {
          pinfect(ppt[k]);
          vertlist->newindex((void **) &parypt);
          *parypt = ppt[k];
        }
      }
      sinfect(subloop);
      caveshlist->newindex((void **) &parysh);
      *parysh = subloop;
      for (i = 0; i < caveshlist->objects; i++) {
        parysh = (face *) fastlookup(caveshlist, i);
        setfacetindex(*parysh, facetindex);
        for (j = 0; j < 3; j++) {
          if (!isshsubseg(*parysh)) {
            spivot(*parysh, neighsh);
            assert(neighsh.sh != NULL);
            if (!sinfected(neighsh)) {
              pa = sapex(neighsh);
              if (!pinfected(pa)) {
                vt = pointtype(pa);
                if ((vt != FREESEGVERTEX) && (vt != FREEFACETVERTEX)) {
                  pinfect(pa);
                  vertlist->newindex((void **) &parypt);
                  *parypt = pa;
                }
              }
              sinfect(neighsh);
              caveshlist->newindex((void **) &parysh1);
              *parysh1 = neighsh;
            }
          }
          senextself(*parysh);
        }
      } // i
      totalvertices += (int) vertlist->objects;
      // Uninfect facet vertices.
      for (k = 0; k < vertlist->objects; k++) {
        parypt = (point *) fastlookup(vertlist, k);
        puninfect(*parypt);
      }
      caveshlist->restart();
      // Save this vertex list.
      facetvertexlist->newindex((void **) &paryvertlist);
      *paryvertlist = vertlist;
      facetindex++;
    } 
    subloop.sh = shellfacetraverse(subfaces);
  }

  // All subfaces are infected. Uninfect them.
  subfaces->traversalinit();
  subloop.sh = shellfacetraverse(subfaces);
  while (subloop.sh != NULL) {
    assert(sinfected(subloop));
    suninfect(subloop);
    subloop.sh = shellfacetraverse(subfaces);
  }

  if (b->verbose) {
    printf("  Found %ld facets.\n", facetvertexlist->objects);
  }

  idx2facetlist = new int[facetindex + 1];
  facetverticeslist = new point[totalvertices];

  totalworkmemory += ((facetindex + 1) * sizeof(int) + 
                      totalvertices * sizeof(point *));

  idx2facetlist[0] = 0;
  for (i = 0, k = 0; i < facetindex; i++) {
    paryvertlist = (arraypool **) fastlookup(facetvertexlist, i);
    vertlist = *paryvertlist;
    idx2facetlist[i + 1] = (idx2facetlist[i] + (int) vertlist->objects);
    for (j = 0; j < vertlist->objects; j++) {
      parypt = (point *) fastlookup(vertlist, j);
      facetverticeslist[k] = *parypt;
      k++;
    }
  }
  assert(k == totalvertices);

  // Free the lists.
  for (i = 0; i < facetvertexlist->objects; i++) {
    paryvertlist = (arraypool **) fastlookup(facetvertexlist, i);
    vertlist = *paryvertlist;
    delete vertlist;
  }
  delete facetvertexlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Check whether two segments, or a segment and a facet, or two facets are   //
// adjacent to each other.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::segsegadjacent(face *seg1, face *seg2)
{
  int segidx1 = getfacetindex(*seg1);
  int segidx2 = getfacetindex(*seg2);

  if (segidx1 == segidx2) return 0;

  point pa1 = segmentendpointslist[segidx1 * 2];
  point pb1 = segmentendpointslist[segidx1 * 2 + 1];
  point pa2 = segmentendpointslist[segidx2 * 2];
  point pb2 = segmentendpointslist[segidx2 * 2 + 1];

  if ((pa1 == pa2) || (pa1 == pb2) || (pb1 == pa2) || (pb1 == pb2)) {
    return 1;
  }
  return 0; 
}

int tetgenmesh::segfacetadjacent(face *subseg, face *subsh)
{
  int segidx = getfacetindex(*subseg);
  point pa = segmentendpointslist[segidx * 2];
  point pb = segmentendpointslist[segidx * 2 + 1];

  pinfect(pa);
  pinfect(pb);

  int fidx = getfacetindex(*subsh);
  int count = 0, i;

  for (i = idx2facetlist[fidx]; i < idx2facetlist[fidx+1]; i++) {
    if (pinfected(facetverticeslist[i])) count++;
  } 

  puninfect(pa);
  puninfect(pb);

  return count == 1;
}

int tetgenmesh::facetfacetadjacent(face *subsh1, face *subsh2)
{
  int count = 0, i;

  int fidx1 = getfacetindex(*subsh1);
  int fidx2 = getfacetindex(*subsh2);

  if (fidx1 == fidx2) return 0;

  for (i = idx2facetlist[fidx1]; i < idx2facetlist[fidx1+1]; i++) {
    pinfect(facetverticeslist[i]);
  }

  for (i = idx2facetlist[fidx2]; i < idx2facetlist[fidx2+1]; i++) {
    if (pinfected(facetverticeslist[i])) count++;
  }

  // Uninfect the vertices.
  for (i = idx2facetlist[fidx1]; i < idx2facetlist[fidx1+1]; i++) {
    puninfect(facetverticeslist[i]);
  }

  return count > 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkseg4encroach()    Check if an edge is encroached upon by a point.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkseg4encroach(point pa, point pb, point checkpt)
{
  // Check if the point lies inside the diametrical sphere of this seg. 
  REAL v1[3], v2[3];

  v1[0] = pa[0] - checkpt[0];
  v1[1] = pa[1] - checkpt[1];
  v1[2] = pa[2] - checkpt[2];
  v2[0] = pb[0] - checkpt[0];
  v2[1] = pb[1] - checkpt[1];
  v2[2] = pb[2] - checkpt[2];

  if (dot(v1, v2) < 0) {
    // Inside.
    if (b->metric) { // -m option.
      if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0)) {
        // The projection of 'checkpt' lies inside the segment [a,b].
        REAL prjpt[3], u, v, t;
        projpt2edge(checkpt, pa, pb, prjpt);
        // Interoplate the mesh size at the location 'prjpt'.
        u = distance(pa, pb);
        v = distance(pa, prjpt);
        t = v / u;
        // 'u' is the mesh size at 'prjpt'
        u = pa[pointmtrindex] + t * (pb[pointmtrindex] - pa[pointmtrindex]);
        v = distance(checkpt, prjpt);
        if (v < u) {
          return 1; // Encroached prot-ball!
        }
      } else {
        return 1; // NO protecting ball. Encroached.
      }
    } else {
      return 1; // Inside! Encroached.
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkseg4split()    Check if we need to split a segment.                  //
//                                                                           //
// A segment needs to be split if it is in the following case:               //
//  (1) It is encroached by an existing vertex.                              //
//  (2) It has bad quality (too long).                                       //
//  (3) Its length is larger than the mesh sizes at its endpoints.           //
//                                                                           //
// Return 1 if it needs to be split, otherwise, return 0.  'pencpt' returns  //
// an encroaching point if there exists. 'qflag' returns '1' if the segment  //
// has a length larger than the desired edge length.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkseg4split(face *chkseg, point& encpt, int& qflag)
{
  REAL ccent[3], len, r;
  int i;

  point forg = sorg(*chkseg);
  point fdest = sdest(*chkseg);

  // Initialize the return values.
  encpt = NULL;
  qflag = 0;

  len = distance(forg, fdest);
  r = 0.5 * len;
  for (i = 0; i < 3; i++) {
    ccent[i] = 0.5 * (forg[i] + fdest[i]);
  }

  // First check its quality.
  if (checkconstraints && (areabound(*chkseg) > 0.0)) {
    if (len > areabound(*chkseg)) {
      qflag = 1;
      return 1;
    }
  }

  if (b->fixedvolume) {
    if ((len * len * len) > b->maxvolume) {
      qflag = 1;
      return 1;
    }
  }

  if (b->metric) { // -m option. Check mesh size. 
    // Check if the ccent lies outside one of the prot.balls at vertices.
    if (((forg[pointmtrindex] > 0) && (r > forg[pointmtrindex])) ||
        ((fdest[pointmtrindex]) > 0 && (r > fdest[pointmtrindex]))) {
      qflag = 1; // Enforce mesh size.
      return 1;
    }
  }


  // Second check if it is encroached.
  // Comment: There may exist more than one encroaching points of this segment. 
  //   The 'encpt' returns the one which is closet to it.
  triface searchtet, spintet;
  point eapex;
  REAL d, diff, smdist = 0;
  int t1ver;

  sstpivot1(*chkseg, searchtet);
  spintet = searchtet;
  while (1) {
    eapex = apex(spintet);
    if (eapex != dummypoint) {
      d = distance(ccent, eapex);
      diff = d - r;
      if (fabs(diff) / r < b->epsilon) diff = 0.0; // Rounding.
      if (diff < 0) {
        // This segment is encroached by eapex.
        if (useinsertradius) {
          if (encpt == NULL) {
            encpt = eapex;
            smdist = d;
          } else {
            // Choose the closet encroaching point.
            if (d < smdist) {
              encpt = eapex;
              smdist = d;
            }
          }
        } else {
          encpt = eapex;
          break;
        }
      }
    }
    fnextself(spintet);
    if (spintet.tet == searchtet.tet) break;
  } // while (1)

  if (encpt != NULL) {
    return 1;
  }

  return 0; // No need to split it.
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// splitsegment()    Split a segment.                                        //
//                                                                           //
// The segment 'splitseg' is intended to be split. It will be split if it    //
// is in one of the following cases:                                         //
//   (1) It is encroached by an existing vertex 'encpt != NULL'; or          //
//   (2) It is in bad quality 'qflag == 1'; or                               //
//   (3) Its length is larger than the mesh sizes at its endpoints.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::splitsegment(face *splitseg, point encpt, REAL rrp, 
                             point encpt1, point encpt2, int qflag, 
                             int chkencflag)
{
  point pa = sorg(*splitseg);
  point pb = sdest(*splitseg);



  if ((encpt == NULL) && (qflag == 0)) {
    if (useinsertradius) {
      // Do not split this segment if the length is smaller than the smaller
      //   insertion radius at its endpoints.
      REAL len = distance(pa, pb);
      REAL smrrv = getpointinsradius(pa);
      REAL rrv = getpointinsradius(pb);
      if (rrv > 0) {
        if (smrrv > 0) {
          if (rrv < smrrv) {
            smrrv = rrv;
          }
        } else {
          smrrv = rrv;
        }
      }
      if (smrrv > 0) {
        if ((fabs(smrrv - len) / len) < b->epsilon) smrrv = len;
        if (len < smrrv) {
          return 0;
        }
      }
    }
  }

  if (b->nobisect) { // With -Y option.
    // Only split this segment if it is allowed to be split.
    if (checkconstraints) {
      // Check if it has a non-zero length bound. 
      if (areabound(*splitseg) == 0) {
        // It is not allowed.  However, if all of facets containing this seg
        //   is allowed to be split, we still split it.
        face parentsh, spinsh;
        //splitseg.shver = 0;
        spivot(*splitseg, parentsh);
        if (parentsh.sh == NULL) {
          return 0; // A dangling segment. Do not split it.
        }
        spinsh = parentsh;
        while (1) {
          if (areabound(spinsh) == 0) break;
          spivotself(spinsh);
          if (spinsh.sh == parentsh.sh) break;
        }
        if (areabound(spinsh) == 0) {
          // All facets at this seg are not allowed to be split.
          return 0;  // Do not split it.
        }
      }
    } else {
      return 0; // Do not split this segment.
    }
  } // if (b->nobisect)

  triface searchtet;
  face searchsh;
  point newpt;
  insertvertexflags ivf;

  makepoint(&newpt, FREESEGVERTEX);
  getsteinerptonsegment(splitseg, encpt, newpt);

  // Split the segment by the Bowyer-Watson algorithm.
  sstpivot1(*splitseg, searchtet);
  ivf.iloc = (int) ONEDGE;
  // Use Bowyer-Watson algorithm. Preserve subsegments and subfaces;
  ivf.bowywat = 3;
  ivf.validflag = 1; // Validate the B-W cavity.
  ivf.lawson = 2; // Do flips to recover Delaunayness.
  ivf.rejflag = 0;     // Do not check encroachment of new segments/facets.
  if (b->metric) {
    ivf.rejflag |= 4;  // Do check encroachment of protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = (int) INSTAR; // ivf.iloc;
  ivf.sbowywat = 3; // ivf.bowywat;  // Surface mesh options.
  ivf.splitbdflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;
  ivf.smlenflag = useinsertradius;


  if (insertpoint(newpt, &searchtet, &searchsh, splitseg, &ivf)) {
    st_segref_count++;
    if (steinerleft > 0) steinerleft--;
    if (useinsertradius) {
      // Update 'rv' (to be the shortest distance).
      REAL rv = ivf.smlen, rp;
      if (pointtype(ivf.parentpt) == FREESEGVERTEX) {
        face parentseg1, parentseg2;
        sdecode(point2sh(newpt), parentseg1);
        sdecode(point2sh(ivf.parentpt), parentseg2);
        if (segsegadjacent(&parentseg1, &parentseg2)) {
          rp = getpointinsradius(ivf.parentpt);
          if (rv < rp) {
            rv = rp; // The relaxed insertion radius of 'newpt'.
          }
        }
      } else if (pointtype(ivf.parentpt) == FREEFACETVERTEX) {
        face parentseg, parentsh;
        sdecode(point2sh(newpt), parentseg);
        sdecode(point2sh(ivf.parentpt), parentsh);
        if (segfacetadjacent(&parentseg, &parentsh)) {
          rp = getpointinsradius(ivf.parentpt);
          if (rv < rp) {
            rv = rp; // The relaxed insertion radius of 'newpt'.
          }            
        }
      }
      setpointinsradius(newpt, rv);
    }
    if (flipstack != NULL) {
      flipconstraints fc;
      fc.chkencflag = chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      unflipqueue->restart();
    }
    return 1;
  } else {
    // Point is not inserted.
    pointdealloc(newpt);
    return 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// repairencsegs()    Repair encroached (sub) segments.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::repairencsegs(int chkencflag)
{
  face *bface;
  point encpt = NULL;
  int qflag = 0;

  // Loop until the pool 'badsubsegs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badsubsegs->items > 0) && (steinerleft != 0)) {
    badsubsegs->traversalinit();
    bface = (face *) badsubsegs->traverse();
    while ((bface != NULL) && (steinerleft != 0)) {
      // Skip a deleleted element.
      if (bface->shver >= 0) {
        // A queued segment may have been deleted (split).
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          // A queued segment may have been processed. 
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
            if (checkseg4split(bface, encpt, qflag)) {
              splitsegment(bface, encpt, 0, NULL, NULL, qflag, chkencflag);
            }
          }
        }
        // Remove this entry from list.
        bface->shver = -1; // Signal it as a deleted element.
        badsubsegs->dealloc((void *) bface);
      }
      bface = (face *) badsubsegs->traverse();
    }
  }

  if (badsubsegs->items > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else {
      assert(0); // Unknown case.
    }
    badsubsegs->traversalinit();
    bface = (face *) badsubsegs->traverse();
    while (bface  != NULL) {
      // Skip a deleleted element.
      if (bface->shver >= 0) {
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
          }
        }
      }
      bface = (face *) badsubsegs->traverse();
    }
    badsubsegs->restart();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// enqueuesubface()    Queue a subface or a subsegment for encroachment chk. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::enqueuesubface(memorypool *pool, face *chkface)
{
  if (!smarktest2ed(*chkface)) {
    smarktest2(*chkface); // Only queue it once.
    face *queface = (face *) pool->alloc();
    *queface = *chkface;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkfac4encroach()    Check if a subface is encroached by a point.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkfac4encroach(point pa, point pb, point pc, point checkpt,
                                  REAL* cent, REAL* r)
{
  REAL rd, len;

  circumsphere(pa, pb, pc, NULL, cent, &rd);
  assert(rd != 0);
  len = distance(cent, checkpt);
  if ((fabs(len - rd) / rd) < b->epsilon) len = rd; // Rounding.
 
  if (len < rd) {
    // The point lies inside the circumsphere of this face.
    if (b->metric) { // -m option.
      if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0) &&
          (pc[pointmtrindex] > 0)) {
        // Get the projection of 'checkpt' in the plane of pa, pb, and pc.
        REAL prjpt[3], n[3];
        REAL a, a1, a2, a3;
        projpt2face(checkpt, pa, pb, pc, prjpt);
        // Get the face area of [a,b,c].
        facenormal(pa, pb, pc, n, 1, NULL);
        a = sqrt(dot(n,n));
        // Get the face areas of [a,b,p], [b,c,p], and [c,a,p].
        facenormal(pa, pb, prjpt, n, 1, NULL);
        a1 = sqrt(dot(n,n));
        facenormal(pb, pc, prjpt, n, 1, NULL);
        a2 = sqrt(dot(n,n));
        facenormal(pc, pa, prjpt, n, 1, NULL);
        a3 = sqrt(dot(n,n));
        if ((fabs(a1 + a2 + a3 - a) / a) < b->epsilon) {
          // This face contains the projection.
          // Get the mesh size at the location of the projection point.
          rd = a1 / a * pc[pointmtrindex]
             + a2 / a * pa[pointmtrindex]
             + a3 / a * pb[pointmtrindex];
          len = distance(prjpt, checkpt);
          if (len < rd) {
            return 1; // Encroached.
          }
        }
      } else {
        return 1;  // No protecting ball. Encroached.
      }
    } else {
      *r = rd;
      return 1;  // Encroached.
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkfac4split()    Check if a subface needs to be split.                 //
//                                                                           //
// A subface needs to be split if it is in the following case:               //
//  (1) It is encroached by an existing vertex.                              //
//  (2) It has bad quality (has a small angle, -q).                          //
//  (3) It's area is larger than a prescribed value (.var).                  //
//                                                                           //
// Return 1 if it needs to be split, otherwise, return 0.                    //
// 'chkfac' represents its longest edge.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkfac4split(face *chkfac, point& encpt, int& qflag, 
                               REAL *cent)
{
  point pa, pb, pc;
  REAL area, rd, len;
  REAL A[4][4], rhs[4], D;
  int indx[4];
  int i;

  encpt = NULL;
  qflag = 0;

  pa = sorg(*chkfac);
  pb = sdest(*chkfac);
  pc = sapex(*chkfac);

  // Compute the coefficient matrix A (3x3).
  A[0][0] = pb[0] - pa[0];
  A[0][1] = pb[1] - pa[1];
  A[0][2] = pb[2] - pa[2]; // vector V1 (pa->pb)
  A[1][0] = pc[0] - pa[0];
  A[1][1] = pc[1] - pa[1];
  A[1][2] = pc[2] - pa[2]; // vector V2 (pa->pc)
  cross(A[0], A[1], A[2]); // vector V3 (V1 X V2)

  area = 0.5 * sqrt(dot(A[2], A[2])); // The area of [a,b,c].

  // Compute the right hand side vector b (3x1).
  rhs[0] = 0.5 * dot(A[0], A[0]); // edge [a,b]
  rhs[1] = 0.5 * dot(A[1], A[1]); // edge [a,c]
  rhs[2] = 0.0;

  // Solve the 3 by 3 equations use LU decomposition with partial 
  //   pivoting and backward and forward substitute.
  if (!lu_decmp(A, 3, indx, &D, 0)) {
    // A degenerate triangle. 
    assert(0);
  }

  lu_solve(A, 3, indx, rhs, 0);
  cent[0] = pa[0] + rhs[0];
  cent[1] = pa[1] + rhs[1];
  cent[2] = pa[2] + rhs[2];
  rd = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);

  if (checkconstraints && (areabound(*chkfac) > 0.0)) {
    // Check if the subface has too big area.
    if (area > areabound(*chkfac)) {
      qflag = 1;
      return 1;
    }
  }

  if (b->fixedvolume) {
    if ((area * sqrt(area)) > b->maxvolume) {
      qflag = 1;
      return 1;
    }
  }

  if (b->varvolume) {
    triface adjtet;
    REAL volbnd;
    int t1ver;

    stpivot(*chkfac, adjtet);
    if (!ishulltet(adjtet)) {
      volbnd = volumebound(adjtet.tet);
      if ((volbnd > 0) && (area * sqrt(area)) > volbnd) {
        qflag = 1;
        return 1;
      }
    }
    fsymself(adjtet);
    if (!ishulltet(adjtet)) {
      volbnd = volumebound(adjtet.tet);
      if ((volbnd > 0) && (area * sqrt(area)) > volbnd) {
        qflag = 1;
        return 1;
      }
    }
  }

  if (b->metric) { // -m option. Check mesh size. 
    // Check if the ccent lies outside one of the prot.balls at vertices.
    if (((pa[pointmtrindex] > 0) && (rd > pa[pointmtrindex])) ||
        ((pb[pointmtrindex] > 0) && (rd > pb[pointmtrindex])) ||
        ((pc[pointmtrindex] > 0) && (rd > pc[pointmtrindex]))) {
      qflag = 1; // Enforce mesh size.
      return 1;
    }
  }

  triface searchtet;
  REAL smlen = 0;

  // Check if this subface is locally encroached.
  for (i = 0; i < 2; i++) {
    stpivot(*chkfac, searchtet);
    if (!ishulltet(searchtet)) {
      len = distance(oppo(searchtet), cent);
      if ((fabs(len - rd) / rd) < b->epsilon) len = rd;// Rounding.
      if (len < rd) {
        if (smlen == 0) {
          smlen = len;
          encpt = oppo(searchtet);
        } else {
          if (len < smlen) {
            smlen = len;
            encpt = oppo(searchtet);
          }
        }
        //return 1;
      }
    }
    sesymself(*chkfac);
  }

  return encpt != NULL; //return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// splitsubface()    Split a subface.                                        //
//                                                                           //
// The subface may be encroached, or in bad-quality. It is split at its cir- //
// cumcenter ('ccent'). Do not split it if 'ccent' encroaches upon any seg-  //
// ment. Instead, one of the encroached segments is split.  It is possible   //
// that none of the encroached segments can be split.                        //
//                                                                           //
// The return value indicates whether a new point is inserted (> 0) or not   //
// (= 0).  Furthermore, it is inserted on an encroached segment (= 1) or     //
// in-side the facet (= 2).                                                  //
//                                                                           //
// 'encpt' is a vertex encroaching upon this subface, i.e., it causes the    //
// split of this subface. If 'encpt' is NULL, then the cause of the split    //
// this subface is a rejected tet circumcenter 'p', and 'encpt1' is the      //
// parent of 'p'.                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::splitsubface(face *splitfac, point encpt, point encpt1, 
                             int qflag, REAL *ccent, int chkencflag)
{
  point pa = sorg(*splitfac);
  point pb = sdest(*splitfac);
  point pc = sapex(*splitfac);



  if (b->nobisect) { // With -Y option.
    if (checkconstraints) {
      // Only split if it is allowed to be split.
      // Check if this facet has a non-zero constraint.
      if (areabound(*splitfac) == 0) {
        return 0; // Do not split it.
      }
    } else {
      return 0;
    }
  } // if (b->nobisect)

  face searchsh;
  insertvertexflags ivf;
  point newpt;
  REAL rv = 0., rp; // Insertion radius of newpt.
  int i;

  // Initialize the inserting point.
  makepoint(&newpt, FREEFACETVERTEX);
  // Split the subface at its circumcenter.
  for (i = 0; i < 3; i++) newpt[i] = ccent[i];

  if (useinsertradius) {
    if (encpt != NULL) {
      rv = distance(newpt, encpt);
      if (pointtype(encpt) == FREESEGVERTEX) {
        face parentseg;
        sdecode(point2sh(encpt), parentseg);
        if (segfacetadjacent(&parentseg, splitfac)) {
          rp = getpointinsradius(encpt);
          if (rv < (sqrt(2.0) * rp)) {
            // This insertion may cause no termination. 
            pointdealloc(newpt);
            return 0; // Reject the insertion of newpt.
          }
        }
      } else if (pointtype(encpt) == FREEFACETVERTEX) {
        face parentsh;
        sdecode(point2sh(encpt), parentsh);
        if (facetfacetadjacent(&parentsh, splitfac)) {
          rp = getpointinsradius(encpt);
          if (rv < rp) {
            pointdealloc(newpt);
            return 0; // Reject the insertion of newpt.
          }
        }
      }
    }
  } // if (useinsertradius)

  // Search a subface which contains 'newpt'.
  searchsh = *splitfac;
  // Calculate an above point. It lies above the plane containing
  //   the subface [a,b,c], and save it in dummypoint. Moreover,
  //   the vector cent->dummypoint is the normal of the plane.
  calculateabovepoint4(newpt, pa, pb, pc);
  //   Parameters: 'aflag' = 1, - above point exists.
  //   'cflag' = 0, - non-convex, check co-planarity of the result.
  //   'rflag' = 0, - no need to round the locating result.
  ivf.iloc = (int) slocate(newpt, &searchsh, 1, 0, 0);

  if (!((ivf.iloc == (int) ONFACE) || (ivf.iloc == (int) ONEDGE))) {
    pointdealloc(newpt);
    return 0;
  }


  triface searchtet;
  face *paryseg;
  int splitflag;

  // Insert the point.
  stpivot(searchsh, searchtet);
  //assert((ivf.iloc == (int) ONFACE) || (ivf.iloc == (int) ONEDGE));
  // Use Bowyer-Watson algorithm. Preserve subsegments and subfaces;
  ivf.bowywat = 3; 
  ivf.lawson = 2;
  ivf.rejflag = 1; // Do check the encroachment of segments.
  if (b->metric) {
    ivf.rejflag |= 4;  // Do check encroachment of protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = (int) INSTAR; // ivf.iloc;
  ivf.sbowywat = 3; // ivf.bowywat;
  ivf.splitbdflag = 1;
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  ivf.refineflag = 2;
  ivf.refinesh = searchsh;
  ivf.smlenflag = useinsertradius; // Update the insertion radius.


  if (insertpoint(newpt, &searchtet, &searchsh, NULL, &ivf)) {
    st_facref_count++;
    if (steinerleft > 0) steinerleft--;
    if (useinsertradius) {
      // Update 'rv' (to be the shortest distance).
      rv = ivf.smlen;
      if (pointtype(ivf.parentpt) == FREESEGVERTEX) {
        face parentseg, parentsh;
        sdecode(point2sh(ivf.parentpt), parentseg);
        sdecode(point2sh(newpt), parentsh);
        if (segfacetadjacent(&parentseg, &parentsh)) {
          rp = getpointinsradius(ivf.parentpt);
          if (rv < (sqrt(2.0) * rp)) {
            rv = sqrt(2.0) * rp; // The relaxed insertion radius of 'newpt'.
          }
        }
      } else if (pointtype(ivf.parentpt) == FREEFACETVERTEX) {
        face parentsh1, parentsh2;
        sdecode(point2sh(ivf.parentpt), parentsh1);
        sdecode(point2sh(newpt), parentsh2);
        if (facetfacetadjacent(&parentsh1, &parentsh2)) {
          rp = getpointinsradius(ivf.parentpt);
          if (rv < rp) {
            rv = rp; // The relaxed insertion radius of 'newpt'.
          }          
        }
      }
      setpointinsradius(newpt, rv);
    } // if (useinsertradius)
    if (flipstack != NULL) {
      flipconstraints fc;
      fc.chkencflag = chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      unflipqueue->restart();
    }
    return 1;
  } else {
    // Point was not inserted.
    pointdealloc(newpt);
    if (ivf.iloc == (int) ENCSEGMENT) {
      // Select an encroached segment and split it.
      splitflag = 0;
      for (i = 0; i < encseglist->objects; i++) {
        paryseg = (face *) fastlookup(encseglist, i);
        if (splitsegment(paryseg, NULL, rv, encpt, encpt1, qflag, 
                         chkencflag | 1)) {
          splitflag = 1; // A point is inserted on a segment.
          break;
        }
      }
      encseglist->restart();
      if (splitflag) {
        // Some segments may need to be repaired.
        repairencsegs(chkencflag | 1);
        // Queue this subface if it is still alive and not queued.
        //if ((splitfac->sh != NULL) && (splitfac->sh[3] != NULL)) {
        //  // Only queue it if 'qflag' is set.
        //  if (qflag) { 
        //    enqueuesubface(badsubfacs, splitfac);
        //  }
        //}
      }
      return splitflag;
    } else {
      return 0;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// repairencfacs()    Repair encroached subfaces.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::repairencfacs(int chkencflag)
{
  face *bface;
  point encpt = NULL;
  int qflag = 0;
  REAL ccent[3];

  // Loop until the pool 'badsubfacs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badsubfacs->items > 0) && (steinerleft != 0)) {
    badsubfacs->traversalinit();
    bface = (face *) badsubfacs->traverse();
    while ((bface != NULL) && (steinerleft != 0)) {
      // Skip a deleted element.
      if (bface->shver >= 0) {
        // A queued subface may have been deleted (split).
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          // A queued subface may have been processed. 
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
            if (checkfac4split(bface, encpt, qflag, ccent)) {
              splitsubface(bface, encpt, NULL, qflag, ccent, chkencflag);
            }
          }
        }
        bface->shver = -1; // Signal it as a deleted element.
        badsubfacs->dealloc((void *) bface); // Remove this entry from list.
      }
      bface = (face *) badsubfacs->traverse();
    }
  }

  if (badsubfacs->items > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else {
      assert(0); // Unknown case.
    }
    badsubfacs->traversalinit();
    bface = (face *) badsubfacs->traverse();
    while (bface  != NULL) {
      // Skip a deleted element.
      if (bface->shver >= 0) {
        if ((bface->sh != NULL) && (bface->sh[3] != NULL)) {
          if (smarktest2ed(*bface)) {
            sunmarktest2(*bface);
          }
        }
      }
      bface = (face *) badsubfacs->traverse();
    }
    badsubfacs->restart();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// enqueuetetrahedron()    Queue a tetrahedron for quality check.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::enqueuetetrahedron(triface *chktet)
{
  if (!marktest2ed(*chktet)) {
    marktest2(*chktet); // Only queue it once.
    triface *quetet = (triface *) badtetrahedrons->alloc();
    *quetet = *chktet;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checktet4split()    Check if the tet needs to be split.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checktet4split(triface *chktet, int &qflag, REAL *ccent) 
{
  point pa, pb, pc, pd, *ppt;
  REAL vda[3], vdb[3], vdc[3];
  REAL vab[3], vbc[3], vca[3];
  REAL N[4][3], L[4], cosd[6], elen[6];
  REAL maxcosd, vol, volbnd, smlen = 0, rd;
  REAL A[4][4], rhs[4], D;
  int indx[4];
  int i, j;

  if (b->convex) { // -c
    // Skip this tet if it lies in the exterior.
    if (elemattribute(chktet->tet, numelemattrib - 1) == -1.0) {
      return 0;
    }
  }

  qflag = 0;

  pd = (point) chktet->tet[7];
  if (pd == dummypoint) {
    return 0; // Do not split a hull tet.
  }

  pa = (point) chktet->tet[4];
  pb = (point) chktet->tet[5];
  pc = (point) chktet->tet[6];

  // Get the edge vectors vda: d->a, vdb: d->b, vdc: d->c.
  // Set the matrix A = [vda, vdb, vdc]^T.
  for (i = 0; i < 3; i++) A[0][i] = vda[i] = pa[i] - pd[i];
  for (i = 0; i < 3; i++) A[1][i] = vdb[i] = pb[i] - pd[i];
  for (i = 0; i < 3; i++) A[2][i] = vdc[i] = pc[i] - pd[i];

  // Get the other edge vectors.
  for (i = 0; i < 3; i++) vab[i] = pb[i] - pa[i];
  for (i = 0; i < 3; i++) vbc[i] = pc[i] - pb[i];
  for (i = 0; i < 3; i++) vca[i] = pa[i] - pc[i];

  if (!lu_decmp(A, 3, indx, &D, 0)) {
    // A degenerated tet (vol = 0).
    // This is possible due to the use of exact arithmetic.  We temporarily
    //   leave this tet. It should be fixed by mesh optimization.
    return 0; 
  }

  // Check volume if '-a#' and '-a' options are used.
  if (b->varvolume || b->fixedvolume) {
    vol = fabs(A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]) / 6.0;
    if (b->fixedvolume) {
      if (vol > b->maxvolume) {
        qflag = 1;
      }
    } 
    if (!qflag && b->varvolume) {
      volbnd = volumebound(chktet->tet);
      if ((volbnd > 0.0) && (vol > volbnd)) {
        qflag = 1;
      }
    }
    if (qflag == 1) {
      // Calculate the circumcenter of this tet.
      rhs[0] = 0.5 * dot(vda, vda);
      rhs[1] = 0.5 * dot(vdb, vdb);
      rhs[2] = 0.5 * dot(vdc, vdc);
      lu_solve(A, 3, indx, rhs, 0);            
      for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
      return 1;
    }
  }

  if (b->metric) { // -m option. Check mesh size. 
    // Calculate the circumradius of this tet.
    rhs[0] = 0.5 * dot(vda, vda);
    rhs[1] = 0.5 * dot(vdb, vdb);
    rhs[2] = 0.5 * dot(vdc, vdc);
    lu_solve(A, 3, indx, rhs, 0);            
    for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
    rd = sqrt(dot(rhs, rhs));
    // Check if the ccent lies outside one of the prot.balls at vertices.
    ppt = (point *) &(chktet->tet[4]);
    for (i = 0; i < 4; i++) {
      if (ppt[i][pointmtrindex] > 0) {
        if (rd > ppt[i][pointmtrindex]) {
          qflag = 1; // Enforce mesh size.
          return 1;
        }
      }
    }
  }

  if (in->tetunsuitable != NULL) {
    // Execute the user-defined meshing sizing evaluation.
    if ((*(in->tetunsuitable))(pa, pb, pc, pd, NULL, 0)) {
      // Calculate the circumcenter of this tet.
      rhs[0] = 0.5 * dot(vda, vda);
      rhs[1] = 0.5 * dot(vdb, vdb);
      rhs[2] = 0.5 * dot(vdc, vdc);
      lu_solve(A, 3, indx, rhs, 0);            
      for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
      return 1;
    }
  }

  if (useinsertradius) {
    // Do not split this tet if the shortest edge is shorter than the
    //   insertion radius of one of its endpoints.
    triface checkedge;
    point e1, e2;
    REAL rrv, smrrv;

    // Get the shortest edge of this tet.
    checkedge.tet = chktet->tet;
    for (i = 0; i < 6; i++) {
      checkedge.ver = edge2ver[i];
      e1 = org(checkedge);
      e2 = dest(checkedge);
      elen[i] = distance(e1, e2);
      if (i == 0) {
        smlen = elen[i];
        j = 0;
      } else {
        if (elen[i] < smlen) {
          smlen = elen[i];
          j = i;
        }
      }
    }
    // Check if the edge is too short.
    checkedge.ver = edge2ver[j];
    // Get the smallest rrv of e1 and e2.
    // Note: if rrv of e1 and e2 is zero. Do not use it.
    e1 = org(checkedge);
    smrrv = getpointinsradius(e1);
    e2 = dest(checkedge);
    rrv = getpointinsradius(e2);
    if (rrv > 0) {
      if (smrrv > 0) {
        if (rrv < smrrv) {
          smrrv = rrv;
        }
      } else {
        smrrv = rrv;
      }
    }
    if (smrrv > 0) {
      // To avoid rounding error, round smrrv before doing comparison.
      if ((fabs(smrrv - smlen) / smlen) < b->epsilon) {
        smrrv = smlen;
      }
      if (smrrv > smlen) {
        return 0;
      }
    }
  } // if (useinsertradius)

  // Check the radius-edge ratio. Set by -q#.
  if (b->minratio > 0) { 
    // Calculate the circumcenter and radius of this tet.
    rhs[0] = 0.5 * dot(vda, vda);
    rhs[1] = 0.5 * dot(vdb, vdb);
    rhs[2] = 0.5 * dot(vdc, vdc);
    lu_solve(A, 3, indx, rhs, 0);            
    for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
    rd = sqrt(dot(rhs, rhs));
    if (!useinsertradius) {
      // Calculate the shortest edge length.
      elen[0] = dot(vda, vda);
      elen[1] = dot(vdb, vdb);
      elen[2] = dot(vdc, vdc);
      elen[3] = dot(vab, vab);
      elen[4] = dot(vbc, vbc);
      elen[5] = dot(vca, vca);
      smlen = elen[0]; //sidx = 0;
      for (i = 1; i < 6; i++) {
        if (smlen > elen[i]) { 
          smlen = elen[i]; //sidx = i; 
        }
      }
      smlen = sqrt(smlen);
    }
    D = rd / smlen;
    if (D > b->minratio) {
      // A bad radius-edge ratio.
      return 1;
    }
  }

  // Check the minimum dihedral angle. Set by -qq#.
  if (b->mindihedral > 0) { 
    // Compute the 4 face normals (N[0], ..., N[3]).
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) N[j][i] = 0.0;
      N[j][j] = 1.0;  // Positive means the inside direction
      lu_solve(A, 3, indx, N[j], 0);
    }
    for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
    // Normalize the normals.
    for (i = 0; i < 4; i++) {
      L[i] = sqrt(dot(N[i], N[i]));
      assert(L[i] > 0);
      //if (L[i] > 0.0) {
        for (j = 0; j < 3; j++) N[i][j] /= L[i];
      //}
    }
    // Calculate the six dihedral angles.
    cosd[0] = -dot(N[0], N[1]); // Edge cd, bd, bc.
    cosd[1] = -dot(N[0], N[2]);
    cosd[2] = -dot(N[0], N[3]);
    cosd[3] = -dot(N[1], N[2]); // Edge ad, ac
    cosd[4] = -dot(N[1], N[3]);
    cosd[5] = -dot(N[2], N[3]); // Edge ab
    // Get the smallest dihedral angle.
    //maxcosd = mincosd = cosd[0];
    maxcosd = cosd[0];
    for (i = 1; i < 6; i++) {
      //if (cosd[i] > maxcosd) maxcosd = cosd[i];
      maxcosd = (cosd[i] > maxcosd ? cosd[i] : maxcosd);
      //mincosd = (cosd[i] < mincosd ? cosd[i] : maxcosd);
    }
    if (maxcosd > cosmindihed) {
      // Calculate the circumcenter of this tet.
      // A bad dihedral angle.
      //if ((b->quality & 1) == 0) {
        rhs[0] = 0.5 * dot(vda, vda);
        rhs[1] = 0.5 * dot(vdb, vdb);
        rhs[2] = 0.5 * dot(vdc, vdc);
        lu_solve(A, 3, indx, rhs, 0);            
        for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
        //*rd = sqrt(dot(rhs, rhs));
      //}
      return 1;
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// splittetrahedron()    Split a tetrahedron.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::splittetrahedron(triface* splittet, int qflag, REAL *ccent, 
                                 int chkencflag)
{
  triface searchtet;
  face *paryseg;
  point newpt;
  badface *bface;
  insertvertexflags ivf;
  int splitflag;
  int i;



  REAL rv = 0.; // Insertion radius of 'newpt'.

  makepoint(&newpt, FREEVOLVERTEX);
  for (i = 0; i < 3; i++) newpt[i] = ccent[i];

  if (useinsertradius) {
    rv = distance(newpt, org(*splittet));
    setpointinsradius(newpt, rv);
  }

  searchtet = *splittet;
  ivf.iloc = (int) OUTSIDE;
  // Use Bowyer-Watson algorithm. Preserve subsegments and subfaces;
  ivf.bowywat = 3;
  ivf.lawson = 2;
  ivf.rejflag = 3;  // Do check for encroached segments and subfaces.
  if (b->metric) {
    ivf.rejflag |= 4; // Reject it if it lies in some protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = ivf.sbowywat = 0; // No use.
  ivf.splitbdflag = 0; // No use.
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  ivf.refineflag = 1;
  ivf.refinetet = *splittet;


  if (insertpoint(newpt, &searchtet, NULL, NULL, &ivf)) {
    // Vertex is inserted.
    st_volref_count++;
    if (steinerleft > 0) steinerleft--;
    if (flipstack != NULL) {
      flipconstraints fc;
      fc.chkencflag = chkencflag;
      fc.enqflag = 2;
      lawsonflip3d(&fc);
      unflipqueue->restart();
    }
    return 1;
  } else {
    // Point is not inserted.
    pointdealloc(newpt);
    // Check if there are encroached segments/subfaces.
    if (ivf.iloc == (int) ENCSEGMENT) {
      splitflag = 0;
      //if (!b->nobisect) { // not -Y option
      if (!b->nobisect || checkconstraints) {  
        // Select an encroached segment and split it.
        for (i = 0; i < encseglist->objects; i++) {
          paryseg = (face *) fastlookup(encseglist, i);
          if (splitsegment(paryseg, NULL, rv, org(*splittet), NULL, qflag, 
                           chkencflag | 3)) {
            splitflag = 1; // A point is inserted on a segment.
            break;
          }
        }
      } // if (!b->nobisect)
      encseglist->restart();
      if (splitflag) {
        // Some segments may need to be repaired.
        repairencsegs(chkencflag | 3);
        // Some subfaces may need to be repaired.
        repairencfacs(chkencflag | 2);
        // Queue the tet if it is still alive and not queued.
        if ((splittet->tet != NULL) && (splittet->tet[4] != NULL)) {
          enqueuetetrahedron(splittet);
        }
      }
      return splitflag;
    } else if (ivf.iloc == (int) ENCSUBFACE) {
      splitflag = 0;
      //if (!b->nobisect) { // not -Y option
      if (!b->nobisect || checkconstraints) {
        // Select an encroached subface and split it.
        for (i = 0; i < encshlist->objects; i++) {
          bface = (badface *) fastlookup(encshlist, i);
          if (splitsubface(&(bface->ss), NULL, org(*splittet), qflag, 
                           bface->cent, chkencflag | 2)){
            splitflag = 1; // A point is inserted on a subface or a segment.
            break;
          }
        }
      } // if (!b->nobisect)
      encshlist->restart();
      if (splitflag) {
        assert(badsubsegs->items == 0l);
        // Some subfaces may need to be repaired.
        repairencfacs(chkencflag | 2);
        // Queue the tet if it is still alive.
        if ((splittet->tet != NULL) && (splittet->tet[4] != NULL)) {
          enqueuetetrahedron(splittet);
        }
      }
      return splitflag;
    }
    return 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// repairbadtets()    Repair bad quality tetrahedra.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::repairbadtets(int chkencflag)
{
  triface *bface;
  REAL ccent[3];
  int qflag = 0;


  // Loop until the pool 'badsubfacs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badtetrahedrons->items > 0) && (steinerleft != 0)) {
    badtetrahedrons->traversalinit();
    bface = (triface *) badtetrahedrons->traverse();
    while ((bface != NULL) && (steinerleft != 0)) {
      // Skip a deleted element.
      if (bface->ver >= 0) {
        // A queued tet may have been deleted.
        if (!isdeadtet(*bface)) {
          // A queued tet may have been processed.
          if (marktest2ed(*bface)) {
            unmarktest2(*bface);
            if (checktet4split(bface, qflag, ccent)) {
              splittetrahedron(bface, qflag, ccent, chkencflag);
            }
          }
        }
        bface->ver = -1; // Signal it as a deleted element.
        badtetrahedrons->dealloc((void *) bface);
      }
      bface = (triface *) badtetrahedrons->traverse();
    }
  }

  if (badtetrahedrons->items > 0) {
    if (steinerleft == 0) {
      if (b->verbose) {
        printf("The desired number of Steiner points is reached.\n");
      }
    } else {
      assert(0); // Unknown case.
    }
    // Unmark all queued tet.
    badtetrahedrons->traversalinit();
    bface = (triface *) badtetrahedrons->traverse();
    while (bface != NULL) {
      // Skip a deleted element.
      if (bface->ver >= 0) {
        if (!isdeadtet(*bface)) {
          if (marktest2ed(*bface)) {
            unmarktest2(*bface);
          }
        }
      }
      bface = (triface *) badtetrahedrons->traverse();
    }
    // Clear the pool.
    badtetrahedrons->restart();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunayrefinement()    Refine the mesh by Delaunay refinement.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunayrefinement()
{
  triface checktet;
  face checksh;
  face checkseg;
  long steinercount;
  int chkencflag;

  long bak_segref_count, bak_facref_count, bak_volref_count;
  long bak_flipcount = flip23count + flip32count + flip44count;

  if (!b->quiet) {
    printf("Refining mesh...\n");
  }

  if (b->verbose) {
    printf("  Min radiu-edge ratio = %g.\n", b->minratio);
    printf("  Min dihedral   angle = %g.\n", b->mindihedral);
    //printf("  Min Edge length = %g.\n", b->minedgelength);
  }

  steinerleft = b->steinerleft;  // Upperbound of # Steiner points (by -S#).
  if (steinerleft > 0) {
    // Check if we've already used up the given number of Steiner points.
    steinercount = st_segref_count + st_facref_count + st_volref_count;
    if (steinercount < steinerleft) {
      steinerleft -= steinercount;
    } else {
      if (!b->quiet) {
        printf("\nWarning:  ");
        printf("The desired number of Steiner points (%d) has reached.\n\n",
               b->steinerleft);
      }
      return; // No more Steiner points.
    }
  }

  if (useinsertradius) {
    if ((b->plc && b->nobisect) || b->refine) { // '-pY' or '-r' option.
      makesegmentendpointsmap();
    }
    makefacetverticesmap();
  }


  encseglist = new arraypool(sizeof(face), 8);
  encshlist = new arraypool(sizeof(badface), 8);


  //if (!b->nobisect) { // if no '-Y' option
  if (!b->nobisect || checkconstraints) {
    if (b->verbose) {
      printf("  Splitting encroached subsegments.\n");
    }

    chkencflag = 1; // Only check encroaching subsegments.
    steinercount = points->items;

    // Initialize the pool of encroached subsegments.
    badsubsegs = new memorypool(sizeof(face), b->shellfaceperblock, 
                                sizeof(void *), 0);

    // Add all segments into the pool.
    subsegs->traversalinit();
    checkseg.sh = shellfacetraverse(subsegs);
    while (checkseg.sh != (shellface *) NULL) {
      enqueuesubface(badsubsegs, &checkseg);
      checkseg.sh = shellfacetraverse(subsegs);
    }

    // Split all encroached segments.
    repairencsegs(chkencflag);

    if (b->verbose) {
      printf("  Added %ld Steiner points.\n", points->items - steinercount);
    }

    if (b->reflevel > 1) { // '-D2' option
      if (b->verbose) {
        printf("  Splitting encroached subfaces.\n");
      }

      chkencflag = 2; // Only check encroaching subfaces.
      steinercount = points->items;
      bak_segref_count = st_segref_count;
      bak_facref_count = st_facref_count;

      // Initialize the pool of encroached subfaces.
      badsubfacs = new memorypool(sizeof(face), b->shellfaceperblock, 
                                  sizeof(void *), 0);

      // Add all subfaces into the pool.
      subfaces->traversalinit();
      checksh.sh = shellfacetraverse(subfaces);
      while (checksh.sh != (shellface *) NULL) {
        enqueuesubface(badsubfacs, &checksh);
        checksh.sh = shellfacetraverse(subfaces);
      }

      // Split all encroached subfaces.
      repairencfacs(chkencflag);

      if (b->verbose) {
        printf("  Added %ld (%ld,%ld) Steiner points.\n",  
               points->items-steinercount, st_segref_count-bak_segref_count,
               st_facref_count-bak_facref_count);
      }
    } // if (b->reflevel > 1)
  } // if (!b->nobisect)

  if (b->reflevel > 2) { // '-D3' option (The default option)
    if (b->verbose) {
      printf("  Splitting bad quality tets.\n");
    }

    chkencflag = 4; // Only check tetrahedra.
    steinercount = points->items;
    bak_segref_count = st_segref_count;
    bak_facref_count = st_facref_count;
    bak_volref_count = st_volref_count;

    // The cosine value of the min dihedral angle (-qq) for tetrahedra.
    cosmindihed = cos(b->mindihedral / 180.0 * PI);

    // Initialize the pool of bad quality tetrahedra.
    badtetrahedrons = new memorypool(sizeof(triface), b->tetrahedraperblock,
                                     sizeof(void *), 0);
    // Add all tetrahedra (no hull tets) into the pool.
    tetrahedrons->traversalinit();
    checktet.tet = tetrahedrontraverse();
    while (checktet.tet != NULL) {
      enqueuetetrahedron(&checktet);
      checktet.tet = tetrahedrontraverse();
    }

    // Split all bad quality tetrahedra.
    repairbadtets(chkencflag);

    if (b->verbose) {
      printf("  Added %ld (%ld,%ld,%ld) Steiner points.\n", 
             points->items - steinercount, 
             st_segref_count - bak_segref_count,
             st_facref_count - bak_facref_count,
             st_volref_count - bak_volref_count);
    }
  } // if (b->reflevel > 2)

  if (b->verbose) {
    if (flip23count + flip32count + flip44count > bak_flipcount) {
      printf("  Performed %ld flips.\n", flip23count + flip32count +
             flip44count - bak_flipcount);
    }
  }

  if (steinerleft == 0) {
    if (!b->quiet) {
      printf("\nWarnning:  ");
      printf("The desired number of Steiner points (%d) is reached.\n\n",
             b->steinerleft);
    }
  }


  delete encseglist;
  delete encshlist;

  //if (!b->nobisect) {
  if (!b->nobisect || checkconstraints) {
    totalworkmemory += (badsubsegs->maxitems * badsubsegs->itembytes);
    delete badsubsegs;
    if (b->reflevel > 1) {
      totalworkmemory += (badsubfacs->maxitems * badsubfacs->itembytes);
      delete badsubfacs;
    }
  }
  if (b->reflevel > 2) {
    totalworkmemory += (badtetrahedrons->maxitems*badtetrahedrons->itembytes);
    delete badtetrahedrons;
  }
}

////                                                                       ////
////                                                                       ////
//// refine_cxx ///////////////////////////////////////////////////////////////

