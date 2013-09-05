#include "../tetgen.h"
//// surface_cxx //////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// calculateabovepoint()    Calculate a point above a facet in 'dummypoint'. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::calculateabovepoint(arraypool *facpoints, point *ppa,
                                     point *ppb, point *ppc)
{
  point *ppt, pa, pb, pc;
  REAL v1[3], v2[3], n[3];
  REAL lab, len, A, area;
  REAL x, y, z;
  int i;

  ppt = (point *) fastlookup(facpoints, 0);
  pa = *ppt; // a is the first point.
  pb = pc = NULL; // Avoid compiler warnings.

  // Get a point b s.t. the length of [a, b] is maximal.
  lab = 0;
  for (i = 1; i < facpoints->objects; i++) {
    ppt = (point *) fastlookup(facpoints, i);
    x = (*ppt)[0] - pa[0];
    y = (*ppt)[1] - pa[1];
    z = (*ppt)[2] - pa[2];
    len = x * x + y * y + z * z;
    if (len > lab) {
      lab = len;
      pb = *ppt;
    }
  }
  lab = sqrt(lab);
  if (lab == 0) {
    if (!b->quiet) {
      printf("Warning:  All points of a facet are coincident with %d.\n",
        pointmark(pa));
    }
    return false;
  }

  // Get a point c s.t. the area of [a, b, c] is maximal.
  v1[0] = pb[0] - pa[0];
  v1[1] = pb[1] - pa[1];
  v1[2] = pb[2] - pa[2];
  A = 0;
  for (i = 1; i < facpoints->objects; i++) {
    ppt = (point *) fastlookup(facpoints, i);
    v2[0] = (*ppt)[0] - pa[0];
    v2[1] = (*ppt)[1] - pa[1];
    v2[2] = (*ppt)[2] - pa[2];
    CROSS(v1, v2, n);
    area = DOT(n, n);
    if (area > A) {
      A = area;
      pc = *ppt;
    }
  }
  if (A == 0) {
    // All points are collinear. No above point.
    if (!b->quiet) {
      printf("Warning:  All points of a facet are collinaer with [%d, %d].\n",
        pointmark(pa), pointmark(pb));
    }
    return false;
  }

  // Calculate an above point of this facet.
  facenormal(pa, pb, pc, n, 1, NULL);
  len = sqrt(DOT(n, n));
  n[0] /= len;
  n[1] /= len;
  n[2] /= len;
  lab /= 2.0; // Half the maximal length.
  dummypoint[0] = pa[0] + lab * n[0];
  dummypoint[1] = pa[1] + lab * n[1];
  dummypoint[2] = pa[2] + lab * n[2];

  if (ppa != NULL) {
    // Return the three points.
    *ppa = pa;
    *ppb = pb;
    *ppc = pc;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Calculate an above point. It lies above the plane containing  the subface //
//   [a,b,c], and save it in dummypoint. Moreover, the vector pa->dummypoint //
//   is the normal of the plane.                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::calculateabovepoint4(point pa, point pb, point pc, point pd)
{
  arraypool *ptarray;
  point *parypt;

  ptarray = new arraypool(sizeof(point), 4);

  ptarray->newindex((void **) &parypt);
  *parypt = pa;
  ptarray->newindex((void **) &parypt);
  *parypt = pb;
  ptarray->newindex((void **) &parypt);
  *parypt = pc;
  ptarray->newindex((void **) &parypt);
  *parypt = pd;

  calculateabovepoint(ptarray, NULL, NULL, NULL);

  delete ptarray;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipshpush()    Push a facet edge into flip stack.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flipshpush(face* flipedge)
{
  badface *newflipface;

  newflipface = (badface *) flippool->alloc();
  newflipface->ss = *flipedge;
  newflipface->forg = sorg(*flipedge);
  newflipface->fdest = sdest(*flipedge);
  newflipface->nextitem = flipstack;
  flipstack = newflipface;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip22()    Remove an edge by transforming 2-to-2 subfaces.               //
//                                                                           //
// 'flipfaces' contains two faces: abc and bad. This routine removes these 2 //
// faces and replaces them by two new faces: cdb and dca.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip22(face* flipfaces, int flipflag, int chkencflag)
{
  face bdedges[4], outfaces[4], infaces[4], bdsegs[4];
  face checkface, checkseg;
  point pa, pb, pc, pd;
  badface *bface;
  int i;

  pa = sorg(flipfaces[0]);
  pb = sdest(flipfaces[0]);
  pc = sapex(flipfaces[0]);
  pd = sapex(flipfaces[1]);

  if (sorg(flipfaces[1]) != pb) {
    sesymself(flipfaces[1]);
  }

  if (b->verbose > 3) {
    printf("        flip 2-to-2: (%d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd));
  }
  flip22count++;

  // Collect the four boundary edges.
  senext(flipfaces[0], bdedges[0]);
  senext2(flipfaces[0], bdedges[1]);
  senext(flipfaces[1], bdedges[2]);
  senext2(flipfaces[1], bdedges[3]);

  // Collect outer boundary faces.
  for (i = 0; i < 4; i++) {
    spivot(bdedges[i], outfaces[i]);
    infaces[i] = outfaces[i];
    sspivot(bdedges[i], bdsegs[i]);
    if (outfaces[i].sh != NULL) {
      sspivot(bdedges[i], checkseg);
      if (checkseg.sh != NULL) {
        spivot(infaces[i], checkface);
        while (checkface.sh != bdedges[i].sh) {
          infaces[i] = checkface;
          spivot(infaces[i], checkface);
        }
      }
    }
  }

  // The flags set in these two subfaces do not change.
  // Shellmark does not change.
  // area constraint does not change.

  // Transform abc -> cdb.
  setshvertices(flipfaces[0], pc, pd, pb);
  // Transform bad -> dca.
  setshvertices(flipfaces[1], pd, pc, pa);

  // Update the point-to-subface map.
  if (pointtype(pa) == FREEFACETVERTEX) {
    setpoint2sh(pa, sencode(flipfaces[1]));
  }
  if (pointtype(pb) == FREEFACETVERTEX) {
    setpoint2sh(pb, sencode(flipfaces[0]));
  }
  if (pointtype(pc) == FREEFACETVERTEX) {
    setpoint2sh(pc, sencode(flipfaces[0]));
  }
  if (pointtype(pd) == FREEFACETVERTEX) {
    setpoint2sh(pd, sencode(flipfaces[0]));
  }

  // Reconnect boundary edges to outer boundary faces.
  for (i = 0; i < 4; i++) {
    if (outfaces[(3 + i) % 4].sh != NULL) {
      // Make sure that the subface has the ori as the segment.
      if (bdsegs[(3 + i) % 4].sh != NULL) {
        bdsegs[(3 + i) % 4].shver = 0;
        if (sorg(bdedges[i]) != sorg(bdsegs[(3 + i) % 4])) {
          sesymself(bdedges[i]);
        }
      }
      sbond1(bdedges[i], outfaces[(3 + i) % 4]);
      sbond1(infaces[(3 + i) % 4], bdedges[i]);
    } else {
      sdissolve(bdedges[i]);
    }
    if (bdsegs[(3 + i) % 4].sh != NULL) {
      ssbond(bdedges[i], bdsegs[(3 + i) % 4]);
      if (chkencflag & 1) {
        // Queue this segment for encroaching check.
        if (!smarktest2ed(bdsegs[(3 + i) % 4])) {
          bface = (badface *) badsubsegs->alloc();
          bface->ss = bdsegs[(3 + i) % 4];
          smarktest2(bface->ss); // Only queue it once.
          bface->forg = sorg(bface->ss); // An alive badface.
        }
      }
    } else {
      ssdissolve(bdedges[i]);
    }
  }

  if (chkencflag & 2) {
    // Queue the flipped subfaces for quality/encroaching checks.
    for (i = 0; i < 2; i++) {
      if (!smarktest2ed(flipfaces[i])) {
        bface = (badface *) badsubfacs->alloc();
        bface->ss = flipfaces[i];
        smarktest2(bface->ss); // Only queue it once.
        bface->forg = sorg(bface->ss); // An alive badface.
      }
    }
  }

  recentsh = flipfaces[0];

  if (flipflag) {
    // Put the boundary edges into flip stack.
    for (i = 0; i < 4; i++) {
      flipshpush(&(bdedges[i]));
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip31()    Remove a vertex by transforming 3-to-1 subfaces.              //
//                                                                           //
// 'flipfaces' is an array of subfaces. Its length is at least 4.  On input, //
// the first three faces are: [p,a,b], [p,b,c], and [p,c,a]. This routine    //
// replaces them by one face [a,b,c], it is returned in flipfaces[3].        //
//                                                                           //
// NOTE: The three old subfaces are not deleted within this routine.  They   //
// still hold pointers to their adjacent subfaces. These informations are    //
// needed by the routine 'sremovevertex()' for recovering a segment.         //
// The caller of this routine must delete the old subfaces after their uses. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip31(face* flipfaces, int flipflag)
{
  face bdedges[3], outfaces[3], infaces[3], bdsegs[3];
  face checkface, checkseg;
  point pa, pb, pc, delpt;
  REAL area;
  int i;

  delpt = sorg(flipfaces[0]);
  pa = sdest(flipfaces[0]);
  pb = sdest(flipfaces[1]);
  pc = sdest(flipfaces[2]);

  if (b->verbose > 3) {
    printf("        flip 3-to-1: (%d, %d, %d) - %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(delpt));
  }
  // flip31count++;

  // Collect all infos at the three boundary edges.
  for (i = 0; i < 3; i++) {
    senext(flipfaces[i], bdedges[i]);
    spivot(bdedges[i], outfaces[i]);
    infaces[i] = outfaces[i];
    sspivot(bdedges[i], bdsegs[i]);
    if (outfaces[i].sh != NULL) {
      sspivot(bdedges[i], checkseg);
      if (checkseg.sh != NULL) {
        spivot(infaces[i], checkface);
        while (checkface.sh != bdedges[i].sh) {
          infaces[i] = checkface;
          spivot(infaces[i], checkface);
        }
      }
    }
  } // i

  // Create a new subface.
  makeshellface(subfaces, &(flipfaces[3]));
  setshvertices(flipfaces[3], pa, pb,pc);
  setshellmark(flipfaces[3], shellmark(flipfaces[0]));
  if (checkconstraints) {
    area = areabound(flipfaces[0]);
    setareabound(flipfaces[3], area);
  }

  // Update the point-to-subface map.
  if (pointtype(pa) == FREEFACETVERTEX) {
    setpoint2sh(pa, sencode(flipfaces[3]));
  }
  if (pointtype(pb) == FREEFACETVERTEX) {
    setpoint2sh(pb, sencode(flipfaces[3]));
  }
  if (pointtype(pc) == FREEFACETVERTEX) {
    setpoint2sh(pc, sencode(flipfaces[3]));
  }

  // Update the three new boundary edges.
  bdedges[0] = flipfaces[3];         // [a,b]
  senext(flipfaces[3], bdedges[1]);  // [b,c]
  senext2(flipfaces[3], bdedges[2]); // [c,a]

  // Reconnect boundary edges to outer boundary faces.
  for (i = 0; i < 3; i++) {
    if (outfaces[i].sh != NULL) {
      // Make sure that the subface has the ori as the segment.
      if (bdsegs[i].sh != NULL) {
        bdsegs[i].shver = 0;
        if (sorg(bdedges[i]) != sorg(bdsegs[i])) {
          sesymself(bdedges[i]);
        }
      }
      sbond1(bdedges[i], outfaces[i]);
      sbond1(infaces[i], bdedges[i]);
    }
    if (bdsegs[i].sh != NULL) {
      ssbond(bdedges[i], bdsegs[i]);
    }
  }

  recentsh = flipfaces[3];

  if (flipflag) {
    // Put the boundary edges into flip stack.
    for (i = 0; i < 3; i++) {
      flipshpush(&(bdedges[i]));
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lawsonflip()    Flip non-locally Delaunay edges.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::lawsonflip()
{
  badface *popface;
  face flipfaces[2];
  face checkseg;
  point pa, pb, pc, pd;
  REAL sign;
  long flipcount;

  if (b->verbose > 2) {
    printf("      Lawson flip %ld edges.\n", flippool->items);
  }
  flipcount = flip22count;

  while (flipstack != (badface *) NULL) {

    // Pop an edge from the stack.
    popface = flipstack;
    flipfaces[0] = popface->ss;
    pa = popface->forg;
    pb = popface->fdest;
    flipstack = popface->nextitem; // The next top item in stack.
    flippool->dealloc((void *) popface);

    // Skip it if it is dead.
    if (flipfaces[0].sh[3] == NULL) continue;
    // Skip it if it is not the same edge as we saved.
    if ((sorg(flipfaces[0]) != pa) || (sdest(flipfaces[0]) != pb)) continue;
    // Skip it if it is a subsegment.
    sspivot(flipfaces[0], checkseg);
    if (checkseg.sh != NULL) continue;

    // Get the adjacent face.
    spivot(flipfaces[0], flipfaces[1]);
    if (flipfaces[1].sh == NULL) continue; // Skip a hull edge.
    pc = sapex(flipfaces[0]);
    pd = sapex(flipfaces[1]);

    sign = incircle3d(pa, pb, pc, pd);

    if (sign < 0) {
      // It is non-locally Delaunay. Flip it.
      flip22(flipfaces, 1, 0);
    }
  }

  if (b->verbose > 2) {
    printf("      %ld edges stacked, %ld flips.\n", flippool->items,
      flip22count - flipcount);
  }
  assert(flippool->items == 0l); // SELF_CHECK

  return flip22count - flipcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// sinsertvertex()    Insert a vertex into a triangulation of a facet.       //
//                                                                           //
// This function uses three global arrays: 'caveshlist', 'caveshbdlist', and //
// 'caveshseglist'. On return, 'caveshlist' contains old subfaces in C(p),   //
// 'caveshbdlist' contains new subfaces in C(p). If the new point lies on a  //
// segment, 'cavesegshlist' returns the two new subsegments.                 //
//                                                                           //
// NOTE: the old subfaces in C(p) are not deleted. Theyare needed in case we //
// want to remove the new point immedately.                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::sinsertvertex(point insertpt, face *searchsh, face *splitseg,
                              int iloc, int bowywat)
{
  triface adjtet;
  face cavesh, neighsh, *parysh;
  face newsh, casout, casin;
  face aseg, bseg, aoutseg, boutseg;
  face checkseg;
  point pa, pb, pc;
  enum locateresult loc;
  REAL sign, ori, area;
  int i, j;

  if (b->verbose > 2) {
    printf("      Insert facet point %d.\n", pointmark(insertpt));
  }

  if (splitseg != NULL) {
    // A segment is going to be split, no point location.
    spivot(*splitseg, *searchsh);
    loc = ONEDGE;
  } else {
    loc = (enum locateresult) iloc;
    if (loc == OUTSIDE) {
      // Do point location in surface mesh.
      if (searchsh->sh == NULL) {
        *searchsh = recentsh;
      }
      // Start searching from 'searchsh'.
      loc = slocate(insertpt, searchsh, 1, 1, 0);
    }
  }

  if (b->verbose > 2) {
    if (searchsh->sh != NULL) {
      pa = sorg(*searchsh);
      pb = sdest(*searchsh);
      pc = sapex(*searchsh);
      printf("      Located subface (%d, %d, %d).\n", pointmark(pa), 
             pointmark(pb), pointmark(pc));
    } else {
      assert(splitseg != NULL);
      pa = sorg(*splitseg);
      pb = sdest(*splitseg);
      printf("      Located segment (%d, %d).\n", pointmark(pa),pointmark(pb));
    }
  }

if (bowywat < 3) { // if (bowywat == 1 || bowywat == 2) {

  // Form the initial sC(p).
  if (loc == ONFACE) {
    if (b->verbose > 2) {
      printf("      Inside face.\n");
    }
    // Add the face into list (in B-W cavity).
    smarktest(*searchsh);
    caveshlist->newindex((void **) &parysh);
    *parysh = *searchsh;
  } else if (loc == ONEDGE) {
    if (b->verbose > 2) {
      printf("      On edge.\n");
    }
    if (splitseg != NULL) {
      splitseg->shver = 0;
      pa = sorg(*splitseg);
    } else {
      pa = sorg(*searchsh);
    }
    if (searchsh->sh != NULL) {
      // Collect all subfaces share at this edge.
      neighsh = *searchsh;
      while (1) {
        // Adjust the origin of its edge to be 'pa'.
        if (sorg(neighsh) != pa) {
          sesymself(neighsh);
        }
        // Add this face into list (in B-W cavity).
        smarktest(neighsh);
        caveshlist->newindex((void **) &parysh);
        *parysh = neighsh;
        // Add this face into face-at-splitedge list.
        cavesegshlist->newindex((void **) &parysh);
        *parysh = neighsh;
        // Go to the next face at the edge.
        spivotself(neighsh);
        // Stop if all faces at the edge have been visited.
        if (neighsh.sh == searchsh->sh) break;
        if (neighsh.sh == NULL) break;
      }
    } // If (not a non-dangling segment).
  } else if (loc == ONVERTEX) {
    if (b->verbose > 2) {
      printf("      On vertex.\n");
    }
    return (int) loc;
  } else if (loc == OUTSIDE) {
    // Comment: This should only happen during the surface meshing step.
    // Enlarge the convex hull of the triangulation by including p.
    // An above point of the facet is set in 'dummypoint' to replace
    // orient2d tests by orient3d tests.
    if (b->verbose > 2) {
      printf("      Outside face.\n");
    }
    // Imagine that the current edge a->b (in 'searchsh') is horizontal in a
    //   plane, and a->b is directed from left to right, p lies above a->b.  
    //   Find the right-most edge of the triangulation which is visible by p.
    neighsh = *searchsh;
    while (1) {
      senext2self(neighsh);
      spivot(neighsh, casout);
      if (casout.sh == NULL) {
        // A convex hull edge. Is it visible by p.
        pa = sorg(neighsh);
        pb = sdest(neighsh);
        ori = orient3d(pa, pb, dummypoint, insertpt);
        if (ori < 0) {
          *searchsh = neighsh; // Visible, update 'searchsh'.
        } else {
          break; // 'searchsh' is the right-most visible edge.
        }
      } else {
        if (sorg(casout) != sdest(neighsh)) sesymself(casout);
        neighsh = casout;
      }
    }
    // Create new triangles for all visible edges of p (from right to left).
    casin.sh = NULL;  // No adjacent face at right.
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    while (1) {
      // Create a new subface on top of the (visible) edge.
      makeshellface(subfaces, &newsh); 
      setshvertices(newsh, pb, pa, insertpt);
      setshellmark(newsh, shellmark(*searchsh));
      if (checkconstraints) {
        area = areabound(*searchsh);
        setareabound(newsh, area);
      }
      // Connect the new subface to the bottom subfaces.
      sbond1(newsh, *searchsh);
      sbond1(*searchsh, newsh);
      // Connect the new subface to its right-adjacent subface.
      if (casin.sh != NULL) {
        senext(newsh, casout);
        sbond1(casout, casin);
        sbond1(casin, casout);
      }
      // The left-adjacent subface has not been created yet.
      senext2(newsh, casin);
      // Add the new face into list (inside the B-W cavity).
      smarktest(newsh);
      caveshlist->newindex((void **) &parysh);
      *parysh = newsh;
      // Move to the convex hull edge at the left of 'searchsh'.
      neighsh = *searchsh;
      while (1) {
        senextself(neighsh);
        spivot(neighsh, casout);
        if (casout.sh == NULL) {
          *searchsh = neighsh;
          break;
        }
        if (sorg(casout) != sdest(neighsh)) sesymself(casout);
        neighsh = casout;
      }
      // A convex hull edge. Is it visible by p.
      pa = sorg(*searchsh);
      pb = sdest(*searchsh);
      ori = orient3d(pa, pb, dummypoint, insertpt);
      // Finish the process if p is not visible by the hull edge.
      if (ori >= 0) break;
    }
  }

} else {

  // Under this case, the sub-cavity sC(p) has already been formed in
  //   insertvertex().  Check it.
  // FOR DEBUG ONLY.
  for (i = 0; i < caveshlist->objects; i++) {
    cavesh = * (face *) fastlookup(caveshlist, i);
    assert(smarktested(cavesh));
  }
  if (splitseg != NULL) {
    assert(smarktested(*splitseg));
  }


}// if (bowywat < 3) 

  // Form the Bowyer-Watson cavity sC(p).
  for (i = 0; i < caveshlist->objects; i++) {
    cavesh = * (face *) fastlookup(caveshlist, i);
    for (j = 0; j < 3; j++) {
      sspivot(cavesh, checkseg);
      if (checkseg.sh == NULL) {
        spivot(cavesh, neighsh);
        if (neighsh.sh != NULL) {
          // The adjacent face exists.
          if (!smarktested(neighsh)) {
            if (bowywat) {
              if (bowywat > 2) {
                // It must be a boundary edge.
                sign = 1;
              } else {
                // Check if this subface is connected to adjacent tet(s).
                stpivot(neighsh, adjtet);
                if (adjtet.tet == NULL) {
                  // Check if the subface is non-Delaunay wrt. the new pt.
                  pa = sorg(neighsh);
                  pb = sdest(neighsh);
                  pc = sapex(neighsh);
                  sign = incircle3d(pa, pb, pc, insertpt);
                } else {
                  // It is connected to an adjacent tet. A boundary edge.
                  sign = 1;
                }
              }
              if (sign < 0) {
                // Add the adjacent face in list (in B-W cavity).
                smarktest(neighsh);
                caveshlist->newindex((void **) &parysh);
                *parysh = neighsh;
              }
            } else {
              sign = 1; // A boundary edge.
            }
          } else {
            sign = -1; // Not a boundary edge.
          }
        } else {
          // No adjacent face. It is a hull edge.
          if (loc == OUTSIDE) {
            // It is a boundary edge if it does not contain p.
            if ((sorg(cavesh) == insertpt) || (sdest(cavesh) == insertpt)) {
              sign = -1; // Not a boundary edge.
            } else {
              sign = 1; // A boundary edge.
            }
          } else {
            sign = 1; // A boundary edge.
          }
        }
      } else {
        // Do not across a segment. It is a boundary edge.
        sign = 1;
      }
      if (sign >= 0) {
        // Add a boundary edge.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = cavesh;
      }
      senextself(cavesh);
    } // j
  } // i

  if (b->verbose > 3) {
    printf("        Size of cavity: %ld faces, %ld bdry edges.\n",
           caveshlist->objects, caveshbdlist->objects);
  }

  // Creating new subfaces.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    sspivot(*parysh, checkseg);
    if ((parysh->shver & 01) != 0) sesymself(*parysh);
    pa = sorg(*parysh);
    pb = sdest(*parysh);
    // Create a new subface.
    makeshellface(subfaces, &newsh); 
    setshvertices(newsh, pa, pb, insertpt);
    setshellmark(newsh, shellmark(*parysh));
    setshelltype(newsh, shelltype(*parysh));
    if (checkconstraints) {
      area = areabound(*parysh);
      setareabound(newsh, area);
    }
    // Update the point-to-subface map.
    if (pointtype(pa) == FREEFACETVERTEX) {
      setpoint2sh(pa, sencode(newsh));
    }
    if (pointtype(pb) == FREEFACETVERTEX) {
      setpoint2sh(pb, sencode(newsh));
    }
    // Connect newsh to outer subfaces.
    spivot(*parysh, casout);
    if (casout.sh != NULL) {
      casin = casout;
      if (checkseg.sh != NULL) {
        // Make sure that newsh has the right ori at this segment.
        checkseg.shver = 0;
        if (sorg(newsh) != sorg(checkseg)) {
          sesymself(newsh);
          sesymself(*parysh); // This side should also be inversed.
        }
        spivot(casin, neighsh);
        while (neighsh.sh != parysh->sh) {
          casin = neighsh;
          spivot(casin, neighsh);
        }
      }
      sbond1(newsh, casout);
      sbond1(casin, newsh);
    }
    if (checkseg.sh != NULL) {
      ssbond(newsh, checkseg);
    }
    // Connect oldsh <== newsh (for connecting adjacent new subfaces).
    //   *parysh and newsh point to the same edge and the same ori.
    sbond1(*parysh, newsh);
  }

  // Set a handle for searching.
  recentsh = newsh;

  // Update the point-to-subface map.
  if (pointtype(insertpt) == FREEFACETVERTEX) {
    setpoint2sh(insertpt, sencode(newsh));
  }

  // Connect adjacent new subfaces together.
  for (i = 0; i < caveshbdlist->objects; i++) {
    // Get an old subface at edge [a, b].
    parysh = (face *) fastlookup(caveshbdlist, i);
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senextself(newsh); // At edge [b, p].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [b, p].
      pb = sdest(*parysh);
      neighsh = *parysh;
      while (1) {
        senextself(neighsh);
        spivotself(neighsh);
        if (neighsh.sh == NULL) break;
        if (!smarktested(neighsh)) break;
        if (sdest(neighsh) != pb) sesymself(neighsh);
      }
      if (neighsh.sh != NULL) {
        // Now 'neighsh' is a new subface at edge [b, #].
        if (sorg(neighsh) != pb) sesymself(neighsh);
        assert(sorg(neighsh) == pb); // SELF_CHECK
        assert(sapex(neighsh) == insertpt); // SELF_CHECK
        senext2self(neighsh); // Go to the open edge [p, b].
        sbond(newsh, neighsh);
      } else {
        // There is no adjacent new face at this side.
        assert(loc == OUTSIDE); // SELF_CHECK
      }
    }
    spivot(*parysh, newsh); // The new subface [a, b, p].
    senext2self(newsh); // At edge [p, a].
    spivot(newsh, neighsh);
    if (neighsh.sh == NULL) {
      // Find the adjacent new subface at edge [p, a].
      pa = sorg(*parysh);
      neighsh = *parysh;
      while (1) {
        senext2self(neighsh);
        spivotself(neighsh);
        if (neighsh.sh == NULL) break;
        if (!smarktested(neighsh)) break;
        if (sorg(neighsh) != pa) sesymself(neighsh);
      }
      if (neighsh.sh != NULL) {
        // Now 'neighsh' is a new subface at edge [#, a].
        if (sdest(neighsh) != pa) sesymself(neighsh);
        assert(sdest(neighsh) == pa); // SELF_CHECK
        assert(sapex(neighsh) == insertpt); // SELF_CHECK
        senextself(neighsh); // Go to the open edge [a, p].
        sbond(newsh, neighsh);
      } else {
        // There is no adjacent new face at this side.
        assert(loc == OUTSIDE); // SELF_CHECK
      }
    }
  }

  if (loc == ONEDGE) {

    // An edge is being split. We distinguish two cases:
    //   (1) the edge is not on the boundary of the cavity;
    //   (2) the edge is on the boundary of the cavity.
    // In case (2), the edge is either a segment or a hull edge. There are
    //   degenerated new faces in the cavity. They must be removed.
    for (i = 0; i < cavesegshlist->objects; i++) {
      // Get the saved old subface.
      parysh = (face *) fastlookup(cavesegshlist, i);
      // Get a possible new degenerated subface.
      spivot(*parysh, cavesh);
      if (sapex(cavesh) == insertpt) {
        // Found a degenerated new subface, i.e., case (2).
        if (cavesegshlist->objects > 1) {
          // There are more than one subface share at this edge.
          j = (i + 1) % (int) cavesegshlist->objects;
          parysh = (face *) fastlookup(cavesegshlist, j);
          spivot(*parysh, neighsh);
          // Adjust cavesh and neighsh both at edge a->b, and has p as apex.
          if (sorg(neighsh) != sorg(cavesh)) {
            sesymself(neighsh);
            assert(sorg(neighsh) == sorg(cavesh)); // SELF_CHECK
          }
          assert(sapex(neighsh) == insertpt); // SELF_CHECK
          // Connect adjacent faces at two other edges of cavesh and neighsh.
          //   As a result, the two degenrated new faces are squessed from the
          //   new triangulation of the cavity. Note that the squeezed faces
          //   still hold the adjacent informations which will be used in 
          //   re-connecting subsegments (if they exist). 
          for (j = 0; j < 2; j++) { 
            senextself(cavesh);
            senextself(neighsh);
            spivot(cavesh, newsh);
            spivot(neighsh, casout);
            sbond1(newsh, casout); // newsh <- casout.
          }
        } else {
          // There is only one subface containing this edge [a,b]. Squeese the
          //   degenerated new face [a,b,c] by disconnecting it from its two 
          //   adjacent subfaces at edges [b,c] and [c,a]. Note that the face
          //   [a,b,c] still hold the connection to them.
          for (j = 0; j < 2; j++) {
            senextself(cavesh);
            spivot(cavesh, newsh);
            sdissolve(newsh);
          }
        }
        recentsh = newsh;
        // Update the point-to-subface map.
        if (pointtype(insertpt) == FREEFACETVERTEX) {
          setpoint2sh(insertpt, sencode(newsh));
        }
      }
    }

    if (splitseg != NULL) {
      if (bowywat < 3) {
        smarktest(*splitseg); // Mark it as being processed.
      }
      
      aseg = *splitseg;
      pa = sorg(*splitseg);
      pb = sdest(*splitseg);
      if (b->verbose > 2) {
        printf("      Split seg (%d, %d) by %d.\n", pointmark(pa), 
               pointmark(pb), pointmark(insertpt));
      }

      // Insert the new point p.
      makeshellface(subsegs, &aseg);
      makeshellface(subsegs, &bseg);

      setshvertices(aseg, pa, insertpt, NULL);
      setshvertices(bseg, insertpt, pb, NULL);
      setshellmark(aseg, shellmark(*splitseg));
      setshellmark(bseg, shellmark(*splitseg));
      setshelltype(aseg, shelltype(*splitseg));
      setshelltype(bseg, shelltype(*splitseg));
      if (checkconstraints) {
        setareabound(aseg, areabound(*splitseg));
        setareabound(bseg, areabound(*splitseg));
      }

      // Connect [#, a]<->[a, p].
      senext2(*splitseg, boutseg); // Temporarily use boutseg.
      spivotself(boutseg);
      if (boutseg.sh != NULL) {
        senext2(aseg, aoutseg);
        sbond(boutseg, aoutseg);
      }
      // Connect [p, b]<->[b, #].
      senext(*splitseg, aoutseg);
      spivotself(aoutseg);
      if (aoutseg.sh != NULL) {
        senext(bseg, boutseg);
        sbond(boutseg, aoutseg);
      }
      // Connect [a, p] <-> [p, b].
      senext(aseg, aoutseg);
      senext2(bseg, boutseg);
      sbond(aoutseg, boutseg);

      // Connect subsegs [a, p] and [p, b] to adjacent new subfaces.
      // Although the degenerated new faces have been squeesed. They still
      //   hold the connections to the actual new faces. 
      for (i = 0; i < cavesegshlist->objects; i++) {        
        parysh = (face *) fastlookup(cavesegshlist, i);
        spivot(*parysh, neighsh);
        // neighsh is a degenerated new face.
        if (sorg(neighsh) != pa) {
          sesymself(neighsh);
        }
        senext2(neighsh, newsh);
        spivotself(newsh); // The edge [p, a] in newsh
        ssbond(newsh, aseg);
        senext(neighsh, newsh);
        spivotself(newsh); // The edge [b, p] in newsh
        ssbond(newsh, bseg);
      }


      // Let the point remember the segment it lies on.
      setpoint2sh(insertpt, sencode(aseg));
      // Update the point-to-seg map.
      setpoint2sh(pa, sencode(aseg));
      setpoint2sh(pb, sencode(bseg));
    } // if (splitseg != NULL)

    // Delete all degenerated new faces.
    for (i = 0; i < cavesegshlist->objects; i++) {
      parysh = (face *) fastlookup(cavesegshlist, i);
      spivotself(*parysh);
      if (sapex(*parysh) == insertpt) {
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
    cavesegshlist->restart();

    if (splitseg != NULL) {
      // Return the two new subsegments (for further process).
      //   Re-use 'cavesegshlist'.
      cavesegshlist->newindex((void **) &parysh);
      *parysh = aseg;
      cavesegshlist->newindex((void **) &parysh);
      *parysh = bseg;
    }

  } // if (loc == ONEDGE)


  return (int) loc;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// sremovevertex()    Remove a vertex from the surface mesh.                 //
//                                                                           //
// 'delpt' (p) is the vertex to be removed. If 'parentseg' is not NULL, p is //
// a segment vertex, and the origin of 'parentseg' is p. Otherwise, p is a   //
// facet vertex, and the origin of 'parentsh' is p.                          //
//                                                                           //
// If 'lawson' > 0, the Lawson flip algorithm is used to recover Delaunay-   //
// ness after p is removed.                                                  //
//                                                                           //
// Within each facet, we first use a sequence of 2-to-2 flips to flip any    //
// edge at p, finally use a 3-to-1 flip to remove p.                         //
//                                                                           //
// All new created subfaces are returned in the global array 'caveshbdlist'. //
// The new segment (when p is on segment) is returned in 'parentseg'.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::sremovevertex(point delpt, face* parentsh, face* parentseg,
                              int lawson)
{
  face flipfaces[4], *parysh;
  face spinsh, startsh, neighsh, nextsh, fakesh;
  face abseg, prevseg, checkseg;
  face adjseg1, adjseg2;
  point pa, pb, pc, pd;
  int it, i, j;

  REAL *norm, n1[3], n2[3];
  REAL len, len1, len2;
  REAL ori1, ori2;

  if (parentseg != NULL) {
    assert(sorg(*parentseg) == delpt);
    assert(parentseg->shver == 0);
    // 'delpt' (p) should be a Steiner point inserted in a segment [a,b],
    //   where 'parentseg' should be [p,b]. Find the segment [a,p].
    senext2(*parentseg, prevseg);
    spivotself(prevseg);
    assert(prevseg.sh != NULL);
    prevseg.shver = 0;
    assert(sdest(prevseg) == delpt);
    // Restore the original segment [a,b].
    pa = sorg(prevseg);
    pb = sdest(*parentseg);
    if (b->verbose > 2) {
      printf("      Remove vertex %d from segment [%d, %d].\n", 
             pointmark(delpt), pointmark(pa), pointmark(pb));
    }
    makeshellface(subsegs, &abseg);
    setshvertices(abseg, pa, pb, NULL);
    setshellmark(abseg, shellmark(*parentseg));
    setshelltype(abseg, shelltype(*parentseg));
    if (checkconstraints) {
      setareabound(abseg, areabound(*parentseg));
    }
    // Connect [#, a]<->[a, b].
    senext2(prevseg, adjseg1);
    spivotself(adjseg1);
    if (adjseg1.sh != NULL) {
      adjseg1.shver = 0;
      assert(sdest(adjseg1) == pa);
      senextself(adjseg1);
      senext2(abseg, adjseg2);
      sbond(adjseg1, adjseg2);
    }
    // Connect [a, b]<->[b, #].
    senext(*parentseg, adjseg1);
    spivotself(adjseg1);
    if (adjseg1.sh != NULL) {
      adjseg1.shver = 0;
      assert(sorg(adjseg1) == pb);
      senext2self(adjseg1);
      senext(abseg, adjseg2);
      sbond(adjseg1, adjseg2);
    }
    // Update the point-to-segment map.
    setpoint2sh(pa, sencode(abseg));
    setpoint2sh(pb, sencode(abseg));

    // Get the faces in face ring at segment [p, b].
    //   Re-use array 'caveshlist'.
    spivot(*parentseg, *parentsh);
    spinsh = *parentsh;
    while (1) {
      // Save this face in list.
      caveshlist->newindex((void **) &parysh);
      *parysh = spinsh;
      // Go to the next face in the ring.
      spivotself(spinsh);
      if (spinsh.sh == NULL) break;
      if (spinsh.sh == parentsh->sh) break;
    } 

    // Create the face ring of the new segment [a,b]. Each face in the ring
    //   is [a,b,p] (degenerated!). It will be removed (automatically).
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      startsh = *parysh;
      if (sorg(startsh) != delpt) {
        sesymself(startsh);
        assert(sorg(startsh) == delpt);
      }      
      // startsh is [p, b, #1], find the subface [a, p, #2].
      neighsh = startsh;
      while (1) {
        senext2self(neighsh);
        sspivot(neighsh, checkseg);
        if (checkseg.sh != NULL) {
          // It must be the segment [a, p].
          assert(checkseg.sh == prevseg.sh);
          break;
        }
        spivotself(neighsh);
        assert(neighsh.sh != NULL);
        if (sorg(neighsh) != delpt) sesymself(neighsh);
      }
      // Now neighsh is [a, p, #2].
      if (neighsh.sh != startsh.sh) {
        // Detach the two subsegments [a,p] and [p,b] from subfaces.
        ssdissolve(startsh);
        ssdissolve(neighsh);
        // Create a degenerated subface [a,b,p]. It is used to: (1) hold the
        //   new segment [a,b]; (2) connect to the two adjacent subfaces
        //   [p,b,#] and [a,p,#].
        makeshellface(subfaces, &fakesh);
        setshvertices(fakesh, pa, pb, delpt);
        setshellmark(fakesh, shellmark(startsh));
        // Connect fakesh to the segment [a,b].
        ssbond(fakesh, abseg);
        // Connect fakesh to adjacent subfaces: [p,b,#1] and [a,p,#2].
        senext(fakesh, nextsh);
        sbond(nextsh, startsh);
        senext2(fakesh, nextsh);
        sbond(nextsh, neighsh);
        smarktest(fakesh); // Mark it as faked.
      } else {
        // Special case. There exists already a degenerated face [a,b,p]!
        //   There is no need to create a faked subface here.
        senext2self(neighsh); // [a,b,p]
        assert(sapex(neighsh) == delpt);
        // Since we will re-connect the face ring using the faked subfaces.
        //   We put the adjacent face of [a,b,p] to the list.
        spivot(neighsh, startsh); // The original adjacent subface.
        if (sorg(startsh) != pa) {
          sesymself(startsh);
        }
        assert(sorg(startsh) == pa);
        assert(sdest(startsh) == pb);
        assert(sapex(startsh) != delpt);
        sdissolve(startsh);
        // Connect fakesh to the segment [a,b].
        ssbond(startsh, abseg);
        fakesh = startsh; // Do not mark it!
        // Delete the degenerated subface.
        shellfacedealloc(subfaces, neighsh.sh);
      }
      // Save the fakesh in list (for re-creating the face ring).
      cavesegshlist->newindex((void **) &parysh);
      *parysh = fakesh;
    } // i
    caveshlist->restart();

    // Re-create the face ring.
    if (cavesegshlist->objects > 1) {
      for (i = 0; i < cavesegshlist->objects; i++) {
        parysh = (face *) fastlookup(cavesegshlist, i);
        fakesh = *parysh;
        // Get the next face in the ring.
        j = (i + 1) % cavesegshlist->objects;
        parysh = (face *) fastlookup(cavesegshlist, j);
        nextsh = *parysh;
        sbond1(fakesh, nextsh);
      }
    }

    // Delete the two subsegments containing p.
    shellfacedealloc(subsegs, parentseg->sh);
    shellfacedealloc(subsegs, prevseg.sh);
    // Return the new segment.
    *parentseg = abseg;
  } else {
    // p is inside the surface.
    if (b->verbose > 2) {
      printf("      Remove vertex %d from surface.\n", pointmark(delpt));
    }
    assert(sorg(*parentsh) == delpt);
    // Let 'delpt' be its apex.
    senextself(*parentsh);
    // For unifying the code, we add parentsh to list.
    cavesegshlist->newindex((void **) &parysh);
    *parysh = *parentsh;
  }

  // Remove the point (p).

  for (it = 0; it < cavesegshlist->objects; it++) {
    parentsh = (face *) fastlookup(cavesegshlist, it); // [a,b,p]
    senextself(*parentsh); // [b,p,a].
    spivotself(*parentsh);
    if (sorg(*parentsh) != delpt) {
      sesymself(*parentsh);
    } 
    // now parentsh is [p,b,#].
    if (sorg(*parentsh) != delpt) {
      // The vertex has already been removed in above special case.
      assert(!smarktested(*parentsh));
      continue;
    }

    while (1) {      
      // Initialize the flip edge list. Re-use 'caveshlist'.
      spinsh = *parentsh; // [p, b, #]
      while (1) {
        caveshlist->newindex((void **) &parysh);
        *parysh = spinsh;
        senext2self(spinsh);
        spivotself(spinsh);
        assert(spinsh.sh != NULL);
        if (spinsh.sh == parentsh->sh) break;
        if (sorg(spinsh) != delpt) {
          sesymself(spinsh);
          assert(sorg(spinsh) == delpt);
        }
      } // while (1)

      if (caveshlist->objects == 3) {
        // Delete the point by a 3-to-1 flip.
        for (i = 0; i < 3; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          flipfaces[i] = *parysh;
        }
        flip31(flipfaces, lawson);
        for (i = 0; i < 3; i++) { 
          shellfacedealloc(subfaces, flipfaces[i].sh);
        }
        caveshlist->restart();
        // Save the new subface.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = flipfaces[3];
        // The vertex is removed.
        break;
      } else {
        // There should be more than 3 subfaces in list.
        assert(caveshlist->objects > 3);
      }

      // Search an edge to flip.
      for (i = 0; i < caveshlist->objects; i++) {
        parysh = (face *) fastlookup(caveshlist, i);
        flipfaces[0] = *parysh;
        spivot(flipfaces[0], flipfaces[1]);
        if (sorg(flipfaces[0]) != sdest(flipfaces[1])) {
          sesymself(flipfaces[1]);
        }
        // Skip this edge if it belongs to a faked subface.
        if (!smarktested(flipfaces[0]) && !smarktested(flipfaces[1])) {
          pa = sorg(flipfaces[0]);
          pb = sdest(flipfaces[0]);
          pc = sapex(flipfaces[0]);
          pd = sapex(flipfaces[1]);
          // Select a base.
          facenormal(pa, pb, pc, n1, 1, NULL);
          len1 = sqrt(DOT(n1, n1));
          facenormal(pa, pb, pd, n2, 1, NULL);
          len2 = sqrt(DOT(n2, n2));
          if (len1 > len2) {
            norm = n1;
            len = len1;
          } else {
            norm = n2;
            len = len2;
          }
          assert(len > 0);
          norm[0] /= len;
          norm[1] /= len;
          norm[2] /= len;
          len = DIST(pa, pb);
          dummypoint[0] = pa[0] + len * norm[0];
          dummypoint[1] = pa[1] + len * norm[1];
          dummypoint[2] = pa[2] + len * norm[2];
          // Check if a 2-to-2 flip is possible.
          ori1 = orient3d(pc, pd, dummypoint, pa);
          ori2 = orient3d(pc, pd, dummypoint, pb);
          if (ori1 * ori2 < 0) {
            // A 2-to-2 flip is found.
            flip22(flipfaces, lawson, 0);
            // The i-th edge is flipped. The i-th and (i-1)-th subfaces are
            //   changed. The 'flipfaces[1]' contains p as its apex.
            senext2(flipfaces[1], *parentsh);
            // Save the new subface.
            caveshbdlist->newindex((void **) &parysh);
            *parysh = flipfaces[0];
            break;
          }
        } //
      } // i
      if (i == caveshlist->objects) {
        // This can happen only if there are 4 edges at p, and they are
        //   orthogonal to each other, see Fig. 2010-11-01.
        assert(caveshlist->objects == 4);
        // Do a flip22 and a flip31 to remove p.
        parysh = (face *) fastlookup(caveshlist, 0);
        flipfaces[0] = *parysh;
        spivot(flipfaces[0], flipfaces[1]);
        if (sorg(flipfaces[0]) != sdest(flipfaces[1])) {
          sesymself(flipfaces[1]);
        }
        flip22(flipfaces, lawson, 0);
        senext2(flipfaces[1], *parentsh);
        // Save the new subface.
        caveshbdlist->newindex((void **) &parysh);
        *parysh = flipfaces[0];
      }
      // The edge list at p are changed.
      caveshlist->restart();
    } // while (1)

  } // it

  cavesegshlist->restart();

  if (b->verbose > 2) {
    printf("      Created %ld new subfaces.\n", caveshbdlist->objects);
  }


  if (lawson) {
    lawsonflip();
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// slocate()    Locate a point in a surface triangulation.                   //
//                                                                           //
// Staring the search from 'searchsh'(it should not be NULL). Perform a line //
// walk search for a subface containing the point (p).                       //
//                                                                           //
// If 'aflag' is set, the 'dummypoint' is pre-calculated so that it lies     //
// above the 'searchsh' in its current orientation. The test if c is CCW to  //
// the line a->b can be done by the test if c is below the oriented plane    //
// a->b->dummypoint.                                                         //
//                                                                           //
// If 'cflag' is not TRUE, the triangulation may not be convex.  Stop search //
// when a segment is met and return OUTSIDE.                                 //
//                                                                           //
// If 'rflag' (rounding) is set, after the location of the point is found,   //
// either ONEDGE or ONFACE, round the result using an epsilon.               //
//                                                                           //
// The returned value inducates the following cases:                         //
//   - ONVERTEX, p is the origin of 'searchsh'.                              //
//   - ONEDGE, p lies on the edge of 'searchsh'.                             //
//   - ONFACE, p lies in the interior of 'searchsh'.                         //
//   - OUTSIDE, p lies outside of the triangulation, p is on the left-hand   //
//     side of the edge 'searchsh'(s), i.e., org(s), dest(s), p are CW.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::locateresult tetgenmesh::slocate(point searchpt, 
  face* searchsh, int aflag, int cflag, int rflag)
{
  face neighsh;
  face checkseg;
  point pa, pb, pc, pd, *parypt;
  enum locateresult loc;
  enum {MOVE_BC, MOVE_CA} nextmove;
  REAL ori, ori_bc, ori_ca;
  REAL dist_bc, dist_ca;
  int i;

  // For finding an approximate location.
  //REAL n[3], len, len3;
  REAL n[3], area_abc, area_abp, area_bcp, area_cap;

  pa = sorg(*searchsh);
  pb = sdest(*searchsh);
  pc = sapex(*searchsh);

  if (!aflag) {
    // No above point is given. Calculate an above point for this facet.
    //   Re-use the 'cavetetvertlist'.
    cavetetvertlist->newindex((void **) &parypt);
    *parypt = pa;
    cavetetvertlist->newindex((void **) &parypt);
    *parypt = pb;
    cavetetvertlist->newindex((void **) &parypt);
    *parypt = pc;
    cavetetvertlist->newindex((void **) &parypt);
    *parypt = searchpt;
    calculateabovepoint(cavetetvertlist, NULL, NULL, NULL);
    cavetetvertlist->restart();
  }

  // 'dummypoint' is given. Make sure it is above [a,b,c]
  ori = orient3d(pa, pb, pc, dummypoint);
  assert(ori != 0); // SELF_CHECK
  if (ori > 0) {
    sesymself(*searchsh); // Reverse the face orientation.
  }

  // Find an edge of the face s.t. p lies on its right-hand side (CCW).
  for (i = 0; i < 3; i++) {
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    ori = orient3d(pa, pb, dummypoint, searchpt);
    if (ori > 0) break;
    senextself(*searchsh);
  }
  assert(i < 3); // SELF_CHECK

  pc = sapex(*searchsh);

  if (pc == searchpt) {
    senext2self(*searchsh);
    return ONVERTEX;
  }

  while (1) {

    ori_bc = orient3d(pb, pc, dummypoint, searchpt);
    ori_ca = orient3d(pc, pa, dummypoint, searchpt);

    if (ori_bc < 0) {
      if (ori_ca < 0) { // (--)
        // Any of the edges is a viable move.
        senext(*searchsh, neighsh); // At edge [b, c].
        spivotself(neighsh);
        if (neighsh.sh != NULL) {
          pd = sapex(neighsh);
          dist_bc = NORM2(searchpt[0] - pd[0], searchpt[1] - pd[1],
            searchpt[2] - pd[2]);
        } else {
          dist_bc = NORM2(xmax - xmin, ymax - ymin, zmax - zmin);
        }
        senext2(*searchsh, neighsh); // At edge [c, a].
        spivotself(neighsh);
        if (neighsh.sh != NULL) {
          pd = sapex(neighsh);
          dist_ca = NORM2(searchpt[0] - pd[0], searchpt[1] - pd[1],
            searchpt[2] - pd[2]);
        } else {
          dist_ca = dist_bc;
        }
        if (dist_ca < dist_bc) {
          nextmove = MOVE_CA;
        } else {
          nextmove = MOVE_BC;
        }
      } else { // (-#)
        // Edge [b, c] is viable.
        nextmove = MOVE_BC;
      }
    } else {
      if (ori_ca < 0) { // (#-)
        // Edge [c, a] is viable.
        nextmove = MOVE_CA;
      } else {
        if (ori_bc > 0) {
          if (ori_ca > 0) { // (++)
            loc = ONFACE;  // Inside [a, b, c].
            break;
          } else { // (+0)
            senext2self(*searchsh); // On edge [c, a].
            loc = ONEDGE;
            break;
          }
        } else { // ori_bc == 0
          if (ori_ca > 0) { // (0+)
            senextself(*searchsh); // On edge [b, c].
            loc = ONEDGE;
            break;
          } else { // (00)
            // p is coincident with vertex c. 
            senext2self(*searchsh);
            return ONVERTEX;
          }
        }
      }
    }

    // Move to the next face.
    if (nextmove == MOVE_BC) {
      senextself(*searchsh);
    } else {
      senext2self(*searchsh);
    }
    if (!cflag) {
      // NON-convex case. Check if we will cross a boundary.
      sspivot(*searchsh, checkseg);
      if (checkseg.sh != NULL) {
        return ENCSEGMENT;
      }
    }
    spivot(*searchsh, neighsh);
    if (neighsh.sh == NULL) {
      return OUTSIDE; // A hull edge.
    }
    // Adjust the edge orientation.
    if (sorg(neighsh) != sdest(*searchsh)) {
      sesymself(neighsh);
    }
    assert(sorg(neighsh) == sdest(*searchsh)); // SELF_CHECK

    // Update the newly discovered face and its endpoints.
    *searchsh = neighsh;
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    pc = sapex(*searchsh);

    if (pc == searchpt) {
      senext2self(*searchsh);
      return ONVERTEX;
    }

  } // while (1)

  // assert(loc == ONFACE || loc == ONEDGE);


  if (rflag) {
    // Round the locate result before return.
    pa = sorg(*searchsh);
    pb = sdest(*searchsh);
    pc = sapex(*searchsh);

    facenormal(pa, pb, pc, n, 1, NULL);
    area_abc = sqrt(dot(n, n));

    facenormal(pb, pc, searchpt, n, 1, NULL);
    area_bcp = sqrt(dot(n, n));
    if ((area_bcp / area_abc) < b->epsilon) {
      area_bcp = 0; // Rounding.
    }

    facenormal(pc, pa, searchpt, n, 1, NULL);
    area_cap = sqrt(dot(n, n));
    if ((area_cap / area_abc) < b->epsilon) {
      area_cap = 0; // Rounding
    }

    if ((loc == ONFACE) || (loc == OUTSIDE)) {
      facenormal(pa, pb, searchpt, n, 1, NULL);
      area_abp = sqrt(dot(n, n));
      if ((area_abp / area_abc) < b->epsilon) {
        area_abp = 0; // Rounding
      }
    } else { // loc == ONEDGE
      area_abp = 0;
    }

    if (area_abp == 0) {
      if (area_bcp == 0) {
        assert(area_cap != 0);
        senextself(*searchsh); 
        loc = ONVERTEX; // p is close to b.
      } else {
        if (area_cap == 0) {
          loc = ONVERTEX; // p is close to a.
        } else {
          loc = ONEDGE; // p is on edge [a,b].
        }
      }
    } else if (area_bcp == 0) {
      if (area_cap == 0) {
        senext2self(*searchsh); 
        loc = ONVERTEX; // p is close to c.
      } else {
        senextself(*searchsh);
        loc = ONEDGE; // p is on edge [b,c].
      }
    } else if (area_cap == 0) {
      senext2self(*searchsh);
      loc = ONEDGE; // p is on edge [c,a].
    } else {
      loc = ONFACE; // p is on face [a,b,c].
    }
  } // if (rflag)

  return loc;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// sscoutsegment()    Look for a segment in surface triangulation.           //
//                                                                           //
// The segment is given by the origin of 'searchsh' and 'endpt'.  Assume the //
// orientation of 'searchsh' is CCW w.r.t. the above point.                  //
//                                                                           //
// If an edge in T is found matching this segment, the segment is "locaked"  //
// in T at the edge.  Otherwise, flip the first edge in T that the segment   //
// crosses. Continue the search from the flipped face.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::interresult 
  tetgenmesh::sscoutsegment(face *searchsh, point endpt)
{
  face flipshs[2], neighsh;
  face newseg, checkseg;
  point startpt, pa, pb, pc, pd;
  enum interresult dir;
  enum {MOVE_AB, MOVE_CA} nextmove;
  REAL ori_ab, ori_ca;
  REAL dist_b, dist_c;
  int shmark = 0;

  // The origin of 'searchsh' is fixed.
  startpt = sorg(*searchsh); // pa = startpt;
  nextmove = MOVE_AB; // Avoid compiler warning.

  if (b->verbose > 2) {
    printf("      Scout segment (%d, %d).\n", pointmark(startpt),
           pointmark(endpt));
  }

  // Search an edge in 'searchsh' on the path of this segment.
  while (1) {

    pb = sdest(*searchsh);
    if (pb == endpt) {
      dir = SHAREEDGE; // Found!
      break;
    }

    pc = sapex(*searchsh);
    if (pc == endpt) {
      senext2self(*searchsh);
      sesymself(*searchsh);
      dir = SHAREEDGE; // Found!
      break;
    }

    ori_ab = orient3d(startpt, pb, dummypoint, endpt);
    ori_ca = orient3d(pc, startpt, dummypoint, endpt);

    if (ori_ab < 0) {
      if (ori_ca < 0) { // (--)
        // Both sides are viable moves.
        spivot(*searchsh, neighsh); // At edge [a, b].
        assert(neighsh.sh != NULL); // SELF_CHECK
        pd = sapex(neighsh);
        dist_b = NORM2(endpt[0] - pd[0], endpt[1] - pd[1], endpt[2] - pd[2]);
        senext2(*searchsh, neighsh); // At edge [c, a].
        spivotself(neighsh);
        assert(neighsh.sh != NULL); // SELF_CHECK
        pd = sapex(neighsh);
        dist_c = NORM2(endpt[0] - pd[0], endpt[1] - pd[1], endpt[2] - pd[2]);
        if (dist_c < dist_b) {
          nextmove = MOVE_CA;
        } else {
          nextmove = MOVE_AB;
        }
      } else { // (-#)
        nextmove = MOVE_AB;
      }
    } else {
      if (ori_ca < 0) { // (#-)
        nextmove = MOVE_CA;
      } else {
        if (ori_ab > 0) {
          if (ori_ca > 0) { // (++)
            // The segment intersects with edge [b, c].
            dir = ACROSSEDGE;
            break;
          } else { // (+0)
            // The segment collinear with edge [c, a].
            senext2self(*searchsh);
            sesymself(*searchsh);
            dir = ACROSSVERT;
            break;
          }
        } else {
          if (ori_ca > 0) { // (0+)
            // The segment collinear with edge [a, b].
            dir = ACROSSVERT;
            break;
          } else { // (00)
            // startpt == endpt. Not possible.
            assert(0); // SELF_CHECK
          }
        }
      }
    }

    // Move 'searchsh' to the next face, keep the origin unchanged.
    if (nextmove == MOVE_AB) {
      spivot(*searchsh, neighsh);
      if (sorg(neighsh) != pb) sesymself(neighsh);
      senext(neighsh, *searchsh);      
    } else {
      senext2(*searchsh, neighsh);
      spivotself(neighsh);
      if (sdest(neighsh) != pc) sesymself(neighsh);
      *searchsh = neighsh;
    }
    assert(sorg(*searchsh) == startpt); // SELF_CHECK

  } // while

  if (dir == SHAREEDGE) {
    // Insert the segment into the triangulation.
    makeshellface(subsegs, &newseg);
    setshvertices(newseg, startpt, endpt, NULL);
    // Set the actual segment marker.
    if (in->facetmarkerlist != NULL) {
      shmark = shellmark(*searchsh);
      setshellmark(newseg, in->facetmarkerlist[shmark - 1]);
    }
    ssbond(*searchsh, newseg);
    spivot(*searchsh, neighsh);
    if (neighsh.sh != NULL) {
      ssbond(neighsh, newseg);
    }
    return dir;
  }

  if (dir == ACROSSVERT) {
    // A point is found collinear with this segment.
    return dir;
  }

  if (dir == ACROSSEDGE) {
    // Edge [b, c] intersects with the segment.
    senext(*searchsh, flipshs[0]);
    sspivot(flipshs[0], checkseg);
    if (checkseg.sh != NULL) {
      printf("Error:  Invalid PLC.\n");
      pb = sorg(flipshs[0]);
      pc = sdest(flipshs[0]);
      printf("  Two segments (%d, %d) and (%d, %d) intersect.\n",
        pointmark(startpt), pointmark(endpt), pointmark(pb), pointmark(pc));
      terminatetetgen(3);
    }
    // Flip edge [b, c], queue unflipped edges (for Delaunay checks).
    spivot(flipshs[0], flipshs[1]);
    assert(flipshs[1].sh != NULL); // SELF_CHECK
    if (sorg(flipshs[1]) != sdest(flipshs[0])) sesymself(flipshs[1]);
    flip22(flipshs, 1, 0);
    // The flip may create an invered triangle, check it.
    pa = sapex(flipshs[1]);
    pb = sapex(flipshs[0]);
    pc = sorg(flipshs[0]);
    pd = sdest(flipshs[0]);
    // Check if pa and pb are on the different sides of [pc, pd]. 
    // Re-use ori_ab, ori_ca for the tests.
    ori_ab = orient3d(pc, pd, dummypoint, pb);
    ori_ca = orient3d(pd, pc, dummypoint, pa);
    //assert(ori_ab * ori_ca != 0); // SELF_CHECK
    if (ori_ab < 0) {
      if (b->verbose > 2) {
        printf("      Queue an inversed triangle (%d, %d, %d) %d\n",
          pointmark(pc), pointmark(pd), pointmark(pb), pointmark(pa));
      }
      flipshpush(&(flipshs[0]));  // push it to 'flipstack'
    } else if (ori_ca < 0) {
      if (b->verbose > 2) {
        printf("      Queue an inversed triangle (%d, %d, %d) %d\n",
          pointmark(pd), pointmark(pc), pointmark(pa), pointmark(pb));
      }
      flipshpush(&(flipshs[1])); // // push it to 'flipstack'
    }
    // Set 'searchsh' s.t. its origin is 'startpt'.
    *searchsh = flipshs[0];
    assert(sorg(*searchsh) == startpt);
  }

  return sscoutsegment(searchsh, endpt);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scarveholes()    Remove triangles not in the facet.                       //
//                                                                           //
// This routine re-uses the two global arrays: caveshlist and caveshbdlist.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::scarveholes(int holes, REAL* holelist)
{
  face *parysh, searchsh, neighsh;
  face checkseg;
  enum locateresult loc;
  int i, j;

  // Get all triangles. Infect unprotected convex hull triangles. 
  smarktest(recentsh);
  caveshlist->newindex((void **) &parysh);
  *parysh = recentsh;
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    searchsh = *parysh;
    searchsh.shver = 0;
    for (j = 0; j < 3; j++) {
      spivot(searchsh, neighsh);
      // Is this side on the convex hull?
      if (neighsh.sh != NULL) {
        if (!smarktested(neighsh)) {
          smarktest(neighsh);
          caveshlist->newindex((void **) &parysh);
          *parysh = neighsh;
        }
      } else {
        // A hull side. Check if it is protected by a segment.
        sspivot(searchsh, checkseg);
        if (checkseg.sh == NULL) {
          // Not protected. Save this face.
          if (!sinfected(searchsh)) {
            sinfect(searchsh);
            caveshbdlist->newindex((void **) &parysh);
            *parysh = searchsh;
          }
        }
      }
      senextself(searchsh);
    }
  }

  // Infect the triangles in the holes.
  for (i = 0; i < 3 * holes; i += 3) {
    searchsh = recentsh;
    loc = slocate(&(holelist[i]), &searchsh, 1, 1, 0);
    if (loc != OUTSIDE) {
      sinfect(searchsh);
      caveshbdlist->newindex((void **) &parysh);
      *parysh = searchsh;
    }
  }

  // Find and infect all exterior triangles.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    searchsh = *parysh;
    searchsh.shver = 0;
    for (j = 0; j < 3; j++) {
      spivot(searchsh, neighsh);
      if (neighsh.sh != NULL) {
        sspivot(searchsh, checkseg);
        if (checkseg.sh == NULL) {
          if (!sinfected(neighsh)) {
            sinfect(neighsh);
            caveshbdlist->newindex((void **) &parysh);
            *parysh = neighsh;
          }
        } else {
          sdissolve(neighsh); // Disconnect a protected face.
        }
      }
      senextself(searchsh);
    }
  }

  // Delete exterior triangles, unmark interior triangles.
  for (i = 0; i < caveshlist->objects; i++) {
    parysh = (face *) fastlookup(caveshlist, i);
    if (sinfected(*parysh)) {
      shellfacedealloc(subfaces, parysh->sh);
    } else {
      sunmarktest(*parysh);
    }
  }

  caveshlist->restart();
  caveshbdlist->restart();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// triangulate()    Create a CDT for the facet.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::triangulate(int shmark, arraypool* ptlist, arraypool* conlist,
                             int holes, REAL* holelist)
{
  face searchsh, newsh, *parysh; 
  face newseg;
  point pa, pb, pc, *ppt, *cons;
  enum locateresult loc;
  int iloc;
  int i, j;

  int idx, fmarker;
  REAL area;

  if (b->verbose > 2) {
    printf("      f%d:  %ld vertices, %ld segments", shmark, ptlist->objects,
           conlist->objects);
    if (holes > 0) {
      printf(", %d holes", holes);
    }
    printf(".\n");
  }

  if (ptlist->objects < 2l) {
    // Not a segment or a facet.
    return;
  } if (ptlist->objects == 2l) {
    pa = * (point *) fastlookup(ptlist, 0);
    pb = * (point *) fastlookup(ptlist, 1);
    if (distance(pa, pb) > 0) {
      // It is a single segment.
      makeshellface(subsegs, &newseg);
      setshvertices(newseg, pa, pb, NULL);
      // Set the actual segment marker.
      if (in->facetmarkerlist != NULL) {
        setshellmark(newseg, in->facetmarkerlist[shmark - 1]);
      }
    }
    if (pointtype(pa) == VOLVERTEX) {
      setpointtype(pa, RIDGEVERTEX);
    }
    if (pointtype(pb) == VOLVERTEX) {
      setpointtype(pb, RIDGEVERTEX);
    }
    return;
  } if (ptlist->objects == 3l) {
    // The facet has only one triangle.
    pa = * (point *) fastlookup(ptlist, 0);
    pb = * (point *) fastlookup(ptlist, 1);
    pc = * (point *) fastlookup(ptlist, 2);
    if (triarea(pa, pb, pc) > 0) {
      makeshellface(subfaces, &newsh);
      setshvertices(newsh, pa, pb, pc);
      setshellmark(newsh, shmark);
      // Create three new segments.
      for (i = 0; i < 3; i++) {
        makeshellface(subsegs, &newseg);
        setshvertices(newseg, sorg(newsh), sdest(newsh), NULL);
        // Set the actual segment marker.
        if (in->facetmarkerlist != NULL) {
          setshellmark(newseg, in->facetmarkerlist[shmark - 1]);
        }
        ssbond(newsh, newseg);
        senextself(newsh);
      }
      if (pointtype(pa) == VOLVERTEX) {
        setpointtype(pa, FACETVERTEX);
      }
      if (pointtype(pb) == VOLVERTEX) {
        setpointtype(pb, FACETVERTEX);
      }
      if (pointtype(pc) == VOLVERTEX) {
        setpointtype(pc, FACETVERTEX);
      }
    }
    return;
  }

  // Calulcate an above point of this facet.
  if (!calculateabovepoint(ptlist, &pa, &pb, &pc)) {
    return; // The point set is degenerate.
  }

  // Create an initial triangulation.
  makeshellface(subfaces, &newsh);
  setshvertices(newsh, pa, pb, pc);
  setshellmark(newsh, shmark);
  recentsh = newsh;

  if (pointtype(pa) == VOLVERTEX) {
    setpointtype(pa, FACETVERTEX);
  }
  if (pointtype(pb) == VOLVERTEX) {
    setpointtype(pb, FACETVERTEX);
  }
  if (pointtype(pc) == VOLVERTEX) {
    setpointtype(pc, FACETVERTEX);
  }

  // Are there area constraints?
  if (b->quality && (in->facetconstraintlist != (REAL *) NULL)) {
    idx = in->facetmarkerlist[shmark - 1]; // The actual facet marker.
    for (i = 0; i < in->numberoffacetconstraints; i++) {
      fmarker = (int) in->facetconstraintlist[i * 2];
      if (fmarker == idx) {
        area = in->facetconstraintlist[i * 2 + 1];
        setareabound(newsh, area);
        break;
      }
    }
  }

  // Incrementally build the triangulation.
  pinfect(pa);
  pinfect(pb);
  pinfect(pc);
  for (i = 0; i < ptlist->objects; i++) {
    ppt = (point *) fastlookup(ptlist, i);
    if (!pinfected(*ppt)) {
      searchsh = recentsh; // Start from 'recentsh'.
      iloc = (int) OUTSIDE;
      if (b->verbose > 2) printf("      # %d", i);
      loc = (enum locateresult) sinsertvertex(*ppt, &searchsh, NULL, iloc, 1);
      assert(loc != ONVERTEX); // SELF_CHECK
      if (pointtype(*ppt) == VOLVERTEX) {
        setpointtype(*ppt, FACETVERTEX);
      }
      // Delete all removed subfaces.
      for (j = 0; j < caveshlist->objects; j++) {
        parysh = (face *) fastlookup(caveshlist, j);
        shellfacedealloc(subfaces, parysh->sh);
      }
      // Clear the global lists.
      caveshbdlist->restart();
      caveshlist->restart();
      cavesegshlist->restart();
    } else {
      puninfect(*ppt); // This point has inserted.
    }
  }

  // Insert the segments.
  for (i = 0; i < conlist->objects; i++) {
    cons = (point *) fastlookup(conlist, i);
    searchsh = recentsh;
    loc = slocate(cons[0], &searchsh, 1, 1, 0);
    assert(loc == ONVERTEX); // SELF_CHECK
    // Recover the segment. Some edges may be flipped.
    sscoutsegment(&searchsh, cons[1]);
    if (flipstack != NULL) {
      // Recover locally Delaunay edges.
      lawsonflip();
    }
  }

  // Remove exterior and hole triangles.
  scarveholes(holes, holelist);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// unifysubfaces()    Unify two identical subfaces.                          //
//                                                                           //
// Two subfaces, f1 [a, b, c] and f2 [a, b, d], share the same edge [a, b].  //
// If c = d, then f1 and f2 are identical. Otherwise, these two subfaces     //
// intersect, and the mesher is stopped.                                     //
//                                                                           //
// If the two subfaces are indentical, we try to replace f2 by f1, i.e, all  //
// neighbors of f2 are re-connected to f1.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::unifysubfaces(face *f1, face *f2)
{
  face casout, casin, neighsh;
  face sseg, checkseg;
  point pa, pb, pc, pd;
  int i;

  pa = sorg(*f1);
  pb = sdest(*f1);
  pc = sapex(*f1);
  pd = sapex(*f2);

  if (pc != pd) {
    printf("Found two facets intersect each other.\n");
    printf("  1st: [%d, %d, %d] #%d\n", 
	   pointmark(pa), pointmark(pb), pointmark(pc), shellmark(*f1));
    printf("  2nd: [%d, %d, %d] #%d\n",
	   pointmark(pa), pointmark(pb), pointmark(pd), shellmark(*f2));
    terminatetetgen(3);
  } else {
    printf("Found two duplicated facets.\n");
    printf("  1st: [%d, %d, %d] #%d\n", 
	   pointmark(pa), pointmark(pb), pointmark(pc), shellmark(*f1));
    printf("  2nd: [%d, %d, %d] #%d\n",
	   pointmark(pa), pointmark(pb), pointmark(pd), shellmark(*f2));
    terminatetetgen(3);
  }

  // f1 and f2 are identical, replace f2 by f1.
  if (!b->quiet) {
    printf("Warning:  Facet #%d is duplicated with Facet #%d. Removed!\n",
           shellmark(*f2), shellmark(*f1));
  }

  // Make possible disconnections/reconnections at neighbors of f2.
  for (i = 0; i < 3; i++) {
    spivot(*f1, casout);
    if (casout.sh == NULL) {
      // f1 has no adjacent subfaces yet.
      spivot(*f2, casout);
      if (casout.sh != NULL) {
        // Re-direct the adjacent connections of f2 to f1.
        casin = casout;
        spivot(casin, neighsh);
        while (neighsh.sh != f2->sh) {
          casin = neighsh;
          spivot(casin, neighsh);
        }
        // Connect casout <= f1 <= casin.
        sbond1(*f1, casout);
        sbond1(casin, *f1);
      }
    }
    sspivot(*f2, sseg); 
    if (sseg.sh != NULL) {
      // f2 has a segment. It must be different to f1's.
      // Disconnect bonds of subfaces to this segment.
      spivot(*f2, casout);
      if (casout.sh != NULL) {
        casin = casout;
        ssdissolve(casin);
        spivot(casin, neighsh);
        while (neighsh.sh != f2->sh) {
          casin = neighsh;
          ssdissolve(casin);
          spivot(casin, neighsh);
        }
      }
      // Delete the segment.
      shellfacedealloc(subsegs, sseg.sh);
    }
    spivot(*f2, casout);
    if (casout.sh != NULL) {
      // Find the subface (casin) pointing to f2.
      casin = casout;
      spivot(casin, neighsh);
      while (neighsh.sh != f2->sh) {
        casin = neighsh;
        spivot(casin, neighsh);
      }
      // Disconnect f2 <= casin.
      sdissolve(casin);
    }
    senextself(*f1);
    senextself(*f2);
  } // i

  // Delete f2.
  shellfacedealloc(subfaces, f2->sh);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// unifysegments()    Remove redundant segments and create face links.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::unifysegments()
{
  badface *facelink = NULL, *newlinkitem, *f1, *f2;
  face *facperverlist, sface;
  face subsegloop, testseg;
  point torg, tdest;
  REAL ori1, ori2, ori3;
  REAL n1[3], n2[3];
  int *idx2faclist;
  int idx, k, m;

  int e1, e2;
  REAL len;

  if (b->verbose > 1) {
    printf("  Unifying segments.\n");
  }

  // Create a mapping from vertices to subfaces.
  makepoint2submap(subfaces, idx2faclist, facperverlist);

  subsegloop.shver = 0;
  subsegs->traversalinit();
  subsegloop.sh = shellfacetraverse(subsegs);
  while (subsegloop.sh != (shellface *) NULL) {
    torg = sorg(subsegloop);
    tdest = sdest(subsegloop);

    idx = pointmark(torg) - in->firstnumber;
    // Loop through the set of subfaces containing 'torg'.  Get all the
    //   subfaces containing the edge (torg, tdest). Save and order them
    //   in 'sfacelist', the ordering is defined by the right-hand rule
    //   with thumb points from torg to tdest.
    for (k = idx2faclist[idx]; k < idx2faclist[idx + 1]; k++) {
      sface = facperverlist[k];
      // The face may be deleted if it is a duplicated face.
      if (sface.sh[3] == NULL) continue;
      // Search the edge torg->tdest.
      assert(sorg(sface) == torg); // SELF_CHECK
      if (sdest(sface) != tdest) {
        senext2self(sface);
        sesymself(sface);
      }
      if (sdest(sface) != tdest) continue;

      // Save the face f in facelink.
      if (flippool->items >= 2) {
        f1 = facelink;
        for (m = 0; m < flippool->items - 1; m++) {
          f2 = f1->nextitem;
          ori1 = orient3d(torg, tdest, sapex(f1->ss), sapex(f2->ss));
          ori2 = orient3d(torg, tdest, sapex(f1->ss), sapex(sface));
          if (ori1 > 0) {
            // apex(f2) is below f1.
            if (ori2 > 0) {
              // apex(f) is below f1 (see Fig.1). 
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break; 
              } else if (ori3 < 0) {
                // apex(f) is above f2, continue.
              } else { // ori3 == 0; 
                // f is coplanar and codirection with f2.
                unifysubfaces(&(f2->ss), &sface);
                break;
              }
            } else if (ori2 < 0) {
              // apex(f) is above f1 below f2, inset it (see Fig. 2).
              break;
            } else { // ori2 == 0;
              // apex(f) is coplanar with f1 (see Fig. 5).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break; 
              } else {
                // f is coplanar and codirection with f1.
                unifysubfaces(&(f1->ss), &sface);
                break;
              }
            }
          } else if (ori1 < 0) {
            // apex(f2) is above f1.
            if (ori2 > 0) {
              // apex(f) is below f1, continue (see Fig. 3).
            } else if (ori2 < 0) {
              // apex(f) is above f1 (see Fig.4).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // apex(f) is below f2, insert it.
                break;
              } else if (ori3 < 0) {
                // apex(f) is above f2, continue.
              } else { // ori3 == 0;
                // f is coplanar and codirection with f2.
                unifysubfaces(&(f2->ss), &sface);
                break;
              }
            } else { // ori2 == 0;
              // f is coplanar and with f1 (see Fig. 6).
              ori3 = orient3d(torg, tdest, sapex(f2->ss), sapex(sface));
              if (ori3 > 0) {
                // f is also codirection with f1.
                unifysubfaces(&(f1->ss), &sface);
                break;
              } else {
                // f is above f2, continue.
              }
            }
          } else { // ori1 == 0;
            // apex(f2) is coplanar with f1. By assumption, f1 is not
            //   coplanar and codirection with f2.
            if (ori2 > 0) {
              // apex(f) is below f1, continue (see Fig. 7).
            } else if (ori2 < 0) {
              // apex(f) is above f1, insert it (see Fig. 7).
              break;
            } else { // ori2 == 0.
              // apex(f) is coplanar with f1 (see Fig. 8).
              // f is either codirection with f1 or is codirection with f2. 
              facenormal(torg, tdest, sapex(f1->ss), n1, 1, NULL);
              facenormal(torg, tdest, sapex(sface), n2, 1, NULL);
              if (DOT(n1, n2) > 0) {
                unifysubfaces(&(f1->ss), &sface);
              } else {
                unifysubfaces(&(f2->ss), &sface);
              }
              break;
            }
          }
          // Go to the next item;
          f1 = f2;
        } // for (m = 0; ...)
        if (sface.sh[3] != NULL) {
          // Insert sface between f1 and f2.
          newlinkitem = (badface *) flippool->alloc();
          newlinkitem->ss = sface;
          newlinkitem->nextitem = f1->nextitem;
          f1->nextitem = newlinkitem;
        }
      } else if (flippool->items == 1) {
        f1 = facelink;
        // Make sure that f is not coplanar and codirection with f1.
        ori1 = orient3d(torg, tdest, sapex(f1->ss), sapex(sface));
        if (ori1 == 0) {
          // f is coplanar with f1 (see Fig. 8).
          facenormal(torg, tdest, sapex(f1->ss), n1, 1, NULL);
          facenormal(torg, tdest, sapex(sface), n2, 1, NULL);
          if (DOT(n1, n2) > 0) {
            // The two faces are codirectional as well.
            unifysubfaces(&(f1->ss), &sface);
          }
        }
        // Add this face to link if it is not deleted.
        if (sface.sh[3] != NULL) {
          // Add this face into link.
          newlinkitem = (badface *) flippool->alloc();
          newlinkitem->ss = sface;
          newlinkitem->nextitem = NULL;
          f1->nextitem = newlinkitem;
        }
      } else {
        // The first face.
        newlinkitem = (badface *) flippool->alloc();
        newlinkitem->ss = sface;
        newlinkitem->nextitem = NULL;
        facelink = newlinkitem;
      }
    } // for (k = idx2faclist[idx]; ...)

    if (b->verbose > 2) {
      printf("      Found %ld segments at (%d  %d).\n", flippool->items,
             pointmark(torg), pointmark(tdest));
    }

    //if (b->nobisect || b->nomerge) { // -Y or -M
      // Set the vertex types of the endpoints of the segment.
      setpointtype(torg, RIDGEVERTEX);
      setpointtype(tdest, RIDGEVERTEX);
    //}

    // Set the connection between this segment and faces containing it,
    //   at the same time, remove redundant segments.
    f1 = facelink;
    for (k = 0; k < flippool->items; k++) {
      sspivot(f1->ss, testseg);
      // If 'testseg' is not 'subsegloop' and is not dead, it is redundant.
      if ((testseg.sh != subsegloop.sh) && (testseg.sh[3] != NULL)) {
        shellfacedealloc(subsegs, testseg.sh);
      }
      // Bonds the subface and the segment together.
      ssbond(f1->ss, subsegloop);
      f1 = f1->nextitem;
    }

    // Create the face ring at the segment.
    if (flippool->items > 1) {
      f1 = facelink;
      for (k = 1; k <= flippool->items; k++) {
        k < flippool->items ? f2 = f1->nextitem : f2 = facelink;
        if (b->verbose > 3) {
          printf("        Bond subfaces (%d, %d, %d) and (%d, %d, %d).\n",
                 pointmark(torg), pointmark(tdest), pointmark(sapex(f1->ss)),
                 pointmark(torg), pointmark(tdest), pointmark(sapex(f2->ss)));
        }
        sbond1(f1->ss, f2->ss);
        f1 = f2;
      }
    }

    // All identified segments has an init marker "0".
    flippool->restart();

    // Are there length constraints?
    if (b->quality && (in->segmentconstraintlist != (REAL *) NULL)) {
      for (k = 0; k < in->numberofsegmentconstraints; k++) {
        e1 = (int) in->segmentconstraintlist[k * 3];
        e2 = (int) in->segmentconstraintlist[k * 3 + 1];
        if (((pointmark(torg) == e1) && (pointmark(tdest) == e2)) ||
            ((pointmark(torg) == e2) && (pointmark(tdest) == e1))) {
          len = in->segmentconstraintlist[k * 3 + 2];
          setareabound(subsegloop, len);
          break;
        }
      }
    }

    subsegloop.sh = shellfacetraverse(subsegs);
  }

  delete [] idx2faclist;
  delete [] facperverlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// mergefacets()    Merge adjacent facets.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::mergefacets()
{
  face parentsh, neighsh, neineish;
  face segloop;
  point pa, pb, pc, pd;
  REAL ang_tol, ang;
  int remsegcount;
  int fidx1, fidx2;
  int fmrk1, fmrk2;

  if (b->verbose > 1) {
    printf("    Merging adjacent facets.\n");
  }

  // The dihedral angle bound for two different facets.
  //   Set by -p option. Default is 179 degree.
  ang_tol = b->facet_ang_tol / 180.0 * PI;
  remsegcount = 0;

  // Loop all segments, merge adjacent coplanar facets.
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != (shellface *) NULL) {
    spivot(segloop, parentsh);
    if (parentsh.sh != NULL) {
      spivot(parentsh, neighsh);
      if (neighsh.sh != NULL) {
        spivot(neighsh, neineish);
        if (neineish.sh == parentsh.sh) {
          // Exactly two subfaces at this segment.
          fidx1 = shellmark(parentsh) - 1;
          fidx2 = shellmark(neighsh) - 1;
          // Only merge them if they are in different facet.
          if (fidx1 != fidx2) {
            // The two subfaces are not in the same facet.
            if (in->facetmarkerlist != NULL) { 
              fmrk1 = in->facetmarkerlist[fidx1];
              fmrk2 = in->facetmarkerlist[fidx2];
            } else {
              fmrk1 = fmrk2 = 0;
            }
            // Only merge them if they have the same boundary marker.
            if (fmrk1 == fmrk2) {
              pa = sorg(segloop);
              pb = sdest(segloop);
              pc = sapex(parentsh);
              pd = sapex(neighsh);
              // Calculate the dihedral angle at the segment [a,b].
              ang = facedihedral(pa, pb, pc, pd);
              if (ang > PI) ang = (2 * PI - ang);
              if (ang > ang_tol) {
                if (b->verbose > 2) {
                  printf("      Merge at segment (%d, %d)-(%d, %d) ang = %g\n",
                         pointmark(pa), pointmark(pb), pointmark(pc), 
                         pointmark(pd), ang / PI * 180.0);
                }
                remsegcount++;
                ssdissolve(parentsh);
                ssdissolve(neighsh);
                shellfacedealloc(subsegs, segloop.sh);
                // Add the edge to flip stack.
                flipshpush(&parentsh);
              } // if (ang > ang_tol)
            } // if (fmrk1 == fmrk2)
          } // if (fidx1 != fidx2)
        } // if (neineish.sh == parentsh.sh)
      }
    }
    segloop.sh = shellfacetraverse(subsegs);
  }

  if (flipstack != NULL) {
    lawsonflip(); // Recover Delaunayness.
  }


  if (b->verbose > 1) {
    printf("    %d segments are removed.\n", remsegcount);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// identifypscedges()    Identify PSC edges.                                 //
//                                                                           //
// The set of PSC edges are provided in the 'in->edgelist'. Each edge should //
// also be an edge in the surface mesh.  We find the corresponding edges in  //
// the surface mesh and make them segments of the mesh.                      //
//                                                                           //
// It is possible to give an edge which is not in any facet, i.e., it is a   //
// dangling edge inside the volume.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::identifypscedges(point *idx2verlist)
{
  face* shperverlist;
  int* idx2shlist;
  face searchsh, neighsh;
  face segloop, checkseg, newseg;
  point checkpt, pa, pb;
  int *endpts;
  int edgemarker;
  int idx, i, j;

  int e1, e2;
  REAL len;

  if (!b->quiet) {
    printf("Inserting edges ...\n");
  }

  // All identified segments have the initial marker '0'.
  // All segments inserted here should have a non-zero marker.

  // Construct a map from points to subfaces.
  makepoint2submap(subfaces, idx2shlist, shperverlist);

  // Process the set of PSC edges.
  for (i = 0; i < in->numberofedges; i++) {
    endpts = &(in->edgelist[(i << 1)]);
    // Find a face contains the edge.
    newseg.sh = NULL;
    searchsh.sh = NULL;
    idx = endpts[0] - in->firstnumber;
    for (j = idx2shlist[idx]; j < idx2shlist[idx + 1]; j++) {
      checkpt = sdest(shperverlist[j]);
      if (pointmark(checkpt) == endpts[1]) {
        searchsh = shperverlist[j];
        break; // Found.
      } else {
        checkpt = sapex(shperverlist[j]);
        if (pointmark(checkpt) == endpts[1]) {
          senext2(shperverlist[j], searchsh);
          sesymself(searchsh);
          break;
        }
      }
    } // j
    edgemarker = 0;
    if (in->edgemarkerlist) {
      edgemarker = in->edgemarkerlist[i];
    }
    if (edgemarker == 0) {
      edgemarker = 1;
    }
    // We should find a subface having this edge.
    if (searchsh.sh != NULL) {
      // Check if this edge is already a segment of the mesh.
      sspivot(searchsh, checkseg);
      if (checkseg.sh != NULL) {
        // There should be no duplicated edges.
        assert(shellmark(checkseg) == 0);
        setshellmark(checkseg, edgemarker);
      } else {
        // Create a new segment at this edge.
        pa = sorg(searchsh);
        pb = sdest(searchsh);
        if (b->verbose > 2) {
          printf("      Create a new segment (%d, %d).\n", 
                 pointmark(pa), pointmark(pb));
        }
        makeshellface(subsegs, &newseg);
        setshvertices(newseg, pa, pb, NULL);
        setshellmark(newseg, edgemarker);
        ssbond(searchsh, newseg);
        spivot(searchsh, neighsh);
        if (neighsh.sh != NULL) {
          ssbond(neighsh, newseg);
          // There should be only two subfaces at this segment.
          spivotself(neighsh); // SELF_CHECK
          assert(neighsh.sh == searchsh.sh);
        }
        if (!b->psc) {
          setpointtype(pa, RIDGEVERTEX);
          setpointtype(pb, RIDGEVERTEX);
        }
      }
    } else {
      // It is a dangling segment (not belong to any facets).
      // Get the two endpoints of this segment.
      pa = idx2verlist[endpts[0]];
      pb = idx2verlist[endpts[1]];
      if (b->verbose > 2) {
        printf("      Create a new segment (%d, %d) - dangling.\n", 
               pointmark(pa), pointmark(pb));
      }
      makeshellface(subsegs, &newseg);
      setshvertices(newseg, pa, pb, NULL);
      setshellmark(newseg, edgemarker);
      //if (!b->psc) {
        setpointtype(pa, RIDGEVERTEX);
        setpointtype(pb, RIDGEVERTEX);
      //}
    }

    if (newseg.sh != NULL) {
      if (b->quality && (in->segmentconstraintlist != (REAL *) NULL)) {
        for (i = 0; i < in->numberofsegmentconstraints; i++) {
          e1 = (int) in->segmentconstraintlist[i * 3];
          e2 = (int) in->segmentconstraintlist[i * 3 + 1];
          if (((pointmark(pa) == e1) && (pointmark(pb) == e2)) ||
              ((pointmark(pa) == e2) && (pointmark(pb) == e1))) {
            len = in->segmentconstraintlist[i * 3 + 2];
            setareabound(newseg, len);
            break;
          }
        }
      }
    }
  } // i

  if (b->psc) {
    // Delete all segments of the mesh with a marker '0'.
    subsegs->traversalinit();
    segloop.sh = shellfacetraverse(subsegs);
    while (segloop.sh != NULL) {
      if (shellmark(segloop) == 0) {
        if (b->verbose > 2) {
          printf("      Remove a segment (%d, %d).\n", 
                 pointmark(sorg(segloop)), pointmark(sdest(segloop)));
        }
        spivot(segloop, searchsh);
        if (searchsh.sh != NULL) {
          ssdissolve(searchsh);
          spivot(searchsh, neighsh);
          if (neighsh.sh != NULL) {
            ssdissolve(neighsh);
            // There should be only two subfaces at this segment.
            spivotself(neighsh); // SELF_CHECK
            assert(neighsh.sh == searchsh.sh);
          }
        }
        shellfacedealloc(subsegs, segloop.sh);
      }
      segloop.sh = shellfacetraverse(subsegs);
    }
  }

  delete [] shperverlist;
  delete [] idx2shlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// meshsurface()    Create a surface mesh of the input PLC.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::meshsurface()
{
  arraypool *ptlist, *conlist;
  point *idx2verlist;
  point tstart, tend, *pnewpt, *cons;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int end1, end2;
  int shmark, i, j;

  if (!b->quiet) {
    printf("Creating surface mesh ...\n");
  }

  // Create a map from indices to points.
  makeindex2pointmap(idx2verlist);

  // Initialize arrays (block size: 2^8 = 256).
  ptlist = new arraypool(sizeof(point *), 8);
  conlist = new arraypool(2 * sizeof(point *), 8);

  // Loop the facet list, triangulate each facet.
  for (shmark = 1; shmark <= in->numberoffacets; shmark++) {

    // Get a facet F.
    f = &in->facetlist[shmark - 1];

    // Process the duplicated points first, they are marked with type
    //   DUPLICATEDVERTEX.  If p and q are duplicated, and p'index > q's,
    //   then p is substituted by q.
    if (dupverts > 0l) {
      // Loop all polygons of this facet.
      for (i = 0; i < f->numberofpolygons; i++) {
        p = &(f->polygonlist[i]);
        // Loop other vertices of this polygon.
        for (j = 0; j < p->numberofvertices; j++) {
          end1 = p->vertexlist[j];
          tstart = idx2verlist[end1];
          if (pointtype(tstart) == DUPLICATEDVERTEX) {
            // Reset the index of vertex-j.
            tend = point2ppt(tstart);
            end2 = pointmark(tend);
            p->vertexlist[j] = end2;
          }
        }
      }
    }

    // Loop polygons of F, get the set of vertices and segments.
    for (i = 0; i < f->numberofpolygons; i++) {
      // Get a polygon.
      p = &(f->polygonlist[i]);
      // Get the first vertex.
      end1 = p->vertexlist[0];
      if ((end1 < in->firstnumber) || 
          (end1 >= in->firstnumber + in->numberofpoints)) {
        if (!b->quiet) {
          printf("Warning:  Invalid the 1st vertex %d of polygon", end1);
          printf(" %d in facet %d.\n", i + 1, shmark);
        }
        continue; // Skip this polygon.
      }
      tstart = idx2verlist[end1];
      // Add tstart to V if it haven't been added yet.
      if (!pinfected(tstart)) {
        pinfect(tstart);
        ptlist->newindex((void **) &pnewpt);
        *pnewpt = tstart;
      }
      // Loop other vertices of this polygon.
      for (j = 1; j <= p->numberofvertices; j++) {
        // get a vertex.
        if (j < p->numberofvertices) {
          end2 = p->vertexlist[j];
        } else {
          end2 = p->vertexlist[0];  // Form a loop from last to first.
        }
        if ((end2 < in->firstnumber) ||
            (end2 >= in->firstnumber + in->numberofpoints)) {
          if (!b->quiet) {
            printf("Warning:  Invalid vertex %d in polygon %d", end2, i + 1);
            printf(" in facet %d.\n", shmark);
          }
        } else {
          if (end1 != end2) {
            // 'end1' and 'end2' form a segment.
            tend = idx2verlist[end2];
            // Add tstart to V if it haven't been added yet.
            if (!pinfected(tend)) {
              pinfect(tend);
              ptlist->newindex((void **) &pnewpt);
              *pnewpt = tend;
            }
            // Save the segment in S (conlist).
            conlist->newindex((void **) &cons);
            cons[0] = tstart;
            cons[1] = tend;
            // Set the start for next continuous segment.
            end1 = end2;
            tstart = tend;
          } else {
            // Two identical vertices mean an isolated vertex of F.
            if (p->numberofvertices > 2) {
              // This may be an error in the input, anyway, we can continue
              //   by simply skipping this segment.
              if (!b->quiet) {
                printf("Warning:  Polygon %d has two identical verts", i + 1);
                printf(" in facet %d.\n", shmark);
              }
            } 
            // Ignore this vertex.
          }
        }
        // Is the polygon degenerate (a segment or a vertex)?
        if (p->numberofvertices == 2) break;
      }
    }
    // Unmark vertices.
    for (i = 0; i < ptlist->objects; i++) {
      pnewpt = (point *) fastlookup(ptlist, i);
      puninfect(*pnewpt);
    }

    // Triangulate F into a CDT.
    triangulate(shmark, ptlist, conlist, f->numberofholes, f->holelist);

    // Clear working lists.
    ptlist->restart();
    conlist->restart();
  }

  if (!b->diagnose) {
    // Remove redundant segments and build the face links.
    unifysegments();
  }

  if (!b->nomerge && !b->nobisect && !b->diagnose) {
    // Merge adjacent coplanar facets.
    mergefacets();
  }

  if (in->numberofedges > 0) { // if (b->psc)
    // There are segments specified by the user. Read and create them.
    identifypscedges(idx2verlist);
  }

  if (b->object == tetgenbehavior::STL) {
    // Remove redundant vertices (for .stl input mesh).
    jettisonnodes();
  }

  if (b->verbose) {
    printf("  %ld (%ld) subfaces (segments).\n", subfaces->items, 
           subsegs->items);
  }

  // The total number of iunput segments.
  insegments = subsegs->items;

  delete [] idx2verlist;
  delete ptlist;
  delete conlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// interecursive()    Recursively do intersection test on a set of triangles.//
//                                                                           //
// Recursively split the set 'subfacearray' of subfaces into two sets using  //
// a cut plane parallel to x-, or, y-, or z-axies.  The split criteria are   //
// follows. Assume the cut plane is H, and H+ denotes the left halfspace of  //
// H, and H- denotes the right halfspace of H; and s be a subface:           //
//                                                                           //
//    (1) If all points of s lie at H+, put it into left array;              //
//    (2) If all points of s lie at H-, put it into right array;             //
//    (3) If some points of s lie at H+ and some of lie at H-, or some       //
//        points lie on H, put it into both arraies.                         //
//                                                                           //
// Partitions by x-axis if axis == '0'; by y-axis if axis == '1'; by z-axis  //
// if axis == '2'. If current cut plane is parallel to the x-axis, the next  //
// one will be parallel to y-axis, and the next one after the next is z-axis,//
// and then alternately return back to x-axis.                               //
//                                                                           //
// Stop splitting when the number of triangles of the input array is not     //
// decreased anymore. Do tests on the current set.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::interecursive(shellface** subfacearray, int arraysize, 
                               int axis, REAL bxmin, REAL bxmax, REAL bymin, 
                               REAL bymax, REAL bzmin, REAL bzmax, 
                               int* internum)
{
  shellface **leftarray, **rightarray;
  face sface1, sface2;
  point p1, p2, p3;
  point p4, p5, p6;
  enum interresult intersect;
  REAL split;
  bool toleft, toright;
  int leftsize, rightsize;
  int i, j;

  if (b->verbose > 2) {
    printf("      Recur %d faces. Bbox (%g, %g, %g),(%g, %g, %g). %s-axis\n",
           arraysize, bxmin, bymin, bzmin, bxmax, bymax, bzmax,
           axis == 0 ? "x" : (axis == 1 ? "y" : "z"));
  }
    
  leftarray = new shellface*[arraysize];
  if (leftarray == NULL) {
    terminatetetgen(1);
  }
  rightarray = new shellface*[arraysize];
  if (rightarray == NULL) {
    terminatetetgen(1);
  }
  leftsize = rightsize = 0;

  if (axis == 0) {
    // Split along x-axis.
    split = 0.5 * (bxmin + bxmax);
  } else if (axis == 1) {
    // Split along y-axis.
    split = 0.5 * (bymin + bymax);
  } else {
    // Split along z-axis.
    split = 0.5 * (bzmin + bzmax);
  }

  for (i = 0; i < arraysize; i++) {
    sface1.sh = subfacearray[i];
    p1 = (point) sface1.sh[3];
    p2 = (point) sface1.sh[4];
    p3 = (point) sface1.sh[5];
    toleft = toright = false;
    if (p1[axis] < split) {
      toleft = true;
      if (p2[axis] >= split || p3[axis] >= split) {
        toright = true;
      } 
    } else if (p1[axis] > split) {
      toright = true;
      if (p2[axis] <= split || p3[axis] <= split) {
        toleft = true;
      } 
    } else {
      // p1[axis] == split;
      toleft = true;
      toright = true;
    }
    // At least one is true;
    assert(!(toleft == false && toright == false));
    if (toleft) {
      leftarray[leftsize] = sface1.sh;
      leftsize++;
    }
    if (toright) {
      rightarray[rightsize] = sface1.sh;
      rightsize++;
    }
  }

  if (leftsize < arraysize && rightsize < arraysize) {
    // Continue to partition the input set. Now 'subfacearray' has been
    //   split into two sets, it's memory can be freed. 'leftarray' and
    //   'rightarray' will be freed in the next recursive (after they're
    //   partitioned again or performing tests).
    delete [] subfacearray;
    // Continue to split these two sets.
    if (axis == 0) {
      interecursive(leftarray, leftsize, 1, bxmin, split, bymin, bymax,
                    bzmin, bzmax, internum);
      interecursive(rightarray, rightsize, 1, split, bxmax, bymin, bymax,
                    bzmin, bzmax, internum);
    } else if (axis == 1) {
      interecursive(leftarray, leftsize, 2, bxmin, bxmax, bymin, split,
                    bzmin, bzmax, internum);
      interecursive(rightarray, rightsize, 2, bxmin, bxmax, split, bymax,
                    bzmin, bzmax, internum);
    } else {
      interecursive(leftarray, leftsize, 0, bxmin, bxmax, bymin, bymax,
                    bzmin, split, internum);
      interecursive(rightarray, rightsize, 0, bxmin, bxmax, bymin, bymax,
                    split, bzmax, internum);
    }
  } else {
    if (b->verbose > 1) {
      printf("  Checking intersecting faces.\n");
    }
    // Perform a brute-force compare on the set.
    for (i = 0; i < arraysize; i++) {
      sface1.sh = subfacearray[i];
      p1 = (point) sface1.sh[3];
      p2 = (point) sface1.sh[4];
      p3 = (point) sface1.sh[5];
      for (j = i + 1; j < arraysize; j++) {
        sface2.sh = subfacearray[j];
        p4 = (point) sface2.sh[3];
        p5 = (point) sface2.sh[4];
        p6 = (point) sface2.sh[5];
        intersect = (enum interresult) tri_tri_inter(p1, p2, p3, p4, p5, p6);
        if (intersect == INTERSECT || intersect == SHAREFACE) {
          if (!b->quiet) {
            if (intersect == INTERSECT) {
              printf("  Facet #%d intersects facet #%d at triangles:\n",
                     shellmark(sface1), shellmark(sface2));
              printf("    (%4d, %4d, %4d) and (%4d, %4d, %4d)\n",
                     pointmark(p1), pointmark(p2), pointmark(p3),
                     pointmark(p4), pointmark(p5), pointmark(p6));
            } else {
              printf("  Facet #%d duplicates facet #%d at triangle:\n",
                     shellmark(sface1), shellmark(sface2));
              printf("    (%4d, %4d, %4d) and (%4d, %4d, %4d)\n",
                     pointmark(p1), pointmark(p2), pointmark(p3),
                     pointmark(p4), pointmark(p5), pointmark(p6));
            }
          }
          // Increase the number of intersecting pairs.
          (*internum)++; 
          // Infect these two faces (although they may already be infected).
          sinfect(sface1);
          sinfect(sface2);
        }
      }
    }
    // Don't forget to free all three arrays. No further partition.
    delete [] leftarray;
    delete [] rightarray;  
    delete [] subfacearray;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// detectinterfaces()    Detect intersecting triangles.                      //
//                                                                           //
// Given a set of triangles,  find the pairs of intersecting triangles from  //
// them.  Here the set of triangles is in 'subfaces' which is a surface mesh //
// of a PLC (.poly or .smesh).                                               //
//                                                                           //
// To detect whether two triangles are intersecting is done by the routine   //
// 'tri_tri_inter()'.  The algorithm for the test is very simple and stable. //
// It is based on geometric orientation test which uses exact arithmetics.   //
//                                                                           //
// Use divide-and-conquer algorithm for reducing the number of intersection  //
// tests.  Start from the bounding box of the input point set, recursively   //
// partition the box into smaller boxes, until the number of triangles in a  //
// box is not decreased anymore. Then perform triangle-triangle tests on the //
// remaining set of triangles.  The memory allocated in the input set is     //
// freed immediately after it has been partitioned into two arrays.  So it   //
// can be re-used for the consequent partitions.                             //
//                                                                           //
// On return, the pool 'subfaces' will be cleared, and only the intersecting //
// triangles remain for output (to a .face file).                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::detectinterfaces()
{
  shellface **subfacearray;
  face shloop;
  int internum;
  int i;

  if (!b->quiet) {
    printf("Detecting self-intersecting facets...\n");
  }

  // Construct a map from indices to subfaces;
  subfacearray = new shellface*[subfaces->items];
  subfaces->traversalinit();
  shloop.sh = shellfacetraverse(subfaces);
  i = 0;
  while (shloop.sh != (shellface *) NULL) {
    subfacearray[i] = shloop.sh;
    shloop.sh = shellfacetraverse(subfaces);
    i++;
  }

  internum = 0;
  // Recursively split the set of triangles into two sets using a cut plane
  //   parallel to x-, or, y-, or z-axies.  Stop splitting when the number
  //   of subfaces is not decreasing anymore. Do tests on the current set.
  interecursive(subfacearray, subfaces->items, 0, xmin, xmax, ymin, ymax,
                zmin, zmax, &internum);

  if (!b->quiet) {
    if (internum > 0) {
      printf("\n!! Found %d pairs of faces are intersecting.\n\n", internum);
    } else {
      printf("\nNo faces are intersecting.\n\n");
    }
  }

  if (internum > 0) {
    // Traverse all subfaces, deallocate those have not been infected (they
    //   are not intersecting faces). Uninfect those have been infected.
    //   After this loop, only intersecting faces remain.
    subfaces->traversalinit();
    shloop.sh = shellfacetraverse(subfaces);
    while (shloop.sh != (shellface *) NULL) {
      if (sinfected(shloop)) {
        suninfect(shloop);
      } else {
        shellfacedealloc(subfaces, shloop.sh);
      }
      shloop.sh = shellfacetraverse(subfaces);
    }
  } else {
    // Deallocate all subfaces.
    subfaces->restart();
  }
}

////                                                                       ////
////                                                                       ////
//// surface_cxx //////////////////////////////////////////////////////////////

