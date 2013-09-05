#include "../tetgen.h"
//// steiner_cxx //////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkflipeligibility()    A call back function for boundary recovery.     //
//                                                                           //
// 'fliptype' indicates which elementary flip will be performed: 1 : 2-to-3, //
// and 2 : 3-to-2, respectively.                                             //
//                                                                           //
// 'pa, ..., pe' are the vertices involved in this flip, where [a,b,c] is    //
// the flip face, and [d,e] is the flip edge. NOTE: 'pc' may be 'dummypoint',//
// other points must not be 'dummypoint'.                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkflipeligibility(int fliptype, point pa, point pb, 
                                     point pc, point pd, point pe,
                                     int level, int edgepivot,
                                     flipconstraints* fc)
{
  int rejflag;
  int i;

  point tmppts[3];
  REAL normal[3], area, len;
  REAL ori1, ori2;
  REAL abovept[3];

  enum interresult dir;
  int types[2], poss[4];
  int intflag;

  rejflag = 0;

  if (fc->seg[0] != NULL) {
    // A constraining edge is given (e.g., for edge recovery).
    if (fliptype == 1) {
      // A 2-to-3 flip: [a,b,c] => [e,d,a], [e,d,b], [e,d,c].
      if (pc != dummypoint) {
        // Do not flip if the newly created faces intersect this edge in 
        //   their interiors.
        tmppts[0] = pa;
        tmppts[1] = pb;
        tmppts[2] = pc;
        for (i = 0; i < 3 && !rejflag; i++) {
          intflag = tri_edge_test(pe, pd, tmppts[i], fc->seg[0], fc->seg[1], 
                                  NULL, 1, types, poss);
          if (intflag == 2) {
            // They intersect at a single point.
            dir = (enum interresult) types[0];
            if (dir == ACROSSFACE) {
              // The interior of [e,d,#] intersect the segment.
              rejflag = 1;
            } else if (dir == ACROSSEDGE) {
              if (poss[0] == 0) {
                // The interior of [e,d] intersect the segment.
                // Since [e,d] is the newly created edge. Reject this flip.
                rejflag = 1; 
              }
            }
          } else if (intflag == 4) {
            // They may intersect at either a point or a line segment.
            dir = (enum interresult) types[0];
            if (dir == ACROSSEDGE) {
              if (poss[0] == 0) {
                // The interior of [e,d] intersect the segment.
                // Since [e,d] is the newly created edge. Reject this flip.
                rejflag = 1;
              }
            }
          } // if (intflag == 4)
        } // i
      } else { // pc == dummypoint
        // Do not flip if the new hull edge [e,d] will intersect this edge
        //   in its interior. 
        // Comment: Here we actually need a 3D edge-edge test.
        //   We only do test if the edge in 'fc' is coplanar with the plane
        //   containing a,b,e,and d.
        // Choose a better triangle [a,b,e] or [a,b,d].
        facenormal(pa, pb, pe, normal, 1, &len);
        area = sqrt(DOT(normal, normal));
        facenormal(pa, pb, pd, normal, 1, &len);
        len = sqrt(DOT(normal, normal)); // Re-use len as area.
        if (area > len) {
          // Choose [a,b,e]
          ori1 = orient3d(pa, pb, pe, fc->seg[0]);
          ori2 = orient3d(pa, pb, pe, fc->seg[1]);
        } else {
          // Choose [a,b,d]
          ori1 = orient3d(pa, pb, pd, fc->seg[0]);
          ori2 = orient3d(pa, pb, pd, fc->seg[1]);
        }
        if ((ori1 == 0) && (ori2 == 0)) {
          calculateabovepoint4(pa, pb, pe, pd);
          for (i = 0; i < 3; i++) {
            abovept[i] = dummypoint[i];
          }
          intflag = tri_edge_test(pe, pd, abovept, fc->seg[0], fc->seg[1], 
                                  NULL, 1, types, poss);
          if (intflag == 2) {
            dir = (enum interresult) types[0];
            assert(dir != ACROSSFACE);
            if (dir == ACROSSEDGE) {
              if (poss[0] == 0) {
                // The interior of [e,d] intersect the segment.
                // Since [e,d] is the newly created edge. Reject this flip.
                rejflag = 1;
              }
            } 
          } else if (intflag == 4) {
            // [e,d,abovept] is coplanar with the constraining edge 'fc'.
            // This is poissible if the edge in 'fc' is just the edge [e,d]
            //   (SHAREEDGE) or they share a common vertex (SHAREVEER)
            dir = (enum interresult) types[0];
            if (dir == ACROSSEDGE) {
              // This case can only happen if [e,d] is coplanar with 'fc'.
              assert(0);  // Not possible.
            }
          }
        }
      } // if (pc == dummypoint)
    } else if (fliptype == 2) {
      // A 3-to-2 flip: [e,d,a], [e,d,b], [e,d,c] => [a,b,c]
      if (pc != dummypoint) {
        if (!rejflag) {
          // Check if the new face [a,b,c] intersect the edge in its interior.
          intflag = tri_edge_test(pa, pb, pc, fc->seg[0], fc->seg[1], NULL, 
                                  1, types, poss);
          if (intflag == 2) {
            // They intersect at a single point.
            dir = (enum interresult) types[0];
            if (dir == ACROSSFACE) {
              // The interior of [a,b,c] intersect the segment.
              rejflag = 1; // Do not flip.
            } else if (dir == ACROSSEDGE) {
              // This case is possible since we allow a previous 2-to-3 flip
              //   even it will create a degenerate tet at edge [a,b].
            }
          } else if (intflag == 4) {
            // [a,b,c] is coplanar with the edge. 
            dir = (enum interresult) types[0];
            if (dir == ACROSSEDGE) {
              // The boundary of [a,b,c] intersect the segment.
              // An example is found in case 'camila.poly', during the recovery
              //   of segment [151, 161] (at linklevel = 2). See: 2011-06-10-a.
              rejflag = 1; // Do not flip.
            }
          }
        } // if (!relflag)
      } else { // pc == dummypoint
        // The flip 3-to-2 will replace [e,d] with a new hull edge [a,b].
        // Only do flip if [a,b] does not intersect the edge of 'fc'.
        // Comment: Here we acutually need a 3D edge-edge intersection test.
        //   We only do test if the edge in 'fc' is coplanar with the plane
        //   containing a,b,e, and d.
        // Choose a better triangle [a,b,e] or [a,b,d].
        facenormal(pa, pb, pe, normal, 1, &len);
        area = sqrt(DOT(normal, normal));
        facenormal(pa, pb, pd, normal, 1, &len);
        len = sqrt(DOT(normal, normal)); // Re-use len as area.
        if (area > len) {
          // Choose [a,b,e]
          ori1 = orient3d(pa, pb, pe, fc->seg[0]);
          ori2 = orient3d(pa, pb, pe, fc->seg[1]);
        } else {
          // Choose [a,b,d]
          ori1 = orient3d(pa, pb, pd, fc->seg[0]);
          ori2 = orient3d(pa, pb, pd, fc->seg[1]);
        }
        if ((ori1 == 0) && (ori2 == 0)) {
          // The edge in 'fc' is coplanar with the plane containing [a,b,e,d].
          calculateabovepoint4(pa, pb, pe, pd);
          for (i = 0; i < 3; i++) {
            abovept[i] = dummypoint[i];
          }
          intflag = tri_edge_test(pa, pb, abovept, fc->seg[0], fc->seg[1], 
                                  NULL, 1, types, poss);
          if (intflag == 2) {
            dir = (enum interresult) types[0];
            assert(dir != ACROSSFACE);
            if (dir == ACROSSEDGE) {
              assert(0); // Check this case.
              rejflag = 1; // Do not flip.
            }
          } else if (intflag == 4) {
            // The edge 'fc' is coplanar with [a,b,abovept]. 
            // This is poissible if the edge in 'fc' is just the edge [a,b]
            //   (SHAREEDGE) or they share a common vertex (SHAREVEER)
            dir = (enum interresult) types[0];
            if (dir == ACROSSEDGE) {
              // This case can only happen if [a,b] is coplanar with 'fc'.
              assert(0);  // Not possible.
            }
          }
        } // if (ori1 == 0 && ori2 == 0)
      }
    } else {
      assert(0); // An unknown flip type.
    }
  } // if (fc->seg[0] != NULL)

  if ((fc->fac[0] != NULL) && !rejflag) {
    // A constraining face is given (e.g., for face recovery).
    if (fliptype == 1) {
      // A 2-to-3 flip.
      // Test if the new edge [e,d] intersects the face.
      intflag = tri_edge_test(fc->fac[0], fc->fac[1], fc->fac[2], pe, pd, 
                              NULL, 1, types, poss);
      if (intflag == 2) {
        // They intersect at a single point.
        dir = (enum interresult) types[0];
        if (dir == ACROSSFACE) {
          rejflag = 1;
        } else if (dir == ACROSSEDGE) {
          rejflag = 1;
        } 
      } else if (intflag == 4) {
        // The edge [e,d] is coplanar with the face.
        // There may be two intersections.
        for (i = 0; i < 2 && !rejflag; i++) {
          dir = (enum interresult) types[i];
          if (dir == ACROSSFACE) {
            rejflag = 1;
          } else if (dir == ACROSSEDGE) {
            rejflag = 1;
          }
        }
      }
    } // if (fliptype == 1)
  } // if (fc->fac[0] != NULL)

  if ((fc->remvert != NULL) && !rejflag) {
    // The vertex is going to be removed. Do not create a new edge which
    //   contains this vertex.
    if (fliptype == 1) {
      // A 2-to-3 flip.
      if ((pd == fc->remvert) || (pe == fc->remvert)) {
        rejflag = 1;
      }
    }
  }

  if (fc->remove_large_angle && !rejflag) {
    // Remove a large dihedral angle. Do not create a new small angle.
    REAL cosmaxd = 0, diff;
    if (fliptype == 1) {
      // We assume that neither 'a' nor 'b' is dummypoint.
      assert((pa != dummypoint) && (pb != dummypoint)); // SELF_CHECK
      // A 2-to-3 flip: [a,b,c] => [e,d,a], [e,d,b], [e,d,c].
      // The new tet [e,d,a,b] will be flipped later. Only two new tets:
      //   [e,d,b,c] and [e,d,c,a] need to be checked.
      if ((pc != dummypoint) && (pe != dummypoint) && (pd != dummypoint)) {
        // Get the largest dihedral angle of [e,d,b,c].
        tetalldihedral(pe, pd, pb, pc, NULL, &cosmaxd, NULL);
        diff = cosmaxd - fc->cosdihed_in;
        if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
        if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
          rejflag = 1;
        } else {
          // Record the largest new angle.
          if (cosmaxd < fc->cosdihed_out) {
            fc->cosdihed_out = cosmaxd; 
          }
          // Get the largest dihedral angle of [e,d,c,a].
          tetalldihedral(pe, pd, pc, pa, NULL, &cosmaxd, NULL);
          diff = cosmaxd - fc->cosdihed_in;
          if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding.
          if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
            rejflag = 1;
          } else {
            // Record the largest new angle.
            if (cosmaxd < fc->cosdihed_out) {
              fc->cosdihed_out = cosmaxd; 
            }
          }
        }
      } // if (pc != dummypoint && ...)
    } else if (fliptype == 2) {
      // A 3-to-2 flip: [e,d,a], [e,d,b], [e,d,c] => [a,b,c]
      // We assume that neither 'e' nor 'd' is dummypoint.
      assert((pe != dummypoint) && (pd != dummypoint)); // SELF_CHECK
      if (level == 0) {
        // Both new tets [a,b,c,d] and [b,a,c,e] are new tets.
        if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
          // Get the largest dihedral angle of [a,b,c,d].
          tetalldihedral(pa, pb, pc, pd, NULL, &cosmaxd, NULL);
          diff = cosmaxd - fc->cosdihed_in;
          if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0; // Rounding
          if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
            rejflag = 1;
          } else {
            // Record the largest new angle.
            if (cosmaxd < fc->cosdihed_out) {
              fc->cosdihed_out = cosmaxd; 
            }
            // Get the largest dihedral angle of [b,a,c,e].
            tetalldihedral(pb, pa, pc, pe, NULL, &cosmaxd, NULL);
            diff = cosmaxd - fc->cosdihed_in;
            if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0;// Rounding
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              if (cosmaxd < fc->cosdihed_out) {
                fc->cosdihed_out = cosmaxd; 
              }
            }
          }
        }
      } else { // level > 0
        assert(edgepivot != 0);
        if (edgepivot == 1) {
          // The new tet [a,b,c,d] will be flipped. Only check [b,a,c,e].
          if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
            // Get the largest dihedral angle of [b,a,c,e].
            tetalldihedral(pb, pa, pc, pe, NULL, &cosmaxd, NULL);
            diff = cosmaxd - fc->cosdihed_in;
            if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0;// Rounding
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              if (cosmaxd < fc->cosdihed_out) {
                fc->cosdihed_out = cosmaxd; 
              }
            }
          }
        } else {
          assert(edgepivot == 2);
          // The new tet [b,a,c,e] will be flipped. Only check [a,b,c,d].
          if ((pa != dummypoint) && (pb != dummypoint) && (pc != dummypoint)) {
            // Get the largest dihedral angle of [b,a,c,e].
            tetalldihedral(pa, pb, pc, pd, NULL, &cosmaxd, NULL);
            diff = cosmaxd - fc->cosdihed_in;
            if (fabs(diff/fc->cosdihed_in) < b->epsilon) diff = 0.0;// Rounding
            if (diff <= 0) { //if (cosmaxd <= fc->cosdihed_in) {
              rejflag = 1;
            } else {
              // Record the largest new angle.
              if (cosmaxd < fc->cosdihed_out) {
                fc->cosdihed_out = cosmaxd; 
              }
            }
          }
        } // edgepivot
      } // level
    }
  }

  return rejflag;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// removeedgebyflips()    Remove an edge by flips.                           //
//                                                                           //
// 'flipedge' is a non-convex or flat edge [a,b,#,#] to be removed.          //
//                                                                           //
// The return value is a positive integer, it indicates whether the edge is  //
// removed or not.  A value "2" means the edge is removed, othereise, the    //
// edge is not removed and the value (must >= 3) is the current number of    //
// tets in the edge star.                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::removeedgebyflips(triface *flipedge, flipconstraints* fc)
{
  triface *abtets, spintet;
  face checkseg, *paryseg;
  int n, nn, i;


  if (b->verbose > 2) {
    printf("      Removing edge (%d, %d)\n", pointmark(org(*flipedge)), 
           pointmark(dest(*flipedge)));
  }

  fc->clearcounters();

  if (checksubsegflag) {
    // Do not flip a segment.
    tsspivot1(*flipedge, checkseg);
    if (checkseg.sh != NULL) {
      if (b->verbose > 2) {
        printf("      Can't flip a segment (%d, %d).\n", 
               pointmark(sorg(checkseg)), pointmark(sdest(checkseg))); 
      }
      fc->encsegcount++;
      if (fc->collectencsegflag) {
        if (!sinfected(checkseg)) {
          // Queue this segment in list.
          sinfect(checkseg);                
          caveencseglist->newindex((void **) &paryseg);
          *paryseg = checkseg;
        }
      }
      return 0;
    }
  }

  // Count the number of tets at edge [a,b].
  n = 0;
  int counter = 0; // Sum of star counters. // SELF_CHECK.
  spintet = *flipedge;
  i = 0;
  while (1) {
    counter += elemcounter(spintet);
    i++;
    fnextself(spintet);
    if (spintet.tet == flipedge->tet) break;
  }
  //assert(i >= 3);
  if (i < 3) {
    // It is only possible when the mesh contains inverted tetrahedra.  
    assert(checkinverttetflag);
    // Since "return 2" means success, we return 0.
    return 0;
  }
  assert(counter == 0); // SELF_CHECK
  n = i;

  flipstarcount++;
  // Record the maximum star size.
  if (n > maxflipstarsize) {
    maxflipstarsize = n;
  }
  if ((b->flipstarsize > 0) && (n > b->flipstarsize)) {
    // The star size exceeds the limit.
    skpflipstarcount++;
    return 0; // Do not flip it.
  }

  // Allocate spaces.
  abtets = new triface[n];
  // Collect the tets at edge [a,b].
  spintet = *flipedge;
  i = 0;
  while (1) {
    abtets[i] = spintet;
    setelemcounter(abtets[i], 1); // Marktest it (in Star(ab)).
    i++;
    fnextself(spintet);
    if (spintet.tet == flipedge->tet) break;
  }


  // Try to flip the edge (level = 0, edgepivot = 0).
  nn = flipnm(abtets, n, 0, 0, fc);


  if (nn == 2) {
    // Edge is flipped.
    if (b->verbose > 2) {
      printf("      Edge is removed.\n");
    }
  } else {
    if (b->verbose > 2) {
      printf("      Edge is not removed. n(%d), nn(%d).\n", n, nn);
    }
    // Edge is not flipped. Unmarktest the remaining tets in Star(ab).
    for (i = 0; i < nn; i++) {
      assert(elemcounter(abtets[i]) == 1);
      setelemcounter(abtets[i], 0);
    }
    // Restore the input edge (needed by Lawson's flip).
    *flipedge = abtets[0];
  }

  // Release the temporary allocated spaces.
  // NOTE: fc->unflip must be 0.
  int bakunflip = fc->unflip;
  fc->unflip = 0;

  flipnm_post(abtets, n, nn, 0, fc);

  fc->unflip = bakunflip;

  delete [] abtets;

  return nn; //return nn == 2;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// removefacebyflips()    Remove a face by flips.                            //
//                                                                           //
// ASSUMPTIONS:                                                              //
//   - 'flipface' must not be a hull face.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::removefacebyflips(triface *flipface, flipconstraints* fc)
{
  triface fliptets[3], flipedge;
  face checksh;
  point pa, pb, pc, pd, pe;
  REAL ori;
  int reducflag, rejflag;

  if (checksubfaceflag) {
    tspivot(*flipface, checksh);
    if (checksh.sh != NULL) {
      if (b->verbose > 2) {
        printf("      Can't flip a subface.\n"); 
      }
      return 0;
    }
  }

  fliptets[0] = *flipface;
  fsym(*flipface, fliptets[1]);
  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

  if (b->verbose > 2) {
    printf("      Removing face (%d, %d, %d) -- %d, %d\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }

  reducflag = 0;

  ori = orient3d(pa, pb, pd, pe);
  if (ori > 0) {
    ori = orient3d(pb, pc, pd, pe);
    if (ori > 0) {
      ori = orient3d(pc, pa, pd, pe);
      if (ori > 0) {
        // Found a 2-to-3 flip.
        reducflag = 1;
      } else {
        eprev(*flipface, flipedge); // [c,a]
      }
    } else {
      enext(*flipface, flipedge); // [b,c]
    }
  } else {
    flipedge = *flipface; // [a,b]
  }

  if (reducflag) {
    // A 2-to-3 flip is found.
    rejflag = 0;
    if (fc != NULL) {
      //rejflag = checkflipeligibility(1, pa, pb, pc, pd, pe, fc);
    }
    if (!rejflag) {
      flip23(fliptets, 0, 0, 0);
      if (b->verbose > 2) {
        printf("      Face is removed by a 2-to-3 flip.\n");
      }
      return 1;
    }
  } else {
    // Try to flip the selected edge of this face.
    if (removeedgebyflips(&flipedge, fc) == 2) {
      if (b->verbose > 2) {
        printf("      Face is removed by removing an edge.\n");
      }
      return 1;
    }
  }

  // Face is not removed.
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoveredge()    Recover an edge in current tetrahedralization.           //
//                                                                           //
// If the edge is recovered, 'searchtet' returns a tet containing the edge.  //
//                                                                           //
// This edge may intersect a set of faces and edges in the mesh.  All these  //
// faces or edges are needed to be flipped.                                  //
//                                                                           //
// If the parameter 'fullsearch' is set, it tries to flip any face or edge   //
// that intersects the recovering edge.  Otherwise, only the face or edge    //
// which is visible by 'startpt' is tried.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::recoveredgebyflips(point startpt, point endpt, 
                                   triface* searchtet, int fullsearch)
{
  triface neightet, spintet;
  point pa, pb, pc, pd;
  badface bakface;
  enum interresult dir, dir1;
  flipconstraints fc;
  int types[2], poss[4], pos = 0;
  int success;
  int i, j;

  if (b->verbose > 2) {
    printf("      Recovering edge (%d, %d)\n", pointmark(startpt), 
           pointmark(endpt));
  }
  fc.seg[0] = startpt;
  fc.seg[1] = endpt;

  // The mainloop of the edge reocvery.
  while (1) { // Loop I

    // Search the edge from 'startpt'.
    point2tetorg(startpt, *searchtet);
    dir = finddirection(searchtet, endpt);
    if (dir == ACROSSVERT) {
      if (dest(*searchtet) == endpt) {
        return 1; // Edge is recovered.
      } else {
        // A PLC problem, or there is a Steiner point.
        terminatetetgen(3); 
      }
    }

    // The edge is missing. 

    // Try to flip the first intersecting face/edge.
    enextesymself(*searchtet); // Go to the opposite face.
    if (dir == ACROSSFACE) {
      // A face is intersected with the segment. Try to flip it.
      if (removefacebyflips(searchtet, &fc)) {
        continue;
      }
    } else if (dir == ACROSSEDGE) {
      // An edge is intersected with the segment. Try to flip it.
      if (removeedgebyflips(searchtet, &fc) == 2) {
        continue;
      }
    } else {
      terminatetetgen(3); //assert(0); // A PLC problem.
    }

    // The edge is missing.

    if (fullsearch) {

      // Try to flip one of the faces/edges which intersects the edge.
      success = 0;

      // Loop through the sequence of intersecting faces/edges from
      //   'startpt' to 'endpt'.
      point2tetorg(startpt, *searchtet);
      dir = finddirection(searchtet, endpt);
      //assert(dir != ACROSSVERT);

      // Go to the face/edge intersecting the searching edge.
      enextesymself(*searchtet); // Go to the opposite face.
      // This face/edge has been tried in previous step.

      while (1) { // Loop I-I

        // Find the next intersecting face/edge.
        fsymself(*searchtet);
        if (dir == ACROSSFACE) {
          neightet = *searchtet;
          j = (neightet.ver & 3); // j is the current face number.
          for (i = j + 1; i < j + 4; i++) {
            neightet.ver = (i % 4);
            pa = org(neightet);
            pb = dest(neightet);
            pc = apex(neightet);
            pd = oppo(neightet); // The above point.
            if (tri_edge_test(pa,pb,pc,startpt,endpt, pd, 1, types, poss)) {
              dir = (enum interresult) types[0];
              pos = poss[0];
              break;
            } else {
              dir = DISJOINT;
              pos = 0;
            }
          } // i
          // There must be an intersection face/edge.
          assert(dir != DISJOINT);  // SELF_CHECK
        } else {
          assert(dir == ACROSSEDGE);
          while (1) { // Loop I-I-I
            // Check the two opposite faces (of the edge) in 'searchtet'.  
            for (i = 0; i < 2; i++) {
              if (i == 0) {
                enextesym(*searchtet, neightet);
              } else {
                eprevesym(*searchtet, neightet);
              }
              pa = org(neightet);
              pb = dest(neightet);
              pc = apex(neightet);
              pd = oppo(neightet); // The above point.
              if (tri_edge_test(pa,pb,pc,startpt,endpt,pd,1, types, poss)) {
                dir = (enum interresult) types[0];
                pos = poss[0];
                break; // for loop
              } else {
                dir = DISJOINT;
                pos = 0;
              }
            } // i
            if (dir != DISJOINT) {
              // Find an intersection face/edge.
              break;  // Loop I-I-I
            }
            // No intersection. Rotate to the next tet at the edge.
            fnextself(*searchtet);
          } // while (1) // Loop I-I-I
        }

        // Adjust to the intersecting edge/vertex.
        for (i = 0; i < pos; i++) {
          enextself(neightet);
        }

        if (dir == SHAREVERT) {
          // Check if we have reached the 'endpt'.
          pd = org(neightet);
          if (pd == endpt) {
            // Failed to recover the edge.
            break; // Loop I-I
          } else {
            // We need to further check this case. It might be a PLC problem
            //   or a Steiner point that was added at a bad location.
            assert(0);
          }
        }

        // The next to be flipped face/edge.
        *searchtet = neightet;

        // Bakup this face (tetrahedron).
        bakface.forg = org(*searchtet);
        bakface.fdest = dest(*searchtet);
        bakface.fapex = apex(*searchtet);
        bakface.foppo = oppo(*searchtet);

        // Try to flip this intersecting face/edge.
        if (dir == ACROSSFACE) {
          if (removefacebyflips(searchtet, &fc)) {
            success = 1;
            break; // Loop I-I 
          }
        } else if (dir == ACROSSEDGE) {
          if (removeedgebyflips(searchtet, &fc) == 2) {
            success = 1;
            break; // Loop I-I
          }
        } else {
          assert(0); // A PLC problem.
        }

        // The face/edge is not flipped.
        if ((searchtet->tet == NULL) ||
            (org(*searchtet) != bakface.forg) ||
            (dest(*searchtet) != bakface.fdest) ||
            (apex(*searchtet) != bakface.fapex) ||
            (oppo(*searchtet) != bakface.foppo)) {
          // 'searchtet' was flipped. We must restore it.
          point2tetorg(bakface.forg, *searchtet);
          dir1 = finddirection(searchtet, bakface.fdest);
          if (dir1 == ACROSSVERT) {
            assert(dest(*searchtet) == bakface.fdest);
            spintet = *searchtet;
            while (1) {
              if (apex(spintet) == bakface.fapex) {
                // Found the face.
                *searchtet = spintet;
                break;
              }
              fnextself(spintet);
              if (spintet.tet == searchtet->tet) {
                searchtet->tet = NULL;
                break; // Not find.
              }
	    } // while (1)
            if (searchtet->tet != NULL) {
              if (oppo(*searchtet) != bakface.foppo) {
                fsymself(*searchtet);
                if (oppo(*searchtet) != bakface.foppo) {
                  assert(0); // Check this case.
                  searchtet->tet = NULL;
                  break; // Not find.
                }
              }
            }
          } else {
            searchtet->tet = NULL; // Not find.
          }
          if (searchtet->tet == NULL) {
            success = 0; // This face/edge has been destroed.
            break; // Loop I-I 
          }
        }
      } // while (1) // Loop I-I

      if (success) {
        // One of intersecting faces/edges is flipped.
        continue;
      }

    } // if (fullsearch)

    // The edge is missing.
    break; // Loop I

  } // while (1) // Loop I

  // The edge is not recovered.
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// add_steinerpt_in_schoenhardtpoly()    Insert a Steiner point in a Schoen- //
//                                       hardt polyhedron.                   //
//                                                                           //
// 'abtets' is an array of n tets which all share at the edge [a,b]. Let the //
// tets are [a,b,p0,p1], [a,b,p1,p2], ..., [a,b,p_(n-2),p_(n-1)].  Moreover, //
// the edge [p0,p_(n-1)] intersects all of the tets in 'abtets'.  A special  //
// case is that the edge [p0,p_(n-1)] is coplanar with the edge [a,b].       //
// Such set of tets arises when we want to recover an edge from 'p0' to 'p_  //
// (n-1)', and the number of tets at [a,b] can not be reduced by any flip.   //
//                                                                           //
// The union of these tets is a polyhedron P. Obviously that P is a star-    //
// shaped polyhedron. The midpoint of [a,b] is visible by all boundary faces //
// of P, push it slightly inside P does not change the visibilty. Indeed     //
// every interior point of [a,b] is visible by the boundary faces of P.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::add_steinerpt_in_schoenhardtpoly(triface *abtets, int n,
                                                 int chkencflag)
{
  triface worktet, *parytet;
  triface faketet1, faketet2;
  point pa, pb, pc, pd;
  point p1, p2, p3;
  point steinerpt;
  insertvertexflags ivf;
  optparameters opm;
  REAL vcd[3], sampt[3], smtpt[3];
  REAL maxminvol = 0.0, minvol = 0.0, ori;
  int success, maxidx = 0;
  int loc;
  int it, i;

  if (b->verbose > 2) {
    printf("      Find a Steiner in Schoenhardt polyhedron (n=%d).\n", n);
  }

  pa = org(abtets[0]);
  pb = dest(abtets[0]);
  pc = apex(abtets[0]);   // pc = p0
  pd = oppo(abtets[n-1]); // pd = p_(n-1)

  // Find an optimial point in edge [c,d]. It is visible by all outer faces
  //   of 'abtets', and it maxmizes the min volume.

  // initialize the list of 2n boundary faces.
  for (i = 0; i < n; i++) {    
    eprev(abtets[i], worktet);
    esymself(worktet); // [a,p_i,p_i+1].
    cavetetlist->newindex((void **) &parytet);
    *parytet = worktet;
    enext(abtets[i], worktet);
    esymself(worktet); // [p_i,b,p_i+1].
    cavetetlist->newindex((void **) &parytet);
    *parytet = worktet;
  }

  // Search the point along the edge [c,d].
  for (i = 0; i < 3; i++) vcd[i] = pd[i] - pc[i];

  // Sample 100 points in edge [c,d].
  for (it = 1; it < 100; it++) {
    for (i = 0; i < 3; i++) {
      sampt[i] = pc[i] + (0.01 * (double) it) * vcd[i];
    }
    for (i = 0; i < cavetetlist->objects; i++) {
      parytet = (triface *) fastlookup(cavetetlist, i);
      p1 = org(*parytet);
      p2 = dest(*parytet);
      p3 = apex(*parytet);
      ori = orient3d(p2, p1, p3, sampt);
      if (i == 0) {
        minvol = ori;
      } else {
        if (minvol > ori) minvol = ori;
      }
    } // i
    if (it == 1) {
      maxminvol = minvol;
      maxidx = it;
    } else {
      if (maxminvol < minvol) {
        maxminvol = minvol;
        maxidx = it;
      } 
    }
  } // it

  if (maxminvol <= 0) {
    if (b->verbose > 2) {
      printf("      Unable to find a initial point: maxminvol = %g\n", 
             maxminvol);
    }
    cavetetlist->restart();
    return 0;
  }

  for (i = 0; i < 3; i++) {
    smtpt[i] = pc[i] + (0.01 * (double) maxidx) * vcd[i];
  }

  // Create two faked tets to hold the two non-existing boundary faces:
  //   [d,c,a] and [c,d,b].
  maketetrahedron(&faketet1);
  setvertices(faketet1, pd, pc, pa, dummypoint);
  cavetetlist->newindex((void **) &parytet);
  *parytet = faketet1;
  maketetrahedron(&faketet2);
  setvertices(faketet2, pc, pd, pb, dummypoint);
  cavetetlist->newindex((void **) &parytet);
  *parytet = faketet2;

  // Point smooth options.
  opm.max_min_volume = 1;
  opm.numofsearchdirs = 20;
  opm.searchstep = 0.001;  
  opm.maxiter = 100; // Limit the maximum iterations.
  opm.initval = 0.0; // Initial volume is zero.

  // Try to relocate the point into the inside of the polyhedron.
  success = smoothpoint(smtpt, cavetetlist, 1, &opm);

  if (success) {
    while (opm.smthiter == 100) {
      // It was relocated and the prescribed maximum iteration reached. 
      // Try to increase the search stepsize.
      opm.searchstep *= 10.0;
      //opm.maxiter = 100; // Limit the maximum iterations.
      opm.initval = opm.imprval;
      opm.smthiter = 0; // Init.
      smoothpoint(smtpt, cavetetlist, 1, &opm);  
    }
  } // if (success)

  // Delete the two faked tets.
  tetrahedrondealloc(faketet1.tet);
  tetrahedrondealloc(faketet2.tet);

  cavetetlist->restart();

  if (!success) {
    if (b->verbose > 2) {
      printf("      Unable to relocate the initial point.\n");
    }
    return 0;
  }


  // Insert the Steiner point.
  makepoint(&steinerpt, FREEVOLVERTEX);
  for (i = 0; i < 3; i++) steinerpt[i] = smtpt[i];

  // Insert the created Steiner point.
  for (i = 0; i < n; i++) {
    infect(abtets[i]);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = abtets[i];
  }
  worktet = abtets[0]; // No need point location.
  ivf.iloc = (int) INSTAR;
  ivf.bowywat = 0; // Do not use Bowyer-Watson algorithm.
  ivf.lawson = 0; //  Do not flip.
  ivf.rejflag = 0;
  ivf.chkencflag = chkencflag;
  ivf.sloc = 0;
  ivf.sbowywat = 0;
  ivf.splitbdflag = 0;
  ivf.validflag = 0;
  ivf.respectbdflag = 0;
  ivf.assignmeshsize = b->metric; 

  // Insert the new point into the tetrahedralization T.
  // Note that T is convex (nonconvex = 0).
  loc = insertvertex(steinerpt, &worktet, NULL, NULL, &ivf);

  if (loc == (int) INSTAR) {
    // The vertex has been inserted.
    st_volref_count++; 
    if (steinerleft > 0) steinerleft--;
    return 1;
  } else {
    // The Steiner point is too close to an existing vertex. Reject it.
    pointdealloc(steinerpt);
    return 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// addsteiner4recoversegment()    Add a Steiner point for recoveing a seg.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::addsteiner4recoversegment(face* misseg, int splitsegflag)
{
  triface *abtets, searchtet, spintet;
  face splitsh;
  face checkseg;
  face *paryseg;
  point startpt, endpt;
  point pa, pb, pd, steinerpt, *parypt;
  enum interresult dir;
  insertvertexflags ivf;
  int types[2], poss[4];
  REAL ip[3], u;
  int n, endi, success;
  int loc;
  int i;

  startpt = sorg(*misseg);
  if (pointtype(startpt) == FREESEGVERTEX) {
    sesymself(*misseg);
    startpt = sorg(*misseg);
  }
  endpt = sdest(*misseg);

  // Try to recover the edge by adding Steiner points.
  point2tetorg(startpt, searchtet);
  dir = finddirection(&searchtet, endpt);
  assert(dir != ACROSSVERT);

  // Get the first intersecting face/edge.
  assert(!ishulltet(searchtet));
  enextself(searchtet); 
  //assert(apex(searchtet) == startpt);

  if (dir == ACROSSFACE) {
    // The segment is crossing at least 3 faces. Find the common edge of 
    //   the first 3 crossing faces.
    esymself(searchtet);
    assert(oppo(searchtet) == startpt);
    fsym(searchtet, spintet);
    pd = oppo(spintet);
    if (pd == endpt) {
      // This should be possible.
      assert(0); // Debug this case.
    }
    for (i = 0; i < 3; i++) {
      pa = org(spintet);
      pb = dest(spintet);
      //pc = apex(neightet);
      if (tri_edge_test(pa, pb, pd, startpt, endpt, NULL, 1, types, poss)) {
        break; // Found the edge.
      }
      enextself(spintet);
      eprevself(searchtet);
    }
    assert(i < 3);
    esymself(searchtet);        
  } else {
    assert(dir == ACROSSEDGE);
    // PLC check.
    tsspivot1(searchtet, checkseg);
    if (checkseg.sh != NULL) {
      printf("Found two segments intersect each other.\n");
      pa = farsorg(*misseg);
      pb = farsdest(*misseg);
      printf("  1st: [%d,%d] %d.\n", pointmark(pa), pointmark(pb), 
             shellmark(*misseg));
      pa = farsorg(checkseg);
      pb = farsdest(checkseg);
      printf("  2nd: [%d,%d] %d.\n", pointmark(pa), pointmark(pb), 
             shellmark(checkseg));
      terminatetetgen(3);
    }
  }
  assert(apex(searchtet) == startpt);

  spintet = searchtet;
  n = 0; endi = -1;
  while (1) {
    // Check if the endpt appears in the star.
    if (apex(spintet) == endpt) {
      endi = n; // Remember the position of endpt.
    }
    n++; // Count a tet in the star.
    fnextself(spintet);
    if (spintet.tet == searchtet.tet) break;
  }
  assert(n >= 3);

  if (endi > 0) {
    // endpt is also in the edge star
    // Get all tets in the edge star.
    abtets = new triface[n];
    spintet = searchtet;
    for (i = 0; i < n; i++) {
      abtets[i] = spintet;
      fnextself(spintet);
    }
    assert(apex(abtets[0]) == startpt);
    assert(apex(abtets[endi]) == endpt);

    success = 0;

    if (dir == ACROSSFACE) {
      // Find a Steiner points inside the polyhedron.
      if (add_steinerpt_in_schoenhardtpoly(abtets, endi, 0)) {
        success = 1;
      }
    } else if (dir == ACROSSEDGE) {
      if (n > 4) {
        // In this case, 'abtets' is separated by the plane (containing the
        //   two intersecting edges) into two parts, P1 and P2, where P1
        //   consists of 'endi' tets: abtets[0], abtets[1], ..., 
        //   abtets[endi-1], and P2 consists of 'n - endi' tets: 
        //   abtets[endi], abtets[endi+1], abtets[n-1].
        if (endi > 2) { // P1
          // There are at least 3 tets in the first part.
          if (add_steinerpt_in_schoenhardtpoly(abtets, endi, 0)) {
            success++;
          }
        }
        if ((n - endi) > 2) { // P2
          // There are at least 3 tets in the first part.
          if (add_steinerpt_in_schoenhardtpoly(&(abtets[endi]), n - endi, 0)) {
            success++;
          }
        }
      } else {
        // In this case, a 4-to-4 flip should be re-cover the edge [c,d].
        //   However, there will be invalid tets (either zero or negtive 
        //   volume). Otherwise, [c,d] should already be recovered by the 
        //   recoveredge() function.
        assert(0); // DEBUG IT
      }
    } else {
      assert(0); // A PLC problem.
    }

    delete [] abtets;

    if (success) {
      // Add the missing segment back to the recovering list.
      subsegstack->newindex((void **) &paryseg);
      *paryseg = *misseg;
      return 1;
    }
  } // if (endi > 0)

  if (!splitsegflag) {
    return 0;
  }

  if (b->verbose > 2) {
    printf("      Splitting segment (%d, %d)\n", pointmark(startpt), 
           pointmark(endpt));
  }

  if (endi == -1) {
    // Let the missing segment be [a,b]. Let the edge [c,d] whose star contains
    // a and intersects [a,b]. We choose the Steiner point at the intersection
    // of the edge star of [c,d] and [a,b] (not a). 
    if (dir == ACROSSFACE) {
      pa = org(searchtet);
      pb = dest(searchtet);

      spintet = searchtet;
      n = 0; endi = -1;
      while (1) {
        n++; // Count a tet in the star.
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) break;
        // Check if the segment leaves the edge star.
        pd = apex(spintet);
        assert(pd != endpt);
        if (!tri_edge_test(pa, pb, pd, startpt, endpt, NULL, 1, types, poss)) {
          if (endi == -1)  endi = (n - 1);
        }
      }
      assert(n >= 3);
      assert(endi != -1); 

      // 'abtets' is only for debug purpose.
      abtets = new triface[endi];
      spintet = searchtet;
      for (i = 0; i < endi; i++) {
        abtets[i] = spintet;
        fnextself(spintet);
      }
      searchtet = abtets[endi - 1]; 
      esymself(searchtet); // The exit face of [startpt, endpt].
      delete [] abtets;
    } else {
      assert(dir == ACROSSEDGE);
      assert(apex(searchtet) == startpt);
      esymself(searchtet); // The exit face of [startpt, endpt].
      //assert(oppo(searchtet) == startpt);
      pa = org(searchtet);
      pb = dest(searchtet);
    }

    pd = apex(searchtet);
    // Get the intersection type (ACROSSFACE or ACROSSEDGE).
    if (tri_edge_test(pa, pb, pd, startpt, endpt, NULL, 1, types, poss)) {
      dir = (enum interresult) types[0];
      assert((dir == ACROSSFACE) || (dir == ACROSSEDGE));
    } else {
      assert(0); // not possible.
    }

    // Calculate the intersection of the face [a,b,d] and the segment.
    planelineint(pa, pb, pd, startpt, endpt, ip, &u);
    assert((u > 0) && (u < 1));

    // Create a Steiner point.
    makepoint(&steinerpt, FREESEGVERTEX);
    for (i = 0; i < 3; i++) steinerpt[i] = ip[i];


    spivot(*misseg, splitsh);
    if (dir == ACROSSFACE) {
      ivf.iloc = (int) ONFACE;
    } else {
      ivf.iloc = (int) ONEDGE;
    }
    ivf.bowywat = 1;
    ivf.lawson = 0;
    ivf.rejflag = 0;
    ivf.chkencflag = 0;
    ivf.sloc = (int) ONEDGE;
    ivf.sbowywat = 1;
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 1;
    ivf.assignmeshsize = b->metric;
    loc = insertvertex(steinerpt, &searchtet, &splitsh, misseg, &ivf);

    if (loc != ivf.iloc) {
      if (loc == (int) NEARVERTEX) {
        // The vertex is rejected. Too close to an existing vertex.
        pointdealloc(steinerpt);
        steinerpt = NULL;
      } else {
        assert(0); // Unknown case. 
      }
    }
  } else { // if (endi > 0)
    steinerpt = NULL;
  }

  if (steinerpt == NULL) {
    // Split the segment at its midpoint.
    makepoint(&steinerpt, FREESEGVERTEX);
    for (i = 0; i < 3; i++) {
      steinerpt[i] = 0.5 * (startpt[i] + endpt[i]);
    }

    // We need to locate the point.
    assert(searchtet.tet != NULL); // Start searching from 'searchtet'.
    spivot(*misseg, splitsh);
    ivf.iloc = (int) OUTSIDE;
    ivf.bowywat = 1;
    ivf.lawson = 0;
    ivf.rejflag = 0;
    ivf.chkencflag = 0;
    ivf.sloc = (int) ONEDGE;
    ivf.sbowywat = 1;
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 1;
    ivf.assignmeshsize = b->metric; 
    loc = insertvertex(steinerpt, &searchtet, &splitsh, misseg, &ivf);

    assert(loc != (int) ONVERTEX);
    assert(loc != (int) NEARVERTEX);
  } // if (endi > 0)

  // Save this Steiner point (for removal).
  //   Re-use the array 'subvertstack'.
  subvertstack->newindex((void **) &parypt);
  *parypt = steinerpt;

  st_segref_count++;
  if (steinerleft > 0) steinerleft--;

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoversegments()    Recover all segments.                                //
//                                                                           //
// All segments need to be recovered are in 'subsegstack'.                   //
//                                                                           //
// This routine first tries to recover each segment by only using flips. If  //
// no flip is possible, and the flag 'steinerflag' is set, it then tries to  //
// insert Steiner points near or in the segment.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::recoversegments(arraypool *misseglist, int fullsearch,
                                int steinerflag)
{
  triface searchtet, spintet;
  face sseg, checkseg, *paryseg;
  point startpt, endpt;
  int success;

  long bak_inpoly_count = st_volref_count; 

  if (b->verbose > 1) {
    printf("    Recover segments [%s level = %2d] #:  %ld.\n",
           (b->fliplinklevel > 0) ? "fixed" : "auto",
           (b->fliplinklevel > 0) ? b->fliplinklevel : autofliplinklevel,
           subsegstack->objects);
  }

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    paryseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *paryseg;

    // Check if this segment has been recovered.
    sstpivot1(sseg, searchtet);
    if (searchtet.tet != NULL) {
      continue; // Not a missing segment.
    }

    startpt = sorg(sseg);
    endpt = sdest(sseg);

    if (b->verbose > 2) {
      printf("      Recover segment (%d, %d).\n", pointmark(startpt), 
             pointmark(endpt));
    }

    success = 0;

    if (recoveredgebyflips(startpt, endpt, &searchtet, 0)) {
      success = 1;
    } else {
      // Try to recover it from the other direction.
      if (recoveredgebyflips(endpt, startpt, &searchtet, 0)) {
        success = 1;
      }
    }

    if (!success && fullsearch) {
      if (recoveredgebyflips(startpt, endpt, &searchtet, fullsearch)) {
        success = 1;
      }
    }

    if (success) {
      // Segment is recovered. Insert it.
      tsspivot1(searchtet, checkseg);  // SELF_CHECK
      assert(checkseg.sh == NULL);
      // Let the segment remember an adjacent tet.
      sstbond1(sseg, searchtet);
      // Bond the segment to all tets containing it.
      spintet = searchtet;
      do {
        tssbond1(spintet, sseg);
        fnextself(spintet);
      } while (spintet.tet != searchtet.tet);
    } else {
      if (steinerflag > 0) {
        // Try to recover the segment but do not split it.
        if (addsteiner4recoversegment(&sseg, 0)) {
          success = 1;
        }
        if (!success && (steinerflag > 1)) {
          // Split the segment.
          addsteiner4recoversegment(&sseg, 1);
          success = 1;
        }
      }
      if (!success) {
        if (misseglist != NULL) {
          // Save this segment.
          misseglist->newindex((void **) &paryseg);
          *paryseg = sseg;
        }
      }
    }

  } // while (subsegstack->objects > 0l)

  if (steinerflag) {
    if (b->verbose > 1) {
      // Report the number of added Steiner points.
      if (st_volref_count > bak_inpoly_count) {
        printf("    Add %ld Steiner points in volume.\n", 
               st_volref_count - bak_inpoly_count);
      }
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoverfacebyflips()    Recover a face by flips.                          //
//                                                                           //
// If 'searchsh' is not NULL, it is a subface to be recovered.  It is only   //
// used for checking self-intersections.                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::recoverfacebyflips(point pa, point pb, point pc, 
                                   face *searchsh, triface* searchtet)
{
  triface spintet, flipedge;
  face checkseg;
  point pd, pe;
  enum interresult dir;
  flipconstraints fc;
  int success, success1;
  int i, j;

  int intflag;
  int types[2], poss[4];

  if (b->verbose > 2) {
    printf("      Recovering face (%d, %d, %d) by flips\n", pointmark(pa), 
           pointmark(pb), pointmark(pc));
  }

  fc.fac[0] = pa;
  fc.fac[1] = pb;
  fc.fac[2] = pc;
  success = 0;

  for (i = 0; i < 3 && !success; i++) {
    while (1) {
      // Get a tet containing the edge [a,b].
      point2tetorg(fc.fac[i], *searchtet);
      dir = finddirection(searchtet, fc.fac[(i+1)%3]);
      //assert(dir == ACROSSVERT);
      assert(dest(*searchtet) == fc.fac[(i+1)%3]);
      // Search the face [a,b,c]
      spintet = *searchtet;
      while (1) {
        if (apex(spintet) == fc.fac[(i+2)%3]) {
          // Found the face.
          *searchtet = spintet;
          // Return the face [a,b,c].
          for (j = i; j > 0; j--) {
            eprevself(*searchtet);
          }
          success = 1;
          break;
        }
        fnextself(spintet);
        if (spintet.tet == searchtet->tet) break;
      } // while (1)
      if (success) break;
      // The face is missing. Try to recover it.
      success1 = 0;
      // Find a crossing edge of this face.
      spintet = *searchtet;
      while (1) {
        pd = apex(spintet);
        pe = oppo(spintet);
        if ((pd != dummypoint) && (pe != dummypoint)) {
          // Check if [d,e] intersects [a,b,c]
          intflag = tri_edge_test(pa, pb, pc, pd, pe, NULL, 1, types, poss);
          if (intflag > 0) {
            // By our assumptions, they can only intersect at a single point.
            if (intflag == 2) {
              // Check the intersection type.
              dir = (enum interresult) types[0];
              if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
                // Go to the edge [d,e].
                eprev(spintet, flipedge);
                esymself(flipedge);
                enextself(flipedge); // [d,e,a,b].
                if (searchsh != NULL) {
                  // Check if [e,d] is a segment.
                  tsspivot1(flipedge, checkseg);
                  if (checkseg.sh != NULL) {
                    if (!b->quiet) {                    
                      printf("Found a segment and a subface intersect.\n");
                      pd = farsorg(checkseg);
                      pe = farsdest(checkseg);
                      printf("  1st: [%d, %d] %d.\n",  pointmark(pd), 
                             pointmark(pe), shellmark(checkseg)); 
                      printf("  2nd: [%d,%d,%d] %d\n", pointmark(pa), 
                        pointmark(pb), pointmark(pc), shellmark(*searchsh));
	            }
                    terminatetetgen(3);
		  }
                }
                // Try to flip the edge [d,e].
                success1 = (removeedgebyflips(&flipedge, &fc) == 2);
              } else {
                if (dir == TOUCHFACE) {
                  point touchpt, *parypt;
                  if (poss[0] == 0) {
                    touchpt = pd; // pd is a coplanar vertex.
                  } else {
                    touchpt = pe; // pe is a coplanar vertex.
                  }
                  if (pointtype(touchpt) == FREEVOLVERTEX) {
                    // A volume Steiner point was added in this subface.
                    // Split this subface by this point.
                    if (b->verbose > 2) {
                      printf("      Shift volume Steiner point %d to facet.\n",
                             pointmark(touchpt));
                    }
                    face checksh, *parysh;
                    int siloc = (int) ONFACE;
                    int sbowat = 0; // Only split this subface.

                    sinsertvertex(touchpt, searchsh, NULL, siloc, sbowat);

                    setpointtype(touchpt, FREEFACETVERTEX);
                    st_volref_count--;
                    st_facref_count++;
                    // Queue this vertex for removal.
                    subvertstack->newindex((void **) &parypt);
                    *parypt = touchpt;
                    // Queue new subfaces for recovery.
                    // Put all new subfaces into stack for recovery.
                    for (i = 0; i < caveshbdlist->objects; i++) {
                      // Get an old subface at edge [a, b].
                      parysh = (face *) fastlookup(caveshbdlist, i);
                      spivot(*parysh, checksh); // The new subface [a, b, p].
                      // Do not recover a deleted new face (degenerated).
                      if (checksh.sh[3] != NULL) {
                        if (b->verbose > 3) {
                          printf("        Queue new subface (%d, %d, %d).\n",
                            pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                            pointmark(sapex(checksh)));
                        }
                        //sdissolve(checksh); // It has not been connected yet.
                        subfacstack->newindex((void **) &parysh);
                        *parysh = checksh;
                      }
                    }
                    // Delete the old subfaces in sC(p).
                    assert(caveshlist->objects == 1);
                    for (i = 0; i < caveshlist->objects; i++) {
                      parysh = (face *) fastlookup(caveshlist, i);
                      shellfacedealloc(subfaces, parysh->sh);
                    }
                    // Clear working lists.
                    caveshlist->restart();
                    caveshbdlist->restart();
                    cavesegshlist->restart();
                    // We can return this function.
                    searchsh->sh = NULL; // It has been split.
                    success1 = 0;
                    success = 1; 
                  } else {
                    // It should be a PLC problem.
                    if (pointtype(touchpt) == FREESEGVERTEX) {
                      // A segment and a subface intersect. 
                    } else if (pointtype(touchpt) == FREEFACETVERTEX) {
                      // Two facets self-intersect.
                    }
                    terminatetetgen(3);
                  }
                } else {
                  assert(0); // Unknown cases. Debug.
                }
              }
              break;
            } else { // intflag == 4. Coplanar case.
              // This may be an input PLC error.
              assert(0);
            }
          } // if (intflag > 0)
        }
        fnextself(spintet);
        assert(spintet.tet != searchtet->tet);
      } // while (1)
      if (!success1) break;
    } // while (1)
  } // i

  return success;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoversubfaces()    Recover all subfaces.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::recoversubfaces(arraypool *misshlist, int steinerflag)
{
  triface searchtet, neightet, spintet;
  face searchsh, neighsh, neineish, *parysh;
  face bdsegs[3], checkseg;
  point startpt, endpt, apexpt, *parypt;
  point steinerpt;
  enum interresult dir;
  insertvertexflags ivf;
  int success;
  int loc;
  int i, j;

  if (b->verbose > 1) {
    printf("    Recover subfaces [%s level = %2d] #:  %ld.\n",
           (b->fliplinklevel > 0) ? "fixed" : "auto",
           (b->fliplinklevel > 0) ? b->fliplinklevel : autofliplinklevel,
           subfacstack->objects);
  }

  // Loop until 'subfacstack' is empty.
  while (subfacstack->objects > 0l) {

    subfacstack->objects--;
    parysh = (face *) fastlookup(subfacstack, subfacstack->objects);
    searchsh = *parysh;

    if (searchsh.sh[3] == NULL) continue; // Skip a dead subface.

    stpivot(searchsh, neightet);
    if (neightet.tet != NULL) continue; // Skip a recovered subface.


    if (b->verbose > 2) {
      printf("      Recover subface (%d, %d, %d).\n", pointmark(sorg(searchsh)),
             pointmark(sdest(searchsh)), pointmark(sapex(searchsh)));
    }

    // The three edges of the face need to be existed first.
    for (i = 0; i < 3; i++) {
      sspivot(searchsh, bdsegs[i]);   
      if (bdsegs[i].sh != NULL) {
        // The segment must exist.
        sstpivot1(bdsegs[i], searchtet);
        if (searchtet.tet == NULL) {
          assert(0);
        }
      } else {
        // This edge is not a segment (due to a Steiner point).
        // Check whether it exists or not.
        success = 0;
        startpt = sorg(searchsh);
        endpt = sdest(searchsh);
        point2tetorg(startpt, searchtet);
        assert(org(searchtet) == startpt); // SELF_CHECK
        dir = finddirection(&searchtet, endpt);
        if (dir == ACROSSVERT) {
          if (dest(searchtet) == endpt) {
            success = 1;  
          } else {
            //assert(0); // A PLC problem.
            terminatetetgen(3);
          }
        } else {
          // The edge is missing. Try to recover it.
          if (recoveredgebyflips(startpt, endpt, &searchtet, 0)) {
            success = 1;
          } else {
            if (recoveredgebyflips(endpt, startpt, &searchtet, 0)) {
              success = 1;
            }
          }
        }
        if (success) {
          // Insert a temporary segment to protect this edge.
          if (b->verbose > 2) {
            printf("      Insert a temp segment to protect edge [%d, %d].\n",
                   pointmark(startpt), pointmark(endpt));
          }
          makeshellface(subsegs, &(bdsegs[i]));
          setshvertices(bdsegs[i], startpt, endpt, NULL);
          //setshellmark(bdsegs[i], -2); // It's a temporary segment.
          smarktest2(bdsegs[i]); // It's a temporary segment.
          // Insert this segment into surface mesh.
          ssbond(searchsh, bdsegs[i]);
          spivot(searchsh, neighsh);
          if (neighsh.sh != NULL) {
            ssbond(neighsh, bdsegs[i]);
          }
          // Insert this segment into tetrahedralization.
          tsspivot1(searchtet, checkseg);  // SELF_CHECK
          assert(checkseg.sh == NULL);
          sstbond1(bdsegs[i], searchtet);
          // Bond the segment to all tets containing it.
          spintet = searchtet;
          do {
            tssbond1(spintet, bdsegs[i]);
            fnextself(spintet);
          } while (spintet.tet != searchtet.tet);
        } else {
          // An edge of this subface is missing. Can't recover this subface.
          // Delete any temporary segment that has been created.
          for (j = (i - 1); j >= 0; j--) {
            if (smarktest2ed(bdsegs[j])) { // if (shellmark(bdsegs[j]) == -2) {
              if (b->verbose > 2) {
                printf("      Remove a temp segment (%d, %d).\n", 
                  pointmark(sorg(bdsegs[j])), pointmark(sdest(bdsegs[j])));
              }
              spivot(bdsegs[j], neineish);
              assert(neineish.sh != NULL);
              //if (neineish.sh != NULL) {
                ssdissolve(neineish);
                spivot(neineish, neighsh);
                if (neighsh.sh != NULL) {
                  ssdissolve(neighsh);
                  // There should be only two subfaces at this segment.
                  spivotself(neighsh); // SELF_CHECK
                  assert(neighsh.sh == neineish.sh);
                }
	      //}
              sstpivot1(bdsegs[j], searchtet);
              assert(searchtet.tet != NULL);
              //if (searchtet.tet != NULL) {
                spintet = searchtet;
                while (1) {
                  tssdissolve1(spintet);
                  fnextself(spintet);
                  if (spintet.tet == searchtet.tet) break;
                }
	      //}
              shellfacedealloc(subsegs, bdsegs[j].sh);
            }
          } // j
          if (steinerflag) {
            // Add a Steiner point at the midpoint of this edge.
            if (b->verbose > 2) {
              printf("      Add a Steiner point in subedge (%d, %d).\n",
                     pointmark(startpt), pointmark(endpt));
            }
            makepoint(&steinerpt, FREEFACETVERTEX);
            for (j = 0; j < 3; j++) {
              steinerpt[j] = 0.5 * (startpt[j] + endpt[j]);
            }

            point2tetorg(startpt, searchtet); // Start from 'searchtet'.
            ivf.iloc = (int) OUTSIDE; // Need point location.
            ivf.bowywat = 1;
            ivf.lawson = 0;
            ivf.rejflag = 0;
            ivf.chkencflag = 0;
            ivf.sloc = (int) ONEDGE;            
            ivf.sbowywat = 1; // Allow flips in facet.
            ivf.splitbdflag = 0;
            ivf.validflag = 1;
            ivf.respectbdflag = 1;
            ivf.assignmeshsize = b->metric;
            loc = insertvertex(steinerpt, &searchtet, &searchsh, NULL, &ivf);
            assert(loc != (int) OUTSIDE);

            // Save this Steiner point (for removal).
            //   Re-use the array 'subvertstack'.
            subvertstack->newindex((void **) &parypt);
            *parypt = steinerpt;

            st_facref_count++;
            if (steinerleft > 0) steinerleft--;
          } // if (steinerflag)
          break;
        }
      }
      senextself(searchsh);
    } // i

    if (i == 3) {
      // Recover the subface.
      startpt = sorg(searchsh);
      endpt   = sdest(searchsh);
      apexpt  = sapex(searchsh);

      success = recoverfacebyflips(startpt,endpt,apexpt,&searchsh,&searchtet);

      // Delete any temporary segment that has been created.
      for (j = 0; j < 3; j++) {
        if (smarktest2ed(bdsegs[j])) { //if (shellmark(bdsegs[j]) == -2) {
          if (b->verbose > 2) {
            printf("      Remove a temp segment (%d, %d).\n", 
                   pointmark(sorg(bdsegs[j])), pointmark(sdest(bdsegs[j])));
          }
          spivot(bdsegs[j], neineish);
          assert(neineish.sh != NULL);
          //if (neineish.sh != NULL) {
            ssdissolve(neineish);
            spivot(neineish, neighsh);
            if (neighsh.sh != NULL) {
              ssdissolve(neighsh);
              // There should be only two subfaces at this segment.
              spivotself(neighsh); // SELF_CHECK
              assert(neighsh.sh == neineish.sh);
            }
	  //}
          sstpivot1(bdsegs[j], neightet);
          assert(neightet.tet != NULL);
          //if (neightet.tet != NULL) {
            spintet = neightet;
            while (1) {
              tssdissolve1(spintet);
              fnextself(spintet);
              if (spintet.tet == neightet.tet) break;
            }
	  //}
          shellfacedealloc(subsegs, bdsegs[j].sh);
        }
      } // j

      if (success) {
        if (searchsh.sh != NULL) {
          // Face is recovered. Insert it.
          tsbond(searchtet, searchsh);
          fsymself(searchtet);
          sesymself(searchsh);
          tsbond(searchtet, searchsh);
        }
      } else {
        if (steinerflag) {
          // Add a Steiner point at the barycenter of this subface.
          if (b->verbose > 2) {
            printf("      Add a Steiner point in subface (%d, %d, %d).\n",
                   pointmark(startpt), pointmark(endpt), pointmark(apexpt));
          }
          makepoint(&steinerpt, FREEFACETVERTEX);
          for (j = 0; j < 3; j++) {
            steinerpt[j] = (startpt[j] + endpt[j] + apexpt[j]) / 3.0;
          }

          point2tetorg(startpt, searchtet); // Start from 'searchtet'.
          ivf.iloc = (int) OUTSIDE; // Need point location.
          ivf.bowywat = 1;
          ivf.lawson = 0;
          ivf.rejflag = 0;
          ivf.chkencflag = 0;
          ivf.sloc = (int) ONFACE;          
          ivf.sbowywat = 1; // Allow flips in facet.
          ivf.splitbdflag = 0;
          ivf.validflag = 1;
          ivf.respectbdflag = 1;
          ivf.assignmeshsize = b->metric; 
          loc = insertvertex(steinerpt, &searchtet, &searchsh, NULL, &ivf);
          assert(loc != (int) OUTSIDE);

          // Save this Steiner point (for removal).
          //   Re-use the array 'subvertstack'.
          subvertstack->newindex((void **) &parypt);
          *parypt = steinerpt;

          st_facref_count++;
          if (steinerleft > 0) steinerleft--;
        } // if (steinerflag)
      }
    } else {
      success = 0;      
    }

    if (!success) {
      if (misshlist != NULL) {
        if (b->verbose > 2) {
          printf("      Subface (%d, %d, %d) is missing.\n", 
                 pointmark(sorg(searchsh)), pointmark(sdest(searchsh)), 
                 pointmark(sapex(searchsh)));
        }
        // Save this subface.
        misshlist->newindex((void **) &parysh);
        *parysh = searchsh;
      }
    }

  } // while (subfacstack->objects > 0l)

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getvertexstar()    Return the star of a vertex.                           //
//                                                                           //
// If the flag 'fullstar' is set, return the complete star of this vertex.   //
// Otherwise, only a part of the star which is bounded by facets is returned.// 
//                                                                           //
// 'tetlist' returns the list of tets in the star of the vertex 'searchpt'.  //
// Every tet in 'tetlist' is at the face oppsiting to 'searchpt'.            //
//                                                                           //
// 'vertlist' returns the list of vertices in the star (exclude 'searchpt'). //
//                                                                           //
// 'shlist' returns the list of subfaces in the star. Each subface must face //
// to the interior of this star.                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::getvertexstar(int fullstar, point searchpt, arraypool* tetlist, 
                              arraypool* vertlist, arraypool* shlist)
{
  triface searchtet, neightet, *parytet;
  face checksh, *parysh;
  //face checkseg;
  point pt, *parypt;
  int collectflag;
  int i, j;

  if (b->verbose > 2) {
    printf("      Form the star of vertex %d.\n", pointmark(searchpt));
  }

  point2tetorg(searchpt, searchtet);

  // Go to the opposite face (the link face) of the vertex.
  enextself(searchtet);
  esymself(searchtet);
  //assert(oppo(searchtet) == searchpt);
  infect(searchtet); // Collect this tet (link face).
  tetlist->newindex((void **) &parytet);
  *parytet = searchtet;
  if (vertlist != NULL) {
    // Collect three (link) vertices.
    for (i = 0; i < 3; i++) {
      pt = org(searchtet);
      pinfect(pt);
      vertlist->newindex((void **) &parypt);
      *parypt = pt;
      enextself(searchtet);
    }
  }

  collectflag = 1;
  esym(searchtet, neightet);
  tspivot(neightet, checksh);
  if (checksh.sh != NULL) {
    if (shlist != NULL) {
      if (!sinfected(checksh)) {
        // Collect this subface (link edge).
        sinfected(checksh);
        shlist->newindex((void **) &parysh);
        *parysh = checksh;
      }
    } // if (checksh.sh != NULL)
    if (!fullstar) {
      collectflag = 0;
    }
  }
  if (collectflag) {
    fsymself(neightet); // Goto the adj tet of this face.
    assert(neightet.tet != NULL);
    esymself(neightet); // Goto the oppo face of this vertex.
    // assert(oppo(neightet) == searchpt);
    infect(neightet); // Collect this tet (link face).
    tetlist->newindex((void **) &parytet);
    *parytet = neightet;
    if (vertlist != NULL) {
      // Collect its apex.
      pt = apex(neightet);
      pinfect(pt);
      vertlist->newindex((void **) &parypt);
      *parypt = pt;
    }
  } // if (collectflag)

  // Continue to collect all tets in the star.
  for (i = 0; i < tetlist->objects; i++) {
    searchtet = * (triface *) fastlookup(tetlist, i);
    // Note that 'searchtet' is a face opposite to 'searchpt', and the neighbor
    //   tet at the current edge is already collected.
    // Check the neighors at the other two edges of this face.
    for (j = 0; j < 2; j++) {
      collectflag = 1;
      enextself(searchtet);
      //fnext(searchtet, neightet);
      esym(searchtet, neightet);
      tspivot(neightet, checksh);
      if (checksh.sh != NULL) {
        if (shlist != NULL) {
          if (!sinfected(checksh)) {
            // Collect this subface (link edge).
            sinfected(checksh);
            shlist->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
        if (!fullstar) {
          collectflag = 0;
        }
      }
      if (collectflag) {
        fsymself(neightet);
        assert(neightet.tet != NULL);
        if (!infected(neightet)) {
          esymself(neightet); // Go to the face opposite to 'searchpt'.
          infect(neightet);
          tetlist->newindex((void **) &parytet);
          *parytet = neightet;
          if (vertlist != NULL) {
            // Check if a vertex is collected.
            pt = apex(neightet);
            if (!pinfected(pt)) {
              pinfect(pt);
              vertlist->newindex((void **) &parypt);
              *parypt = pt;
            }
          }
        } // if (!infected(neightet))
      } // if (collectflag)
    } // j
  } // i

  if (b->verbose > 2) {
    printf("      Collected %ld tets", tetlist->objects);
    if (vertlist != NULL) {
      printf(", %ld vertices", vertlist->objects);
    }
    if (shlist != NULL) {
      printf(", %ld subfaces", shlist->objects);
    }
    printf(".\n");
  }

  // Uninfect the list of tets and vertices.
  for (i = 0; i < tetlist->objects; i++) {
    parytet = (triface *) fastlookup(tetlist, i);
    uninfect(*parytet);
  }

  if (vertlist != NULL) {
    for (i = 0; i < vertlist->objects; i++) {
      parypt = (point *) fastlookup(vertlist, i);
      puninfect(*parypt);
    }
  }

  if (shlist != NULL) {
    for (i = 0; i < shlist->objects; i++) {
      parysh = (face *) fastlookup(shlist, i);
      suninfect(*parysh);
    }
  }

  return (int) tetlist->objects;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getedge()    Get a tetrahedron having the two endpoints.                  //
//                                                                           //
// The method here is to search the second vertex in the link faces of the   //
// first vertex. The global array 'cavetetlist' is re-used for searching.    //
//                                                                           //
// This function is used for the case when the mesh is non-convex. Otherwise,//
// the function finddirection() should be faster than this.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::getedge(point e1, point e2, triface *tedge)
{
  triface searchtet, neightet, *parytet;
  point pt;
  int done;
  int i, j;

  if (b->verbose > 2) {
    printf("      Get edge from %d to %d.\n", pointmark(e1), pointmark(e2));
  }

  // Quickly check if 'tedge' is just this edge.
  if (!isdeadtet(*tedge)) {
    if (org(*tedge) == e1) {
      if (dest(*tedge) == e2) {
        return 1;
      }
    } else if (org(*tedge) == e2) {
      if (dest(*tedge) == e1) {
        esymself(*tedge);
        return 1;
      }
    }
  }

  // Search for the edge [e1, e2].
  point2tetorg(e1, *tedge);
  finddirection(tedge, e2);
  if (dest(*tedge) == e2) {
    return 1;
  } else {
    // Search for the edge [e2, e1].
    point2tetorg(e2, *tedge);
    finddirection(tedge, e1);
    if (dest(*tedge) == e1) {
      esymself(*tedge);
      return 1;
    }
  }


  // Go to the link face of e1.
  point2tetorg(e1, searchtet);
  enextself(searchtet);
  esymself(searchtet);
  //assert(oppo(searchtet) == e1);

  assert(cavetetlist->objects == 0l); // It will re-use this list.

  // Search e2.
  for (i = 0; i < 3; i++) {
    pt = apex(searchtet);
    if (pt == e2) {
      // Found. 'searchtet' is [#,#,e2,e1].
      enext(searchtet, *tedge);
      esymself(*tedge);
      eprevself(*tedge); // [e1,e2,#,#].
      return 1;
    }
    enextself(searchtet);
  }

  // Get the adjacent link face at 'searchtet'.
  fnext(searchtet, neightet);
  esymself(neightet);
  // assert(oppo(neightet) == e1);
  pt = apex(neightet);
  if (pt == e2) {
    // Found. 'neightet' is [#,#,e2,e1].
    enext(neightet, *tedge);
    esymself(*tedge);
    eprevself(*tedge); // [e1,e2,#,#].
    return 1;
  }

  // Continue searching in the link face of e1.
  infect(searchtet);
  cavetetlist->newindex((void **) &parytet);
  *parytet = searchtet;
  infect(neightet);
  cavetetlist->newindex((void **) &parytet);
  *parytet = neightet;

  done = 0;

  for (i = 0; (i < cavetetlist->objects) && !done; i++) {
    parytet = (triface *) fastlookup(cavetetlist, i);
    searchtet = *parytet;
    for (j = 0; (j < 2) && !done; j++) {
      enextself(searchtet);
      fnext(searchtet, neightet);
      if (!infected(neightet)) {        
        esymself(neightet);
        pt = apex(neightet);
        if (pt == e2) {
          // Found. 'neightet' is [#,#,e2,e1].
          enext(neightet, *tedge);
          esymself(*tedge);
          eprevself(*tedge); // [e1,e2,#,#].
          done = 1;
        } else {
          infect(neightet);
          cavetetlist->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
    } // j
  } // i 

  // Uninfect the list of visited tets.
  for (i = 0; i < cavetetlist->objects; i++) {
    parytet = (triface *) fastlookup(cavetetlist, i);
    uninfect(*parytet);
  }
  cavetetlist->restart();

  return done;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// reduceedgesatvertex()    Reduce the number of edges at a given vertex.    //
//                                                                           //
// 'endptlist' contains the endpoints of edges connecting at the vertex.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::reduceedgesatvertex(point startpt, arraypool* endptlist)
{
  triface searchtet;
  face checkseg;
  point *pendpt, *parypt;
  enum interresult dir;
  flipconstraints fc;
  int reduceflag;
  int count;
  int n, i, j;

  if (b->verbose > 2) {
    printf("      Initial edge degree = %ld.\n", endptlist->objects);
  }
  assert(endptlist->objects >= 4l);

  // Reduce the number of edges.
  fc.remvert = startpt;

  while (1) {

    count = 0;

    for (i = 0; i < endptlist->objects; i++) {
      pendpt = (point *) fastlookup(endptlist, i);
      if (*pendpt == dummypoint) {
        continue; // Do not reduce a virtual edge.
      }
      reduceflag = 0;
      // Find the edge.
      if (nonconvex) {
        if (getedge(startpt, *pendpt, &searchtet)) {
          dir = ACROSSVERT;
        } else {
          // The edge does not exist (was flipped).
          dir = INTERSECT;
        }
      } else {
        point2tetorg(startpt, searchtet);
        dir = finddirection(&searchtet, *pendpt);
      }
      if (dir == ACROSSVERT) {
        if (dest(searchtet) == *pendpt) {
          // Do not flip a segment.
          tsspivot1(searchtet, checkseg);
          if (checkseg.sh == NULL) {
            n = removeedgebyflips(&searchtet, &fc);
            if (n == 2) {
              reduceflag = 1;
            }
          }
        } else {
          assert(0); // A plc problem.
        }
      } else {
        // The edge has been flipped.
        reduceflag = 1;
      }
      if (reduceflag) {
        count++;
        // Move the last vertex into this slot.
        j = endptlist->objects - 1;
        parypt = (point *) fastlookup(endptlist, j);
        *pendpt = *parypt;
        endptlist->objects--;
        i--;
      }
    } // i

    if (count == 0) {
      // No edge is reduced.
      break;
    }

  } // while (1)

  if (b->verbose > 2) {
    printf("      Final edge degree = %ld.\n", endptlist->objects);
  }

  return (int) endptlist->objects;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// removevertexbyflips()    Remove a vertex by flips.                        //
//                                                                           //
// This routine attempts to remove the given vertex 'rempt' (p) from the cur-//
// rent tetrahedralization (T) by a sequence of elementary flips.            //
//                                                                           //
// The algorithm used here is a simple edge reduce method. Suppose there are //
// n edges connected at p. We try to reduce the number of edges by flipping  //
// any edge (not a segment) that is connecting at p.                         //
//                                                                           //
// Unless T is a Delaunay tetrahedralization, there is no guarantee that 'p' //
// can be successfully removed.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::removevertexbyflips(point steinerpt)
{
  triface *fliptets = NULL, wrktets[4];
  triface searchtet, spintet, neightet;
  face parentsh, spinsh, checksh;
  face leftseg, rightseg, checkseg;
  point lpt = NULL, rpt = NULL, apexpt, *parypt;
  enum verttype vt;
  enum locateresult loc;
  int valence, removeflag;
  int slawson;
  int n, i;

  vt = pointtype(steinerpt);

  if (vt == FREESEGVERTEX) {
    sdecode(point2sh(steinerpt), leftseg);
    assert(leftseg.sh != NULL);
    leftseg.shver = 0;
    if (sdest(leftseg) == steinerpt) {
      senext(leftseg, rightseg);
      spivotself(rightseg);
      assert(rightseg.sh != NULL);
      rightseg.shver = 0;
      assert(sorg(rightseg) == steinerpt);
    } else {
      assert(sorg(leftseg) == steinerpt);
      rightseg = leftseg;
      senext2(rightseg, leftseg);
      spivotself(leftseg);
      assert(leftseg.sh != NULL);
      leftseg.shver = 0;
      assert(sdest(leftseg) == steinerpt);
    }
    lpt = sorg(leftseg);
    rpt = sdest(rightseg);
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in segment (%d, %d).\n",
             pointmark(steinerpt), pointmark(lpt), pointmark(rpt));

    }
  } else if (vt == FREEFACETVERTEX) {
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in facet.\n",
             pointmark(steinerpt));
    }
  } else if (vt == FREEVOLVERTEX) {
    if (b->verbose > 2) {
      printf("      Removing Steiner point %d in volume.\n",
             pointmark(steinerpt));
    }
  } else {
    // It is not a Steiner point.
    return 0;
  }

  // Try to reduce the number of edges at 'p' by flips.
  getvertexstar(1, steinerpt, cavetetlist, cavetetvertlist, NULL);
  cavetetlist->restart(); // This list may be re-used.
  if (cavetetvertlist->objects > 3l) {
    valence = reduceedgesatvertex(steinerpt, cavetetvertlist);
  } else {
    valence = cavetetvertlist->objects;
  }
  assert(cavetetlist->objects == 0l);
  cavetetvertlist->restart();

  removeflag = 0;

  if (valence < 3) {
    assert(0); // Unknown cases.
  }

  if (valence == 3) {
    // Only three edges at this vertex. This is only possible when there are
    //   Inverted elements.
    getvertexstar(1, steinerpt, cavetetlist, NULL, NULL);
    if (cavetetlist->objects == 2) {
      printf("to be continued...\n");
      assert(0);
    } else {
      assert(0); // Unknown cases.
    }
    cavetetlist->restart();
    loc = OUTSIDE;
    removeflag = 1;
  } else if (valence == 4) {
    // Only 4 vertices (4 tets) left! 'p' is inside the convex hull of the 4
    //   vertices. This case is due to that 'p' is not exactly on the segment.
    point2tetorg(steinerpt, searchtet);
    loc = INTETRAHEDRON;
    removeflag = 1;
  } else if (valence == 5) {
    // There are 5 edges.
    if (vt == FREESEGVERTEX) {
      sstpivot1(leftseg, searchtet);
      if (org(searchtet) != steinerpt) {
        esymself(searchtet);
      }
      assert(org(searchtet) == steinerpt);
      assert(dest(searchtet) == lpt);
      i = 0; // Count the numbe of tet at the edge [p,lpt].
      neightet.tet = NULL; // Init the face.
      spintet = searchtet;
      while (1) {
        i++;
        if (apex(spintet) == rpt) {
          // Remember the face containing the edge [lpt, rpt].
          neightet = spintet;
        }
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) break;
      }
      if (i == 3) {
        // This case has been checked below.
      } else if (i == 4) {
        // There are 4 tets sharing at [p,lpt]. There must be 4 tets sharing
        //   at [p,rpt].  There must be a face [p, lpt, rpt].  
        if (apex(neightet) == rpt) {
          // The edge (segment) has been already recovered!  At first, this is 
          //   due to the same reason as the case 'valence == 4'.  Second, 
          //   there are 4 vertices (including p, lpt, rpt) exactly coplanar. 
          // We can do a 6-to-2 flip to remove p and recover a face 
          //  [lpt, rpt, c] = [a,b,c].
          // Let 'searchtet' be [p,d,a,b]
          esym(neightet, searchtet);
          enextself(searchtet);
          loc = ONFACE;
          removeflag = 1;
        }
      }
    } else if (vt == FREEFACETVERTEX) {
      point2tetorg(steinerpt, searchtet);
      // Get the three faces of 'searchtet' which share at p.
      //    All faces has p as origin.
      wrktets[0] = searchtet;
      wrktets[1] = searchtet;
      esymself(wrktets[1]);
      enextself(wrktets[1]);
      wrktets[2] = searchtet;
      eprevself(wrktets[2]);
      esymself(wrktets[2]);
      // Get the one which has a subface (should be only 1).
      n = -1;
      valence = 0; // Re-use it as a counter.
      for (i = 0; i < 3; i++) {
        tspivot(wrktets[i], checksh);
        if (checksh.sh != NULL) {
          n = i;
          valence++; 
        }
      }
      assert(valence == 1);
      searchtet = wrktets[n];
      esymself(searchtet);
      enextself(searchtet);
      loc = ONFACE;
      removeflag = 1;
    } else {
      // assert(0); DEBUG IT
    }
    //removeflag = 1;
  } else { // valence > 5.
    
  } // if (valence > 5)

  if (!removeflag) {
    if (vt == FREESEGVERTEX) { 
      // Check is it possible to recover the edge [lpt,rpt].
      // The condition to check is:  Whether each tet containing 'leftseg' is
      //   adjacent to a tet containing 'rightseg'.
      sstpivot1(leftseg, searchtet);
      if (org(searchtet) != steinerpt) {
        esymself(searchtet);
      }
      assert(org(searchtet) == steinerpt);
      assert(dest(searchtet) == lpt);
      spintet = searchtet;
      while (1) {
        // Go to the bottom face of this tet.
        eprev(spintet, neightet);
        esymself(neightet);  // [steinerpt, p1, p2, lpt]
        // Get the adjacent tet.
        fsymself(neightet);  // [p1, steinerpt, p2, rpt]
        if (oppo(neightet) != rpt) {
          // Found a non-matching adjacent tet.
          break;
        }
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) {
          // 'searchtet' is [p,d,p1,p2].
          loc = ONEDGE;
          removeflag = 1;
          break;
        }
      }
    } // if (vt == FREESEGVERTEX)
  }

  if (!removeflag) {
    if (vt == FREESEGVERTEX) {
      // Check if the edge [lpt, rpt] exists.
      if (getedge(lpt, rpt, &searchtet)) {
        // We have recovered this edge. Shift the vertex into the volume.
        // We can recover this edge if the subfaces are not recovered yet.
        if (!checksubfaceflag) {
          // Remove the vertex from the surface mesh.
          //   This will re-create the segment [lpt, rpt] and re-triangulate
          //   all the facets at the segment.
          // Detach the subsegments from their surronding tets.
          for (i = 0; i < 2; i++) {
            checkseg = (i == 0) ? leftseg : rightseg;
            sstpivot1(checkseg, neightet);
            spintet = neightet;
            while (1) {
              tssdissolve1(spintet);
              fnextself(spintet);
              if (spintet.tet == neightet.tet) break;
            }
            sstdissolve1(checkseg);
          } // i
          slawson = 1; // Do lawson flip after removal.
          spivot(rightseg, parentsh); // 'rightseg' has p as its origin.
          sremovevertex(steinerpt, &parentsh, &rightseg, slawson);
          // Clear the list for new subfaces.
          caveshbdlist->restart();
          // Insert the new segment.
          assert(org(searchtet) == lpt);
          assert(dest(searchtet) == rpt);
          sstbond1(rightseg, searchtet);
          spintet = searchtet;
          while (1) {
            tsspivot1(spintet, checkseg); // FOR DEBUG ONLY
            assert(checkseg.sh == NULL);  // FOR DEBUG ONLY
            tssbond1(spintet, rightseg);
            fnextself(spintet);
            if (spintet.tet == searchtet.tet) break;
          }
          // The Steiner point has been shifted into the volume.
          setpointtype(steinerpt, FREEVOLVERTEX);          
          st_segref_count--;
          st_volref_count++;
          // Save this Steiner points in (global) list.
          suppsteinerptlist->newindex((void **) &parypt);
          *parypt = steinerpt;
          return 1;
        } // if (!checksubfaceflag)
      } // if (getedge(...))
    } // if (vt == FREESEGVERTEX)
  } // if (!removeflag)

  if (!removeflag) {
    if (b->verbose > 2) {
      printf("      Unable to remove Steiner point %d val(%d).\n",
             pointmark(steinerpt), valence);
    }
    return 0;
  }

  assert(org(searchtet) == steinerpt);

  if (vt == FREESEGVERTEX) {
    // Detach the subsegments from their surronding tets.
    for (i = 0; i < 2; i++) {
      checkseg = (i == 0) ? leftseg : rightseg;
      sstpivot1(checkseg, neightet);
      spintet = neightet;
      while (1) {
        tssdissolve1(spintet);
        fnextself(spintet);
        if (spintet.tet == neightet.tet) break;
      }
      sstdissolve1(checkseg);
    } // i
    if (checksubfaceflag) {
      // Detach the subfaces at the subsegments from their attached tets.
      for (i = 0; i < 2; i++) {
        checkseg = (i == 0) ? leftseg : rightseg;
        spivot(checkseg, parentsh);
        if (parentsh.sh != NULL) {
          spinsh = parentsh;
          while (1) {
            stpivot(spinsh, neightet);
            if (neightet.tet != NULL) {
              tsdissolve(neightet);
            }
            sesymself(spinsh);
            stpivot(spinsh, neightet);
            if (neightet.tet != NULL) {
              tsdissolve(neightet);
            }
            stdissolve(spinsh);
            spivotself(spinsh); // Go to the next subface.
            if (spinsh.sh == parentsh.sh) break;
          }
        }
      } // i
    } // if (checksubfaceflag)
  }

  if (loc == INTETRAHEDRON) {
    // Collect the four tets containing 'p'.
    fliptets = new triface[4];
    fliptets[0] = searchtet; // [p,d,a,b]
    for (i = 0; i < 2; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,b,c], [p,d,c,a]
    }
    eprev(fliptets[0], fliptets[3]);
    fnextself(fliptets[3]); // it is [a,p,b,c]
    eprevself(fliptets[3]);
    esymself(fliptets[3]); // [a,b,c,p].
    // Remove p by a 4-to-1 flip.
    flip41(fliptets, 1, 0, 0);
    //recenttet = fliptets[0];
  } else if (loc == ONFACE) {
    // Let the original two tets be [a,b,c,d] and [b,a,c,e]. And p is in
    //   face [a,b,c].  Let 'searchtet' be the tet [p,d,a,b].
    // Collect the six tets containing 'p'.
    fliptets = new triface[6];
    fliptets[0] = searchtet; // [p,d,a,b]
    for (i = 0; i < 2; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,b,c], [p,d,c,a]
    }
    eprev(fliptets[0], fliptets[3]);
    fnextself(fliptets[3]); // [a,p,b,e]
    esymself(fliptets[3]);  // [p,a,e,b]
    eprevself(fliptets[3]); // [e,p,a,b]
    for (i = 3; i < 5; i++) {
      fnext(fliptets[i], fliptets[i+1]); // [e,p,b,c], [e,p,c,a]
    }
    // Remove p by a 6-to-2 flip, which is a combination of two flips:
    //   a 3-to-2 (deletes the edge [e,p]), and
    //   a 4-to-1 (deletes the vertex p).
    // First do a 3-to-2 flip on [e,p,a,b],[e,p,b,c],[e,p,c,a]. It creates
    //   two new tets: [a,b,c,p] and [b,a,c,e].  The new tet [a,b,c,p] is
    //   degenerate (has zero volume). It will be deleted in the followed
    //   4-to-1 flip.
    flip32(&(fliptets[3]), 1, 0, 0);
    // Second do a 4-to-1 flip on [p,d,a,b],[p,d,b,c],[p,d,c,a],[a,b,c,p].
    //   This creates a new tet [a,b,c,d].
    flip41(fliptets, 1, 0, 0);
    //recenttet = fliptets[0];
  } else if (loc == ONEDGE) {
    // Let the original edge be [e,d] and p is in [e,d]. Assume there are n
    //   tets sharing at edge [e,d] originally.  We number the link vertices
    //   of [e,d]: p_0, p_1, ..., p_n-1. 'searchtet' is [p,d,p_0,p_1].
    // Count the number of tets at edge [e,p] and [p,d] (this is n).
    n = 0;
    spintet = searchtet;
    while (1) {
      n++;
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }
    assert(n >= 3);
    // Collect the 2n tets containing 'p'.
    fliptets = new triface[2 * n];
    fliptets[0] = searchtet; // [p,b,p_0,p_1] 
    for (i = 0; i < (n - 1); i++) {
      fnext(fliptets[i], fliptets[i+1]); // [p,d,p_i,p_i+1].
    }
    eprev(fliptets[0], fliptets[n]);
    fnextself(fliptets[n]); // [p_0,p,p_1,e]
    esymself(fliptets[n]);  // [p,p_0,e,p_1]
    eprevself(fliptets[n]); // [e,p,p_0,p_1]
    for (i = n; i <  (2 * n - 1); i++) {
      fnext(fliptets[i], fliptets[i+1]); // [e,p,p_i,p_i+1].
    }
    // Remove p by a 2n-to-n flip, it is a sequence of n flips:
    // - Do a 2-to-3 flip on 
    //     [p_0,p_1,p,d] and 
    //     [p,p_1,p_0,e].
    //   This produces: 
    //     [e,d,p_0,p_1], 
    //     [e,d,p_1,p] (degenerated), and 
    //     [e,d,p,p_0] (degenerated).
    wrktets[0] = fliptets[0]; // [p,d,p_0,p_1]
    eprevself(wrktets[0]);    // [p_0,p,d,p_1]
    esymself(wrktets[0]);     // [p,p_0,p_1,d]
    enextself(wrktets[0]);    // [p_0,p_1,p,d] [0]
    wrktets[1] = fliptets[n]; // [e,p,p_0,p_1]
    enextself(wrktets[1]);    // [p,p_0,e,p_1]
    esymself(wrktets[1]);     // [p_0,p,p_1,e]
    eprevself(wrktets[1]);    // [p_1,p_0,p,e] [1]
    flip23(wrktets, 1, 0, 0);
    // Save the new tet [e,d,p,p_0] (degenerated).
    fliptets[n] = wrktets[2];
    // Save the new tet [e,d,p_0,p_1].
    fliptets[0] = wrktets[0];
    // - Repeat from i = 1 to n-2: (n - 2) flips
    //   - Do a 3-to-2 flip on 
    //       [p,p_i,d,e], 
    //       [p,p_i,e,p_i+1], and 
    //       [p,p_i,p_i+1,d]. 
    //     This produces: 
    //       [d,e,p_i+1,p_i], and
    //       [e,d,p_i+1,p] (degenerated).
    for (i = 1; i < (n - 1); i++) {
      wrktets[0] = wrktets[1]; // [e,d,p_i,p] (degenerated).
      enextself(wrktets[0]);   // [d,p_i,e,p] (...)
      esymself(wrktets[0]);    // [p_i,d,p,e] (...) 
      eprevself(wrktets[0]);   // [p,p_i,d,e] (degenerated) [0].
      wrktets[1] = fliptets[n+i];  // [e,p,p_i,p_i+1]
      enextself(wrktets[1]);       // [p,p_i,e,p_i+1] [1]
      wrktets[2] = fliptets[i]; // [p,d,p_i,p_i+1]
      eprevself(wrktets[2]);    // [p_i,p,d,p_i+1]
      esymself(wrktets[2]);     // [p,p_i,p_i+1,d] [2]
      flip32(wrktets, 1, 0, 0);
      // Save the new tet [e,d,p_i,p_i+1].         // FOR DEBUG ONLY
      fliptets[i] = wrktets[0]; // [d,e,p_i+1,p_i] // FOR DEBUG ONLY
      esymself(fliptets[i]);    // [e,d,p_i,p_i+1] // FOR DEBUG ONLY
    }
    // - Do a 4-to-1 flip on 
    //     [p,p_0,e,d],     [d,e,p_0,p],
    //     [p,p_0,d,p_n-1], [e,p_n-1,p_0,p], 
    //     [p,p_0,p_n-1,e], [p_0,p_n-1,d,p], and
    //     [e,d,p_n-1,p]. 
    //   This produces 
    //     [e,d,p_n-1,p_0] and 
    //     deletes p.
    wrktets[3] = wrktets[1];  // [e,d,p_n-1,p] (degenerated) [3]
    wrktets[0] = fliptets[n]; // [e,d,p,p_0] (degenerated)
    eprevself(wrktets[0]);    // [p,e,d,p_0] (...)
    esymself(wrktets[0]);     // [e,p,p_0,d] (...)
    enextself(wrktets[0]);    // [p,p_0,e,d] (degenerated) [0]
    wrktets[1] = fliptets[n-1];   // [p,d,p_n-1,p_0]
    esymself(wrktets[1]);         // [d,p,p_0,p_n-1]
    enextself(wrktets[1]);        // [p,p_0,d,p_n-1] [1]
    wrktets[2] = fliptets[2*n-1]; // [e,p,p_n-1,p_0]
    enextself(wrktets[2]);        // [p_p_n-1,e,p_0]
    esymself(wrktets[2]);         // [p_n-1,p,p_0,e]
    enextself(wrktets[2]);        // [p,p_0,p_n-1,e] [2]
    flip41(wrktets, 1, 0, 0);
    // Save the new tet [e,d,p_n-1,p_0]             // FOR DEBUG ONLY
    fliptets[n-1] = wrktets[0];  // [e,d,p_n-1,p_0] // FOR DEBUG ONLY
    //recenttet = fliptets[0];
  } else {
    assert(0); // Unknown location.
  } // if (iloc == ...)

  delete [] fliptets;

  if (vt == FREESEGVERTEX) {
    // Remove the vertex from the surface mesh.
    //   This will re-create the segment [lpt, rpt] and re-triangulate
    //   all the facets at the segment.
    // Only do lawson flip when subfaces are not recovery yet.
    slawson = (checksubfaceflag ? 0 : 1);
    spivot(rightseg, parentsh); // 'rightseg' has p as its origin.
    sremovevertex(steinerpt, &parentsh, &rightseg, slawson);

    // The original segment is returned in 'rightseg'. 
    rightseg.shver = 0;
    assert(sorg(rightseg) == lpt);
    assert(sdest(rightseg) == rpt);

    // Insert the new segment.
    point2tetorg(lpt, searchtet);
    finddirection(&searchtet, rpt);
    assert(dest(searchtet) == rpt);
    sstbond1(rightseg, searchtet);
    spintet = searchtet;
    while (1) {
      tsspivot1(spintet, checkseg); // FOR DEBUG ONLY
      assert(checkseg.sh == NULL);  // FOR DEBUG ONLY
      tssbond1(spintet, rightseg);
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    }

    if (checksubfaceflag) {
      // Insert subfaces at segment [lpt,rpt] into the tetrahedralization.
      spivot(rightseg, parentsh);
      if (parentsh.sh != NULL) {
        spinsh = parentsh;
        while (1) {
          if (sorg(spinsh) != lpt) {
            sesymself(spinsh);
            assert(sorg(spinsh) == lpt);
          }
          assert(sdest(spinsh) == rpt);
          apexpt = sapex(spinsh);
          // Find the adjacent tet of [lpt,rpt,apexpt];
          spintet = searchtet;
          while (1) {
            if (apex(spintet) == apexpt) {
              tsbond(spintet, spinsh);
              sesymself(spinsh); // Get to another side of this face.
              fsym(spintet, neightet);
              tsbond(neightet, spinsh);
              sesymself(spinsh); // Get back to the original side.
              break;
            }
            fnextself(spintet);
            assert(spintet.tet != searchtet.tet);
            //if (spintet.tet == searchtet.tet) break;
          }
          spivotself(spinsh);
          if (spinsh.sh == parentsh.sh) break;
        }
      }
    } // if (checksubfaceflag)

    // Clear the set of new subfaces.
    caveshbdlist->restart();
  } // if (vt == FREESEGVERTEX)

  // The point has been removed.
  setpointtype(steinerpt, UNUSEDVERTEX);
  unuverts++;
  // Update the correspinding counters.
  if (vt == FREESEGVERTEX) {
    st_segref_count--;
  } else if (vt == FREEFACETVERTEX) {
    st_facref_count--;
  } else if (vt == FREEVOLVERTEX) {
    st_volref_count--;
  }
  if (steinerleft > 0) steinerleft++;

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// suppresssteinerpoint()    Suppress a Steiner point.                       //
//                                                                           //
// Remove a Steiner point 'p' from the segment it lies on. It is replaced by //
// a set of volume Steiner points in each sector at the segment.             //
//                                                                           //
// The list of volume Steiner points is returned in 'suppsteinerptlist'.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::suppressssteinerpoint(point steinerpt)
{
  triface searchtet, neightet, spintet, *parytet;
  triface newtet, newface;
  face parentsh, spinsh, *parysh;
  face newsh, neighsh;
  face leftseg, rightseg, checkseg, *splitseg;
  point lpt = NULL, rpt = NULL, newpt, *parypt;
  point pa, pb, pc;
  verttype vt;
  long bak_supp_steiners;
  int slawson;
  int i, j, k;

  vt = pointtype(steinerpt);

  if (vt == FREESEGVERTEX) {
    sdecode(point2sh(steinerpt), leftseg);
    assert(leftseg.sh != NULL);
    leftseg.shver = 0;
    if (sdest(leftseg) == steinerpt) {
      senext(leftseg, rightseg);
      spivotself(rightseg);
      assert(rightseg.sh != NULL);
      rightseg.shver = 0;
      assert(sorg(rightseg) == steinerpt);
    } else {
      assert(sorg(leftseg) == steinerpt);
      rightseg = leftseg;
      senext2(rightseg, leftseg);
      spivotself(leftseg);
      assert(leftseg.sh != NULL);
      leftseg.shver = 0;
      assert(sdest(leftseg) == steinerpt);
    }
    lpt = sorg(leftseg);
    rpt = sdest(rightseg);
    if (b->verbose > 2) {
      printf("      Suppressing point %d from segment (%d, %d).\n",
             pointmark(steinerpt), pointmark(lpt), pointmark(rpt));
    }
  } else if (vt == FREEFACETVERTEX) {
    if (b->verbose > 2) {
      printf("      Suppressing point %d from facet.\n",
             pointmark(steinerpt));
    }
    //point2shorg(steinerpt, parentsh);
    getvertexstar(0, steinerpt, cavetetlist, NULL, caveshlist);
    parysh = (face *) fastlookup(caveshlist, 0);
    parentsh = *parysh;
    //assert(sapex(parentsh) == steinerpt);
    senext2self(parentsh);
    assert(sorg(parentsh) == steinerpt);
    cavetetlist->restart();
    caveshlist->restart();
  } else {
    // Do nothing.
    return 0;
  }

  if (vt == FREESEGVERTEX) {
    // Check if this edge [lpt, rpt] already exists.
    if (getedge(lpt, rpt, &searchtet)) {
      tsspivot1(searchtet, checkseg);  // SELF_CHECK
      assert(checkseg.sh == NULL);
      return 0;
    }
  }

  bak_supp_steiners = suppsteinerptlist->objects;

  if (vt == FREESEGVERTEX) {
    // Get all subfaces at the left segment [lpt, steinerpt].
    spivot(leftseg, parentsh);
    spinsh = parentsh;
    while (1) {
      cavesegshlist->newindex((void **) &parysh);
      *parysh = spinsh;
      spivotself(spinsh);
      if (spinsh.sh == NULL) break;
      if (spinsh.sh == parentsh.sh) break;
    }
    if (cavesegshlist->objects < 2) {
      // It is a single segment. Not handle it yet.
      cavesegshlist->restart();
      return 0;
    }
    // Detach 'leftseg' and 'rightseg' from their adjacent tets.
    //   These two subsegments will be deleted. 
    sstpivot1(leftseg, neightet);
    spintet = neightet;
    while (1) {
      tssdissolve1(spintet);
      fnextself(spintet);
      if (spintet.tet == neightet.tet) break;
    }
    sstpivot1(rightseg, neightet);
    spintet = neightet;
    while (1) {
      tssdissolve1(spintet);
      fnextself(spintet);
      if (spintet.tet == neightet.tet) break;
    }
  } else { // vt == FREEFACETVERTEX
    // A facet Steiner point. There are exactly two sectors.
    for (i = 0; i < 2; i++) {
      cavesegshlist->newindex((void **) &parysh);
      *parysh = parentsh;
      sesymself(parentsh);
    }
  }

  // Loop through all sectors bounded by facets at this segment.
  //   Within each sector, create a new Steiner point 'np', and replace 'p'
  //   by 'np' for all tets in this sector.
  for (i = 0; i < cavesegshlist->objects; i++) {
    parysh = (face *) fastlookup(cavesegshlist, i);
    // 'parysh' is the face [lpt, steinerpt, #].
    stpivot(*parysh, neightet);
    // Get all tets in this sector.
    setpoint2tet(steinerpt, encode(neightet));
    getvertexstar(0, steinerpt, cavetetlist, NULL, caveshlist);
    assert(caveshlist->objects > 0);
    // Create a new vertex 'np'. 
    makepoint(&newpt, FREEVOLVERTEX);
    st_volref_count++;
    // Init 'np' at the same location of 'p'.
    for (j = 0; j < 3; j++) newpt[j] = steinerpt[j];
    // Within the tet, replace 'p' by 'np'.
    for (j = 0; j < cavetetlist->objects; j++) {
      parytet = (triface *) fastlookup(cavetetlist, j);
      setoppo(*parytet, newpt);
    } // j
    // Save the new Steiner point in list.
    suppsteinerptlist->newindex((void **) &parypt);
    *parypt = newpt;
    // Disconnect the set of boundary faces. They're temporarily open faces.
    //   They will be connected to the new tets after 'p' is removed.
    for (j = 0; j < caveshlist->objects; j++) {
      // Get a boundary face.
      parysh = (face *) fastlookup(caveshlist, j);
      stpivot(*parysh, neightet);
      assert(apex(neightet) == newpt);
      // Clear the connection at this face.
      dissolve(neightet);
      tsdissolve(neightet);
    }
    // Clear the working lists.
    cavetetlist->restart();
    caveshlist->restart();
  } // i
  cavesegshlist->restart();

  // Remove p from the segment.
  slawson = 0; // Do not do flip afterword.
  if (vt == FREESEGVERTEX) { 
    spivot(rightseg, parentsh); // 'rightseg' has p as its origin.
    splitseg = &rightseg;
  } else {
    assert(sorg(parentsh) == steinerpt);
    splitseg = NULL;
  }
  sremovevertex(steinerpt, &parentsh, splitseg, slawson);

  if (vt == FREESEGVERTEX) {
    // The original segment is returned in 'rightseg'. 
    rightseg.shver = 0;
    assert(sorg(rightseg) == lpt);
    assert(sdest(rightseg) == rpt);
  }
  // The set of new subfaces are found in 'caveshbdlist'.
  assert(caveshbdlist->objects > 0);


  // For each new subface, create two new tets at each side of it.
  //   Both of the two new tets have its opposite be dummypoint. 
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    sinfect(*parysh); // Mark it for connecting new tets.
    newsh = *parysh;
    pa = sorg(newsh);
    pb = sdest(newsh);
    pc = sapex(newsh);
    maketetrahedron(&newtet);
    maketetrahedron(&neightet);
    setvertices(newtet, pa, pb, pc, dummypoint);
    setvertices(neightet, pb, pa, pc, dummypoint);
    bond(newtet, neightet);
    tsbond(newtet, newsh);
    sesymself(newsh);
    tsbond(neightet, newsh);
  }

  if (vt == FREESEGVERTEX) {
    // Connecting new tets at the recovered segment.
    spivot(rightseg, parentsh);
    assert(parentsh.sh != NULL);
    spinsh = parentsh;
    while (1) {
      assert(sinfected(spinsh));
      if (sorg(spinsh) != lpt) sesymself(spinsh);
      assert(sorg(spinsh) == lpt);
      assert(sdest(spinsh) == rpt);
      // Get the new tet at this subface.
      stpivot(spinsh, newtet);
      assert(oppo(newtet) == dummypoint);
      tssbond1(newtet, rightseg);
      // Go to the other face at this segment.
      esymself(newtet);
      assert(org(newtet) == rpt);
      assert(newtet.tet[newtet.ver & 3] == NULL);
      // Get the adjacent tet at this segment.
      spivot(spinsh, neighsh);
      if (sorg(neighsh) != lpt) sesymself(neighsh);
      sesymself(neighsh);
      stpivot(neighsh, neightet);
      assert(oppo(neightet) == dummypoint);
      tssbond1(neightet, rightseg);
      sstbond1(rightseg, neightet); 
      // Go to the other face at this segment.
      esymself(neightet);
      assert(org(neightet) == lpt);
      assert(neightet.tet[neightet.ver & 3] == NULL);
      // Connect the two tets (at rightseg) together.
      bond(newtet, neightet);
      // Go to the next subface.
      spivotself(spinsh);
      if (spinsh.sh == parentsh.sh) break;
    }
  }

  // Connecting new tets at new subfaces together.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    newsh = *parysh;
    //assert(sinfected(newsh));
    // Each new subface contains two new tets.
    for (k = 0; k < 2; k++) {
      stpivot(newsh, newtet);
      for (j = 0; j < 3; j++) {
        // Check if this side is open.
        esym(newtet, newface);
        if (newface.tet[newface.ver & 3] == NULL) {
          // An open face. Connect it to its adjacent tet.
          sspivot(newsh, checkseg);
          if (checkseg.sh != NULL) {
            // A segment. It must not be the recovered segment.
            assert(checkseg.sh != rightseg.sh);
            tssbond1(newtet, checkseg);
            //sstbond1(checkseg, newtet);
          }
          spivot(newsh, neighsh);
          if (neighsh.sh != NULL) {
            // The adjacent subface exists. It's not a dangling segment.
            if (sorg(neighsh) != sdest(newsh)) sesymself(neighsh);
            stpivot(neighsh, neightet);
            if (sinfected(neighsh)) {
              esymself(neightet);
              assert(neightet.tet[neightet.ver & 3] == NULL);          
            } else {
              // Search for an open face at this edge.
              spintet = neightet;
              while (1) {
                esym(spintet, searchtet);
                fsym(searchtet, spintet);
                if (spintet.tet == NULL) break;
                assert(spintet.tet != neightet.tet);
              }
              // Found an open face at 'searchtet'.
              neightet = searchtet;
            }
          } else {
            // The edge (at 'newsh') is a dangling segment.
            assert(checkseg.sh != NULL);
            // Get an adjacent tet at this segment.
            sstpivot1(checkseg, neightet);
            assert(!isdeadtet(neightet));
            if (org(neightet) != sdest(newsh)) esymself(neightet);
            assert((org(neightet) == sdest(newsh)) &&
                   (dest(neightet) == sorg(newsh)));
            // Search for an open face at this edge.
            spintet = neightet;
            while (1) {
              esym(spintet, searchtet);
              fsym(searchtet, spintet);
              if (spintet.tet == NULL) break;
              assert(spintet.tet != neightet.tet);
            }
            // Found an open face at 'searchtet'.
            neightet = searchtet;
          }
          pc = apex(newface);
          if (pc == dummypoint) {
            setapex(newface, apex(neightet));
          } else {
            assert(pc == apex(neightet));
          }
          bond(newface, neightet);
        } // if (newface.tet[newface.ver & 3] == NULL)
        enextself(newtet);
        senextself(newsh);
      } // j
      sesymself(newsh);
    } // k
  } // i

  // Unmark all new subfaces.
  for (i = 0; i < caveshbdlist->objects; i++) {
    parysh = (face *) fastlookup(caveshbdlist, i);
    suninfect(*parysh);
  }
  caveshbdlist->restart();


  setpointtype(steinerpt, UNUSEDVERTEX);
  unuverts++;
  if (vt == FREESEGVERTEX) {
    st_segref_count--;
  } else { // vt == FREEFACETVERTEX
    st_facref_count--;
  }
  if (steinerleft > 0) steinerleft++;

  if (b->verbose > 2) {
    printf("      Duplicated %ld Steiner points.\n", 
           suppsteinerptlist->objects - bak_supp_steiners);
  }

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// suppresssteinerpoints()    Suppress Steiner points.                       //
//                                                                           //
// Each Steiner point is either removed or shifted into the interior.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::suppresssteinerpoints()
{
  triface *parytet;
  point rempt, *parypt, *plastpt, *ppt;
  optparameters opm;
  REAL ori;
  int bak_fliplinklevel;
  int remcount, smtcount;
  int count, nt;
  int i, j;

  if (!b->quiet) {
    printf("Suppressing Steiner points ...\n");
  }

  bak_fliplinklevel = b->fliplinklevel;
  b->fliplinklevel = 100000; // Unlimited flip level.
  remcount = 0;

  if (b->nobisect_param > 1) { // -Y2
    // Try to remove all the Steiner points.
    for (i = 0; i < subvertstack->objects; i++) {
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (pointtype(rempt) != UNUSEDVERTEX) {
        if (removevertexbyflips(rempt)) {
          remcount++;
        }
      }
    }
    if (b->verbose) {
      if (remcount > 0) {
        printf("  Removed %d Steiner points.\n", remcount);
      }
    }
    subvertstack->restart();
  }

  remcount = smtcount = 0;

  // Try to remove the suppressed Steiner points.
  for (i = 0; i < suppsteinerptlist->objects; i++) {
    // Get the Steiner point.
    parypt = (point *) fastlookup(suppsteinerptlist, i);
    rempt = *parypt;
    if (pointtype(rempt) != UNUSEDVERTEX) {
      assert((pointtype(rempt) == FREESEGVERTEX) ||
             (pointtype(rempt) == FREEFACETVERTEX) ||
             (pointtype(rempt) == FREEVOLVERTEX));
      if (removevertexbyflips(rempt)) {
        // Move the last entry to fill the current one.
        j = (int) (suppsteinerptlist->objects - 1);
        plastpt = (point *) fastlookup(suppsteinerptlist, j);
        *parypt = *plastpt;
        suppsteinerptlist->objects--;
        i--;
        remcount++;
      }
    } else {
      // The point has been removed.
      // Move the last entry to fill the current one.
      j = (int) (suppsteinerptlist->objects - 1);
      plastpt = (point *) fastlookup(suppsteinerptlist, j);
      *parypt = *plastpt;
      suppsteinerptlist->objects--;
      i--;
    }
  } // i

  if (b->verbose) {
    if (remcount > 0) {
      printf("  Removed %d suppressed Steiner points.\n", remcount);
    }
  }

  if (suppsteinerptlist->objects == 0l) {
    b->fliplinklevel = bak_fliplinklevel;
    return remcount;
  }

  // Point smooth options.
  opm.max_min_volume = 1;
  opm.numofsearchdirs = 20;
  opm.searchstep = 0.001;

  nt = 0;

  while (1) {
    // Try to smooth volume Steiner points.
    count = 0;

    for (i = 0; i < suppsteinerptlist->objects; i++) {
      parypt = (point *) fastlookup(suppsteinerptlist, i);
      rempt = *parypt;
      if (pointtype(rempt) == FREEVOLVERTEX) {
        getvertexstar(1, rempt, cavetetlist, NULL, NULL);
        // Calculate the initial smallest volume (maybe zero or negative).
        for (j = 0; j < cavetetlist->objects; j++) {
          parytet = (triface *) fastlookup(cavetetlist, j);
          ppt = (point *) &(parytet->tet[4]);
          ori = orient3d(ppt[1], ppt[0], ppt[2], ppt[3]);
          if (j == 0) {
            opm.initval = ori;
          } else {
            if (opm.initval > ori) opm.initval = ori; 
          }
        }
        if (smoothpoint(rempt, cavetetlist, 1, &opm)) {
          count++;
        }
        cavetetlist->restart();
      }
    } // i

    smtcount += count;

    if (count == 0) {
      // No point has been smoothed.
      break;
    }

    nt++;
    if (nt > 2) {
      break; // Already three iterations.
    }
  } // while

  // The mesh should not contain inverted (or degenrrated) tets now.
  checkinverttetflag = 0;

  if (b->verbose) {
    if (smtcount > 0) {
      printf("  Smoothed %d Steiner points.\n", smtcount); 
    }
  }

  b->fliplinklevel = bak_fliplinklevel;

  return smtcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoverboundary()    Recover segments and facets.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::recoverboundary(clock_t& tv)
{
  arraypool *misseglist, *misshlist;
  arraypool *bdrysteinerptlist;
  face searchsh, *parysh;
  face searchseg, *paryseg;
  point rempt, *parypt;
  long ms; // The number of missing segments/subfaces.
  int nit; // The number of iterations.
  int s, i;

  // Counters.
  long bak_segref_count, bak_facref_count, bak_volref_count;
  long bak_supp_count;

  if (!b->quiet) {
    printf("Recovering boundaries...\n");
  }

  if (b->verbose) {
    printf("  Flip link level = %d\n", b->fliplinklevel);
  }


  if (b->verbose) {
    printf("  Recovering segments.\n");
  }

  // Segments will be introduced.
  checksubsegflag = 1;

  misseglist = new arraypool(sizeof(face), 8);
  bdrysteinerptlist = new arraypool(sizeof(point), 8);

    // In random order.
    subsegs->traversalinit();
    for (i = 0; i < subsegs->items; i++) {
      s = randomnation(i + 1);
      // Move the s-th seg to the i-th.
      subsegstack->newindex((void **) &paryseg);
      *paryseg = * (face *) fastlookup(subsegstack, s);
      // Put i-th seg to be the s-th.
      searchseg.sh = shellfacetraverse(subsegs);
      //sinfect(searchseg);  // Only save it once.
      paryseg = (face *) fastlookup(subsegstack, s);
      *paryseg = searchseg;
    }

  // The init number of missing segments.
  ms = subsegs->items;
  nit = 0; 
  if (b->fliplinklevel < 0) {
    autofliplinklevel = 1; // Init value.
  }

  while (1) {
    recoversegments(misseglist, 0, 0);

    if (misseglist->objects > 0) {
      if (b->fliplinklevel >= 0) {
        break;
      } else {
        if (misseglist->objects >= ms) {
          nit++;
          if (nit >= 3) {
            //break;
            // Do the last round with unbounded flip link level.
            b->fliplinklevel = 100000;
          }
        } else {
          ms = misseglist->objects;
          if (nit > 0) {
            nit--;
          }
        }
        for (i = 0; i < misseglist->objects; i++) {
          subsegstack->newindex((void **) &paryseg);
          *paryseg = * (face *) fastlookup(misseglist, i);
        }
        misseglist->restart();
        autofliplinklevel+=b->fliplinklevelinc;
      }
    } else {
      // All segments are recovered.
      break;
    }
  } // while (1)

  if (b->verbose) {
    printf("  %ld (%ld) segments are recovered (missing).\n", 
           subsegs->items - misseglist->objects, misseglist->objects);
  }

  if (misseglist->objects > 0) {
    // There are missing segments. Increase the fliplevel.
    nit = 0;
    while (misseglist->objects > 0) {
      ms = misseglist->objects;
      for (i = 0; i < misseglist->objects; i++) {
        subsegstack->newindex((void **) &paryseg);
        *paryseg = * (face *) fastlookup(misseglist, i);
      }
      misseglist->restart();

      // Recover the missing segments by doing more flips.
      recoversegments(misseglist, 1, 0);

      if (misseglist->objects < ms) {
        // The number of missing segments is reduced.
        continue;
      } else {
        nit++;
        if (nit >= 3) {
          break;
        }
      }
    }
    if (b->verbose) {
      printf("  %ld (%ld) segments are recovered (missing).\n", 
             subsegs->items - misseglist->objects, misseglist->objects);
    }
  }

  if (misseglist->objects > 0) {
    // There are missing segments. Add Steiner points in volume.
    nit = 0;
    while (misseglist->objects > 0) {
      ms = misseglist->objects;
      for (i = 0; i < misseglist->objects; i++) {
        subsegstack->newindex((void **) &paryseg);
        *paryseg = * (face *) fastlookup(misseglist, i);
      }
      misseglist->restart();

      // Recover the missing segments (with Steiner points).
      recoversegments(misseglist, 1, 1);

      if (misseglist->objects < ms) {
        // The number of missing segments is reduced.
        continue;
      } else {
        nit++;
        if (nit >= 3) {
          break;
        }
      }
    }
    if (b->verbose) {
      printf("  Added %ld Steiner points in volume.\n", st_volref_count);
    }
  }

  if (misseglist->objects > 0) {
    // There are missing segments. Add Steiner points to split them.
    long bak_inpoly_count = st_volref_count; //st_inpoly_count;
    for (i = 0; i < misseglist->objects; i++) {
      subsegstack->newindex((void **) &paryseg);
      *paryseg = * (face *) fastlookup(misseglist, i);
    }
    misseglist->restart();

    // Recover the missing segments (with Steiner points).
    recoversegments(misseglist, 1, 2);

    if (b->verbose) {
      printf("  Added %ld Steiner points in segments.\n", st_segref_count);
      if (st_volref_count > bak_inpoly_count) {
        printf("  Added another %ld Steiner points in volume.\n", 
               st_volref_count - bak_inpoly_count);
      }
    }
    assert(misseglist->objects == 0l);
  }


  if (st_segref_count > 0) {
    // Try to remove the Steiner points added in segments.
    bak_segref_count = st_segref_count;
    bak_volref_count = st_volref_count;
    for (i = 0; i < subvertstack->objects; i++) {
      // Get the Steiner point.
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (!removevertexbyflips(rempt)) {
        // Save it in list.
        bdrysteinerptlist->newindex((void **) &parypt);
        *parypt = rempt;
      }
    }
    if (b->verbose) {
      if (st_segref_count < bak_segref_count) {
        if (bak_volref_count < st_volref_count) {
          printf("  Suppressed %ld Steiner points in segments.\n", 
                 st_volref_count - bak_volref_count);
        }
        if ((st_segref_count + (st_volref_count - bak_volref_count)) <
            bak_segref_count) {
          printf("  Removed %ld Steiner points in segments.\n", 
                 bak_segref_count - 
                   (st_segref_count + (st_volref_count - bak_volref_count)));
        }
      }
    }
    subvertstack->restart();
  }


  tv = clock();

  if (b->verbose) {
    printf("  Recovering facets.\n");
  }

  // Subfaces will be introduced.
  checksubfaceflag = 1;

  misshlist = new arraypool(sizeof(face), 8);

  // Randomly order the subfaces.
  subfaces->traversalinit();
  for (i = 0; i < subfaces->items; i++) {
    s = randomnation(i + 1);
    // Move the s-th subface to the i-th.
    subfacstack->newindex((void **) &parysh);
    *parysh = * (face *) fastlookup(subfacstack, s);
    // Put i-th subface to be the s-th.
    searchsh.sh = shellfacetraverse(subfaces);
    parysh = (face *) fastlookup(subfacstack, s);
    *parysh = searchsh;
  }

  ms = subfaces->items;
  nit = 0; 
  b->fliplinklevel = -1; // Init.
  if (b->fliplinklevel < 0) {
    autofliplinklevel = 1; // Init value.
  }

  while (1) {
    recoversubfaces(misshlist, 0);
    if (misshlist->objects > 0) {
      if (b->fliplinklevel >= 0) {
        break;
      } else {
        if (misshlist->objects >= ms) {
          nit++;
          if (nit >= 3) {
            //break;
            // Do the last round with unbounded flip link level.
            b->fliplinklevel = 100000;
          }
        } else {
          ms = misshlist->objects;
          if (nit > 0) {
            nit--;
          }
        }
        for (i = 0; i < misshlist->objects; i++) {
          subfacstack->newindex((void **) &parysh);
          *parysh = * (face *) fastlookup(misshlist, i);
        }
        misshlist->restart();
        autofliplinklevel+=b->fliplinklevelinc;
      }
    } else {
      // All subfaces are recovered.
      break;
    }
  } // while (1)

  if (b->verbose) {
    printf("  %ld (%ld) subfaces are recovered (missing).\n", 
           subfaces->items - misshlist->objects, misshlist->objects);
  }

  if (misshlist->objects > 0) {
    // There are missing subfaces. Add Steiner points.
    for (i = 0; i < misshlist->objects; i++) {
      subfacstack->newindex((void **) &parysh);
      *parysh = * (face *) fastlookup(misshlist, i);
    }
    misshlist->restart();

    recoversubfaces(NULL, 1);

    if (b->verbose) {
      printf("  Added %ld Steiner points in facets.\n", st_facref_count);
    }
  }


  if (st_facref_count > 0) {
    // Try to remove the Steiner points added in facets.
    bak_facref_count = st_facref_count;
    for (i = 0; i < subvertstack->objects; i++) {
      // Get the Steiner point.
      parypt = (point *) fastlookup(subvertstack, i);
      rempt = *parypt;
      if (!removevertexbyflips(*parypt)) {
        // Save it in list.
        bdrysteinerptlist->newindex((void **) &parypt);
        *parypt = rempt;
      }
    }
    if (b->verbose) {
      if (st_facref_count < bak_facref_count) {
        printf("  Removed %ld Steiner points in facets.\n", 
               bak_facref_count - st_facref_count);
      }
    }
    subvertstack->restart();
  }


  if ((bdrysteinerptlist->objects > 0) && (b->nobisect_param > 0)) { // -Y1
    bak_supp_count = 0;
    b->fliplinklevel = 100000; // Unlimited flip levels.
    do {
      // Suppress boundary Steiner points.
      for (i = 0; i < bdrysteinerptlist->objects; i++) {
        parypt = (point *) fastlookup(bdrysteinerptlist, i);
        rempt = *parypt;
        suppressssteinerpoint(rempt);
        bak_supp_count++;
      }
      bdrysteinerptlist->restart();
      // There may be subfaces need to be recover.
      if (subfacstack->objects > 0l) {
        recoversubfaces(NULL, 1);
      }
    } while (bdrysteinerptlist->objects > 0);
    if (b->verbose) {
      printf("  Suppressed %ld Steiner points from boundary.\n", 
             bak_supp_count);
    }
    // The mesh contains inverted (or degenrrated) tets now.
    checkinverttetflag = 1;
  } // if


  delete bdrysteinerptlist;
  delete misseglist;
  delete misshlist;
}

////                                                                       ////
////                                                                       ////
//// steiner_cxx //////////////////////////////////////////////////////////////


