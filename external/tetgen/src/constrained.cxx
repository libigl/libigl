#include "../tetgen.h"
//// constrained_cxx //////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// markacutevertices()    Classify vertices as ACUTEVERTEXs or RIDGEVERTEXs. //
//                                                                           //
// Initially all segment vertices have type RIDGEVERTEX.  A segment is acute //
// if there are at least two segments incident at it form an angle less than //
// theta (= 60 degree).                                                      //
//                                                                           //
// The minimum segment-segment angle (minfaceang) is calculated.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::markacutevertices()
{
  face* segperverlist;
  int* idx2seglist;
  point pa, pb, pc;
  REAL anglimit, ang;
  bool acuteflag;
  int acutecount;
  int idx, i, j;

  REAL sharpanglimit;
  int sharpsegcount;

  if (b->verbose) {
    printf("  Marking acute vertices.\n");
  }
  anglimit = PI / 3.0;  // 60 degree.
  sharpanglimit = 5.0 / 180.0 * PI; // 5 degree. 
  minfaceang = PI; // 180 degree.
  acutecount = sharpsegcount = 0;

  // Construct a map from points to segments.
  makepoint2submap(subsegs, idx2seglist, segperverlist);

  // Loop over the set of vertices.
  points->traversalinit();
  pa = pointtraverse();
  while (pa != NULL) {
    idx = pointmark(pa) - in->firstnumber;
    // Mark it if it is an endpoint of some segments.
    if (idx2seglist[idx + 1] > idx2seglist[idx]) {
      if (b->psc) {
        // Only test it if it is an input vertex.
        if (pointtype(pa) == FREESEGVERTEX) {
          pa = pointtraverse();
          continue;
        }
      }
      acuteflag = false;
      // Do a brute-force pair-pair check.
      for (i=idx2seglist[idx]; i<idx2seglist[idx + 1]; i++) {
        pb = sdest(segperverlist[i]);
        for (j = i + 1; j < idx2seglist[idx + 1]; j++) {
          pc = sdest(segperverlist[j]);
          ang = interiorangle(pa, pb, pc, NULL); 
          if (!acuteflag) {
            acuteflag = ang < anglimit;
          }
          // Remember the smallest angle.
          if (ang < minfaceang) minfaceang = ang;
          // Mark segments at extremely small angle.
          if (ang < sharpanglimit) {
            if (shelltype(segperverlist[i]) != SHARP) {
              setshelltype(segperverlist[i], SHARP);
              sharpsegcount++;
            }
            if (shelltype(segperverlist[j]) != SHARP) {
              setshelltype(segperverlist[j], SHARP);
              sharpsegcount++;
            }
          }
        } // j
      } // i
      if (!acuteflag) {
        if ((idx2seglist[idx + 1] - idx2seglist[idx]) > 4) {
          // There are at least 5 segments shared at this vertices.
          acuteflag = true;
        }
      }
      if (acuteflag) {
        if (b->verbose > 2) {
          printf("      Mark %d as ACUTEVERTEX.\n", pointmark(pa));
        }
        setpointtype(pa, ACUTEVERTEX);
        acutecount++;
      }
    }
    pa = pointtraverse();
  }

  if (b->verbose) {
    if (acutecount > 0) {
      printf("  Found %d acute vertices.\n", acutecount);
    }
    if (sharpsegcount > 0) {
      printf("  Found %d sharp segments.\n", sharpsegcount);
    }
    printf("  Minimum seg-seg angle = %g.\n", minfaceang / PI * 180.0);
  }

  delete [] idx2seglist;
  delete [] segperverlist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// reportselfintersect()    Report a self-intersection.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::reportselfintersect(face *checkseg, face *checksh)
{
  face parentsh;
  point pa, pb, pc, pd, pe;
  point fa, fb;

  pa = sorg(*checkseg);
  pb = sdest(*checkseg);
  fa = farsorg(*checkseg);
  fb = farsdest(*checkseg);

  pc = sorg(*checksh);
  pd = sdest(*checksh);
  pe = sapex(*checksh);

  printf("  !! Detected a self-intersection between:\n");
  printf("     A segment [%d,%d] < [%d,%d], \n", pointmark(pa), pointmark(pb),
         pointmark(fa), pointmark(fb));
  printf("     a subface [%d,%d,%d] in facet #%d.\n", pointmark(pc), 
         pointmark(pd), pointmark(pe), shellmark(*checksh));

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// finddirection()    Find the tet on the path from one point to another.    //
//                                                                           //
// The path starts from 'searchtet''s origin and ends at 'endpt'. On finish, //
// 'searchtet' contains a tet on the path, its origin does not change.       //
//                                                                           //
// The return value indicates one of the following cases (let 'searchtet' be //
// abcd, a is the origin of the path):                                       //
//   - ACROSSVERT, edge ab is collinear with the path;                       //
//   - ACROSSEDGE, edge bc intersects with the path;                         //
//   - ACROSSFACE, face bcd intersects with the path.                        //
//                                                                           //
// WARNING: This routine is designed for convex triangulations, and will not //
// generally work after the holes and concavities have been carved.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::interresult 
  tetgenmesh::finddirection(triface* searchtet, point endpt)
{
  triface neightet;
  point pa, pb, pc, pd;
  enum {HMOVE, RMOVE, LMOVE} nextmove;
  REAL hori, rori, lori;
  int s;

  // The origin is fixed.
  pa = org(*searchtet);
  if ((point) searchtet->tet[7] == dummypoint) {
    // A hull tet. Choose the neighbor of its base face.
    searchtet->ver = 11;
    fsymself(*searchtet);
    // Reset the origin to be pa.
    if ((point) searchtet->tet[4] == pa) {
      searchtet->ver = 11;
    } else if ((point) searchtet->tet[5] == pa) {
      searchtet->ver = 3;
    } else if ((point) searchtet->tet[6] == pa) {
      searchtet->ver = 7;
    } else {
      assert((point) searchtet->tet[7] == pa); 
      searchtet->ver = 0;
    }
  }

  pb = dest(*searchtet);
  // Check whether the destination or apex is 'endpt'.
  if (pb == endpt) {
    // pa->pb is the search edge.
    return ACROSSVERT;
  }

  pc = apex(*searchtet);
  if (pc == endpt) {
    // pa->pc is the search edge.
    eprevself(*searchtet);
    esymself(*searchtet);
    return ACROSSVERT;
  }

  // Walk through tets around pa until the right one is found.
  while (1) {

    pd = oppo(*searchtet);

    if (b->verbose > 3) {
      printf("        From tet (%d, %d, %d, %d) to %d.\n", pointmark(pa),
        pointmark(pb), pointmark(pc), pointmark(pd), pointmark(endpt));
    }

    // Check whether the opposite vertex is 'endpt'.
    if (pd == endpt) {
      // pa->pd is the search edge.
      esymself(*searchtet);
      enextself(*searchtet);
      return ACROSSVERT;
    }
    // Check if we have entered outside of the domain.
    if (pd == dummypoint) {
      // This is possible when the mesh is non-convex.
      assert(nonconvex);
      return ACROSSSUB; // Hit a bounday.
    }

    // Now assume that the base face abc coincides with the horizon plane,
    //   and d lies above the horizon.  The search point 'endpt' may lie
    //   above or below the horizon.  We test the orientations of 'endpt'
    //   with respect to three planes: abc (horizon), bad (right plane),
    //   and acd (left plane). 
    hori = orient3d(pa, pb, pc, endpt);
    rori = orient3d(pb, pa, pd, endpt);
    lori = orient3d(pa, pc, pd, endpt);

    // Now decide the tet to move.  It is possible there are more than one
    //   tets are viable moves. Is so, randomly choose one. 
    if (hori > 0) {
      if (rori > 0) {
        if (lori > 0) {
          // Any of the three neighbors is a viable move.
          s = randomnation(3); 
          if (s == 0) {
            nextmove = HMOVE;
          } else if (s == 1) {
            nextmove = RMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // Two tets, below horizon and below right, are viable.
          s = randomnation(2); 
          if (s == 0) {
            nextmove = HMOVE;
          } else {
            nextmove = RMOVE;
          }
        }
      } else {
        if (lori > 0) {
          // Two tets, below horizon and below left, are viable.
          s = randomnation(2); 
          if (s == 0) {
            nextmove = HMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // The tet below horizon is chosen.
          nextmove = HMOVE;
        }
      }
    } else {
      if (rori > 0) {
        if (lori > 0) {
          // Two tets, below right and below left, are viable.
          s = randomnation(2); 
          if (s == 0) {
            nextmove = RMOVE;
          } else {
            nextmove = LMOVE;
          }
        } else {
          // The tet below right is chosen.
          nextmove = RMOVE;
        }
      } else {
        if (lori > 0) {
          // The tet below left is chosen.
          nextmove = LMOVE;
        } else {
          // 'endpt' lies either on the plane(s) or across face bcd.
          if (hori == 0) {
            if (rori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pb.
              return ACROSSVERT;
            }
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pc.
              eprevself(*searchtet);
              esymself(*searchtet); // [a,c,d]
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pc.
            return ACROSSEDGE;
          }
          if (rori == 0) {
            if (lori == 0) {
              // pa->'endpt' is COLLINEAR with pa->pd.
              esymself(*searchtet); // face bad.
              enextself(*searchtet); // face [a,d,b]
              return ACROSSVERT;
            }
            // pa->'endpt' crosses the edge pb->pd.
            esymself(*searchtet); // face bad.
            enextself(*searchtet); // face adb
            return ACROSSEDGE;
          }
          if (lori == 0) {
            // pa->'endpt' crosses the edge pc->pd.
            eprevself(*searchtet);
            esymself(*searchtet); // face acd
            return ACROSSEDGE;
          }
          // pa->'endpt' crosses the face bcd.
          return ACROSSFACE;
        }
      }
    }

    // Move to the next tet, fix pa as its origin.
    if (nextmove == RMOVE) {
      fnextself(*searchtet);
    } else if (nextmove == LMOVE) {
      eprevself(*searchtet);
      fnextself(*searchtet);
      enextself(*searchtet);
    } else { // HMOVE
      fsymself(*searchtet);
      enextself(*searchtet);
    }
    assert(org(*searchtet) == pa); 
    pb = dest(*searchtet);
    pc = apex(*searchtet);

  } // while (1)

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsegment()    Look for a given segment in the tetrahedralization T.   //
//                                                                           //
// Search an edge in the tetrahedralization that matches the given segmment. //
// If such an edge exists, the segment is 'locked' at the edge. 'searchtet'  //
// returns this (constrained) edge. Otherwise, the segment is missing.       //
//                                                                           //
// The returned value indicates one of the following cases:                  //
//   - SHAREEDGE, the segment exists and is inserted in T;                   //
//   - ACROSSEDGE, the segment intersects an edge (in 'searchtet').          //
//   - ACROSSFACE, the segment crosses a face (in 'searchtet').              //
//                                                                           //
// The following cases can happen when the input PLC is not valid.           //
//   - ACROSSVERT, the segment intersects a vertex ('refpt').                //
//   - ACROSSSEG, the segment intersects a segment(returned by 'searchtet'). //
//   - ACROSSSUB, the segment intersects a subface(returned by 'searchtet'). //
//                                                                           //
// If the returned value is ACROSSEDGE or ACROSSFACE, i.e., the segment is   //
// missing, 'refpt' returns the reference point for splitting thus segment,  //
// 'searchtet' returns a tet containing the 'refpt'.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::interresult 
  tetgenmesh::scoutsegment(point startpt, point endpt, triface* searchtet, 
                           point* refpt, arraypool* intfacelist)
{
  triface neightet, reftet;
  face checkseg, checksh;
  point pa, pb, pc, pd;
  badface *bface;
  enum interresult dir;
  REAL angmax, ang;
  long facecount;
  int types[2], poss[4];
  int pos, i, j;

  if (b->verbose > 2) {
    printf("      Scout seg (%d, %d).\n", pointmark(startpt), pointmark(endpt));
  }

  point2tetorg(startpt, *searchtet);
  dir = finddirection(searchtet, endpt);

  if (dir == ACROSSVERT) {
    pd = dest(*searchtet);
    if (pd == endpt) {
      // The job is done. 
      return SHAREEDGE;
    } else {
      // A point is on the path.
      *refpt = pd;
      return ACROSSVERT;
    }
  } // if (dir == ACROSSVERT)

  if (b->verbose > 2) {
    printf("      Seg is missing.\n");
  }
  // dir is either ACROSSEDGE or ACROSSFACE.

  enextesymself(*searchtet); // Go to the opposite face.
  fsymself(*searchtet); // Enter the adjacent tet.

  if (dir == ACROSSEDGE) {
    // Check whether two segments are intersecting.
    tsspivot1(*searchtet, checkseg);
    if (checkseg.sh != NULL) {
      return ACROSSSEG;
    }
    across_edge_count++;
  } else if (dir == ACROSSFACE) {
    if (checksubfaceflag) {
      // Check whether a segment and a subface are intersecting.
      tspivot(*searchtet, checksh);
      if (checksh.sh != NULL) {
        return ACROSSSUB;
      }
    }
  }

  if (refpt == NULL) {
    return dir;
  }

  if (b->verbose > 2) {
    printf("      Scout a ref-point for it.\n");
  }
  facecount = across_face_count;

  pa = org(*searchtet);
  angmax = interiorangle(pa, startpt, endpt, NULL);
  *refpt = pa;
  pb = dest(*searchtet);
  ang = interiorangle(pb, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pb;
  }
  pc = apex(*searchtet);
  ang = interiorangle(pc, startpt, endpt, NULL);
  if (ang > angmax) {
    angmax = ang;
    *refpt = pc;
  }
  reftet = *searchtet; // Save the tet containing the refpt.

  // Search intersecting faces along the segment.
  while (1) {

    if (intfacelist != NULL) {
      if (dir == ACROSSFACE) { 
        // Save the intersecting face.
        intfacelist->newindex((void **) &bface);
        bface->tt = *searchtet;
        bface->forg = org(*searchtet);
        bface->fdest = dest(*searchtet);
        bface->fapex = apex(*searchtet);
        // Save the intersection type (ACROSSFACE or ACROSSEDGE).
        bface->key = (REAL) dir;
      } else { // dir == ACROSSEDGE
        i = 0;      
        if (intfacelist->objects > 0l) {
          // Get the last saved one.
          bface = (badface *) fastlookup(intfacelist, intfacelist->objects - 1);
          if (((enum interresult) (int) bface->key) == ACROSSEDGE) {
            // Skip this edge if it is the same as the last saved one.
            if (((bface->forg == org(*searchtet)) &&
                 (bface->fdest == dest(*searchtet))) ||
                ((bface->forg == dest(*searchtet)) &&
                 (bface->fdest == org(*searchtet)))) {
              i = 1; // Skip this edge.
            }
          }
        }
        if (i == 0) {
          // Save this crossing edge.
          intfacelist->newindex((void **) &bface);
          bface->tt = *searchtet;
          bface->forg = org(*searchtet);
          bface->fdest = dest(*searchtet);
          // bface->fapex = apex(*searchtet);
          // Save the intersection type (ACROSSFACE or ACROSSEDGE).
          bface->key = (REAL) dir;
        }
      }
    }

    pd = oppo(*searchtet);
    assert(pd != dummypoint);  // SELF_CHECK

    if (b->verbose > 3) {
      printf("        Passing face (%d, %d, %d, %d), dir(%d).\n", 
             pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd), 
             (int) dir);
    }
    across_face_count++;

    // Stop if we meet 'endpt'.
    if (pd == endpt) break;

    ang = interiorangle(pd, startpt, endpt, NULL);
    if (ang > angmax) {
      angmax = ang;
      *refpt = pd;
      reftet = *searchtet;
    }

    // Find a face intersecting the segment.
    if (dir == ACROSSFACE) {
      // One of the three oppo faces in 'searchtet' intersects the segment.
      neightet = *searchtet;
      j = (neightet.ver & 3); // j is the current face number.
      for (i = j + 1; i < j + 4; i++) {
        neightet.ver = (i % 4);
        pa = org(neightet);
        pb = dest(neightet);
        pc = apex(neightet);
        pd = oppo(neightet); // The above point.
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum interresult) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
      assert(dir != DISJOINT);  // SELF_CHECK
    } else { // dir == ACROSSEDGE
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
        if (tri_edge_test(pa, pb, pc, startpt, endpt, pd, 1, types, poss)) {
          dir = (enum interresult) types[0];
          pos = poss[0];
          break;
        } else {
          dir = DISJOINT;
          pos = 0;
        }
      }
      if (dir == DISJOINT) {
        // No intersection. Rotate to the next tet at the edge.
        dir = ACROSSEDGE;
        fnextself(*searchtet);
        continue;
      }
    }

    if (dir == ACROSSVERT) {
      // This segment passing a vertex. Choose it and return.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
      pd = org(neightet);
      if (b->verbose > 2) {
        angmax = interiorangle(pd, startpt, endpt, NULL);
      }
      *refpt = pd;
      // break;
      return ACROSSVERT;
    } else if (dir == ACROSSEDGE) {
      // Get the edge intersects with the segment.
      for (i = 0; i < pos; i++) {
        enextself(neightet);
      }
    }
    // Go to the next tet.
    fsym(neightet, *searchtet);

    if (dir == ACROSSEDGE) {
      // Check whether two segments are intersecting.
      tsspivot1(*searchtet, checkseg);
      if (checkseg.sh != NULL) {
        return ACROSSSEG;
      }
      across_edge_count++;
    } else if (dir == ACROSSFACE) {
      if (checksubfaceflag) {
        // Check whether a segment and a subface are intersecting.
        tspivot(*searchtet, checksh);
        if (checksh.sh != NULL) {
          return ACROSSSUB;
        }
      }
    }

  } // while (1)

  // A valid reference point should inside the diametrial circumsphere of
  //   the missing segment, i.e., it encroaches upon it.
  if (2.0 * angmax < PI) {
    *refpt = NULL;
  }

  // dir is either ACROSSVERT, or ACROSSEDGE, or ACROSSFACE.
  if (b->verbose > 2) {
    if (*refpt != NULL) {
      printf("      Refpt %d (%g), visited %ld faces.\n", pointmark(*refpt),
             angmax / PI * 180.0, across_face_count - facecount);
    } else {
      printf("      No refpt (%g) is found, visited %ld faces.\n", 
             angmax / PI * 180.0, across_face_count - facecount);
    }
  }
  if (across_face_count - facecount > across_max_count) {
    across_max_count = across_face_count - facecount;
  }

  *searchtet = reftet;
  return dir;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getsteinerpointonsegment()    Get a Steiner point on a segment.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::getsteinerptonsegment(face* seg, point refpt, point steinpt)
{
  point ei, ej;
  REAL Li, Lj, L, dj, dr;
  REAL ti = 0.0, tj = 0.0, t;
  int type, eid = 0, i;

  REAL diff, stept = 0.0, L1;
  int iter;

  ei = sorg(*seg);
  ej = sdest(*seg);


  if (b->verbose > 2) {
    printf("      Get Steiner point on seg [%d (%c), %d (%c)].\n", 
           pointmark(ei), pointtype(ei) == ACUTEVERTEX ? 'A' : 'N',  
           pointmark(ej), pointtype(ej) == ACUTEVERTEX ? 'A' : 'N');
  }

  if (b->psc) {
    eid = shellmark(*seg);
    if (pointtype(ei) != FREESEGVERTEX) {
      ti = in->getvertexparamonedge(in->geomhandle, pointmark(ei), eid);
    } else {
      ti = pointgeomuv(ei, 0);
    }
    if (pointtype(ej) != FREESEGVERTEX) {
      tj = in->getvertexparamonedge(in->geomhandle, pointmark(ej), eid);
    } else {
      tj = pointgeomuv(ej, 0);
    }
  }

  if (refpt != NULL) {
    if (pointtype(ei) == ACUTEVERTEX) {
      if (pointtype(ej) == ACUTEVERTEX) {
        // Choose the vertex which is closer to refpt.
        Li = distance(ei, refpt);
        Lj = distance(ej, refpt);
        if (Li > Lj) {
          // Swap ei and ej;
          sesymself(*seg);
          ei = sorg(*seg);
          ej = sdest(*seg);
          t = ti;
          ti = tj;
          tj = t;
        }
        type = 1;
      } else {
        type = 1;
      }
    } else {
      if (pointtype(ej) == ACUTEVERTEX) {
        type = 1;
        // Swap ei and ej;
        sesymself(*seg);
        ei = sorg(*seg);
        ej = sdest(*seg);
        t = ti;
        ti = tj;
        tj = t;
      } else {
        type = 0;
      }
    }
  } else {
    type = 0;
  }

  if (type == 1) {
    L = distance(ei, ej);
    Li = distance(ei, refpt);
    // Cut the segment by a sphere centered at ei with radius Li.
    if (b->psc) {
      stept = (tj - ti) / 100.0;
      iter = 0;
      t = ti + (Li / L) * (tj - ti);
      while (1) {
        in->getsteineronedge(in->geomhandle, eid, t, steinpt);
        L1 = distance(steinpt, ei);
        diff = L1 - Li;
        if ((fabs(diff) / L) < 1e-3) {
          break;
        }
        if (diff > 0) {
          t -= stept; // Move it towards ei.
        } else {
          t += stept; // Move it towards ej.
        }
        iter++;
        if (iter > 10) {
          printf("Warning:  Get the right Steiner point failed.\n");
          break;
        }
      } // while (1)
    } else {
      t = Li / L;
      for (i = 0; i < 3; i++) {
        steinpt[i] = ei[i] + t * (ej[i] - ei[i]);
      }
    }
    // Avoid creating a too short edge.
    dj = distance(steinpt, ej);
    dr = distance(steinpt, refpt);
    if (dj < dr) {
      // Cut the segment by the radius equal to Li / 2.
      if (b->psc) {
        iter = 0;
        t = ti + ((Li / 2.0) / L) * (tj - ti);
        while (1) {
          in->getsteineronedge(in->geomhandle, eid, t, steinpt);
          L1 = distance(steinpt, ei);
          diff = L1 - (Li / 2.0);
          if ((fabs(diff) / L) < 1e-3) {
            break;
          }
          if (diff > 0) {
            t -= stept; // Move it towards ei.
          } else {
            t += stept; // Move it towards ej.
          }
          iter++;
          if (iter > 10) {
            printf("Warning:  Get the right Steiner point failed.\n");
            break;
          }
        } // while (1)
      } else {
        t = (Li / 2.0) / L;
        for (i = 0; i < 3; i++) {
          steinpt[i] = ei[i] + t * (ej[i] - ei[i]);
        }
      }
      r3count++;
    } else {
      r2count++;
    }
  } else {
    // Split the point at the middle.
    if (b->psc) {
      t = 0.5 * (ti + tj);
      in->getsteineronedge(in->geomhandle, eid, t, steinpt);
    } else {
      t = 0.5;
      for (i = 0; i < 3; i++) {
        steinpt[i] = ei[i] + t * (ej[i] - ei[i]);
      }
    }
    r1count++;
  }

  if (b->psc) {
    setpointgeomuv(steinpt, 0, t);
    setpointgeomtag(steinpt, eid);
  }

  if (pointtype(steinpt) == UNUSEDVERTEX) {
    setpointtype(steinpt, FREESEGVERTEX);
  }

  if (b->verbose > 2) {
    printf("      Split at t(%g)", t);
    if (b->psc) {
      printf(", ti(%g), tj(%g)", ti, tj);
    }
    printf(".\n");
  }
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizesegments()    Recover segments in a DT.                          //
//                                                                           //
// All segments need to be recovered are in 'subsegstack' (Q).  They will be //
// be recovered one by one (in a random order).                              //
//                                                                           //
// Given a segment s in the Q, this routine first queries s in the DT, if s  //
// matches an edge in DT, it is 'locked' at the edge. Otherwise, s is split  //
// by inserting a new point p in both the DT and itself. The two new subseg- //
// ments of s are queued in Q.  The process continues until Q is empty.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunizesegments()
{
  triface searchtet, spintet;
  face searchsh, checksh;
  face sseg, checkseg, *psseg;
  point refpt, newpt;
  enum interresult dir;
  insertvertexflags ivf;
  int loc;

  // For reporting PLC problems.
  point forg1, fdest1;   // The 1st segment.
  point forg2, fdest2, fapex2;   // The 2nd segment.

  // Does this mesh containing subfaces? 
  if (checksubfaceflag) {
    ivf.bowywat = 2;   // The mesh is a CDT. 
    ivf.lawson = 2;    // Do flip to recover Delaunayness.
    ivf.validflag = 1; // Validation is needed.
  } else {
    ivf.bowywat = 1;   // The mesh is a DT.
    ivf.lawson = 0;    // No need to do flip.
    ivf.validflag = 0; // No need to valid the B-W cavity.
  }

  searchsh.sh = NULL;

  // Loop until 'subsegstack' is empty.
  while (subsegstack->objects > 0l) {
    // seglist is used as a stack.
    subsegstack->objects--;
    psseg = (face *) fastlookup(subsegstack, subsegstack->objects);
    sseg = *psseg;

    assert(!sinfected(sseg)); 
    // Check if this segment has been recovered.
    sstpivot1(sseg, searchtet);
    if (searchtet.tet != NULL) {
      // Check if the tet contains the same segment.
      tsspivot1(searchtet, checkseg); // SELF_CHECK
      assert(checkseg.sh == sseg.sh);
      continue; // Not a missing segment.
    }

    // Search the segment.
    dir = scoutsegment(sorg(sseg), sdest(sseg), &searchtet, &refpt, NULL);

    if (dir == SHAREEDGE) {
      // Found this segment, insert it.
      tsspivot1(searchtet, checkseg);  // SELF_CHECK
      if (checkseg.sh == NULL) {
        // Let the segment remember an adjacent tet.
        sstbond1(sseg, searchtet);
        // Bond the segment to all tets containing it.
        spintet = searchtet;
        do {
          tssbond1(spintet, sseg);
          fnextself(spintet);
        } while (spintet.tet != searchtet.tet);
      } else {
        // Collision! Should not happen.
        assert(0);
      }
    } else {
      if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
        // The segment is missing. Split it.
        // Create a new point.
        makepoint(&newpt, FREESEGVERTEX);
        //setpointtype(newpt, FREESEGVERTEX);
        getsteinerptonsegment(&sseg, refpt, newpt);

        // Start searching from the 'searchtet'.
        ivf.iloc = (int) OUTSIDE;
        //ivf.bowywat;
        //ivf.lawson;
        ivf.rejflag = 0;
        ivf.chkencflag = 0;
        ivf.sloc = ivf.iloc;
        ivf.sbowywat = ivf.bowywat;
        ivf.splitbdflag = 0;
        // ivf.validflag
        ivf.respectbdflag = 0;
        ivf.assignmeshsize = b->metric;
        // Insert the new point into the tetrahedralization T.
        //   Missing segments and subfaces are queued for recovery.
        //   Note that T is convex (nonconvex = 0).
        loc = insertvertex(newpt, &searchtet, &searchsh, &sseg, &ivf);

        assert(loc != (int) ONVERTEX);
        if (loc != (int) NEARVERTEX) {
          // The new point has been inserted.
          if (ivf.lawson > 0) {
            // For CDT, use flips to reocver Delaunayness.
            lawsonflip3d(newpt, ivf.lawson, 0, 0, 0);
          }
          st_segref_count++;
          if (steinerleft > 0) steinerleft--;
        } else {
          // The new point is either ON or VERY CLOSE to an existing point.
          refpt = point2ppt(newpt);
          printf("  !! Avoid to create a short edge (length = %g)\n",
                 distance(newpt, refpt));

          // It is probably an input problem. Two possible cases are:
          //   (1) An input vertex is very close an input segment; or
          //   (2) Two input segments are nearly intersect each other.
          forg1 = farsorg(sseg);
          fdest1 = farsdest(sseg);

          if ((pointtype(refpt) == RIDGEVERTEX) || 
	      (pointtype(refpt) == ACUTEVERTEX) ||
              (pointtype(refpt) == VOLVERTEX)) {
            // Case (1)
            printf("  !! Point %d is very close to segment (%d, %d).\n", 
                   pointmark(refpt), pointmark(forg1), pointmark(fdest1));
          } else if (pointtype(refpt) == FREESEGVERTEX) {
            // Case (2). Find a subsegment contain 'refpt'.
            subsegs->traversalinit();
            checkseg.sh = shellfacetraverse(subsegs);
            while (checkseg.sh != NULL) {
              if (((point) checkseg.sh[3] == refpt) || 
                  ((point) checkseg.sh[4] == refpt)) break;
              checkseg.sh = shellfacetraverse(subsegs);
            }
            assert(checkseg.sh != NULL);
            checkseg.shver = 0;
            forg2 = farsorg(checkseg);
            fdest2 = farsdest(checkseg);
            printf("  !! Two segments are very close to each other.\n");
            printf("  1st: (%d, %d), 2nd: (%d, %d)\n", pointmark(forg1), 
                   pointmark(fdest1), pointmark(forg2), pointmark(fdest2));
          } else {
            // Unknown case
            assert(0);
          }
          // Indicate it may be an input problem.
          printf("  Short edge length bound is: %g. Tolerance is %g.\n", 
                 b->minedgelength, b->epsilon);
          terminatetetgen(4);
        }
      } else {
        // The input PLC contains self-intersections.
        if (dir == ACROSSVERT) {
          // refpt is the vertex intersecting the segment.
          forg1 = farsorg(sseg);
          fdest1 = farsdest(sseg);
          if ((pointtype(refpt) == RIDGEVERTEX) || 
	      (pointtype(refpt) == ACUTEVERTEX) ||
              (pointtype(refpt) == FACETVERTEX) ||
              (pointtype(refpt) == VOLVERTEX)) {
            printf("Point %d is on segment (%d, %d).\n", 
                   pointmark(refpt), pointmark(forg1), pointmark(fdest1));
          } else if (pointtype(refpt) == FREESEGVERTEX) {
            // Case (2). Find a subsegment contain 'refpt'.
            subsegs->traversalinit();
            checkseg.sh = shellfacetraverse(subsegs);
            while (checkseg.sh != NULL) {
              if (((point) checkseg.sh[3] == refpt) || 
                  ((point) checkseg.sh[4] == refpt)) break;
              checkseg.sh = shellfacetraverse(subsegs);
            }
            assert(checkseg.sh != NULL);
            checkseg.shver = 0;
            forg2 = farsorg(checkseg);
            fdest2 = farsdest(checkseg);
            printf("Two segments intersect.\n");
            printf("    1st: (%d, %d), 2nd: (%d, %d)", pointmark(forg1), 
                   pointmark(fdest1), pointmark(forg2), pointmark(fdest2));
          } else if (pointtype(refpt) == FREEFACETVERTEX) {
            assert(0); // Report this case.
          }
        } else if (dir == ACROSSSEG) {
          tsspivot1(searchtet, checkseg);
          if (!b->quiet) {
            printf("Two segments intersect.\n");
            forg1 = farsorg(sseg);
            fdest1 = farsdest(sseg);
            forg2 = farsorg(checkseg);
            fdest2 = farsdest(checkseg);
            printf("  1st: (%d, %d), 2nd: (%d, %d).\n", pointmark(forg1), 
                   pointmark(fdest1), pointmark(forg2), pointmark(fdest2));
          }
        } else if (dir == ACROSSSUB) {
          tspivot(searchtet, checksh);
          if (!b->quiet) {
            printf("A segment and a subface intersect.\n");
            forg1 = farsorg(sseg);
            fdest1 = farsdest(sseg);
            forg2 = sorg(checksh);
            fdest2 = sdest(checksh);
            fapex2 = sapex(checksh);
            printf("  Seg: (%d, %d), Sub: (%d, %d, %d).\n", 
                   pointmark(forg1), pointmark(fdest1), 
                   pointmark(forg2), pointmark(fdest2), pointmark(fapex2));
          }
        } else {
          // Unknown cases.
          assert(0);
        }
        // Indicate it is an input problem.
        terminatetetgen(3);
      }
    }
  } // while

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutsubface()    Look for a given subface in the tetrahedralization T.   //
//                                                                           //
// 'searchsh' is searched in T. If it exists, it is 'locked' at the face in  //
// T. 'searchtet' refers to the face. Otherwise, it is missing.              //
//                                                                           //
// The return value indicates one of the following cases:                    //
//   - SHAREFACE, 'searchsh' exists and is inserted in T.                    //
//   - COLLISIONFACE, 'searchsh' exists in T, but it conflicts with another  //
//     subface which was inserted earlier. It is not inserted.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

enum tetgenmesh::interresult 
  tetgenmesh::scoutsubface(face* searchsh, triface* searchtet)
{
  triface spintet;
  face checksh;
  point pa, pb, pc;
  enum interresult dir;

  pa = sorg(*searchsh);
  pb = sdest(*searchsh);

  if (b->verbose > 2) {
    printf("      Scout subface (%d, %d, %d).\n", pointmark(pa), pointmark(pb),
           pointmark(sapex(*searchsh)));
  }

  // Get a tet whose origin is a.
  point2tetorg(pa, *searchtet);
  // Search the edge [a,b].
  dir = finddirection(searchtet, pb);
  if (dir == ACROSSVERT) {
    // Check validity of a PLC.
    if (dest(*searchtet) != pb) {
      // A vertex lies on the search edge. Return it.
      enextself(*searchtet);
      return TOUCHEDGE;
    }
    // The edge exists. Check if the face exists.
    pc = sapex(*searchsh);
    // Searchtet holds edge [a,b]. Search a face with apex c.
    spintet = *searchtet;
    while (1) {
      if (apex(spintet) == pc) {
        // Found a face matching to 'searchsh'!
        tspivot(spintet, checksh);
        if (checksh.sh == NULL) {
          // Insert 'searchsh'.
          tsbond(spintet, *searchsh);
          fsymself(spintet);
          sesymself(*searchsh);
          tsbond(spintet, *searchsh);
          *searchtet = spintet;
          return SHAREFACE;
        } else {
          // Another subface is already inserted.
          assert(checksh.sh != searchsh->sh); // SELF_CHECK
          // This is possibly an input problem, i.e., two facets overlap.
          // Report this problem and exit.
          printf("Warning:  Found two facets nearly overlap.\n");
          terminatetetgen(5);
          // unifysubfaces(&checksh, searchsh);
          *searchtet = spintet;
          return COLLISIONFACE;
        }
      }
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    }
  }

  // dir is either ACROSSEDGE or ACROSSFACE.
  return dir; //ACROSSTET;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// formmissingregion()    Form the missing region of a missing subface.      //
//                                                                           //
// 'missh' is a missing subface. From it we form a missing region R which is //
// a collection of missing subfaces connected through adjacent edges.        //
//                                                                           //
// The missing region R is returned in the array 'missingshs'.  All subfaces //
// in R are oriented as 'missh'. The array 'missingshverts' returns all ver- //
// tices of R. All subfaces and vertices of R are marktested.                //
//                                                                           //
// 'adjtets' returns a list of tetrahedra adjacent to R.  They are used to   //
// search a crossing tetrahedron of R.                                       //
//                                                                           //
// Many ways are possible to form the missing region.  The method used here  //
// is to search missing edges in R. Starting from 'missh', its three edges   //
// are checked. If one of the edges is missing, then the adjacent subface at //
// this edge is also missing. It is added to the array. By an incrementally  //
// broad-first searching, we can find all subfaces of R.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::formmissingregion(face* missh, arraypool* missingshs,  
                                   arraypool* missingshbds, 
                                   arraypool* missingshverts, 
                                   arraypool* adjtets)
{
  triface searchtet, *parytet;
  face neighsh, *parysh;
  point pa, pb, *parypt;
  enum interresult dir;
  int i, j;

  if (b->verbose > 2) {
    printf("      Form missing region from subface (%d, %d, %d)\n", 
           pointmark(sorg(*missh)), pointmark(sdest(*missh)), 
           pointmark(sapex(*missh)));
  }
  smarktest(*missh);
  missingshs->newindex((void **) &parysh);
  *parysh = *missh;

  // Incrementally find other missing subfaces.
  for (i = 0; i < missingshs->objects; i++) {
    missh = (face *) fastlookup(missingshs, i);
    for (j = 0; j < 3; j++) {
      pa = sorg(*missh);
      pb = sdest(*missh);
      // Get a tet whose origin is a.
      point2tetorg(pa, searchtet);
      // Search the edge [a,b].
      dir = finddirection(&searchtet, pb);
      if (dir != ACROSSVERT) {
        // This edge is missing. Its neighbor is a missing subface.
        spivot(*missh, neighsh);
        assert(neighsh.sh != NULL);
        if (!smarktested(neighsh)) {
          // Adjust the face orientation.
          if (sorg(neighsh) != pb) {
            sesymself(neighsh);
          }
          if (b->verbose > 3) {
            printf("      Add a missing subface (%d, %d, %d)\n", 
                   pointmark(pb), pointmark(pa), pointmark(sapex(neighsh)));
          }
          smarktest(neighsh);
          missingshs->newindex((void **) &parysh);
          *parysh = neighsh;
        }
      } else {
        if (dest(searchtet) == pb) {
          // Remember an existing edge for searching the first crossing tet.
          adjtets->newindex((void **) &parytet);
          *parytet = searchtet;
          // Found an existing edge, it must be a boundary edge of R.
          if (b->verbose > 3) {
            printf("      -- A boundary edge (%d, %d)\n", pointmark(pa), 
                   pointmark(pb));
          }
          missingshbds->newindex((void **) &parysh);
          *parysh = *missh; // It is only queued once.
        } else {
          // The input PLC has problem.
          //assert(0);
          terminatetetgen(3);
        }
      }
      // Collect the vertices of R.
      if (!pmarktested(pa)) {
        pmarktest(pa);
        missingshverts->newindex((void **) &parypt);
        *parypt = pa;
      }
      senextself(*missh);
    } // j
  } // i

  if (b->verbose > 2) {
    printf("      Region has: %ld subfaces, %ld vertices\n", 
           missingshs->objects, missingshverts->objects);
  }

  if (missingshs->objects > maxregionsize) {
    maxregionsize = missingshs->objects;
  }

  // Unmarktest collected missing subfaces.
  for (i = 0; i < missingshs->objects; i++) {
    missh = (face *) fastlookup(missingshs, i);
    sunmarktest(*missh);
  }

  // Comment: All vertices in R are pmarktested.
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutcrossedge()    Search an edge that crosses the missing region.       //
//                                                                           //
// Assumption: All vertices of the missing region are marktested.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::scoutcrossedge(triface& crosstet, arraypool* adjtets, 
                               arraypool* missingshs)
{
  triface *searchtet, spintet;
  face *parysh;
  face checkseg;
  point pa, pb, pc, pd, pe;
  enum interresult dir;
  REAL ori;
  int types[2], poss[4];
  int searchflag, interflag;
  int i, j;

  if (b->verbose > 2) {
    printf("      Search a crossing edge.\n");
  }
  searchflag = 0;

  for (j = 0; j < adjtets->objects && !searchflag; j++) {
    searchtet = (triface *) fastlookup(adjtets, j);
    interflag = 0;
    // Let 'spintet' be [#,#,d,e] where [#,#] is the boundary edge of R.
    spintet = *searchtet;
    while (1) {
      pd = apex(spintet);
      pe = oppo(spintet);
      // Skip a hull edge.
      if ((pd != dummypoint) && (pe != dummypoint)) {
        // Skip an edge containing a vertex of R.
        if (!pmarktested(pd) && !pmarktested(pe)) {
          // Check if [d,e] intersects R.
          for (i = 0; i < missingshs->objects && !interflag; i++) {
            parysh = (face *) fastlookup(missingshs, i);
            pa = sorg(*parysh);
            pb = sdest(*parysh);
            pc = sapex(*parysh);
            interflag = tri_edge_test(pa, pb, pc, pd, pe, NULL, 1, types, poss);
            if (interflag > 0) { 
              if (interflag == 2) {
                // They intersect at a single point.
                dir = (enum interresult) types[0];
                if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
                  //pos = poss[0];
                  // Go to the crossing edge [d,e,#,#].
                  eprev(spintet, crosstet);
                  esymself(crosstet);
                  enextself(crosstet); // [d,e,#,#].
                  // Check if it is a segment.
                  tsspivot1(crosstet, checkseg);
                  if (checkseg.sh != NULL) {
                    reportselfintersect(&checkseg, parysh);
                    terminatetetgen(3);
                  }
                  // Adjust the edge such that d lies below [a,b,c].
                  ori = orient3d(pa, pb, pc, pd);
                  assert(ori != 0);
                  if (ori < 0) {
                    esymself(crosstet);
                  }
                  if (b->verbose > 2) {
                    printf("      Found edge (%d, %d) intersect", pointmark(pd),
                           pointmark(pe));
                    printf(" face (%d, %d, %d)\n", pointmark(pa), pointmark(pb),
                           pointmark(pc));
                  }
                  // Save the corners of this subface.
                  plane_pa = pa;
                  plane_pb = pb;
                  plane_pc = pc;              
                  searchflag = 1;                  
                } else {
                  // An improper intersection type.
                  // Maybe it is a PLC problem.
                  // At the moment, just ignore it.
                }
              }
              break;
            } // if (interflag > 0)
          }
        } 
      }
      // Leave search at this bdry edge if an intersection is found.
      if (interflag > 0) break;
      // Go to the next tetrahedron.
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break; 
    } // while (1)
  } // j

  adjtets->restart();
  return searchflag;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// formcavity()    Form the cavity of a missing region.                      //
//                                                                           //
// The missing region R is formed by a set of missing subfaces 'missingshs'. //
// In the following, we assume R is horizontal and oriented. (All subfaces   //
// of R are oriented in the same way.)  'searchtet' is a tetrahedron [d,e,#, //
// #] which intersects R in its interior, where the edge [d,e] intersects R, //
// and d lies below R.                                                       //
//                                                                           //
//                                                                           //
// 'crosstets' returns the set of crossing tets. Every tet in it has the     //
// form [d,e,#,#] where [d,e] is a crossing edge, and d lies below R.  The   //
// set of tets form the cavity C, which is divided into two parts by R, one  //
// at top and one at bottom. 'topfaces' and 'botfaces' return the upper and  //
// lower boundary faces of C. 'toppoints' contains vertices of 'crosstets'   //
// in the top part of C, and so does 'botpoints'. Both 'toppoints' and       //
// 'botpoints' contain vertices of R.                                        //
//                                                                           //
// NOTE: 'toppoints' may contain points which are not vertices of any top    //
// faces, and so may 'botpoints'. Such points may belong to other facets and //
// need to be present after the recovery of this cavity (P1029.poly).        //
//                                                                           //
// A pair of boundary faces: 'firsttopface' and 'firstbotface', are saved.   //
// They share the same edge in the boundary of the missing region.           //
//                                                                           //
// Important: This routine assumes all vertices of the facet containing this //
// subface are marked, i.e., pmarktested(p) returns true.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::formcavity(triface* searchtet, arraypool* missingshs,
                            arraypool* crosstets, arraypool* topfaces, 
                            arraypool* botfaces, arraypool* toppoints, 
                            arraypool* botpoints)
{
  arraypool *crossedges, *testededges;
  triface spintet, neightet, *parytet;
  face checksh, *parysh = NULL;
  face checkseg; // *paryseg;
  point pa, pd, pe, *parypt;
  enum interresult dir; 
  bool testflag, invalidflag;
  int types[2], poss[4];
  int i, j, k;

  // Temporarily re-use 'topfaces' for all crossing edges.
  crossedges = topfaces;
  // Temporarily re-use 'botfaces' for all tested edges.
  testededges = botfaces; // Only used by 'b->psc'.

  if (b->verbose > 2) {
    printf("      Form the cavity of missing region.\n"); 
  }
  missingsubfacecount += missingshs->objects;
  // Mark this edge to avoid testing it later.
  markedge(*searchtet);
  crossedges->newindex((void **) &parytet);
  *parytet = *searchtet;

  invalidflag = 0; 

  // Collect all crossing tets.  Each cross tet is saved in the standard
  //   form [d,e,#,#], where [d,e] is a corossing edge, d lies below R.
  //   NEITHER d NOR e is a vertex of R (!pmarktested). 
  for (i = 0; i < crossedges->objects; i++) {
    // Get a crossing edge [d,e,#,#].
    searchtet = (triface *) fastlookup(crossedges, i);

    // Sort vertices into the bottom and top arrays.
    pd = org(*searchtet);
    assert(!pmarktested(pd)); // pd is not on R.
    if (!pinfected(pd)) {
      pinfect(pd);
      botpoints->newindex((void **) &parypt);
      *parypt = pd;
    }
    pe = dest(*searchtet);
    assert(!pmarktested(pe)); // pe is not on R.
    if (!pinfected(pe)) {
      pinfect(pe);
      toppoints->newindex((void **) &parypt);
      *parypt = pe;
    }

    // All tets sharing this edge are crossing tets.
    spintet = *searchtet;
    while (1) {
      if (!infected(spintet)) {
        if (b->verbose > 3) {
          printf("      Add a crossing tet (%d, %d, %d, %d)\n",
                 pointmark(org(spintet)), pointmark(dest(spintet)),
                 pointmark(apex(spintet)), pointmark(oppo(spintet)));
        }
        infect(spintet);
        crosstets->newindex((void **) &parytet);
        *parytet = spintet;
      }
      // Go to the next crossing tet.
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)

    // Detect new crossing edges.
    spintet = *searchtet;
    while (1) {
      // spintet is [d,e,a,#], where d lies below R, and e lies above R. 
      pa = apex(spintet);
      if (pa != dummypoint) {
        if (!pmarktested(pa) || b->psc) {
	  // There exists a crossing edge, either [e,a] or [a,d]. First check
          //   if the crossing edge has already be added. This is to check if
          //   a tetrahedron at this edge is marked.
          testflag = true;
          for (j = 0; j < 2 && testflag; j++) {
            if (j == 0) {
              enext(spintet, neightet);
            } else {
              eprev(spintet, neightet);
            }
            while (1) {
              if (edgemarked(neightet)) {
                // This crossing edge has already been tested. Skip it.
                testflag = false;
                break;
              }
              fnextself(neightet);
              if (neightet.tet == spintet.tet) break;
            }
          } // j
          if (testflag) {
            // Test if [e,a] or [a,d] intersects R.
            // Do a brute-force search in the set of subfaces of R. Slow!
            //   Need to be improved!
            pd = org(spintet);
            pe = dest(spintet);
            for (k = 0; k < missingshs->objects; k++) {
              parysh = (face *) fastlookup(missingshs, k);
              plane_pa = sorg(*parysh);
              plane_pb = sdest(*parysh);
              plane_pc = sapex(*parysh);
              // Test if this face intersects [e,a].
              if (tri_edge_test(plane_pa, plane_pb, plane_pc, pe, pa, 
                                NULL, 1, types, poss)) {
                // Found intersection. 'a' lies below R.
                enext(spintet, neightet);
                dir = (enum interresult) types[0];
                if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
                  // A valid intersection.
                } else {
                  // A non-valid intersection. Maybe a PLC problem.
                  invalidflag = 1;
                }
                break;
              }
              // Test if this face intersects [a,d].
              if (tri_edge_test(plane_pa, plane_pb, plane_pc, pa, pd, 
                                NULL, 1, types, poss)) {
                // Found intersection. 'a' lies above R.
                eprev(spintet, neightet);
                dir = (enum interresult) types[0];
                if ((dir == ACROSSFACE) || (dir == ACROSSEDGE)) {
                  // A valid intersection.
                } else {
                  // A non-valid intersection. Maybe a PLC problem.
                  invalidflag = 1;
                }
                break;
              }
            } // k
            if (k < missingshs->objects) {
              // Found a pair of triangle - edge interseciton.
              if (invalidflag) {
                if (b->verbose > 2) {
                  printf("      A non-valid subface - edge intersection\n");
                  printf("      subface: (%d, %d, %d) edge: (%d, %d)\n",
                         pointmark(plane_pa), pointmark(plane_pb), 
                         pointmark(plane_pc), pointmark(org(neightet)),
                         pointmark(dest(neightet)));
                }
                // It may be a PLC problem.
                terminatetetgen(3);
              } else if (b->psc) {
                if (pmarktested(pa)) {
                  // The intersection is invalid.
                  if (b->verbose > 2) {
                    printf("    A non-valid subface - edge intersection\n");
                    printf("      subface: (%d, %d, %d) edge: (%d, %d)\n",
                           pointmark(plane_pa), pointmark(plane_pb), 
                           pointmark(plane_pc), pointmark(org(neightet)),
                           pointmark(dest(neightet)));
                  }
                  // Split the subface intersecting this edge.
                  recentsh = *parysh;
                  recenttet = neightet; // For point location.
                  invalidflag = 1;
                  break;
                } // if (pmarktested(pa))
              } // if (b->psc)
              // Adjust the edge direction, so that its origin lies below R,
              //   and its destination lies above R.
              esymself(neightet);
              // Check if this edge is a segment.
              tsspivot1(neightet, checkseg);
              if (checkseg.sh != NULL) {
                // Invalid PLC!
                reportselfintersect(&checkseg, parysh);
                terminatetetgen(3);
              }
              if (b->verbose > 3) {
                printf("      Add a crossing edge (%d, %d)\n", 
                       pointmark(org(neightet)), pointmark(dest(neightet)));
              }
              // Mark this edge to avoid testing it again.
              markedge(neightet);
              crossedges->newindex((void **) &parytet);
              *parytet = neightet;            
            } else {
              // No intersection is found. It may be a PLC problem.
              //assert(b->psc);              
              // Mark this edge to avoid testing it again.
              //markedge(neightet);
              //testededges->newindex((void **) &parytet);
              //*parytet = neightet;
              invalidflag = 1;
              // Split the subface intersecting [d,e].
              for (k = 0; k < missingshs->objects; k++) {
                parysh = (face *) fastlookup(missingshs, k);
                plane_pa = sorg(*parysh);
                plane_pb = sdest(*parysh);
                plane_pc = sapex(*parysh);
                // Test if this face intersects [e,a].
                if (tri_edge_test(plane_pa, plane_pb, plane_pc, pd, pe, 
                                  NULL, 1, types, poss)) {
                  break;
                }
              } // k
              assert(k < missingshs->objects);
              recentsh = *parysh;
              recenttet = spintet; // For point location.
              break; // the while (1) loop
            } // if (k == missingshs->objects)
          } // if (testflag)
	} // if (!pmarktested(pa) || b->psc)
      }
      // Go to the next crossing tet.
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)

    //if (b->psc) {
      if (invalidflag) break;
    //}
  } // i

  if (b->verbose > 2) {
    printf("      Formed cavity: %ld (%ld) cross tets (edges).\n", 
           crosstets->objects, crossedges->objects);
  }
  crossingtetcount += crosstets->objects;

  // Unmark all marked edges.
  for (i = 0; i < crossedges->objects; i++) {
    searchtet = (triface *) fastlookup(crossedges, i);
    assert(edgemarked(*searchtet)); // SELF_CHECK
    unmarkedge(*searchtet);
  }
  crossedges->restart();

  if (b->psc) {
    // Unmark all marked edges.
    for (i = 0; i < testededges->objects; i++) {
      searchtet = (triface *) fastlookup(testededges, i);
      assert(edgemarked(*searchtet)); // SELF_CHECK
      unmarkedge(*searchtet);
    }
    testededges->restart();
  } else { // only p->plc
    assert(testededges->objects == 0l);
  }

  if (invalidflag) {
    // Unmark all collected tets.
    for (i = 0; i < crosstets->objects; i++) {
      searchtet = (triface *) fastlookup(crosstets, i);
      uninfect(*searchtet);
    }
    // Unmark all collected vertices.
    for (i = 0; i < botpoints->objects; i++) {
      parypt = (point *) fastlookup(botpoints, i);
      puninfect(*parypt);
    }
    for (i = 0; i < toppoints->objects; i++) {
      parypt = (point *) fastlookup(toppoints, i);
      puninfect(*parypt);
    }
    crosstets->restart();
    botpoints->restart();
    toppoints->restart();
    return false;
  }

  // Find a pair of cavity boundary faces from the top and bottom sides of
  //   the facet each, and they share the same edge. Save them in the
  //   global variables: firsttopface, firstbotface. They will be used in
  //   fillcavity() for gluing top and bottom new tets.
  for (i = 0; i < crosstets->objects; i++) {
    searchtet = (triface *) fastlookup(crosstets, i);
    // Crosstet is [d,e,a,b].
    enextesym(*searchtet, spintet);
    eprevself(spintet); // spintet is [b,a,e,d]
    fsym(spintet, neightet); // neightet is [a,b,e,#]
    if (!infected(neightet)) {
      // A top face.
      firsttopface = neightet;
    } else {
      continue; // Go to the next cross tet.
    }
    eprevesym(*searchtet, spintet);
    enextself(spintet); // spintet is [a,b,d,e]
    fsym(spintet, neightet); // neightet is [b,a,d,#]
    if (!infected(neightet)) {
      // A bottom face.
      firstbotface = neightet;
    } else {
      continue;
    }
    break;
  } // i
  assert(i < crosstets->objects); // SELF_CHECK

  // Collect the top and bottom faces and the middle vertices. Since all top
  //   and bottom vertices have been infected. Uninfected vertices must be
  //   middle vertices (i.e., the vertices of R).
  // NOTE 1: Hull tets may be collected. Process them as a normal one.
  // NOTE 2: Some previously recovered subfaces may be completely inside the
  //   cavity. In such case, we remove these subfaces from the cavity and put     //   them into 'subfacstack'. They will be recovered later.
  // NOTE 3: Some segments may be completely inside the cavity, e.g., they
  //   attached to a subface which is inside the cavity. Such segments are
  //   put in 'subsegstack'. They will be recovered later. 
  // NOTE4 : The interior subfaces and segments mentioned in NOTE 2 and 3
  //   are identified in the routine "carvecavity()". 

  for (i = 0; i < crosstets->objects; i++) {
    searchtet = (triface *) fastlookup(crosstets, i);
    // searchtet is [d,e,a,b].
    enextesym(*searchtet, spintet);
    eprevself(spintet); // spintet is [b,a,e,d]
    fsym(spintet, neightet); // neightet is [a,b,e,#]
    if (!infected(neightet)) {
      // A top face.
      topfaces->newindex((void **) &parytet);
      *parytet = neightet;
    } 
    eprevesym(*searchtet, spintet);
    enextself(spintet); // spintet is [a,b,d,e]
    fsym(spintet, neightet); // neightet is [b,a,d,#]
    if (!infected(neightet)) {
      // A bottom face.
      botfaces->newindex((void **) &parytet);
      *parytet = neightet;
    }
    // Add middle vertices if there are (skip dummypoint).
    pa = org(neightet);
    if (!pinfected(pa)) {
      if (pa != dummypoint) {
        pinfect(pa);
        botpoints->newindex((void **) &parypt);
        *parypt = pa;
        toppoints->newindex((void **) &parypt);
        *parypt = pa;
      }
    }
    pa = dest(neightet);
    if (!pinfected(pa)) {
      if (pa != dummypoint) {
        pinfect(pa);
        botpoints->newindex((void **) &parypt);
        *parypt = pa;
        toppoints->newindex((void **) &parypt);
        *parypt = pa;
      }
    }
  } // i

  // Uninfect all collected top, bottom, and middle vertices.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    puninfect(*parypt);
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    puninfect(*parypt);
  }
  cavitycount++;

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// delaunizecavity()    Fill a cavity by Delaunay tetrahedra.                //
//                                                                           //
// The cavity C to be tetrahedralized is the top or bottom part of a whole   //
// cavity. 'cavfaces' contains the boundary faces of C. NOTE: faces in 'cav- //
// faces' do not form a closed polyhedron.  The "open" side are subfaces of  //
// the missing facet. These faces will be recovered later in fillcavity().   //
//                                                                           //
// This routine first constructs the DT of the vertices. Then it identifies  //
// the half boundary faces of the cavity in DT. Possiblely the cavity C will //
// be enlarged.                                                              //
//                                                                           //
// The DT is returned in 'newtets'.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunizecavity(arraypool *cavpoints, arraypool *cavfaces, 
                                 arraypool *cavshells, arraypool *newtets, 
                                 arraypool *crosstets, arraypool *misfaces)
{
  triface searchtet, neightet, spintet, *parytet, *parytet1;
  face checksh, tmpsh, *parysh;
  face checkseg;
  point pa, pb, pc, pd, pt[3], *parypt;
  enum interresult dir;
  insertvertexflags ivf;
  REAL ori; //, ang, len;
  long baknum, bakhullsize;
  int bakchecksubsegflag, bakchecksubfaceflag;
  //int iloc;
  int i, j;


  if (b->verbose > 2) {
    printf("      Delaunizing cavity: %ld points, %ld faces.\n", 
           cavpoints->objects, cavfaces->objects);
  }
  // Remember the current number of crossing tets. It may be enlarged later.
  baknum = crosstets->objects;
  bakhullsize = hullsize;
  bakchecksubsegflag = checksubsegflag;
  bakchecksubfaceflag = checksubfaceflag;
  hullsize = 0l;
  checksubsegflag = 0;
  checksubfaceflag = 0;
  b->verbose--;  // Suppress informations for creating Delaunay tetra.
  b->plc = 0; // Do not do unifypoint();

  // Get four non-coplanar points (no dummypoint).
  parytet = (triface *) fastlookup(cavfaces, 0);
  pa = org(*parytet);
  pb = dest(*parytet);
  pc = apex(*parytet);
  pd = NULL;
  for (i = 1; i < cavfaces->objects; i++) {
    parytet = (triface *) fastlookup(cavfaces, i);
    pt[0] = org(*parytet);
    pt[1] = dest(*parytet);
    pt[2] = apex(*parytet);
    for (j = 0; j < 3; j++) {
      if (pt[j] != dummypoint) { // Do not include a hull point.
        // if (!pinfected(pt[j])) {
          ori = orient3d(pa, pb, pc, pt[j]);
          if (ori != 0) {
            pd = pt[j];
            if (ori > 0) {  // Swap pa and pb.
              pt[j] = pa; pa = pb; pb = pt[j]; 
            }
            break;
          }
	// }
      }
    }
    if (pd != NULL) break;
  }
  assert(i < cavfaces->objects); // SELF_CHECK

  // Create an init DT.
  initialdelaunay(pa, pb, pc, pd);

  // Incrementally insert the vertices (duplicated vertices are ignored).
  for (i = 0; i < cavpoints->objects; i++) {
    pt[0] = * (point *) fastlookup(cavpoints, i);
    assert(pt[0] != dummypoint); // SELF_CHECK
    searchtet = recenttet;
    ivf.iloc = (int) OUTSIDE;
    ivf.bowywat = 1;
    insertvertex(pt[0], &searchtet, NULL, NULL, &ivf);
  }

  if (b->verbose > 2) {
    printf("      Identfying %ld boundary faces of the cavity.\n", 
           cavfaces->objects);
  }

  while (1) {

    // Identify boundary faces. Mark interior tets. Save missing faces.
    for (i = 0; i < cavfaces->objects; i++) {
      parytet = (triface *) fastlookup(cavfaces, i);
      // Skip an interior face (due to the enlargement of the cavity).
      if (infected(*parytet)) continue;
      // This face may contain dummypoint (See fig/dum-cavity-case2).
      //   If so, dummypoint must be its apex.
      j = (parytet->ver & 3); // j is the face number.
      parytet->ver = epivot[j]; // [4,5,2,11]
      pt[0] = org(*parytet);
      pt[1] = dest(*parytet);
      pt[2] = apex(*parytet);
      // Create a temp subface.
      makeshellface(subfaces, &tmpsh);
      setshvertices(tmpsh, pt[0], pt[1], pt[2]);
      // Insert tmpsh in DT.
      searchtet.tet = NULL; 
      dir = scoutsubface(&tmpsh, &searchtet);
      if (dir == SHAREFACE) {
        // Inserted. Make sure that tmpsh connects an interior tet of C.
        stpivot(tmpsh, neightet);
        // neightet and tmpsh refer to the same edge [pt[0], pt[1]].
        //   If the origin of neightet is pt[1], it is inside.
        if (org(neightet) != pt[1]) {
          fsymself(neightet);
          assert(org(neightet) == pt[1]); // SELF_CHECK
          // Make sure that tmpsh is connected with an interior tet.
          sesymself(tmpsh);
          tsbond(neightet, tmpsh);
        }
        assert(dest(neightet) == pt[0]); // SELF_CHECK
      } else if (dir == COLLISIONFACE) {
        // This case is not possible anymore. 2010-02-01
        assert(0);
      } else {
        if (b->verbose > 2) {
          printf("        bdry face (%d, %d, %d) -- %d is missing\n",
                 pointmark(pt[0]), pointmark(pt[1]), pointmark(pt[2]), i);
        }
        shellfacedealloc(subfaces, tmpsh.sh);
        // Save this face in list.
        misfaces->newindex((void **) &parytet1);
        *parytet1 = *parytet;
        continue;
      }
      // Remember the boundary tet (outside the cavity) in tmpsh 
      //   (use the adjacent tet slot). 
      tmpsh.sh[0] = (shellface) encode(*parytet);
      // Save this subface.
      cavshells->newindex((void **) &parysh);
      *parysh = tmpsh;
    } // i

    if (misfaces->objects > 0) {
      if (b->verbose > 2) {
        printf("      Enlarging the cavity. %ld missing bdry faces\n", 
               misfaces->objects);
      }

      // Removing all tempoaray subfaces.
      for (i = 0; i < cavshells->objects; i++) {
        parysh = (face *) fastlookup(cavshells, i);
        stpivot(*parysh, neightet);
        tsdissolve(neightet); // Detach it from adj. tets.
        fsymself(neightet);
        tsdissolve(neightet);
        shellfacedealloc(subfaces, parysh->sh);
      }
      cavshells->restart();

      // Infect the points which are of the cavity.
      for (i = 0; i < cavpoints->objects; i++) {
        pt[0] = * (point *) fastlookup(cavpoints, i);
        pinfect(pt[0]); // Mark it as inserted.
      }

      // Enlarge the cavity.
      for (i = 0; i < misfaces->objects; i++) {
        // Get a missing face.
        parytet = (triface *) fastlookup(misfaces, i);
        if (!infected(*parytet)) {
          // Put it into crossing tet list.
          infect(*parytet);
          crosstets->newindex((void **) &parytet1);
          *parytet1 = *parytet;
          // Insert the opposite point if it is not in DT.
          pd = oppo(*parytet);
          if (!pinfected(pd)) {
            searchtet = recenttet;
            ivf.iloc = (int) OUTSIDE;
            ivf.bowywat = 1;
            insertvertex(pd, &searchtet, NULL, NULL, &ivf);
            if (b->verbose > 2) {
              printf("      Add point %d into list.\n", pointmark(pd));
            }
            pinfect(pd);
            cavpoints->newindex((void **) &parypt);
            *parypt = pd;
          }
          // Add three opposite faces into the boundary list.
          for (j = 0; j < 3; j++) {
            esym(*parytet, neightet);
            fsymself(neightet);
            if (!infected(neightet)) {
              if (b->verbose > 2) {
                printf("      Add a cavface (%d, %d, %d).\n",
                       pointmark(org(neightet)), pointmark(dest(neightet)),
                       pointmark(apex(neightet)));
              }
              cavfaces->newindex((void **) &parytet1);
              *parytet1 = neightet;
            } else {
            }
            enextself(*parytet);
          } // j
        } // if (!infected(parytet))
      } // i

      // Uninfect the points which are of the cavity.
      for (i = 0; i < cavpoints->objects; i++) {
        pt[0] = * (point *) fastlookup(cavpoints, i);
        puninfect(pt[0]);
      }

      misfaces->restart();
      continue;
    } // if (misfaces->objects > 0)

    break;

  } // while (1)

  // Collect all tets of the DT. All new tets are marktested.
  marktest(recenttet);
  newtets->newindex((void **) &parytet);
  *parytet = recenttet;
  for (i = 0; i < newtets->objects; i++) {
    searchtet = * (triface *) fastlookup(newtets, i);
    for (searchtet.ver = 0; searchtet.ver < 4; searchtet.ver++) {
      fsym(searchtet, neightet);
      if (!marktested(neightet)) {
        marktest(neightet);
        newtets->newindex((void **) &parytet);
        *parytet = neightet;
      }
    }
  }

  cavpoints->restart();
  cavfaces->restart();

  if (cavshells->objects > maxcavsize) {
    maxcavsize = cavshells->objects;
  }
  if (crosstets->objects > baknum) {
    // The cavity has been enlarged.
    cavityexpcount++;
  }

  // Restore the original values.
  hullsize = bakhullsize;
  checksubsegflag = bakchecksubsegflag;
  checksubfaceflag = bakchecksubfaceflag;
  b->verbose++;
  b->plc = 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// fillcavity()    Fill new tets into the cavity.                            //
//                                                                           //
// The new tets are stored in two disjoint sets(which share the same facet). //
// 'topfaces' and 'botfaces' are the boundaries of these two sets, respect-  //
// ively. 'midfaces' is empty on input, and will store faces in the facet.   //
//                                                                           //
// Important: This routine assumes all vertices of the missing region R are  //
// marktested, i.e., pmarktested(p) returns true.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::fillcavity(arraypool* topshells, arraypool* botshells,
                            arraypool* midfaces, arraypool* missingshs)
{
  arraypool *cavshells;
  triface *parytet, bdrytet, toptet, bottet, midface;
  triface neightet, spintet;
  face checksh, *parysh;
  face checkseg;
  point pa, pb, pc, pf, pg; //, *pts;
  int types[2], poss[4];
  //REAL elen[3]; //ori, len, n[3];
  bool mflag, bflag;
  int i, j, k;

  // Connect newtets to tets outside the cavity.  These connections are needed
  //   for identifying the middle faces (which belong to R).
  for (k = 0; k < 2; k++) {
    cavshells = (k == 0 ? topshells : botshells);
    if (cavshells != NULL) {
      for (i = 0; i < cavshells->objects; i++) {
        // Get a temp subface.
        parysh = (face *) fastlookup(cavshells, i);
        // Get the boundary tet outside the cavity.
        decode(parysh->sh[0], bdrytet);
        pa = org(bdrytet);
        pb = dest(bdrytet);
        pc = apex(bdrytet);
        // Get the adjacent new tet.
        stpivot(*parysh, neightet);
        assert(org(neightet) == pb); // SELF_CHECK
        assert(dest(neightet) == pa); // SELF_CHECK
        // Mark neightet as an interior tet of this cavity, 2009-04-24.
        // Comment: We know neightet is an interior tet.
        if (!infected(neightet)) {
          infect(neightet);
        }
        assert(oppo(bdrytet) != NULL); // No faked tet.
        // if (oppo(bdrytet) != NULL) {
          // Bond the two tets.
          bond(bdrytet, neightet); // Also cleared the pointer to tmpsh.
	// }
        tsdissolve(neightet); // Clear the pointer to tmpsh.
        // Update the point-to-tets map.
        setpoint2tet(pa, encode(neightet));
        setpoint2tet(pb, encode(neightet));
        setpoint2tet(pc, encode(neightet));
        // Delete the temp subface.
        // shellfacedealloc(subfacepool, parysh->sh);
      } // i
    } // if (cavshells != NULL)
  } // k

  mflag = true;  // Initialize it.

  if (midfaces != NULL) {

    // The first pair of top and bottom tets share the same edge [a, b].
    // toptet = * (triface *) fastlookup(topfaces, 0);
    if (infected(firsttopface)) {
      // This is due to he enlargement of the cavity. Find the updated top
      //   boundary face at edge [a,b].
      // Comment: An uninfected tet at [a,b] should be found since [a,b] is a
      //   boundary edge of the missing region R. It should not be enclosed
      //   by the enlarged cavity.
      pa = apex(firsttopface); // SELF_CHECK
      while (1) {
        fnextself(firsttopface);
        if (!infected(firsttopface)) break;
        assert(apex(firsttopface) != pa); // SELF_CHECK
      }
    }
    toptet = firsttopface;
    pa = apex(toptet);
    fsymself(toptet);
    // Search a subface from the top mesh.
    while (1) {
      esymself(toptet); // The next face in the same tet.
      pc = apex(toptet);
      assert(pc != pa);  // We should not return to the starting point.
      if (pmarktested(pc)) break; // [a,b,c] is a subface.
      fsymself(toptet); // Go to the adjacent tet.
    }
    // Search the subface [a,b,c] in the bottom mesh.
    // bottet = * (triface *) fastlookup(botfaces, 0);
    if (infected(firstbotface)) {
      pa = apex(firstbotface); // SELF_CHECK
      while (1) {
        fnextself(firstbotface);
        if (!infected(firstbotface)) break;
        assert(apex(firstbotface) != pa); // SELF_CHECK
      }
    }
    bottet = firstbotface;
    pa = apex(bottet);
    fsymself(bottet);
    while (1) {
      esymself(bottet); // The next face in the same tet.
      pf = apex(bottet);
      assert(pf != pa); // We should not return to the starting point.
      if (pf == pc) break; // Face matched.
      if (pmarktested(pf)) {
        mflag = false; break; // Not matched.
      }
      fsymself(bottet);
    }
    if (mflag) {
      // Connect the two tets together.
      bond(toptet, bottet);
      // Both are interior tets.
      infect(toptet);
      infect(bottet);
      // Add this face into search list.
      markface(toptet);
      midfaces->newindex((void **) &parytet);
      *parytet = toptet;
    }

    // Match pairs of subfaces (middle faces), connect top and bottom tets.
    for (i = 0; i < midfaces->objects && mflag; i++) {
      // Get a matched middle face [a, b, c]
      midface = * (triface *) fastlookup(midfaces, i);
      // The tet must be a new created tet (marktested).
      assert(marktested(midface)); // SELF_CHECK

      // Check the neighbors at edges [b, c] and [c, a].
      for (j = 0; j < 2 && mflag; j++) {
        enextself(midface); // [b, c] or [c, a].
        pg = apex(midface);
        toptet = midface;
        bflag = false;
        while (1) {
          // Go to the next face in the same tet.
          esymself(toptet);
          pc = apex(toptet);
          if (pmarktested(pc)) {
            break; // Find a subface.
          }
          if (pc == dummypoint) {
            break; // Find a subface.
          }
          // Go to the adjacent tet.
          fsymself(toptet);
          // Do we walk outside the cavity? 
          if (!marktested(toptet)) {
            // Yes, the adjacent face is not a middle face.
            bflag = true; break; 
          }
        }
        if (!bflag) {
          // assert(marktested(toptet)); // SELF_CHECK
          if (!facemarked(toptet)) {
            fsym(midface, bottet);
            while (1) {
              esymself(bottet);
              pf = apex(bottet);
              if (pf == pc) break; // Face matched.
              if (pmarktested(pf)) {
                mflag = false; break; // Not matched
              }
              fsymself(bottet);
            }
            if (mflag) {
              if (marktested(bottet)) {
                // Connect two tets together.
                bond(toptet, bottet);
                // Both are interior tets.
                infect(toptet);
                infect(bottet);
                // Add this face into list.
                markface(toptet);
                midfaces->newindex((void **) &parytet);
                *parytet = toptet;
              } else {
                // The 'bottet' is not inside the cavity! 
                // This case can happen when the cavity was enlarged, and the
                //   'toptet' is a co-facet (sub)face adjacent to the missing
                //   region, and it is a boundary face of the top cavity.
                // So the toptet and bottet should be bonded already through
                //   a temp subface. See fig/dump-cavity-case18. Check it.
                fsym(toptet, neightet);
                assert(neightet.tet == bottet.tet); // SELF_CHECK
                assert(neightet.ver == bottet.ver); // SELF_CHECK
                // Do not add this face into 'midfaces'.
              }
            }
          } // if (!facemarked(toptet))
        }
      } // j
    } // i

  } // if (midfaces != NULL)

  if (mflag) {
    if (midfaces != NULL) {
      if (b->verbose > 2) {
        printf("      Found %ld middle subfaces.\n", midfaces->objects);
      }
      if (midfaces->objects > maxregionsize) {
        maxregionsize = midfaces->objects;
      }
      // Unmark middle faces.
      for (i = 0; i < midfaces->objects; i++) {
        // Get a matched middle face [a, b, c]
        midface = * (triface *) fastlookup(midfaces, i);
        assert(facemarked(midface)); // SELF_CHECK
        unmarkface(midface);
      }
    }
  } else {
    // Faces at top and bottom are not matched. There exists non-Delaunay
    //   subedges. See fig/dump-cavity-case5.lua. 
    pa = org(toptet);
    pb = dest(toptet);
    pc = apex(toptet);
    pf = apex(bottet);

    pf = oppo(toptet);
    pg = oppo(bottet);
    // Find a subface in R which intersects the edge [f,g].
    for (i = 0; i < missingshs->objects; i++) {
      parysh = (face *) fastlookup(missingshs, i);
      pa = sorg(*parysh);
      pb = sdest(*parysh);
      pc = sapex(*parysh);
      if (tri_edge_test(pa, pb, pc, pf, pg, NULL, 1, types, poss)) {
        // Found a subface.
        break;
      }
    }

    if (i < missingshs->objects) { 
      // Such subface exist.
      recentsh = *parysh;
    } else {
      assert(0); // Debug this case.
    }


    // Set a tet for searching the new point.
    recenttet = firsttopface;
  }

  // Delete the temp subfaces.
  for (k = 0; k < 2; k++) {
    cavshells = (k == 0 ? topshells : botshells);
    if (cavshells != NULL) {
      for (i = 0; i < cavshells->objects; i++) {
        parysh = (face *) fastlookup(cavshells, i);
        shellfacedealloc(subfaces, parysh->sh);
      }
    }
  }

  topshells->restart();
  if (botshells != NULL) {
    botshells->restart();
  }
  if (midfaces != NULL) {
    midfaces->restart();
  }

  return mflag;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// carvecavity()    Delete old tets and outer new tets of the cavity.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::carvecavity(arraypool *crosstets, arraypool *topnewtets,
                             arraypool *botnewtets)
{
  arraypool *newtets;
  triface *parytet, *pnewtet, newtet, neightet, spintet;
  face checksh, *parysh;
  face checkseg, *paryseg;
  int i, j, k;

  if (b->verbose > 2) {
    printf("      Carve cavity: %ld old tets.\n", crosstets->objects);
  }

  // First process subfaces and segments which are adjacent to the cavity.
  //   They must be re-connected to new tets in the cavity.
  // Comment: It is possible that some subfaces and segments are completely
  //   inside the cavity. This can happen even if the cavity is not enlarged. 
  //   Before deleting the old tets, find and queue all interior subfaces
  //   and segments. They will be recovered later. 2010-05-06.

  // Collect all subfaces and segments which attached to the old tets.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    assert(infected(*parytet)); // SELF_CHECK
    for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
      tspivot(*parytet, checksh);
      if (checksh.sh != NULL) {
        if (!sinfected(checksh)) {
          sinfect(checksh);
          cavetetshlist->newindex((void **) &parysh);
          *parysh = checksh;
        }
      }
    }
    for (j = 0; j < 6; j++) {
      parytet->ver = edge2ver[j];
      tsspivot1(*parytet, checkseg);
      if (checkseg.sh != NULL) {
        if (!sinfected(checkseg)) {
          sinfect(checkseg);
          cavetetseglist->newindex((void **) &paryseg);
          *paryseg = checkseg;
        }
      }
    }
  } // i
  // Uninfect collected subfaces.
  for (i = 0; i < cavetetshlist->objects; i++) {
    checksh = * (face *) fastlookup(cavetetshlist, i);
    suninfect(checksh);
  }
  // Uninfect collected segments.
  for (i = 0; i < cavetetseglist->objects; i++) {
    checkseg = * (face *) fastlookup(cavetetseglist, i);
    suninfect(checkseg);
  }

  // Connect subfaces to new tets.
  for (i = 0; i < cavetetshlist->objects; i++) {
    parysh = (face *) fastlookup(cavetetshlist, i);
    // Get an adjacent tet at this subface.
    stpivot(*parysh, neightet);
    // Does this tet lie inside the cavity.
    if (infected(neightet)) {
      // Yes. Get the other adjacent tet at this subface.
      sesymself(*parysh);
      stpivot(*parysh, neightet);
      // Does this tet lie inside the cavity.
      if (infected(neightet)) {
        checksh = *parysh;
        if (b->verbose > 2) {
          printf("      Found an interior subface (%d, %d, %d)\n", 
                 pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                 pointmark(sapex(checksh)));
        }
        stdissolve(checksh);
        caveencshlist->newindex((void **) &parysh);
        *parysh = checksh;
      }
    }
    if (!infected(neightet)) {
      // Found an outside tet. Re-connect this subface to a new tet.
      fsym(neightet, newtet);
      assert(marktested(newtet)); // It's a new tet.
      sesymself(*parysh);
      tsbond(newtet, *parysh);
    }
  } // i
  if (b->verbose > 2) {
    printf("      %ld (%ld) cavity (interior) subfaces.\n", 
           cavetetshlist->objects, caveencshlist->objects);
  }

  for (i = 0; i < cavetetseglist->objects; i++) {
    checkseg = * (face *) fastlookup(cavetetseglist, i);
    // Check if the segment is inside the cavity.
    sstpivot1(checkseg, neightet);
    spintet = neightet;
    while (1) {
      if (!infected(spintet)) {
        // This segment is on the boundary of the cavity.
        break;
      }
      fnextself(spintet);
      if (spintet.tet == neightet.tet) {
        if (b->verbose > 2) {
          printf("      Found an interior seg (%d, %d)\n", 
	         pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
        }
        sstdissolve1(checkseg);
        caveencseglist->newindex((void **) &paryseg);
        *paryseg = checkseg;
        break;
      }
    }
    if (!infected(spintet)) {
      // A boundary segment. Connect this segment to the new tets.
      sstbond1(checkseg, spintet);
      neightet = spintet;
      while (1) {
        tssbond1(spintet, checkseg);
        fnextself(spintet);
        if (spintet.tet == neightet.tet) break;
      }
    }
  } // i
  if (b->verbose > 2) {
    printf("      %ld (%ld) cavity (interior) segments.\n", 
           cavetetseglist->objects, caveencseglist->objects);
  }

  cavetetshlist->restart();
  cavetetseglist->restart();

  // Delete the old tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    tetrahedrondealloc(parytet->tet);
  }

  crosstets->restart(); // crosstets will be re-used.

  // Collect new tets in cavity.  Some new tets have already been found 
  //   (and infected) in the fillcavity(). We first collect them.
  for (k = 0; k < 2; k++) {
    newtets = (k == 0 ? topnewtets : botnewtets);
    if (newtets != NULL) {
      for (i = 0; i < newtets->objects; i++) {
        parytet = (triface *) fastlookup(newtets, i);
        if (infected(*parytet)) {
          crosstets->newindex((void **) &pnewtet);
          *pnewtet = *parytet;
        }
      } // i
    }
  } // k

  // Now we collect all new tets in cavity.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    if (i == 0) {
      recenttet = *parytet; // Remember a live handle.
    }
    for (j = 0; j < 4; j++) {
      decode(parytet->tet[j], neightet);
      if (marktested(neightet)) { // Is it a new tet?
        if (!infected(neightet)) {
          // Find an interior tet.
          assert((point) neightet.tet[7] != dummypoint); // SELF_CHECK
          infect(neightet);
          crosstets->newindex((void **) &pnewtet);
          *pnewtet = neightet;
        }
      }
    } // j
  } // i

  // Delete outer new tets.
  for (k = 0; k < 2; k++) {
    newtets = (k == 0 ? topnewtets : botnewtets);
    if (newtets != NULL) {
      for (i = 0; i < newtets->objects; i++) {
        parytet = (triface *) fastlookup(newtets, i);
        if (infected(*parytet)) {
          // This is an interior tet.
          uninfect(*parytet);
          unmarktest(*parytet);
        } else {
          // An outer tet. Delete it.
          tetrahedrondealloc(parytet->tet);
        }
      }
    }
  }

  crosstets->restart();
  topnewtets->restart();
  if (botnewtets != NULL) {
    botnewtets->restart();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// restorecavity()    Reconnect old tets and delete new tets of the cavity.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::restorecavity(arraypool *crosstets, arraypool *topnewtets,
                               arraypool *botnewtets)
{
  triface *parytet, neightet;
  face checksh;
  face checkseg;
  point *ppt;
  int i, j;

  // Reconnect crossing tets to cavity boundary.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    assert(infected(*parytet)); // SELF_CHECK
    if (i == 0) {
      recenttet = *parytet; // Remember a live handle.
    }
    parytet->ver = 0;
    for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
      fsym(*parytet, neightet);
      if (!infected(neightet)) {
        // Restore the old connections of tets.
        bond(*parytet, neightet);
      }
    }
    // Update the point-to-tet map.
    parytet->ver = 0;
    ppt = (point *) &(parytet->tet[4]);
    for (j = 0; j < 4; j++) {
      setpoint2tet(ppt[j], encode(*parytet));
    }
  }

  // Uninfect all crossing tets.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    uninfect(*parytet);
  }

  // Delete new tets.
  for (i = 0; i < topnewtets->objects; i++) {
    parytet = (triface *) fastlookup(topnewtets, i);
    tetrahedrondealloc(parytet->tet);
  }

  if (botnewtets != NULL) {
    for (i = 0; i < botnewtets->objects; i++) {
      parytet = (triface *) fastlookup(botnewtets, i);
      tetrahedrondealloc(parytet->tet);
    }
  }

  crosstets->restart();
  topnewtets->restart();
  if (botnewtets != NULL) {
    botnewtets->restart();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipcertify()    Insert a crossing face into priority queue.              //
//                                                                           //
// A crossing face of a facet must have at least one top and one bottom ver- //
// tex of the facet.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flipcertify(triface *chkface, badface **pqueue)
{
  badface *parybf, *prevbf, *nextbf;
  triface neightet;
  face checksh;
  point p[5];
  REAL w[5];
  REAL insph, ori4;
  int topi, boti;
  int i;

  // Compute the flip time \tau.
  fsym(*chkface, neightet);

  p[0] = org(*chkface);
  p[1] = dest(*chkface);
  p[2] = apex(*chkface);
  p[3] = oppo(*chkface);
  p[4] = oppo(neightet);

  // Check if the face is a crossing face.
  topi = boti = 0;
  for (i = 0; i < 3; i++) {
    if (pmarktest2ed(p[i])) topi++;
    if (pmarktest3ed(p[i])) boti++;
  }
  if ((topi == 0) || (boti == 0)) {
    // It is not a crossing face.
    // return;
    for (i = 3; i < 5; i++) {
      if (pmarktest2ed(p[i])) topi++;
      if (pmarktest3ed(p[i])) boti++;
    }
    if ((topi == 0) || (boti == 0)) {
      // The two tets sharing at this face are on one side of the facet.
      // Check if this face is locally Delaunay (due to rounding error).
      if ((p[3] != dummypoint) && (p[4] != dummypoint)) {
        // Do not check it if it is a subface.
        tspivot(*chkface, checksh);
        if (checksh.sh == NULL) {
          insph = insphere_s(p[1], p[0], p[2], p[3], p[4]);
          assert(insph != 0);
          if (insph > 0) {
            // Add the face into queue.
            if (b->verbose > 2) {
              printf("      A locally non-Delanay face (%d, %d, %d)-%d,%d\n", 
                     pointmark(p[0]), pointmark(p[1]), pointmark(p[2]), 
                     pointmark(p[3]), pointmark(p[4]));
            }
            parybf = (badface *) flippool->alloc();
            parybf->key = 0.;  // tau = 0, do immediately.
            parybf->tt = *chkface;
            parybf->forg = p[0];
            parybf->fdest = p[1];
            parybf->fapex = p[2];
            parybf->foppo = p[3];
            parybf->noppo = p[4];
            // Add it at the top of the priority queue.
            if (*pqueue == NULL) {
              *pqueue = parybf;
              parybf->nextitem = NULL;
            } else {
              parybf->nextitem = *pqueue;
              *pqueue = parybf;
            }
          } // if (insph > 0)
        } // if (checksh.sh == NULL)
      }
      //return;
    }
    return; // Test: omit this face.
  }

  // Decide the "height" for each point.
  for (i = 0; i < 5; i++) {
    if (pmarktest2ed(p[i])) {
      // A top point has a positive weight.
      w[i] = orient3d(plane_pa, plane_pb, plane_pc, p[i]);      
      if (w[i] < 0) w[i] = -w[i];
      assert(w[i] != 0);
    } else {
      w[i] = 0;
    }
  }

  // Make sure orient3d(p[1], p[0], p[2], p[3]) > 0;
  //   Hence if (insphere(p[1], p[0], p[2], p[3], p[4]) > 0) means that
  //     p[4] lies inside the circumsphere of p[1], p[0], p[2], p[3].
  //   The same if orient4d(p[1], p[0], p[2], p[3], p[4]) > 0 means that
  //     p[4] lies below the oriented hyperplane passing through 
  //     p[1], p[0], p[2], p[3].

  insph = insphere(p[1], p[0], p[2], p[3], p[4]);
  ori4 = orient4d(p[1], p[0], p[2], p[3], p[4], w[1], w[0], w[2], w[3], w[4]);

  if (b->verbose > 2) {
    printf("      Heights: (%g, %g, %g, %g, %g)\n", w[0],w[1],w[2],w[3],w[4]);
    printf("      Insph: %g, ori4: %g, tau = %g\n", insph, ori4, -insph/ori4);
  }

  if (ori4 > 0) {
    // Add the face into queue.
    if (b->verbose > 2) {
      printf("      Insert face (%d, %d, %d) - %d, %d\n", pointmark(p[0]),
        pointmark(p[1]), pointmark(p[2]), pointmark(p[3]), pointmark(p[4]));
    }
    
    parybf = (badface *) flippool->alloc();

    parybf->key = -insph / ori4;
    parybf->tt = *chkface;
    parybf->forg = p[0];
    parybf->fdest = p[1];
    parybf->fapex = p[2];
    parybf->foppo = p[3];
    parybf->noppo = p[4];

    // Push the face into priority queue.
    //pq.push(bface);
    if (*pqueue == NULL) {
      *pqueue = parybf;
      parybf->nextitem = NULL;
    } else {
      // Search an item whose key is larger or equal to current key.
      prevbf = NULL;
      nextbf = *pqueue;
      //if (!b->flipinsert_random) { // Default use a priority queue.
        // Insert the item into priority queue.
        while (nextbf != NULL) {
          if (nextbf->key < parybf->key) {
            prevbf = nextbf;
            nextbf = nextbf->nextitem;
          } else {
            break;
          }
        }
      //} // if (!b->flipinsert_random)
      // Insert the new item between prev and next items.
      if (prevbf == NULL) {
        *pqueue = parybf;
      } else {
        prevbf->nextitem = parybf;
      }
      parybf->nextitem = nextbf;
    }
  } else if (ori4 == 0) {
    
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipinsertfacet()    Insert a facet into a CDT by flips.                  //
//                                                                           //
// The algorithm is described in Shewchuk's paper "Updating and Constructing //
// Constrained Delaunay and Constrained Regular Triangulations by Flips", in //
// Proc. 19th Ann. Symp. on Comput. Geom., 86--95, 2003.                     //
//                                                                           //
// 'crosstets' contains the set of crossing tetrahedra (infected) of the     //
// facet.  'toppoints' and 'botpoints' are points lies above and below the   //
// facet, not on the facet.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flipinsertfacet(arraypool *crosstets, arraypool *toppoints,
                                 arraypool *botpoints, arraypool *midpoints)
{
  arraypool *crossfaces, *bfacearray;
  triface fliptets[5], baktets[2], fliptet, newface;
  triface neightet, *parytet;
  face checksh;
  face checkseg;
  badface *pqueue;
  badface *popbf, bface;
  point p1, p2, pd, pe;
  point *parypt;
  REAL ori[3];
  int convcount, copcount;
  int flipflag, fcount;
  int n, i;

  long f23count, f32count, f44count;
  long totalfcount;

  f23count = flip23count;
  f32count = flip32count;
  f44count = flip44count;

  // Get three affinely independent vertices in the missing region R.
  calculateabovepoint(midpoints, &plane_pa, &plane_pb, &plane_pc);

  // Mark top and bottom points. Do not mark midpoints.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    if (!pmarktested(*parypt)) {
      pmarktest2(*parypt);
    }
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    if (!pmarktested(*parypt)) {
      pmarktest3(*parypt);
    }
  }

  // Collect crossing faces. 
  crossfaces = cavetetlist;  // Re-use array 'cavetetlist'.

  // Each crossing face contains at least one bottom vertex and
  //   one top vertex.
  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    fliptet = *parytet;
    for (fliptet.ver = 0; fliptet.ver < 4; fliptet.ver++) {
      fsym(fliptet, neightet);
      if (infected(neightet)) { // It is an interior face.
        if (!marktested(neightet)) { // It is an unprocessed face.
          crossfaces->newindex((void **) &parytet);
          *parytet = fliptet;
        }
      }
    }
    marktest(fliptet);
  }

  if (b->verbose > 1) {
    printf("    Found %ld crossing faces.\n", crossfaces->objects);
  }
  if (crossfaces->objects > maxcrossfacecount) {
    maxcrossfacecount = crossfaces->objects;
  }

  for (i = 0; i < crosstets->objects; i++) {
    parytet = (triface *) fastlookup(crosstets, i);
    unmarktest(*parytet);
    uninfect(*parytet);
  }

  // Initialize the priority queue.
  pqueue = NULL;

  for (i = 0; i < crossfaces->objects; i++) {
    parytet = (triface *) fastlookup(crossfaces, i);
    flipcertify(parytet, &pqueue);
  }
  crossfaces->restart();

  // The list for temporarily storing unflipable faces.
  bfacearray = new arraypool(sizeof(triface), 4);


  fcount = 0;  // Count the number of flips.

  // Flip insert the facet.
  while (pqueue != NULL) {

    // Pop a face from the priotity queue.
    popbf = pqueue;
    bface = *popbf;

    // Update the queue.
    pqueue = pqueue->nextitem;

    // Delete the popped item from the pool.
    flippool->dealloc((void *) popbf);

    if (!isdeadtet(bface.tt)) {
      if ((org(bface.tt) == bface.forg) && (dest(bface.tt) == bface.fdest) &&
          (apex(bface.tt) == bface.fapex) && (oppo(bface.tt) == bface.foppo)) {
        // It is still a crossing face of R.
        fliptet = bface.tt;
        fsym(fliptet, neightet);
        assert(!isdeadtet(neightet));
        if (oppo(neightet) == bface.noppo) {
          pd = oppo(fliptet);
          pe = oppo(neightet);

          if (b->verbose > 2) {
            printf("      Get face (%d, %d, %d) - %d, %d, tau = %.17g\n",
                   pointmark(bface.forg), pointmark(bface.fdest),
                   pointmark(bface.fapex), pointmark(bface.foppo),
                   pointmark(bface.noppo), bface.key);
          }
          flipflag = 0;

          // Check for which type of flip can we do.
          convcount = 3;
          copcount = 0;
          for (i = 0; i < 3; i++) {
            p1 = org(fliptet);
            p2 = dest(fliptet);
            ori[i] = orient3d(p1, p2, pd, pe);
            if (ori[i] < 0) {
              convcount--;
              //break;
            } else if (ori[i] == 0) {
              convcount--; // Possible 4-to-4 flip.
              copcount++;
              //break;
            }
            enextself(fliptet);
          }

          if (convcount == 3) {
            // A 2-to-3 flip is found.
            // The face should not be a subface.
            tspivot(fliptet, checksh);
            assert(checksh.sh == NULL);

            fliptets[0] = fliptet; // abcd, d may be the new vertex.
            fliptets[1] = neightet; // bace.
            flip23(fliptets, 1, 0, 0);  
            // Put the link faces into check list.
            for (i = 0; i < 3; i++) {
              eprevesym(fliptets[i], newface);
              crossfaces->newindex((void **) &parytet);
              *parytet = newface;
            }
            for (i = 0; i < 3; i++) {
              enextesym(fliptets[i], newface);
              crossfaces->newindex((void **) &parytet);
              *parytet = newface;
            }
            flipflag = 1;
          } else if (convcount == 2) {
            assert(copcount <= 1);
            //if (copcount <= 1) {
            // A 3-to-2 or 4-to-4 may be possible.
            // Get the edge which is locally non-convex or flat. 
            for (i = 0; i < 3; i++) {
              if (ori[i] <= 0) break;
              enextself(fliptet);
            }
            // The edge should not be a segment.
            tsspivot1(fliptet, checkseg);
            assert(checkseg.sh == NULL);

            // Collect tets sharing at this edge.
            // NOTE: This operation may collect tets which lie outside the
            //   cavity, e.g., when the edge lies on the boundary of the
            //   cavity. Do not flip if there are outside tets at this edge.
            //   2012-07-27.
            esym(fliptet, fliptets[0]); // [b,a,d,c]
            n = 0;
            do {
              p1 = apex(fliptets[n]);
              if (!(pmarktested(p1) || pmarktest2ed(p1) || pmarktest3ed(p1))) {
                // This apex is not on the cavity. Hence the face does not
                //   lie inside the cavity. Do not flip this edge.
                n = 1000; break;
              }
              fnext(fliptets[n], fliptets[n + 1]);
              n++;
            } while ((fliptets[n].tet != fliptet.tet) && (n < 5));

            if (n == 3) {
              // Found a 3-to-2 flip.
              flip32(fliptets, 1, 0, 0);
              // Put the link faces into check list.
              for (i = 0; i < 3; i++) {
                esym(fliptets[0], newface);
                crossfaces->newindex((void **) &parytet);
                *parytet = newface;
                enextself(fliptets[0]);
              }
              for (i = 0; i < 3; i++) {
                esym(fliptets[1], newface);
                crossfaces->newindex((void **) &parytet);
                *parytet = newface;
                enextself(fliptets[1]);
              }
              flipflag = 1;
            } else if (n == 4) {
              if (copcount == 1) {                
                // Found a 4-to-4 flip. 
                // Let the six vertices are: a,b,c,d,e,f, where
                //   fliptets[0] = [b,a,d,c]
                //           [1] = [b,a,c,e]
                //           [2] = [b,a,e,f]
                //           [3] = [b,a,f,d]
                // After the 4-to-4 flip, edge [a,b] is flipped, edge [e,d]
                //   is created.
                // First do a 2-to-3 flip.
                // Comment: This flip temporarily creates a degenerated
                //   tet (whose volume is zero). It will be removed by the 
                //   followed 3-to-2 flip.
                fliptets[0] = fliptet; // = [a,b,c,d], d is the new vertex.
                // fliptets[1];        // = [b,a,c,e].
                baktets[0] = fliptets[2]; // = [b,a,e,f]
                baktets[1] = fliptets[3]; // = [b,a,f,d]
                // The flip may involve hull tets.
                flip23(fliptets, 1, 0, 0);
                // Put the "outer" link faces into check list.
                //   fliptets[0] = [e,d,a,b] => will be flipped, so 
                //   [a,b,d] and [a,b,e] are not "outer" link faces.
                for (i = 1; i < 3; i++) {
                  eprevesym(fliptets[i], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                }
                for (i = 1; i < 3; i++) {
                  enextesym(fliptets[i], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                }
                // Then do a 3-to-2 flip.
                enextesymself(fliptets[0]);  // fliptets[0] is [e,d,a,b].
                eprevself(fliptets[0]); // = [b,a,d,c], d is the new vertex.
                fliptets[1] = baktets[0]; // = [b,a,e,f]
                fliptets[2] = baktets[1]; // = [b,a,f,d]
                flip32(fliptets, 1, 0, 0);
                // Put the "outer" link faces into check list.
                //   fliptets[0] = [d,e,f,a]
                //   fliptets[1] = [e,d,f,b]
                //   Faces [a,b,d] and [a,b,e] are not "outer" link faces.
                enextself(fliptets[0]);
                for (i = 1; i < 3; i++) {
                  esym(fliptets[0], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                  enextself(fliptets[0]);
                }
                enextself(fliptets[1]);
                for (i = 1; i < 3; i++) {
                  esym(fliptets[1], newface);
                  crossfaces->newindex((void **) &parytet);
                  *parytet = newface;
                  enextself(fliptets[1]);
                }
                flip23count--;
                flip32count--;
                flip44count++;
                flipflag = 1;
              } else {
                //n == 4, convflag != 0; assert(0);
              }
            } else { 
              // n > 4 => unflipable. //assert(0);
            }
          } else {
            // There are more than 1 non-convex or coplanar cases.
            flipflag = -1; // Ignore this face.
            if (b->verbose > 2) {
              printf("        Ignore face (%d, %d, %d) - %d, %d, tau = %.17g\n",
                     pointmark(bface.forg), pointmark(bface.fdest),
                     pointmark(bface.fapex), pointmark(bface.foppo),
                     pointmark(bface.noppo), bface.key);
            }
            dbg_ignore_facecount++;
          } // if (convcount == 1)

          if (flipflag == 1) {
            // Update the priority queue.
            for (i = 0; i < crossfaces->objects; i++) {
              parytet = (triface *) fastlookup(crossfaces, i);
              flipcertify(parytet, &pqueue);
            }
            crossfaces->restart();
            if (1) { // if (!b->flipinsert_random) {
              // Insert all queued unflipped faces.
              for (i = 0; i < bfacearray->objects; i++) {
                parytet = (triface *) fastlookup(bfacearray, i);
                // This face may be changed.
                if (!isdeadtet(*parytet)) {
                  flipcertify(parytet, &pqueue);
                }
              }
              bfacearray->restart();
            }
            fcount++;
          } else if (flipflag == 0) {
            // Queue an unflippable face. To process it later.
            bfacearray->newindex((void **) &parytet);
            *parytet = fliptet;
          }
        } // if (pe == bface.noppo)  
      } // if ((pa == bface.forg) && ...)
    } // if (bface.tt != NULL)

  } // while (pqueue != NULL)

  if (bfacearray->objects > 0) {
    if (fcount == 0) {
      printf("!! No flip is found in %ld faces.\n", bfacearray->objects);
      assert(0);
    }
  }

  // 'bfacearray' may be not empty (for what reason ??).
  dbg_unflip_facecount += bfacearray->objects;

  assert(flippool->items == 0l);
  delete bfacearray;

  // Un-mark top and bottom points.
  for (i = 0; i < toppoints->objects; i++) {
    parypt = (point *) fastlookup(toppoints, i);
    punmarktest2(*parypt);
  }
  for (i = 0; i < botpoints->objects; i++) {
    parypt = (point *) fastlookup(botpoints, i);
    punmarktest3(*parypt);
  }

  f23count = flip23count - f23count;
  f32count = flip32count - f32count;
  f44count = flip44count - f44count;
  totalfcount = f23count + f32count + f44count;

  if (totalfcount > maxflipsequence) {
    maxflipsequence = totalfcount;
  }

  if (b->verbose > 2) {
    printf("      Total %ld flips. f23(%ld), f32(%ld), f44(%ld).\n",
           totalfcount, f23count, f32count, f44count);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// fillregion()    Fill the missing region by a set of new subfaces.         //
//                                                                           //
// 'missingshs' contains the list of subfaces in R.  Moreover, each subface  //
// (except the first one) in this list represents an interior edge of R.     //
// Note: All subfaces in R are smarktested.                                  //
//                                                                           //
// Note: We assume that all vertices of R are marktested so we can detect    //
// new subface by checking the flag in apexes.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenmesh::fillregion(arraypool* missingshs, arraypool* missingshbds, 
                            arraypool* newshs)
{
  badface *newflipface, *popface;
  triface searchtet, spintet;
  face oldsh, newsh, opensh, *parysh;
  face casout, casin, neighsh, checksh;
  face checkseg, fakeseg;
  point pc, pd, pe, pf, ppt[2];
  enum interresult dir;
  REAL n[3], len; // elen[3];
  bool insideflag;
  int types[2], poss[4];
  int i, j, k;

  if (b->verbose > 2) {
    printf("      Fill region: %ld old subfaces (%ld).\n", missingshs->objects,
           fillregioncount);
  }

  // Search the first constrained face of R. It is found from the set of
  //   faces sharing at a boundary edge [a,b]. Such face must be found.
  //   The search takes the following two steps:
  //   - First, finds a candidate face [a,b,c] where c is also a vertex of R;
  //     Note that [a,b,c] may not be the right face to fill R. For instance,
  //     when R is concave at b.
  //   - Second, check if [a,b,c] can fill R. This can be checked if an
  //     adjacent tet of [a,b,c] intersects R. This is a tetrahedron-triangle
  //     intersection test. It can be reduced to two triangle-edge intersect
  //     tests, i.e., intersect the two faces not containing the edge [a,b] in
  //     this tet with all interior edges of R.

  // We start from the first boundary edge of R.
  oldsh = * (face *) fastlookup(missingshbds, 0);
  ppt[0] = sorg(oldsh);
  ppt[1] = sdest(oldsh);
  point2tetorg(ppt[0], searchtet);
  dir = finddirection(&searchtet, ppt[1]);
  assert(dir == ACROSSVERT); // SELF_CHECK

  insideflag = false;

  // Each face has two adjacent tets.
  for (k = 0; k < 2; k++) {
    if (b->verbose > 2) {
      printf("      Search an interior face from edge (%d, %d).\n",
             pointmark(ppt[0]), pointmark(ppt[1]));
    }
    spintet = searchtet;
    while (1) {
      pc = apex(spintet);
      if (pmarktested(pc)) {
        // Found a candidate face. Check if it is inside R.
        if (missingshs->objects > 2l) {
          // pd = oppo(spintet);
          // if (pd == dummypoint) {
            // Calculate an above point for this subface.
	    facenormal(ppt[0], ppt[1], pc, n, 1, NULL);
            len = sqrt(DOT(n, n));
            n[0] /= len;
            n[1] /= len;
            n[2] /= len;
            len = DIST(ppt[0], ppt[1]);
            len += DIST(ppt[1], pc);
            len += DIST(pc, ppt[0]);
            len /= 3.0;
            dummypoint[0] = ppt[0][0] + len * n[0];
            dummypoint[1] = ppt[0][1] + len * n[1];
            dummypoint[2] = ppt[0][2] + len * n[2];
            pd = dummypoint;
	  // }
          //if (pd != dummypoint) {
            for (j = 0; j < 2 && !insideflag; j++) {
              for (i = 1; i < missingshs->objects && !insideflag; i++) {
                parysh = (face *) fastlookup(missingshs, i);
                // Get an interior edge of R.
                pe = sorg(*parysh);
                pf = sdest(*parysh);
                if (tri_edge_test(ppt[j],pc,pd,pe,pf,NULL,1,types,poss)) {
                  dir = (enum interresult) types[0];
                  if (dir == ACROSSFACE) {
                    searchtet = spintet;
                    insideflag = true;
                  } else if (dir == ACROSSEDGE) {
                    searchtet = spintet;
                    insideflag = true;
                  }
                }
              } // i
            } // j
	  // }
          // if (pd == dummypoint) {
            dummypoint[0] = 0;
            dummypoint[1] = 0;
            dummypoint[2] = 0;
	  // }
        } else {
          // It is a simple 2-to-2 flip.
          searchtet = spintet;
          insideflag = true;
        }
      } // if (pmarktested(pc))
      if (insideflag) break;
      fnextself(spintet);
      if (spintet.tet == searchtet.tet) break;
    } // while (1)
    if (insideflag) break;
    esymself(searchtet);
    ppt[0] = org(searchtet);
    ppt[1] = dest(searchtet);
  } // k

  if (!insideflag) {
    // Something strange is happening.
    // Refine the missing region by adding a Steiner point.
    recentsh = oldsh;
    recenttet = searchtet; // For point location.
    return false;
  }

  // Create a new subface at the boundary edge.
  if (b->verbose > 2) {
    printf("      Create a new subface (%d, %d, %d)\n", pointmark(ppt[0]),
           pointmark(ppt[1]), pointmark(pc));
  }
  makeshellface(subfaces, &newsh);
  setsorg(newsh, ppt[0]);
  setsdest(newsh, ppt[1]);
  setsapex(newsh, pc);
  // The new subface gets its markers from the old one.
  setshellmark(newsh, shellmark(oldsh));
  if (checkconstraints) {
    setareabound(newsh, areabound(oldsh));
  }
  // Connect the new subface to adjacent tets.
  tspivot(searchtet, checksh); // SELF_CHECK
  assert(checksh.sh == NULL); // SELF_CHECK 
  tsbond(searchtet, newsh);
  fsymself(searchtet);
  sesymself(newsh);
  tsbond(searchtet, newsh);
  // Connect newsh to outer subfaces.
  sspivot(oldsh, checkseg);
  spivot(oldsh, casout);
  if (casout.sh != NULL) {
    casin = casout;
    if (checkseg.sh != NULL) {
      // Make sure that the subface has the right ori at the segment.
      checkseg.shver = 0;
      if (sorg(newsh) != sorg(checkseg)) {
        sesymself(newsh);
      }
      spivot(casin, neighsh);
      while (neighsh.sh != oldsh.sh) {
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
  // Add this new subface into list.
  sinfect(newsh);
  newshs->newindex((void **) &parysh);
  *parysh = newsh;

  // Push two "open" side of the new subface into stack.
  for (i = 0; i < 2; i++) {
    senextself(newsh);
    newflipface = (badface *) flippool->alloc();
    newflipface->ss = newsh;
    newflipface->nextitem = flipstack;
    flipstack = newflipface;
  }

  // Every other boundary edge of R is identified as a segment. Insert a faked
  //   segments at the place if it is not a segment.
  for (i = 1; i < missingshbds->objects; i++) {
    parysh = (face *) fastlookup(missingshbds, i);
    ppt[0] = sorg(*parysh);
    ppt[1] = sdest(*parysh);
    point2tetorg(ppt[0], searchtet);
    dir = finddirection(&searchtet, ppt[1]);
    assert(dir == ACROSSVERT); // SELF_CHECK
    tsspivot1(searchtet, checkseg);
    if (checkseg.sh == NULL) {
      // Insert a fake segment at this tet.
      if (b->verbose > 2) {
        printf("      Insert a fake segment (%d, %d)\n", pointmark(ppt[0]),
               pointmark(ppt[1]));
      }
      makeshellface(subsegs, &fakeseg);
      setsorg(fakeseg, ppt[0]);
      setsdest(fakeseg, ppt[1]);
      sinfect(fakeseg); // Mark it as faked.
      // Connect it to all tets at this edge.
      spintet = searchtet;
      while (1) {
        tssbond1(spintet, fakeseg);
        fnextself(spintet);
        if (spintet.tet == searchtet.tet) break;
      }
      checkseg = fakeseg;
    }
    // Let the segment hold the old subface.
    checkseg.shver = 0;
    sbond1(checkseg, *parysh);
    // Remember it to free it later.
    *parysh = checkseg;
  }

  // Loop until 'flipstack' is empty.
  while (flipstack != NULL) {

    // Pop an "open" side from the stack.
    popface = flipstack;
    opensh = popface->ss;
    flipstack = popface->nextitem; // The next top item in stack.
    flippool->dealloc((void *) popface);

    // Process it if it is still open.
    spivot(opensh, casout);
    if (casout.sh == NULL) {
      if (b->verbose > 2) {
        printf("      Get an open side (%d, %d) - %d.\n",
               pointmark(sorg(opensh)), pointmark(sdest(opensh)),
               pointmark(sapex(opensh)));
      }
      // Search a neighbor to close this side.
      stpivot(opensh, searchtet);
      tsspivot1(searchtet, checkseg);
      if (checkseg.sh == NULL) {
        // No segment. It is inside R. Search for a new face to fill in R.
        //   Note that the face may not be found (see fig 2010-05-25-c).
        spintet = searchtet;
        fnextself(spintet); // Skip the current face.
        while (1) {
          pc = apex(spintet);
          if (pmarktested(pc)) {
            // Found a place for a new subface inside R -- Case (i).
            tspivot(spintet, checksh);
            if (checksh.sh == NULL) {
              // Create a new subface.
              if (b->verbose > 2) {
                printf("      Create a new subface (%d, %d, %d)\n", 
                       pointmark(org(spintet)), pointmark(dest(spintet)),
                       pointmark(pc));
              }
              makeshellface(subfaces, &newsh);
              setsorg(newsh, org(spintet));
              setsdest(newsh, dest(spintet));
              setsapex(newsh, pc);
              // The new subface gets its markers from its neighbor.
              setshellmark(newsh, shellmark(opensh));
              if (checkconstraints) {
                setareabound(newsh, areabound(opensh));
              }
              // Connect the new subface to adjacent tets.
              tsbond(spintet, newsh);
              fsymself(spintet);
              sesymself(newsh);
              tsbond(spintet, newsh);
              // Connect newsh to its adjacent subface.
              sbond(newsh, opensh);
              // Add this new subface into list.
              sinfect(newsh);
              newshs->newindex((void **) &parysh);
              *parysh = newsh;
              // Push two "open" side of the new subface into stack.
              for (i = 0; i < 2; i++) {
                senextself(newsh);
                newflipface = (badface *) flippool->alloc();
                newflipface->ss = newsh;
                newflipface->nextitem = flipstack;
                flipstack = newflipface;
              }
            } else {
              // A new subface has already been created.
              assert(sinfected(checksh)); // It must be in stack.
              spivot(checksh, neighsh); // SELF_CHECK
              assert(neighsh.sh == NULL); // Its side must be open.
              if (b->verbose > 2) {
                printf("      Connect to another open side (%d, %d, %d)\n", 
                       pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                       pointmark(sapex(checksh)));
              }
              sbond(opensh, checksh); // Simply connect them.
            }
            break; // -- Case (i)
          }
          fnextself(spintet);
          if (spintet.tet == searchtet.tet) {
            // Not find any face to fill in R at this side.
            // TO DO: suggest a point to split the edge.
            assert(0);
          }
        } // while (1)
      } else {
        // This side coincident with a boundary edge of R.
        checkseg.shver = 0;
        spivot(checkseg, oldsh);
        if (sinfected(checkseg)) {
          // It's a faked segment. Delete it.
          if (b->verbose > 2) {
            printf("      Delete a fake segment (%d, %d)\n", 
                   pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
          }
          spintet = searchtet;
          while (1) {
            tssdissolve1(spintet);
            fnextself(spintet);
            if (spintet.tet == searchtet.tet) break;
          }
          shellfacedealloc(subsegs, checkseg.sh);
        }
        if (b->verbose > 2) {
          printf("      Connect to a boundary edge (%d, %d, %d)\n", 
                 pointmark(sorg(oldsh)), pointmark(sdest(oldsh)),
                 pointmark(sapex(oldsh)));
        }
        sspivot(oldsh, checkseg);
        spivot(oldsh, casout);
        if (casout.sh != NULL) {
          casin = casout;
          if (checkseg.sh != NULL) {
            // Make sure that the subface has the right ori at the segment.
            checkseg.shver = 0;
            if (sorg(opensh) != sorg(checkseg)) {
              sesymself(opensh);
	    }
            spivot(casin, neighsh);
            while (neighsh.sh != oldsh.sh) {
              casin = neighsh;
              spivot(casin, neighsh);
            }
          }
          sbond1(opensh, casout);
          sbond1(casin, opensh);
        }
        if (checkseg.sh != NULL) {
          ssbond(opensh, checkseg);
        }
      }

    } // if (casout.sh == NULL)

  } // while (flipstack != NULL)

  // Uninfect all new subfaces.
  for (i = 0; i < newshs->objects; i++) {
    parysh = (face *) fastlookup(newshs, i);
    suninfect(*parysh);
  }

  if (b->verbose > 2) {
    printf("      Created %ld new subfaces.\n", newshs->objects);
  }
  fillregioncount++;

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// refineregion()    Refine a missing region by inserting points.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::refineregion()
{
  triface searchtet;
  face splitsh;
  face *paryseg, sseg;
  point steinpt, pa, pb, pc;
  insertvertexflags ivf;
  REAL auv[2], buv[2], newuv[2], t;
  int fmark, fid, eid;
  int loc; // iloc, sloc;
  int s, i;

  // The mesh is a CDT.
  assert(subsegstack->objects == 0l); // SELF_CHECK

  // Create a new point.
  makepoint(&steinpt, FREEFACETVERTEX);

  // The 'recentsh' saved an edge to be split.
  splitsh = recentsh;
  // Add the Steiner point at the barycenter of the face.
  pa = sorg(splitsh);
  pb = sdest(splitsh);
  pc = sapex(splitsh);

  if (b->psc) {
    assert(in->facetmarkerlist != NULL);
    fmark = shellmark(splitsh) - 1;
    fid = in->facetmarkerlist[fmark];
    if (pointtype(pa) == RIDGEVERTEX) {
      in->getvertexparamonface(in->geomhandle, pointmark(pa), fid, auv);
    } else if (pointtype(pa) == FREESEGVERTEX) {
      eid = pointgeomtag(pa); // The Edge containing this Steiner point.
      t = pointgeomuv(pa, 0);  // The Steiner point's parameter on Edge.
      in->getedgesteinerparamonface(in->geomhandle, eid, t, fid, auv);
    } else if (pointtype(pa) == FREEFACETVERTEX) {
      auv[0] = pointgeomuv(pa, 0);
      auv[1] = pointgeomuv(pa, 1);
    } else {
      assert(0);
    }
    if (pointtype(pb) == RIDGEVERTEX) {
      in->getvertexparamonface(in->geomhandle, pointmark(pb), fid, buv);
    } else if (pointtype(pb) == FREESEGVERTEX) {
      eid = pointgeomtag(pb); // The Edge containing this Steiner point.
      t = pointgeomuv(pb, 0);  // The Steiner point's parameter on Edge.
      in->getedgesteinerparamonface(in->geomhandle, eid, t, fid, buv);
    } else if (pointtype(pb) == FREEFACETVERTEX) {
      buv[0] = pointgeomuv(pb, 0);
      buv[1] = pointgeomuv(pb, 1);
    } else {
      assert(0);
    }
    newuv[0] = 0.5 * (auv[0] + buv[0]);
    newuv[1] = 0.5 * (auv[1] + buv[1]);
    in->getsteineronface(in->geomhandle, fid, newuv, steinpt);
    setpointgeomuv(steinpt, 0, newuv[0]);
    setpointgeomuv(steinpt, 1, newuv[1]);
    setpointgeomtag(steinpt, fid);
  } else {
    for (i = 0; i < 3; i++) {
      steinpt[i] = (pa[i] + pb[i] + pc[i]) / 3.0;
    }
  }

  // Start searching it from 'recentet'.
  searchtet = recenttet;
  // Now insert the point p. The flags are chosen as follows: 
  //   - boywat  = 2, the current T is a CDT, 
  //   - lawson  = 2, do flip after inserting p, some existing segments
  //                  and subfaces may be flipped, they are queued and
  //                  and will be recovered.
  //   - rejflag = 1, reject p if it encroaches upon at least one segment,
  //                  queue encroached segments.
  ivf.iloc = (int) OUTSIDE;
  ivf.bowywat = 2;
  ivf.lawson = 2;
  ivf.rejflag = 1;
  ivf.chkencflag = 0;
  ivf.sloc = (int) ONFACE;
  ivf.sbowywat = 2;
  ivf.splitbdflag = 0;
  ivf.validflag = 1;
  ivf.respectbdflag = 0;
  ivf.assignmeshsize = b->metric;
  loc = insertvertex(steinpt, &searchtet, &splitsh, NULL, &ivf);

  assert((loc != OUTSIDE) && (loc != ONVERTEX));
  if (loc == NEARVERTEX) {
    // The new point is either ON or VERY CLOSE to an existing point.
    pa = point2ppt(steinpt);
    printf("  !! Avoid to create a short edge (length = %g)\n",
           distance(steinpt, pa));
    // Indicate it may be an input problem.
    printf("  Short edge length bound is: %g. Tolerance is %g.\n", 
           b->minedgelength, b->epsilon);
    terminatetetgen(4);
  }

  if (loc == ENCSEGMENT) {
    // Some segments are encroached and queued.
    assert(encseglist->objects > 0l);
    // Randomly pick one encroached segment to split.
    s = randomnation(encseglist->objects);
    paryseg = (face *) fastlookup(encseglist, s);
    sseg = *paryseg;
    // The new point p is the midpoint of this segment.
    getsteinerptonsegment(&sseg, NULL, steinpt);
    setpointtype(steinpt, FREESEGVERTEX);
    encseglist->restart();  // Clear the queue.

    // Start searching from an adjacent tetrahedron (containing the segment).
    sstpivot1(sseg, searchtet);
    spivot(sseg, splitsh);
    // Insert the point p. The flags are chosen as follows: 
    //   - boywat  = 2, the current T is a CDT, 
    //   - lawson  = 2, do flip after inserting p, some existing segments
    //                  and subfaces may be flipped, they are queued and
    //                  and will be recovered.
    //   - rejflag = 0, always insert p, even it will cause some segments
    //                  or subfaces missing, queue missing boundaries.
    ivf.iloc = (int) ONEDGE;
    ivf.bowywat = 2;
    ivf.lawson = 2;
    ivf.rejflag = 0;
    ivf.chkencflag = 0;
    ivf.sloc = (int) ONEDGE;
    ivf.sbowywat = 2;
    ivf.splitbdflag = 0;
    ivf.validflag = 1;
    ivf.respectbdflag = 0;
    ivf.assignmeshsize = b->metric;
    loc = insertvertex(steinpt, &searchtet, &splitsh, &sseg, &ivf);

    if (loc == NEARVERTEX) {
      // The new point is either ON or VERY CLOSE to an existing point.
      pa = point2ppt(steinpt);
      printf("  !! Avoid to create a short edge (length = %g)\n",
             distance(steinpt, pa));
      // Indicate it may be an input problem.
      printf("  Short edge length bound is: %g. Tolerance is %g.\n", 
             b->minedgelength, b->epsilon);
      terminatetetgen(4);
    }

    st_segref_count++;
  } else {
    st_facref_count++;
  }
  if (steinerleft > 0) steinerleft--;

  // Do flip to recover Delaunayniess.
  lawsonflip3d(steinpt, 2, 0, 0, 0);

  // Some vertices may be queued, recover them.
  if (subvertstack->objects > 0l) {
    assert(0); //delaunizevertices();
  }

  // Some subsegments may be queued, recover them.
  if (subsegstack->objects > 0l) {
    delaunizesegments();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// constrainedfacets()    Recover subfaces saved in 'subfacestack'.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::constrainedfacets()
{
  arraypool *tg_crosstets, *tg_topnewtets, *tg_botnewtets;
  arraypool *tg_topfaces, *tg_botfaces, *tg_midfaces;
  arraypool *tg_topshells, *tg_botshells, *tg_facfaces; 
  arraypool *tg_toppoints, *tg_botpoints;
  arraypool *tg_missingshs, *tg_missingshbds, *tg_missingshverts;

  triface searchtet, neightet;
  face searchsh, neighsh, *parysh;
  face checkseg, *paryseg;
  point refpt, *parypt;
  enum interresult dir;
  bool success;
  int facetcount;
  //int bakhullsize;
  int crossflag;
  int i, j;

  // Initialize arrays.
  tg_crosstets      = new arraypool(sizeof(triface), 10);
  tg_topnewtets     = new arraypool(sizeof(triface), 10);
  tg_botnewtets     = new arraypool(sizeof(triface), 10);
  tg_topfaces       = new arraypool(sizeof(triface), 10);
  tg_botfaces       = new arraypool(sizeof(triface), 10);
  tg_midfaces       = new arraypool(sizeof(triface), 10);
  tg_toppoints      = new arraypool(sizeof(point), 8);
  tg_botpoints      = new arraypool(sizeof(point), 8);
  tg_facfaces       = new arraypool(sizeof(face), 10);
  tg_topshells      = new arraypool(sizeof(face), 10);
  tg_botshells      = new arraypool(sizeof(face), 10);
  tg_missingshs     = new arraypool(sizeof(face), 10);
  tg_missingshbds   = new arraypool(sizeof(face), 10);
  tg_missingshverts = new arraypool(sizeof(point), 8);

  // This is a global array used by refineregion().
  encseglist        = new arraypool(sizeof(face), 4); 

  facetcount = 0;

  // Loop until 'subfacstack' is empty.
  while (subfacstack->objects > 0l) {
    subfacstack->objects--;
    parysh = (face *) fastlookup(subfacstack, subfacstack->objects);
    searchsh = *parysh;

    if (searchsh.sh[3] == NULL) continue; // Skip a dead subface.

    stpivot(searchsh, neightet);
    if (neightet.tet == NULL) {
      // Find an unrecovered subface.
      smarktest(searchsh);
      tg_facfaces->newindex((void **) &parysh);
      *parysh = searchsh;
      // Collect all non-recovered subfaces of the same facet.
      for (i = 0; i < tg_facfaces->objects; i++) {
        searchsh = * (face *) fastlookup(tg_facfaces, i);
        for (j = 0; j < 3; j++) {
          sspivot(searchsh, checkseg);
          if (checkseg.sh == NULL) {
            spivot(searchsh, neighsh);
            assert(neighsh.sh != NULL); // SELF_CHECK
            if (!smarktested(neighsh)) {
              // It may be already recovered.
              stpivot(neighsh, neightet);
              if (neightet.tet == NULL) {
                smarktest(neighsh);
                tg_facfaces->newindex((void **) &parysh);
                *parysh = neighsh;
              }
            }
          }
          senextself(searchsh);
        } // j
      } // i
      // Have found all facet subfaces (vertices). Uninfect them.
      for (i = 0; i < tg_facfaces->objects; i++) {
        parysh = (face *) fastlookup(tg_facfaces, i);
        sunmarktest(*parysh);
      }

      if (b->verbose > 2) {
        printf("    Recover facet #%d: %ld subfaces.\n", facetcount + 1, 
               tg_facfaces->objects);
      }
      facetcount++;

      // Loop until 'tg_facfaces' is empty.
      while (tg_facfaces->objects > 0l) {
        // Get the last subface of this array.
        tg_facfaces->objects--;
        parysh = (face *) fastlookup(tg_facfaces, tg_facfaces->objects);
        searchsh = *parysh;

        if (searchsh.sh[3] == NULL) continue; // Skip a dead subface.

        stpivot(searchsh, neightet);
        if (neightet.tet != NULL) continue; // Not a missing subface.

        // Insert the subface.
        searchtet.tet = NULL;
        dir = scoutsubface(&searchsh, &searchtet);
        if (dir == SHAREFACE) continue; // The subface is inserted.
        if (dir == COLLISIONFACE) continue; // The subface is removed.

        // The subface is missing. Form the missing region.
        //   Re-use 'tg_crosstets' for 'adjtets'.
        formmissingregion(&searchsh, tg_missingshs, tg_missingshbds,
                          tg_missingshverts, tg_crosstets);

        // Search for a crossing edge (tg_crosstets is cleared).
        crossflag = scoutcrossedge(searchtet, tg_crosstets, tg_missingshs);

        if (crossflag == 1) {
          // Recover subfaces by local retetrahedralization.
          // Form a cavity of crossing tets.
          if (formcavity(&searchtet, tg_missingshs, tg_crosstets, tg_topfaces, 
                         tg_botfaces, tg_toppoints, tg_botpoints)) {
            if (!b->flipinsert) {
              // Tetrahedralize the top part. Re-use 'tg_midfaces'.
              delaunizecavity(tg_toppoints, tg_topfaces, tg_topshells,
                              tg_topnewtets, tg_crosstets, tg_midfaces);
              // Tetrahedralize the bottom part. Re-use 'tg_midfaces'.
              delaunizecavity(tg_botpoints, tg_botfaces, tg_botshells,
                              tg_botnewtets, tg_crosstets, tg_midfaces);
              // Fill the cavity with new tets.
              success = fillcavity(tg_topshells, tg_botshells, tg_midfaces,
                                   tg_missingshs);
              if (success) {
                // Cavity is remeshed. Delete old tets and outer new tets.
                carvecavity(tg_crosstets, tg_topnewtets, tg_botnewtets);
                // Insert the missing region into cavity.
                j = 0; // FOR DEBUG! Count the number of non-recovered faces. 
                for (i = 0; i < tg_missingshs->objects; i++) {
                  searchsh = * (face *) fastlookup(tg_missingshs, i);
                  searchtet.tet = NULL;
                  dir = scoutsubface(&searchsh, &searchtet);
                  assert(dir != COLLISIONFACE); // SELF_CHECK
                  if (dir != SHAREFACE) {
                    // A subface is missing. This is possible that the subface
                    //   is not actually a constrained Delaunay face in T. 
                    // Add this face at the end of the list, so it will be
                    //   processed immediately. This is necessary because we
                    //   have created some non-locally Delaunay face (by the
                    //   remesh of the cavity). We have to insert the subfaces
                    //   to make these face constrained Delaunay.
                    tg_facfaces->newindex((void **) &parysh);
                    *parysh = searchsh;
                    j++; // FOR DEBUG!
                  }
                } // i
                // Recover interior subfaces.
                for (i = 0; i < caveencshlist->objects; i++) {
                  searchsh = * (face *) fastlookup(caveencshlist, i);
                  searchtet.tet = NULL;
                  dir = scoutsubface(&searchsh, &searchtet);
                  assert(dir != COLLISIONFACE); // SELF_CHECK
                  if (dir != SHAREFACE) {
                    // The subface is missing. This is possible that the subface
                    //   is removed by the enlargement of the cavity. It has to
                    //   be recovered. 
                    // Add this face at the end of the list, so it will be
                    //   processed immediately. We have to insert the subfaces
                    //   to make these face constrained Delaunay.
                    tg_facfaces->newindex((void **) &parysh);
                    *parysh = searchsh;
                    j++; // FOR DEBUG!
                  }
                } // i
                // Recover interior segments. This should always be recovered.
                for (i = 0; i < caveencseglist->objects; i++) {
                  paryseg = (face *) fastlookup(caveencseglist, i);
                  searchtet.tet = NULL;
                  refpt = NULL;
                  dir = scoutsegment(sorg(*paryseg),sdest(*paryseg),&searchtet,
                                     &refpt, NULL);
                  assert(dir == SHAREEDGE);
                  // Insert this segment.
                  tsspivot1(searchtet, checkseg);  // SELF_CHECK
                  if (checkseg.sh == NULL) {
                    // Let the segment remember an adjacent tet.
                    sstbond1(*paryseg, searchtet);
                    // Bond the segment to all tets containing it.
                    neightet = searchtet;
                    do {
                      tssbond1(neightet, *paryseg);
                      fnextself(neightet);
                    } while (neightet.tet != searchtet.tet);
                  } else {
                    // Collision! Should not happen.
                    assert(0);
                  }
                } // i
                caveencshlist->restart();
                caveencseglist->restart();
              } else {
                // Restore old tets and delete new tets.
                restorecavity(tg_crosstets, tg_topnewtets, tg_botnewtets);
                // Set a handle for searching subface.
                //recentsh = searchsh;
              }
            } else {
              // Use the flip algorithm of Shewchuk to recover the subfaces.
              flipinsertfacet(tg_crosstets, tg_toppoints, tg_botpoints, 
                              tg_missingshverts);
              // Check the missing subfaces again.
              j = 0; // FOR DEBUG! Count the number of non-recovered faces. 
              for (i = 0; i < tg_missingshs->objects; i++) {
                searchsh = * (face *) fastlookup(tg_missingshs, i);
                searchtet.tet = NULL;
                dir = scoutsubface(&searchsh, &searchtet);
                assert(dir != COLLISIONFACE); // SELF_CHECK
                if (dir != SHAREFACE) {
                  // A subface is missing. This is possible that the subface
                  //   is not actually a constrained Delaunay face in T. 
                  // Add this face at the end of the list, so it will be
                  //   processed immediately. This is necessary because we
                  //   have created some non-locally Delaunay face (by the
                  //   remesh of the cavity). We have to insert the subfaces
                  //   to make these face constrained Delaunay.
                  tg_facfaces->newindex((void **) &parysh);
                  *parysh = searchsh;
                  j++; // FOR DEBUG!
                }
              } // i
              // Clear working lists.
              tg_crosstets->restart();
              tg_topfaces->restart();
              tg_botfaces->restart();
              tg_toppoints->restart();
              tg_botpoints->restart();
              success = true;
            } // if (b->flipinsert)
          } else {
            // Formcavity failed.
            success = false;
          }
        } else { //if (crossflag == 0) {
          // Recover subfaces by retriangulate the surface mesh.
          //   Re-use tg_topshells for newshs.
          success = fillregion(tg_missingshs, tg_missingshbds, tg_topshells);
          if (success) {
            // Region is remeshed. Delete old subfaces (in tg_missingshs).
            for (i = 0; i < tg_missingshs->objects; i++) {
              parysh = (face *) fastlookup(tg_missingshs, i);
              shellfacedealloc(subfaces, parysh->sh);
            }
            tg_topshells->restart();
          } else {
            // Search a handle for searching tetrahedron.
            recenttet = searchtet;
          }
        }

        // Unmarktest all points of the missing region.
        for (i = 0; i < tg_missingshverts->objects; i++) {
          parypt = (point *) fastlookup(tg_missingshverts, i);
          punmarktest(*parypt);
        }
        tg_missingshverts->restart();
        tg_missingshbds->restart();
        tg_missingshs->restart();

        if (!success) {
          // The missing region can not be recovered. Refine it.
          refineregion();
          // Clean the current list of facet subfaces.
          //tg_facfaces->restart();
        }
      } // while (tg_facfaces->objects > 0l)

    } // if (neightet.tet == NULL)
  } // while (subfacstack->objects > 0l)

  // Delete arrays.
  delete tg_crosstets;
  delete tg_topnewtets;
  delete tg_botnewtets;
  delete tg_topfaces;
  delete tg_botfaces;
  delete tg_midfaces;
  delete tg_toppoints;
  delete tg_botpoints;
  delete tg_facfaces;
  delete tg_topshells;
  delete tg_botshells;
  delete tg_missingshs;
  delete tg_missingshbds;
  delete tg_missingshverts;
  delete encseglist;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// constraineddelaunay()    Create a constrained Delaunay tetrahedralization.//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::constraineddelaunay(clock_t& tv)
{
  face searchsh, *parysh;
  face searchseg, *paryseg;
  int s, i;

  // Statistics.
  long bakfillregioncount;
  long bakcavitycount, bakcavityexpcount;

  if (!b->quiet) {
    printf("Constrained Delaunay...\n");
  }

  // Identify acute vertex for PLC inputs.
  markacutevertices();

  if (b->verbose) {
    printf("  Delaunizing segments.\n");
  }

  checksubsegflag = 1;

  // Put all segments into the list.
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

  // Recover non-Delaunay segments.
  delaunizesegments();

  if (b->verbose) {
    printf("  %ld Steiner points.\n", st_segref_count); 
  }

  tv = clock();

  if (b->verbose) {
    printf("  Constraining facets.\n");
  }

  if (b->flipinsert) {
    // Clear the counters.
    flip23count = flip32count = flip44count = 0l;
  }

  // Subfaces will be introduced.
  checksubfaceflag = 1;

  bakfillregioncount = fillregioncount;
  bakcavitycount = cavitycount;
  bakcavityexpcount = cavityexpcount;

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

  // Recover facets.
  constrainedfacets();

  if (b->verbose) {
    if (fillregioncount > bakfillregioncount) {
      printf("  Remeshed %ld regions.\n", fillregioncount-bakfillregioncount);
    }
    if (cavitycount > bakcavitycount) {
      printf("  Remeshed %ld cavities", cavitycount - bakcavitycount);
      if (cavityexpcount - bakcavityexpcount) {
        printf(" (%ld enlarged)", cavityexpcount - bakcavityexpcount);
      }
      printf(".\n");
    }
    if (st_segref_count + st_facref_count > 0) {
      printf("  Inserted %ld (%ld, %ld) refine points.\n", 
             st_segref_count + st_facref_count, st_segref_count,
             st_facref_count);
    }
  }
}

////                                                                       ////
////                                                                       ////
//// constrained_cxx //////////////////////////////////////////////////////////

