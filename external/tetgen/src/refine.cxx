#include "../tetgen.h"
//// refine_cxx ///////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// marksharpsegments()    Mark sharp segments.                               //
//                                                                           //
// A segment is SHARP if there are two facets intersecting at it with an     //
// internal dihedral angle (*) less than an angle \theta.                    //
//                                                                           //
// A theoretical value of \theta is arccos(1/3) \approx 70.54 degree.  It is //
// possible to relax it in practice. Here we choose \theta = 65 degree.      //
//                                                                           //
// The minimum dihedral angle between facets (minfacetdihed) is calulcated.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::marksharpsegments()
{
  triface adjtet;
  face startsh, spinsh, neighsh;
  face segloop, nextseg, prevseg;
  point eorg, edest;
  REAL ang, smallang;
  bool issharp;
  int sharpcount;

  // For storing extremely small dihedral angle.
  face *parysh, *parysh1;
  REAL exsmallang;
  int exsharpcount;
  int i, j, k;

  if (b->verbose > 0) {
    printf("  Marking sharp segments.\n");
  }

  minfacetdihed = PI;
  smallang = 65.0 * PI / 180.0; // 65 degree.
  exsmallang = 5.0 * PI / 180.0; // 5 degree.
  sharpcount = exsharpcount = 0;

  // A segment s may have been split into many subsegments. Operate the one
  //   which contains the origin of s. Then mark the rest of subsegments.
  subsegs->traversalinit();
  segloop.sh = shellfacetraverse(subsegs);
  while (segloop.sh != (shellface *) NULL) {
    segloop.shver = 0;
    senext2(segloop, prevseg);
    spivotself(prevseg);
    if (prevseg.sh == NULL) {
      // Operate on this seg s.
      issharp = false;
      spivot(segloop, startsh);
      if (startsh.sh != NULL) {
        // First check if two facets form an acute dihedral angle at s.
        eorg = sorg(segloop);
        edest = sdest(segloop);
        spinsh = startsh;
        while (1) {
          if (sorg(spinsh) != eorg) sesymself(spinsh);
          // Only do test when the spinsh is faceing inward.
          stpivot(spinsh, adjtet);
          if (adjtet.tet != NULL) {
            if (!ishulltet(adjtet)) {
              // Get the subface on the adjacent facet.
              spivot(spinsh, neighsh);
              // Do not calculate if it is self-bonded.
              if ((neighsh.sh != NULL) && (neighsh.sh != spinsh.sh)) {
                // Calculate the dihedral angle between the two subfaces.
                ang = facedihedral(eorg, edest, sapex(spinsh), sapex(neighsh));
                // Only do check if a sharp angle has not been found.
                if (!issharp) issharp = (ang < smallang);
                // Remember the smallest facet dihedral angle.
                minfacetdihed = minfacetdihed < ang ? minfacetdihed : ang;
                if (ang < exsmallang) {
                  // It's an extremely small dihedral angle.
                  // Mark the two facets. 
                  // To avoid too many Steiner points, do not refine them.
                  if (shelltype(spinsh) != SHARP) {
                    setshelltype(spinsh, SHARP);
                    cavesegshlist->newindex((void **) &parysh);
                    *parysh = spinsh;
                  }
                  if (shelltype(neighsh) != SHARP) {                 
                    setshelltype(neighsh, SHARP);
                    cavesegshlist->newindex((void **) &parysh);
                    *parysh = neighsh;
                  }
                  exsharpcount++;
                }
              }
            }
          }
          // Go to the next facet.
          spivotself(spinsh);
          if (spinsh.sh == NULL) break; // A single subface case.
          if (spinsh.sh == startsh.sh) break;
        }
      } // if (startsh.sh != NULL)
      if (issharp) {
        if (b->verbose > 2) {
          printf("      Mark a sharp segment (%d, %d).\n",
                 pointmark(eorg), pointmark(edest));
        }
        setshelltype(segloop, SHARP);
        // The endpoint of this segment is acute.
        if (pointtype(eorg) == RIDGEVERTEX) {
          setpointtype(eorg, ACUTEVERTEX);
        } else {
          assert(pointtype(eorg) == ACUTEVERTEX); // SELF_CHECK
        }
        // Set the type for all subsegments at forwards.
        edest = sdest(segloop);
        senext(segloop, nextseg);
        spivotself(nextseg);
        while (nextseg.sh != NULL) {
          setshelltype(nextseg, SHARP);
          // Adjust the direction of nextseg.
          nextseg.shver = 0;
          if (sorg(nextseg) != edest) {
            sesymself(nextseg);
          }
          assert(sorg(nextseg) == edest);
          edest = sdest(nextseg);
          // Go the next connected subsegment at edest.
          senextself(nextseg);
          spivotself(nextseg);
        }
        // The endpoint of this segment is acute.
        if (pointtype(edest) == RIDGEVERTEX) {
          setpointtype(edest, ACUTEVERTEX);
        } else {
          assert(pointtype(edest) == ACUTEVERTEX); // SELF_CHECK
        }
        sharpcount++;
      } // if (issharp)
    } // if (prevseg.sh == NULL)
    segloop.sh = shellfacetraverse(subsegs);
  }

  // Mark all facets at extremely small dihedral angles.
  if (cavesegshlist->objects > 0) {
    for (i = 0; i < cavesegshlist->objects; i++) {
      parysh = (face *) fastlookup(cavesegshlist, i);
      caveshlist->newindex((void **) &parysh1);
      *parysh1 = *parysh;
      for (j = 0; j < caveshlist->objects; j++) {
        parysh1 = (face *) fastlookup(caveshlist, j);
        spinsh = *parysh1;
        for (k = 0; k < 3; k++) {
          sspivot(spinsh, nextseg);
          if (nextseg.sh == NULL) {
            spivot(spinsh, neighsh);
            if (shelltype(neighsh) != SHARP) {                 
              setshelltype(neighsh, SHARP);
              caveshlist->newindex((void **) &parysh1);
              *parysh1 = neighsh;
            }
          }
          senextself(spinsh);
        } // k
      } // j
      caveshlist->restart();
    } // i
    cavesegshlist->restart();
  } // if (cavesegshlist->objects > 0)

  if (b->verbose) {
    if (sharpcount > 0) {
      printf("  Found %d (%d) sharp segments.\n", sharpcount, exsharpcount);
    }
    printf("  Minimum fac-fac angle = %g.\n", minfacetdihed / PI * 180.0);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// decidefeaturepointsizes()    Calculate sizes for all feature points.      //
//                                                                           //
// A feature point is either an acute vertex or a Steiner point on a sharp   //
// segment.  Each feature point p will be protected by a ball whose radius   //
// is called its "feature size".                                             //
//                                                                           //
// NOTE: we should have already marked all features points in the two func-  //
// tions: markacutevertices() and marksharpsegments().  Each feature point   //
// has the type ACUTEVERTEX or FREESEGVERTEX.                                //
//                                                                           //
// The feature size of a vertex is the minimum of the following sizes:       //
//   (0) the (approximated) local feature size (the distance to the second   //
//       nearest boundary) of the vertex;
//   (1) the value specified in .mtr file (-m option);                       //
//   (2) the cubic root of a fixed maximal volume constraint ('-a__');       //
//   (3) the cubic root of a maximal volume constraint in a region ('-a');   //
//   (4) the square root of a maximal area constraint in a .var file;        //
//   (5) a maximal length constraint in a .var file;                         //
//                                                                           //
// If 'b->nobisect' ('-Y' option) is set, every input vertex has a size. It  //
// is used to prevent creating too close Steiner points.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::decidefeaturepointsizes()
{
  arraypool *tetlist, *verlist;
  triface starttet, *parytet;
  face checksh, parentsh, shloop;
  face checkseg, prevseg, nextseg, testseg;
  point ploop, adjpt, e1, e2, *parypt;
  REAL lfs_0, lfs_1, lfs_2;
  REAL len, vol, maxlen = 0.0, varlen;
  REAL ang, a, a1, a2, a3, prjpt[3], n[3];
  int featureflag, featurecount;
  int i, j;

  if (b->verbose > 0) {
    printf("  Deciding feature-point sizes.\n");
  }

  // Initialize working lists.
  tetlist = cavetetlist;
  verlist = cavetetvertlist;

  if (b->fixedvolume) {
    // A fixed volume constraint is imposed. This gives an upper bound of
    //   the maximal radius of the protect ball of a vertex.
    maxlen = pow(6.0 * b->maxvolume, 1.0 / 3.0);
  }

  // First, assign a size of p if p is a feature point or an input point and 
  //   the -Y option is used.
  featurecount = 0;
  points->traversalinit();
  ploop = pointtraverse();
  while (ploop != (point) NULL) {
    // Check if it is a feature point.
    featureflag = 0;
    // Only calculate the size if it has a size zero.
    // The point may already has a positive size (-m option).
    if (ploop[pointmtrindex] == 0) {
      if (pointtype(ploop) == ACUTEVERTEX) {
        featureflag = 1;
      } else {
        if (b->nobisect) { // '-Y' option
          if ((pointtype(ploop) == RIDGEVERTEX) ||
              (pointtype(ploop) == FACETVERTEX) ||
              (pointtype(ploop) == VOLVERTEX)) {
            featureflag = 1;  // It is an input vertex.
          }
        }
      }
    }
    if (featureflag) {
      // Form star(p).
      getvertexstar(1, ploop, tetlist, verlist, NULL);
      // Calculate lfs_0(p), i.e., the smallest distance from p to a vertex.
      // We approximate it by taking the distance of p to its nearest
      //   vertex in Link(p).
      lfs_0 = longest;
      for (i = 0; i < verlist->objects; i++) {
        parypt = (point *) fastlookup(verlist, i);
        adjpt = * parypt;
        if (adjpt == dummypoint) {
          continue; // Skip a dummypoint.
        }
        if (pointtype(adjpt) == FREESEGVERTEX) {
          // A Steiner point. Get the subsegment.
          sdecode(point2sh(adjpt), checkseg);
          assert(checkseg.sh != NULL);
          checkseg.shver = 0;
          if (sdest(checkseg) != adjpt) {
            sesymself(checkseg);
          }
          assert(sdest(checkseg) == adjpt);
          // It is possible that the original segment of 'adjpt' does not
          //   have 'ploop' as an endpoint.
          if (sorg(checkseg) == ploop) {
            // Find the other end point of the original segment.
            nextseg = checkseg;
            while (1) {
              senext(nextseg, testseg);
              spivotself(testseg);
              if (testseg.sh == NULL) break;
              // Go to the next subseg.
              nextseg = testseg;
              // Adjust the direction of the nextseg.
              nextseg.shver = 0;
              if (sorg(nextseg) != adjpt) {
                sesymself(nextseg);
              }
              assert(sorg(nextseg) == adjpt);
              adjpt = sdest(nextseg);
            }
          }
	} else if (pointtype(adjpt) == FREEFACETVERTEX) {
          // Ignore a Steiner point on facet.
          continue;
        } else if (pointtype(adjpt) == FREEVOLVERTEX) {
          // Ignore a Steiner point in volume.
          continue;
        }  
        len = distance(ploop, adjpt);
        if (lfs_0 > len) lfs_0 = len;
      } // i
      assert(lfs_0 < longest); // SELF_CHECK
      ploop[pointmtrindex] = lfs_0;
      // Calculate lfs_1(p), i.e., the smallest distance from p to a segment.
      //   We approximate it by restricting the segments in Link(p).
      lfs_1 = lfs_0;
      for (i = 0; i < tetlist->objects; i++) {
        parytet = (triface *) fastlookup(tetlist, i);
        for (j = 0; j < 3; j++) {
          tsspivot1(*parytet, checkseg);
          if (checkseg.sh != NULL) {
            e1 = sorg(checkseg);
            e2 = sdest(checkseg);
            // Only do calculation if the projeciton of 'p' lies inside the
            //   segment [e1, e2].
            ang = interiorangle(ploop, e1, e2, NULL);
            ang *= 2.0;
            if (ang > PI) { 
              len = shortdistance(ploop, e1, e2);
              if (lfs_1 > len) {
                lfs_1 = len;
              }
            }
          }
          enextself(*parytet);
        } // j
      } // i
      if (ploop[pointmtrindex] > lfs_1) {
        ploop[pointmtrindex] = lfs_1;
      }
      // Calculate lfs_2(p), i.e., the smallest distance from p to a facet.
      //   We approximate it by restricting the facets in Link(p).
      lfs_2 = lfs_0; 
      for (i = 0; i < tetlist->objects; i++) {
        parytet = (triface *) fastlookup(tetlist, i);
        tspivot(*parytet, checksh);
        if (checksh.sh != NULL) {
          adjpt = sorg(checksh);
          e1 = sdest(checksh);
          e2 = sapex(checksh);
          // Only do calculation if the projeciton of 'p' lies inside the
          //   subface [adjpt, e1, e2].
          projpt2face(ploop, adjpt, e1, e2, prjpt);
          facenormal(adjpt, e1, e2, n, 1, NULL);
          a = sqrt(dot(n, n)); // area of [adjpt, e1, e2].
          if (a > 0) {
            facenormal(adjpt, e1, prjpt, n, 1, NULL);
            a1 = sqrt(dot(n, n));
            facenormal(e1, e2, prjpt, n, 1, NULL);
            a2 = sqrt(dot(n, n));
            facenormal(e2, adjpt, prjpt, n, 1, NULL);
            a3 = sqrt(dot(n, n));
            if ((fabs(a1 + a2 + a3 - a) / a) < b->epsilon) {
              len = distance(ploop, prjpt);
              if (lfs_2 > len) {
                lfs_2 = len;
              }
            }
          } else {
            assert(0); // a degenerate triangle.
          } // if (a > 0)
        }
      }
      if (ploop[pointmtrindex] > lfs_2) {
        ploop[pointmtrindex] = lfs_2;
      }
      if (b->fixedvolume) {
        // A fixed volume constraint is imposed. Adjust H(p) <= maxlen.
        if (ploop[pointmtrindex] > maxlen) {
          ploop[pointmtrindex] = maxlen;
        }
      }
      if (b->varvolume) {
        // Variant volume constraints are imposed. Adjust H(p) <= varlen.
        for (i = 0; i < tetlist->objects; i++) {
          parytet = (triface *) fastlookup(tetlist, i);
          starttet = *parytet;
          vol = volumebound(starttet.tet);
          if (vol > 0.0) {
            varlen = pow(6 * vol, 1.0 / 3.0);
            if (ploop[pointmtrindex] > varlen) {
              ploop[pointmtrindex] = varlen;
            }
          }
        }
      }
      // The size is calculated.
      assert(ploop[pointmtrindex] > 0); // SELF_CHECK
      // Clear working lists.
      tetlist->restart();
      verlist->restart();
      featurecount++;
    } // if (featureflag)
    ploop = pointtraverse();
  }

  if (b->verbose) {
    printf("  %d feature points.\n", featurecount);
  }

  // Second only assign sizes for all Steiner points which were inserted on
  //   sharp segments. The sizes are interpolated from the endpoints of
  //   the segments.
  featurecount = 0;
  points->traversalinit();
  ploop = pointtraverse();
  while (ploop != (point) NULL) {
    if (ploop[pointmtrindex] == 0.0) {
      if (pointtype(ploop) == FREESEGVERTEX) {
        // A Steiner point on segment.
        featureflag = 0;
        sdecode(point2sh(ploop), checkseg);
        assert(checkseg.sh != NULL);
        checkseg.shver = 0;
        e1 = farsorg(checkseg);  // The origin of this seg.        
        e2 = farsdest(checkseg); // The dest of this seg.      
        if (b->nobisect) { // '-Y' option.
          assert(e1[pointmtrindex] > 0); // SELF_CHECK
          assert(e2[pointmtrindex] > 0); // SELF_CHECK
          featureflag = 1;
        } else {
          if ((e1[pointmtrindex] > 0) && (e2[pointmtrindex] > 0)) {
            featureflag = 1;
          }
        }
        if (featureflag) {
          len = distance(e1, e2);
          lfs_0 = distance(e1, ploop); // Re-use lfs_0.
          ploop[pointmtrindex] = e1[pointmtrindex]
            + (lfs_0 / len) * (e2[pointmtrindex] - e1[pointmtrindex]);
          featurecount++;
        } // if (featureflag)
      } 
    } // if (ploop[pointmtrindex] == 0.0)
    ploop = pointtraverse();
  }

  if (b->verbose && (featurecount > 0)) {
    printf("  %d Steiner feature points.\n", featurecount);
  }

  if (checkconstraints) {
    // A .var file exists. Adjust feature sizes. And make sure that every
    //   corner of a constraining facet get a size.
    if (in->facetconstraintlist) {
      // Have facet area constrains.
      subfaces->traversalinit();
      shloop.sh = shellfacetraverse(subfaces);
      while (shloop.sh != (shellface *) NULL) {
        varlen = areabound(shloop);
        if (varlen > 0.0) {
          // Check if the three corners are feature points.
          varlen = sqrt(varlen);
          for (j = 0; j < 3; j++) {
            ploop = (point) shloop.sh[3 + j];
            if (ploop[pointmtrindex] > 0) {
              if (ploop[pointmtrindex] > varlen) {
                ploop[pointmtrindex] = varlen;
              }
            } else {
              // This corner has no size yet. Set it.
              ploop[pointmtrindex] = varlen;
            }
          } // j
        }
        shloop.sh = shellfacetraverse(subfaces);
      }
    }
    if (in->segmentconstraintlist) {
      // Have facet area constrains.
      subsegs->traversalinit();
      shloop.sh = shellfacetraverse(subsegs);
      while (shloop.sh != (shellface *) NULL) {
        varlen = areabound(shloop);
        if (varlen > 0.0) {
          // Check if the two endpoints are feature points.
          for (j = 0; j < 2; j++) {
            ploop = (point) shloop.sh[3 + j];
            if (ploop[pointmtrindex] > 0.0) {
              if (ploop[pointmtrindex] > varlen) {
                ploop[pointmtrindex] = varlen;
              }
            } else {
              ploop[pointmtrindex] = varlen;
            }
          } // j
        }
        shloop.sh = shellfacetraverse(subsegs);
      }
    }
  } // if (checkconstraints)
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkseg4encroach()    Check if an edge is encroached upon by a point.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkseg4encroach(point pa, point pb, point checkpt)
{
  REAL ang;
  REAL prjpt[3], u, v, t;

  // Check if the point lies inside the diametrical sphere of this seg. 
  ang = interiorangle(checkpt, pa, pb, NULL);
  ang *= 2.0; // Compare it to PI/2 (90 degree).

  if (ang > PI) {
    // Inside.
    if (b->metric || b->nobisect) { // -m or -Y option.
      if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0)) {
        // In this case, we're sure that the projection of 'checkpt' lies
        //   inside the segment [a,b]. Check if 'checkpt' lies inside the
        //   protecting region of this seg.
        projpt2edge(checkpt, pa, pb, prjpt);
        // Get the mesh size at the location 'prjpt'.
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
  triface searchtet, spintet;
  point forg, fdest, eapex;
  REAL ccent[3], len, r, d, diff;
  int i;

  REAL ti, tj, t, midpt[3];
  REAL ang;
  int eid;

  forg = sorg(*chkseg);
  fdest = sdest(*chkseg);

  if (b->verbose > 2) {
    printf("      Check segment (%d, %d)\n", pointmark(forg), pointmark(fdest));
  }

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
      if (b->verbose > 2) {
        printf("      has too large size, len = %g (> %g)\n", len, 
               areabound(*chkseg));
      }
      qflag = 1;
      return 1;
    }
  }

  if (b->fixedvolume) { // if (b->varvolume || b->fixedvolume) {
    if ((len * len * len) > b->maxvolume) {
      if (b->verbose > 2) {
        printf("      has too large size, len^3 = %g (> %g)\n", len*len*len, 
               b->maxvolume);
      }
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

  if (b->psc) {
    // Check if it satisfies the approximation requirement.
    eid = shellmark(*chkseg);
    if ((pointtype(forg) == ACUTEVERTEX)||(pointtype(forg) == RIDGEVERTEX)) {
      ti = in->getvertexparamonedge(in->geomhandle, pointmark(forg), eid);
    } else {
      ti = pointgeomuv(forg, 0);
    }
    if ((pointtype(fdest) == ACUTEVERTEX)||(pointtype(fdest) == RIDGEVERTEX)) {
      tj = in->getvertexparamonedge(in->geomhandle, pointmark(fdest), eid);
    } else {
      tj = pointgeomuv(fdest, 0);
    }
    t = 0.5 * (ti + tj);
    in->getsteineronedge(in->geomhandle, eid, t, midpt);
    ang = interiorangle(midpt, forg, fdest, NULL) / PI * 180.0;
    if (ang < b->facet_ang_tol) {
      // Refine this segment.
      if (b->verbose > 2) {
        printf("      has bad approx, ang = %g\n", ang);
      }
      qflag = 1;
      return 1;
    }
  } // if (b->psc)

  // Second check if it is encroached.
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
        encpt = eapex;
        break;
      }
    }
    fnextself(spintet);
    if (spintet.tet == searchtet.tet) break;
  } // while (1)

  if (encpt != NULL) {
    if (b->verbose > 2) {
      printf("      is encroached by %d\n", pointmark(encpt));
    }
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

int tetgenmesh::splitsegment(face *splitseg, point encpt, int qflag, 
                             int chkencflag)
{
  triface searchtet;
  face searchsh;
  point newpt, pa, pb;
  insertvertexflags ivf;
  REAL len; //, len1;
  int loc;
  //int i;

  pa = sorg(*splitseg);
  pb = sdest(*splitseg);
  len = distance(pa, pb);

  if (b->verbose > 2) {
    printf("      Split segment (%d, %d).\n", pointmark(pa), pointmark(pb));
  }

  if (qflag == 0) {
    if (shelltype(*splitseg) == SHARP) {
      // Do not split it (due to a very small angle) even it is encroached.
      // Avoid creating too many Steiner points.
      return 0;
    }
    // Quickly check if we CAN split this segment.
    if (encpt == NULL) {
      // Do not split this segment if the length is smaller than the mesh
      //   size at one of its endpoints.    
      if ((len < pa[pointmtrindex]) || (len < pb[pointmtrindex])) {
        return 0;
      }
    }
  }

  makepoint(&newpt, FREESEGVERTEX);
  getsteinerptonsegment(splitseg, encpt, newpt);


  // Split the segment by the Bowyer-Watson algorithm.
  sstpivot1(*splitseg, searchtet);
  ivf.iloc = (int) ONEDGE;
  if (b->psc) {
    ivf.bowywat = 0;   // Do not enlarge the initial cavity.
    ivf.validflag = 0; // Do not validate the initial cavity.
  } else {
    ivf.bowywat = 3;   // Preserve subsegments and subfaces;
    ivf.validflag = 1; // Validate the B-W cavity.
  }
  ivf.lawson = b->conforming ? 3 : 1; // Check flip for internal new faces?.
  ivf.rejflag = 0;     // Do not check encroachment of new segments/facets.
  if ((encpt == NULL) && (qflag == 0)) {
    ivf.rejflag |= 4;  // Do check encroachment of protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = ivf.iloc;
  ivf.sbowywat = ivf.bowywat;  // Surface mesh options.
  ivf.splitbdflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = 1;

  loc = insertvertex(newpt, &searchtet, &searchsh, splitseg, &ivf);

  if (loc == (int) ONEDGE) {
    if (b->verbose > 2) {
      printf("      Point inserted successfully on segment.\n");
    }
    // Flip non-locally Delaunay faces at the link of its star.
    lawsonflip3d(newpt, 4, 0, chkencflag, 0);
    st_segref_count++;
    if (steinerleft > 0) steinerleft--;
    return 1;
  } else {
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
  badface *bface;
  point encpt = NULL;
  int qflag = 0;

  // Loop until the pool 'badsubsegs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badsubsegs->items > 0) && (steinerleft != 0)) {
    badsubsegs->traversalinit();
    bface = badfacetraverse(badsubsegs);
    while ((bface != NULL) && (steinerleft != 0)) {
      // A queued segment may have been deleted (split).
      if (bface->ss.sh[3] != NULL) {
        // A queued segment may have been processed. 
        if (smarktest2ed(bface->ss)) {
          sunmarktest2(bface->ss);
          if (checkseg4split(&(bface->ss), encpt, qflag)) {
            splitsegment(&(bface->ss), encpt, qflag, chkencflag);
          }
        }
      }
      badfacedealloc(badsubsegs, bface); // Remove this entry from list.
      bface = badfacetraverse(badsubsegs);
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
    bface = badfacetraverse(badsubsegs);
    while (bface  != NULL) {
      if (bface->ss.sh[3] != NULL) {
        if (smarktest2ed(bface->ss)) {
          sunmarktest2(bface->ss);
        }
      }
      bface = badfacetraverse(badsubsegs);
    }
    badsubsegs->restart();
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
  REAL prjpt[3], n[3];
  REAL a, a1, a2, a3;

  circumsphere(pa, pb, pc, NULL, cent, &rd);
  assert(rd != 0);
  len = distance(cent, checkpt);
  if ((fabs(len - rd) / rd) < b->epsilon) len = rd; // Rounding.
 
  if (len < rd) {
    // The point lies inside the circumsphere of this face.
    if (b->metric || b->nobisect) { // -m or -Y option.
      if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0) &&
          (pc[pointmtrindex] > 0)) {
        // Get the projection of 'checkpt' in the plane of pa, pb, and pc.
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
        } else {
          // The projection lies outside the face.
          // In this case, 'p' must close to another face or a segment than
          //   to this one. We ignore this boundary face. 
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
  triface searchtet;
  face checksh; // *parysh;
  face checkseg;
  point pa, pb, pc;
  REAL area, rd, len, sintheta;
  REAL A[4][4], rhs[4], D;
  int indx[4];
  REAL elen[3];
  int i;

  encpt = NULL;
  qflag = 0;

  pa = sorg(*chkfac);
  pb = sdest(*chkfac);
  pc = sapex(*chkfac);

  if (b->verbose > 2) {
    printf("      Check subface (%d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc));
  }

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
  elen[0] = dot(A[0], A[0]); // edge [a,b]
  elen[1] = dot(A[1], A[1]); // edge [a,c]
  rhs[0] = 0.5 * elen[0];
  rhs[1] = 0.5 * elen[1];
  rhs[2] = 0.0;

  // Solve the 3 by 3 equations use LU decomposition with partial 
  //   pivoting and backward and forward substitute..
  if (lu_decmp(A, 3, indx, &D, 0)) {
    lu_solve(A, 3, indx, rhs, 0);
    cent[0] = pa[0] + rhs[0];
    cent[1] = pa[1] + rhs[1];
    cent[2] = pa[2] + rhs[2];
    rd = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);

    if (b->verbose > 3) {
      printf("        circent: (%g, %g, %g)\n", cent[0], cent[1], cent[2]);
      printf("        cirradi: %g\n", rd);
    }

    // Check the quality (radius-edge ratio) of this subface.
    //   Re-use variables 'A', 'rhs', and 'D'.
    A[2][0] = pb[0] - pc[0];
    A[2][1] = pb[1] - pc[1];
    A[2][2] = pb[2] - pc[2];
    elen[2] = dot(A[2], A[2]); // edge [b,c]
    // Get the shortest edge length in 'D'.
    D = elen[0]; // edge [a,b]
    for (i = 1; i < 3; i++) {
      if (D > elen[i]) D = elen[i];
    }


    D = sqrt(D);
    if (b->verbose > 3) {
      printf("        shortest edge length = %g\n", D);
    }

    rhs[3] = rd / D; // The radius-edge ratio.

    // Check if this subface is nearly degenerate.
    sintheta = 1.0 / (2.0 * rhs[3]);
    if (sintheta < sintheta_tol) {
      // Do not split this subface. Save it in list.
      if (b->verbose > 1) {
        printf("  !! A degenerated subface, theta = %g (deg)\n",
               asin(sintheta) / PI * 180.0);
      }
      return 0; // Do not split a degenerated subface.
    }

    if (checkconstraints && (areabound(*chkfac) > 0.0)) {
      // Check if the subface has too big area.
      if (area > areabound(*chkfac)) {
        if (b->verbose > 2) {
          printf("      has too big area: %g (> %g)\n", area, 
                 areabound(*chkfac));
        }
        qflag = 1;
        return 1;
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


    // Check if this subface is locally encroached.
    for (i = 0; i < 2; i++) {
      stpivot(*chkfac, searchtet);
      if (!ishulltet(searchtet)) {
        len = distance(oppo(searchtet), cent);
        if ((fabs(len - rd) / rd) < b->epsilon) len = rd;// Rounding.
        if (len < rd) {
          if (b->verbose > 2) {
            printf("      is encroached by point %d\n", 
                   pointmark(oppo(searchtet)));
          }
          encpt = oppo(searchtet);
          return 1;
        }
      }
      sesymself(*chkfac);
    }
  } 

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// splitsubface()    Split a subface.                                        //
//                                                                           //
// The subface may be encroached, or in bad-quality. It is split at its cir- //
// cumcenter ('ccent'). Do not split it if 'ccent' encroaches upon any seg-  //
// ments. Instead, one of the encroached segments is split.  It is possible  //
// that none of the encorached segments can be split.                        //
//                                                                           //
// The return value indicates whether a new point is inserted (> 0) or not   //
// (= 0). Furthermore, it is inserted on an encorached segment (= 1) or in-  //
// side the facet (= 2).                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::splitsubface(face *splitfac, point encpt, int qflag,
                             REAL *ccent, int chkencflag)
{
  badface *bface;
  triface searchtet;
  face searchsh;
  face checkseg, *paryseg;
  point newpt, pa, pb, pc;
  insertvertexflags ivf;
  REAL rd;
  int splitflag;
  int loc;
  int i;


  pa = sorg(*splitfac);
  pb = sdest(*splitfac);
  pc = sapex(*splitfac);

  if (b->verbose > 2) {
    printf("      Split subface (%d, %d, %d).\n", pointmark(pa), pointmark(pb),
           pointmark(pc));
  }


  // Quickly check if we CAN split this subface.
  if (qflag == 0) {
    // Do not split this subface if it forms a very small dihedral with
    //   another facet. Avoid creating too many Steiner points.
    if (shelltype(*splitfac) == SHARP) {
      return 0;
    }
    // Do not split this subface if the 'ccent' lies inside the protect balls
    //   of one of its vertices.
    rd = distance(ccent, pa);
    if ((rd <= pa[pointmtrindex]) || (rd <= pb[pointmtrindex]) ||
        (rd <= pc[pointmtrindex])) {
      return 0;
    }
  }

  // Initialize the inserting point.
  makepoint(&newpt, FREEFACETVERTEX);

    // Split the subface at its circumcenter.
    for (i = 0; i < 3; i++) newpt[i] = ccent[i];
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
    if ((ivf.iloc == (int) ONFACE) || (ivf.iloc == (int) ONEDGE)) {
      // Insert this point.
    } else {
      pointdealloc(newpt);
      return 0;
    }


  // Insert the point.
  stpivot(searchsh, searchtet);
  //assert((ivf.iloc == (int) ONFACE) || (ivf.iloc == (int) ONEDGE));
  // Split the subface by the Bowyer-Watson algorithm.
  ivf.bowywat = 3; // Preserve segments and subfaces.
  ivf.lawson = b->conforming ? 3 : 1;
  ivf.rejflag = 1; // Do check the encroachment of segments.
  if (qflag == 0) {
    ivf.rejflag |= 4; // Reject it if it encroached upon any vertex.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = ivf.iloc;
  ivf.sbowywat = ivf.bowywat;
  ivf.splitbdflag = 1;
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = 1;

  ivf.refineflag = 2;
  ivf.refinesh = searchsh;

  loc = insertvertex(newpt, &searchtet, &searchsh, NULL, &ivf);

  if (loc == (int) ivf.iloc) {
    if (b->verbose > 2) {
      printf("      Point inserted successfully on facet.\n");
    }
    // Flip not locally Delaunay link facets.
    lawsonflip3d(newpt, 4, 0, chkencflag, 0);
    st_facref_count++;
    if (steinerleft > 0) steinerleft--;
    return 1;
  } else {
    // Point was not inserted.
    if (loc == (int) ENCSEGMENT) {
      if (b->verbose > 2) {
        printf("      Point encroached upon %ld segments.\n", 
               encseglist->objects);
      }
      assert(encseglist->objects > 0);
      pointdealloc(newpt);
      // Select an encroached segment and split it.
      splitflag = 0;
      for (i = 0; i < encseglist->objects; i++) {
        paryseg = (face *) fastlookup(encseglist, i);
        if (splitsegment(paryseg, NULL, qflag, chkencflag | 1)) {
          splitflag = 1; // A point is inserted on a segment.
          break;
        }
      }
      encseglist->restart();
      if (splitflag) {
        // Some segments may need to be repaired.
        repairencsegs(chkencflag | 1);
        // Queue this subface if it is still alive and not queued.
        if (splitfac->sh[3] != NULL) {
          if (!smarktest2ed(*splitfac)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = *splitfac;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(*splitfac); // An alive badface.
          }
        }
      }
      return splitflag;
    } else {
      pointdealloc(newpt);
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
  badface *bface;
  point encpt = NULL;
  int qflag = 0;
  REAL ccent[3];

  // Loop until the pool 'badsubfacs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badsubfacs->items > 0) && (steinerleft != 0)) {
    badsubfacs->traversalinit();
    bface = badfacetraverse(badsubfacs);
    while ((bface != NULL) && (steinerleft != 0)) {
      // A queued subface may have been deleted (split).
      if (bface->ss.sh[3] != NULL) {
        // A queued subface may have been processed. 
        if (smarktest2ed(bface->ss)) {
          sunmarktest2(bface->ss);
          if (checkfac4split(&(bface->ss), encpt, qflag, ccent)) {
            splitsubface(&(bface->ss), encpt, qflag, ccent, chkencflag);
          }
        }
      }
      badfacedealloc(badsubfacs, bface); // Remove this entry from list.
      bface = badfacetraverse(badsubfacs);
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
    bface = badfacetraverse(badsubfacs);
    while (bface  != NULL) {
      if (bface->ss.sh[3] != NULL) {
        if (smarktest2ed(bface->ss)) {
          sunmarktest2(bface->ss);
        }
      }
      bface = badfacetraverse(badsubfacs);
    }
    badsubfacs->restart();
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
  REAL maxcosd, vol, volbnd, smlen, rd;
  REAL A[4][4], rhs[4], D;
  int indx[4];
  int i, j;

  qflag = 0;

  pd = (point) chktet->tet[7];
  if (pd == dummypoint) {
    return 0; // Do not split a hull tet.
  }

  pa = (point) chktet->tet[4];
  pb = (point) chktet->tet[5];
  pc = (point) chktet->tet[6];

  if (b->verbose > 2) {
    printf("      Check tet (%d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd));
  }

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
    if (b->verbose > 3) {
      printf("        Min dihed = 0 (degree)\n");
    }
    // Return its barycenter.
    for (i = 0; i < 3; i++) {
      ccent[i] = 0.25 * (pa[i] + pb[i] + pc[i] + pd[i]);
    }
    return 1;
  }

  // Check volume if '-a#' and '-a' options are used.
  if (b->varvolume || b->fixedvolume) {
    vol = fabs(A[indx[0]][0] * A[indx[1]][1] * A[indx[2]][2]) / 6.0;
    if (b->verbose > 3) {
      printf("        volume = %g.\n", vol);
    }
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
    } else {
      return 0; // Do not split this tet.
    }
  }

  // Check the radius-edge ratio. Set by -q#.
  if (b->minratio > 0) { 
    // Calculate the circumcenter and radius of this tet.
    rhs[0] = 0.5 * dot(vda, vda);
    rhs[1] = 0.5 * dot(vdb, vdb);
    rhs[2] = 0.5 * dot(vdc, vdc);
    lu_solve(A, 3, indx, rhs, 0);            
    for (i = 0; i < 3; i++) ccent[i] = pd[i] + rhs[i];
    rd = sqrt(dot(rhs, rhs));
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
    D = rd / smlen;
    if (b->verbose > 3) {
      printf("        Ratio-edge ratio = %g, smlen = %g\n", D, smlen);
    }
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
    // Get the smallest diehedral angle.
    //maxcosd = mincosd = cosd[0];
    maxcosd = cosd[0];
    for (i = 1; i < 6; i++) {
      //if (cosd[i] > maxcosd) maxcosd = cosd[i];
      maxcosd = (cosd[i] > maxcosd ? cosd[i] : maxcosd);
      //mincosd = (cosd[i] < mincosd ? cosd[i] : maxcosd);
    }
    if (b->verbose > 3) {
      printf("        Min dihed = %g (degree)\n", acos(maxcosd) / PI * 180.0);
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
  badface *bface;
  triface searchtet;
  face checkseg, *paryseg;
  point newpt, pa, *ppt = NULL;
  insertvertexflags ivf;
  REAL rd;
  int splitflag;
  int loc;
  int i;

  if (b->verbose > 2) {
    ppt = (point *) &(splittet->tet[4]);
    printf("      Split tet (%d, %d, %d, %d).\n", pointmark(ppt[0]), 
           pointmark(ppt[1]), pointmark(ppt[2]), pointmark(ppt[3]));
  }


  if (qflag == 0) {
    // It is a bad quality tet (not due to mesh size).
    // It can be split if 'ccent' does not encroach upon any prot. balls.
    //   Do a quick check if the 'ccent' lies inside the protect balls
    //   of one of the vertices of this tet.
    ppt = (point *) &(splittet->tet[4]);
    rd = distance(ccent, ppt[0]);
    if ((rd <= ppt[0][pointmtrindex]) || (rd <= ppt[1][pointmtrindex]) ||
        (rd <= ppt[2][pointmtrindex]) || (rd <= ppt[3][pointmtrindex])) {
      if (b->verbose > 2) {
        printf("      Encroaching a protecting ball. Rejected.\n");
      }
      return 0;
    }
  }

  makepoint(&newpt, FREEVOLVERTEX);
  for (i = 0; i < 3; i++) newpt[i] = ccent[i];


  searchtet = *splittet;
  ivf.iloc = (int) OUTSIDE;
  ivf.bowywat = 3;  // Preserve subsegments and subfaces;
  ivf.lawson = b->conforming ? 3 : 1;
  ivf.rejflag = 3;  // Do check for encroached segments and subfaces.
  if (qflag == 0) {
    ivf.rejflag |= 4; // Reject it if it lies in some protecting balls.
  }
  ivf.chkencflag = chkencflag;
  ivf.sloc = ivf.sbowywat = 0; // No use.
  ivf.splitbdflag = 0; // No use.
  ivf.validflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = 1;

  ivf.refineflag = 1;
  ivf.refinetet = *splittet;

  loc = insertvertex(newpt, &searchtet, NULL, NULL, &ivf);

  if (loc == (int) ENCSEGMENT) {
    if (b->verbose > 2) {
      printf("      Point encroached upon %ld segments.\n", 
             encseglist->objects);
    }
    pointdealloc(newpt);
    assert(encseglist->objects > 0);
    splitflag = 0;
    if (!b->nobisect) { // not -Y option
      // Select an encroached segment and split it.
      for (i = 0; i < encseglist->objects; i++) {
        paryseg = (face *) fastlookup(encseglist, i);
        if (splitsegment(paryseg, NULL, qflag, chkencflag | 3)) {
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
      if (splittet->tet[4] != NULL) {
        if (!marktest2ed(*splittet)) {
          bface = (badface *) badtetrahedrons->alloc();
          bface->tt = *splittet;
          marktest2(bface->tt); // Only queue it once.
          bface->forg = org(*splittet); // An alive badface.
        }
      }
    }
    return splitflag;
  } else if (loc == (int) ENCSUBFACE) {
    if (b->verbose > 2) {
      printf("      Point encroached upon %ld subfaces.\n", 
             encshlist->objects);
    }
    pointdealloc(newpt);
    assert(encshlist->objects > 0);
    splitflag = 0;
    if (!b->nobisect) { // not -Y option
      // Select an encroached subface and split it.
      for (i = 0; i < encshlist->objects; i++) {
        bface = (badface *) fastlookup(encshlist, i);
        if (splitsubface(&(bface->ss),NULL,qflag,bface->cent,chkencflag | 2)) {
          splitflag = 1; // A point is inserted on a subface or a segment.
          break;
        }
      }
    } // if (!b->nobisect)
    encshlist->restart();
    if (splitflag) {
      assert(badsubsegs->items == 0l); // repairencsegs(chkencflag | 3);
      // Some subfaces may need to be repaired.
      repairencfacs(chkencflag | 2);
      // Queue the tet if it is still alive.
      if (splittet->tet[4] != NULL) {
        if (!marktest2ed(*splittet)) {
          bface = (badface *) badtetrahedrons->alloc();
          bface->tt = *splittet;
          marktest2(bface->tt); // Only queue it once.
          bface->forg = org(*splittet); // An alive badface.
        }
      }
    }
    return splitflag;
  } else if (loc == (int) OUTSIDE) {
    // There exists not boundary conforming segments/subfaces.
    pointdealloc(newpt);
  } else if (loc == (int) ONVERTEX) {
    // Found a coincident vertex. It should be a Steiner point.
    pa = org(searchtet);
    assert(pointtype(pa) == FREEVOLVERTEX);
    // Delete this new point.
    pointdealloc(newpt);
  } else if (loc == (int) NEARVERTEX) {
    // The point lies very close to an existing point.
    pa = point2ppt(newpt);
    assert(pointtype(pa) == FREEVOLVERTEX);
    // Delete this new point.
    pointdealloc(newpt);
  } else if (loc == (int) ENCVERTEX) {
    // The new point encoraches upon some protecting balls. Rejected.
    pointdealloc(newpt);
  } else if (loc == (int) BADELEMENT) {
    pointdealloc(newpt);
  } else {
    if (b->verbose > 2) {
      printf("      Point inserted successfully.\n");
    }
    // Recover Delaunayness.
    lawsonflip3d(newpt, 4, 0, chkencflag, 0);
    // Vertex is inserted.
    st_volref_count++;
    if (steinerleft > 0) steinerleft--;
    return 1;
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// repairbadtets()    Repair bad quality tetrahedra.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::repairbadtets(int chkencflag)
{
  badface *bface;
  REAL ccent[3];
  int qflag = 0;

  // Loop until the pool 'badsubfacs' is empty. Note that steinerleft == -1
  //   if an unlimited number of Steiner points is allowed.
  while ((badtetrahedrons->items > 0) && (steinerleft != 0)) {
    badtetrahedrons->traversalinit();
    bface = badfacetraverse(badtetrahedrons);
    while ((bface != NULL) && (steinerleft != 0)) {
      // A queued tet may have been deleted.
      if (!isdeadtet(bface->tt)) {
        // A queued tet may have been processed.
        if (marktest2ed(bface->tt)) {
          unmarktest2(bface->tt);
          if (checktet4split(&(bface->tt), qflag, ccent)) {
            splittetrahedron(&(bface->tt), qflag, ccent, chkencflag);
          }
        }
      }
      badfacedealloc(badtetrahedrons, bface);
      bface = badfacetraverse(badtetrahedrons);
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
    bface = badfacetraverse(badtetrahedrons);
    while (bface != NULL) {
      if (!isdeadtet(bface->tt)) {
        if (marktest2ed(bface->tt)) {
          unmarktest2(bface->tt);
        }
      }
      bface = badfacetraverse(badtetrahedrons);
    }
    // Clear the pool.
    badtetrahedrons->restart();
  }
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// enforcequality()    Refine the mesh.                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::delaunayrefinement()
{
  badface *bface;
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

  if (b->refine || b->nobisect) { // '-r' or '-Y' option.
    markacutevertices();
  }

  marksharpsegments();

  decidefeaturepointsizes();


  encseglist = new arraypool(sizeof(face), 8);
  encshlist = new arraypool(sizeof(badface), 8);

  if (!b->nobisect) { // if no '-Y' option
    if (b->verbose) {
      printf("  Splitting encroached subsegments.\n");
    }

    chkencflag = 1; // Only check encroaching subsegments.
    steinercount = points->items;

    // Initialize the pool of encroached subsegments.
    badsubsegs = new memorypool(sizeof(badface), b->shellfaceperblock, 
                                memorypool::POINTER, 0);

    // Add all segments into the pool.
    subsegs->traversalinit();
    checkseg.sh = shellfacetraverse(subsegs);
    while (checkseg.sh != (shellface *) NULL) {
      bface = (badface *) badsubsegs->alloc();
      bface->ss = checkseg;
      smarktest2(bface->ss); // Only queue it once.
      bface->forg = sorg(checkseg); // An alive badface.
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
      badsubfacs = new memorypool(sizeof(badface), b->shellfaceperblock, 
                                  memorypool::POINTER, 0);

      // Add all subfaces into the pool.
      subfaces->traversalinit();
      checksh.sh = shellfacetraverse(subfaces);
      while (checksh.sh != (shellface *) NULL) {
        bface = (badface *) badsubfacs->alloc();
        bface->ss = checksh;
        smarktest2(bface->ss); // Only queue it once.
        bface->forg = sorg(checksh); // An alive badface.
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
    badtetrahedrons = new memorypool(sizeof(badface), b->tetrahedraperblock,
                                     memorypool::POINTER, 0);

    // Add all tetrahedra (no hull tets) into the pool.
    tetrahedrons->traversalinit();
    checktet.tet = tetrahedrontraverse();
    while (checktet.tet != NULL) {
      bface = (badface *) badtetrahedrons->alloc();
      bface->tt = checktet;
      marktest2(bface->tt); // Only queue it once.
      bface->forg = org(checktet); // An alive badface.
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
      printf("  Performed %ld flips.\n", flip23count + flip32count +
             flip44count - bak_flipcount);
    }
  } // if (b->reflevel > 2)

  if (steinerleft == 0) {
    if (!b->quiet) {
      printf("\nWarnning:  ");
      printf("The desired number of Steiner points (%d) is reached.\n\n",
             b->steinerleft);
    }
  }

  delete encseglist;
  delete encshlist;

  if (!b->nobisect) {
    delete badsubsegs;
    if (b->reflevel > 1) {
      delete badsubfacs;
    }
  }
  if (b->reflevel > 2) {
    delete badtetrahedrons;
  }
}

////                                                                       ////
////                                                                       ////
//// refine_cxx ///////////////////////////////////////////////////////////////

