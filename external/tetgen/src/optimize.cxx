#include "../tetgen.h"
//// optimize_cxx /////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lawsonflip3d()    A three-dimensional Lawson's algorithm.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::lawsonflip3d(flipconstraints *fc)
{
  triface fliptets[5], neightet, hulltet;
  face checksh, casingout;
  badface *popface, *bface;
  point pd, pe, *pts;
  REAL sign, ori;
  long flipcount, totalcount = 0l;
  long sliver_peels = 0l;
  int t1ver;
  int i;


  while (1) {

    if (b->verbose > 2) {
      printf("      Lawson flip %ld faces.\n", flippool->items);
    }
    flipcount = 0l;

    while (flipstack != (badface *) NULL) {
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

      if (ishulltet(fliptets[0])) continue;

      fsym(fliptets[0], fliptets[1]);
      if (ishulltet(fliptets[1])) {
        if (nonconvex) {
          // Check if 'fliptets[0]' it is a hull sliver.
          tspivot(fliptets[0], checksh);
          for (i = 0; i < 3; i++) {
            if (!isshsubseg(checksh)) {
              spivot(checksh, casingout);
              //assert(casingout.sh != NULL);
              if (sorg(checksh) != sdest(casingout)) sesymself(casingout);
              stpivot(casingout, neightet);
              if (neightet.tet == fliptets[0].tet) {
                // Found a hull sliver 'neightet'. Let it be [e,d,a,b], where 
                //   [e,d,a] and [d,e,b] are hull faces.
                edestoppo(neightet, hulltet); // [a,b,e,d]
                fsymself(hulltet); // [b,a,e,#]
                if (oppo(hulltet) == dummypoint) {
                  pe = org(neightet);
                  if ((pointtype(pe) == FREEFACETVERTEX) ||
                      (pointtype(pe) == FREESEGVERTEX)) {
                    removevertexbyflips(pe);
                  }
                } else {
                  eorgoppo(neightet, hulltet); // [b,a,d,e]
                  fsymself(hulltet); // [a,b,d,#]
                  if (oppo(hulltet) == dummypoint) {
                    pd = dest(neightet);
                    if ((pointtype(pd) == FREEFACETVERTEX) ||
                        (pointtype(pd) == FREESEGVERTEX)) {
                      removevertexbyflips(pd);
                    }
                  } else {
                    // Perform a 3-to-2 flip to remove the sliver.
                    fliptets[0] = neightet;          // [e,d,a,b]
                    fnext(fliptets[0], fliptets[1]); // [e,d,b,c]
                    fnext(fliptets[1], fliptets[2]); // [e,d,c,a]
                    flip32(fliptets, 1, fc);
                    // Update counters.
                    flip32count--;
                    flip22count--;
                    sliver_peels++;
                    if (fc->remove_ndelaunay_edge) {
                      // Update the volume (must be decreased).
                      //assert(fc->tetprism_vol_sum <= 0);
                      tetprism_vol_sum += fc->tetprism_vol_sum;
                      fc->tetprism_vol_sum = 0.0; // Clear it.
                    }
                  }
                }
                break;
              } // if (neightet.tet == fliptets[0].tet)
            } // if (!isshsubseg(checksh))
            senextself(checksh);
          } // i
        } // if (nonconvex)
        continue;
      }

      if (checksubfaceflag) {
        // Do not flip if it is a subface.
        if (issubface(fliptets[0])) continue;
      }

      // Test whether the face is locally Delaunay or not.
      pts = (point *) fliptets[1].tet; 
      sign = insphere_s(pts[4], pts[5], pts[6], pts[7], oppo(fliptets[0]));

      if (sign < 0) {
        // A non-Delaunay face. Try to flip it.
        pd = oppo(fliptets[0]);
        pe = oppo(fliptets[1]);

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
          if (fc->remove_ndelaunay_edge) {
            // Update the volume (must be decreased).
            //assert(fc->tetprism_vol_sum <= 0);
            tetprism_vol_sum += fc->tetprism_vol_sum;
            fc->tetprism_vol_sum = 0.0; // Clear it.
          }
          continue;
        } else { // ori <= 0
          // The edge ('fliptets[0]' = [a',b',c',d]) is non-convex or flat,
          //   where the edge [a',b'] is one of [a,b], [b,c], and [c,a].
          if (checksubsegflag) {
            // Do not flip if it is a segment.
            if (issubseg(fliptets[0])) continue;
          }
          // Check if there are three or four tets sharing at this edge.        
          esymself(fliptets[0]); // [b,a,d,c]
          for (i = 0; i < 3; i++) {
            fnext(fliptets[i], fliptets[i+1]);
          }
          if (fliptets[3].tet == fliptets[0].tet) {
            // A 3-to-2 flip is found. (No hull tet.)
            flip32(fliptets, 0, fc); 
            flipcount++;
            if (fc->remove_ndelaunay_edge) {
              // Update the volume (must be decreased).
              //assert(fc->tetprism_vol_sum <= 0);
              tetprism_vol_sum += fc->tetprism_vol_sum;
              fc->tetprism_vol_sum = 0.0; // Clear it.
            }
            continue;
          } else {
            // There are more than 3 tets at this edge.
            fnext(fliptets[3], fliptets[4]);
            if (fliptets[4].tet == fliptets[0].tet) {
              // There are exactly 4 tets at this edge.
              if (nonconvex) {
                if (apex(fliptets[3]) == dummypoint) {
                  // This edge is locally non-convex on the hull.
                  // It can be removed by a 4-to-4 flip.                  
                  ori = 0;
                }
              } // if (nonconvex)
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
                flip23count--;
                flip32count--;
                flip44count++;
                if (fc->remove_ndelaunay_edge) {
                  // Update the volume (must be decreased).
                  //assert(fc->tetprism_vol_sum <= 0);
                  tetprism_vol_sum += fc->tetprism_vol_sum;
                  fc->tetprism_vol_sum = 0.0; // Clear it.
                }
                continue;
              } // if (ori == 0)
            }
          }
        } // if (ori <= 0)

        // This non-Delaunay face is unflippable. Save it.
        unflipqueue->newindex((void **) &bface);
        bface->tt = fliptets[0];
        bface->forg  = org(fliptets[0]);
        bface->fdest = dest(fliptets[0]);
        bface->fapex = apex(fliptets[0]);
      } // if (sign < 0)
    } // while (flipstack)

    if (b->verbose > 2) {
      if (flipcount > 0) {
        printf("      Performed %ld flips.\n", flipcount);
      }
    }
    // Accumulate the counter of flips.
    totalcount += flipcount;

    assert(flippool->items == 0l);
    // Return if no unflippable faces left.
    if (unflipqueue->objects == 0l) break; 
    // Return if no flip has been performed.
    if (flipcount == 0l) break;

    // Try to flip the unflippable faces.
    for (i = 0; i < unflipqueue->objects; i++) {
      bface = (badface *) fastlookup(unflipqueue, i);
      if (!isdeadtet(bface->tt) && 
          (org(bface->tt) == bface->forg) &&
          (dest(bface->tt) == bface->fdest) &&
          (apex(bface->tt) == bface->fapex)) {
        flippush(flipstack, &(bface->tt));
      }
    }
    unflipqueue->restart();

  } // while (1)

  if (b->verbose > 2) {
    if (totalcount > 0) {
      printf("      Performed %ld flips.\n", totalcount);
    }
    if (sliver_peels > 0) {
      printf("      Removed %ld hull slivers.\n", sliver_peels);
    }
    if (unflipqueue->objects > 0l) {
      printf("      %ld unflippable edges remained.\n", unflipqueue->objects);
    }
  }

  return totalcount + sliver_peels;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoverdelaunay()    Recovery the locally Delaunay property.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::recoverdelaunay()
{
  arraypool *flipqueue, *nextflipqueue, *swapqueue;
  triface tetloop, neightet, *parytet;
  badface *bface, *parybface;
  point *ppt;
  flipconstraints fc;
  int i, j;

  if (!b->quiet) {
    printf("Recovering Delaunayness...\n");
  }

  tetprism_vol_sum = 0.0; // Initialize it.

  // Put all interior faces of the mesh into 'flipstack'.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != NULL) {
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      decode(tetloop.tet[tetloop.ver], neightet);
      if (!facemarked(neightet)) {
        flippush(flipstack, &tetloop);
      }
    }
    ppt = (point *) &(tetloop.tet[4]);
    tetprism_vol_sum += tetprismvol(ppt[0], ppt[1], ppt[2], ppt[3]);
    tetloop.tet = tetrahedrontraverse();
  }

  // Calulate a relatively lower bound for small improvement. 
  //   Used to avoid rounding error in volume calculation.
  fc.bak_tetprism_vol = tetprism_vol_sum * b->epsilon * 1e-3;

  if (b->verbose) {
    printf("  Initial obj = %.17g\n", tetprism_vol_sum);
  }

  if (b->verbose > 1) {
    printf("    Recover Delaunay [Lawson] : %ld\n", flippool->items);
  }

  // First only use the basic Lawson's flip.
  fc.remove_ndelaunay_edge = 1;
  fc.enqflag = 2;

  lawsonflip3d(&fc);

  if (b->verbose > 1) {
    printf("    obj (after Lawson) = %.17g\n", tetprism_vol_sum);
  }

  if (unflipqueue->objects == 0l) {
    return; // The mesh is Delaunay.
  }

  fc.unflip = 1; // Unflip if the edge is not flipped.
  fc.collectnewtets = 1; // new tets are returned in 'cavetetlist'.
  fc.enqflag = 0;

  autofliplinklevel = 1; // Init level.
  b->fliplinklevel = -1; // No fixed level.

  // For efficiency reason, we limit the maximium size of the edge star.
  int bakmaxflipstarsize = b->flipstarsize;
  b->flipstarsize = 10; // default

  flipqueue = new arraypool(sizeof(badface), 10);
  nextflipqueue = new arraypool(sizeof(badface), 10);
  
  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  while (flipqueue->objects > 0l) {

    if (b->verbose > 1) {
      printf("    Recover Delaunay [level = %2d] #:  %ld.\n",
             autofliplinklevel, flipqueue->objects);
    }

    for (i = 0; i < flipqueue->objects; i++) {
      bface  = (badface *) fastlookup(flipqueue, i);
      if (getedge(bface->forg, bface->fdest, &bface->tt)) {
        if (removeedgebyflips(&(bface->tt), &fc) == 2) {
          tetprism_vol_sum += fc.tetprism_vol_sum;
          fc.tetprism_vol_sum = 0.0; // Clear it.
          // Queue new faces for flips.
          for (j = 0; j < cavetetlist->objects; j++) {
            parytet = (triface *) fastlookup(cavetetlist, j);
            // A queued new tet may be dead.
            if (!isdeadtet(*parytet)) {
              for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
                // Avoid queue a face twice.
                decode(parytet->tet[parytet->ver], neightet);
                if (!facemarked(neightet)) {
                  flippush(flipstack, parytet);
                }
              } // parytet->ver
            }
          } // j
          cavetetlist->restart();
          // Remove locally non-Delaunay faces. New non-Delaunay edges
          //   may be found. They are saved in 'unflipqueue'.
          fc.enqflag = 2;
          lawsonflip3d(&fc);
          fc.enqflag = 0;
          // There may be unflipable faces. Add them in flipqueue.
          for (j = 0; j < unflipqueue->objects; j++) {
            bface  = (badface *) fastlookup(unflipqueue, j);
            flipqueue->newindex((void **) &parybface);
            *parybface = *bface;
          }
          unflipqueue->restart();
        } else {
          // Unable to remove this edge. Save it.
          nextflipqueue->newindex((void **) &parybface);
          *parybface = *bface;
          // Normally, it should be zero. 
          //assert(fc.tetprism_vol_sum == 0.0);
          // However, due to rounding errors, a tiny value may appear.
          fc.tetprism_vol_sum = 0.0;
        }
      }
    } // i

    if (b->verbose > 1) {
      printf("    obj (after level %d) = %.17g.\n", autofliplinklevel,
             tetprism_vol_sum);
    }
    flipqueue->restart();

    // Swap the two flip queues.
    swapqueue = flipqueue;
    flipqueue = nextflipqueue;
    nextflipqueue = swapqueue;

    if (flipqueue->objects > 0l) {
      // default 'b->delmaxfliplevel' is 1.
      if (autofliplinklevel >= b->delmaxfliplevel) {
        // For efficiency reason, we do not search too far.
        break;
      }
      autofliplinklevel+=b->fliplinklevelinc;
    }
  } // while (flipqueue->objects > 0l)

  if (flipqueue->objects > 0l) {
    if (b->verbose > 1) {
      printf("    %ld non-Delaunay edges remained.\n", flipqueue->objects);
    }
  }

  if (b->verbose) {
    printf("  Final obj  = %.17g\n", tetprism_vol_sum);
  }

  b->flipstarsize = bakmaxflipstarsize;
  delete flipqueue;
  delete nextflipqueue;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// gettetrahedron()    Get a tetrahedron which have the given vertices.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::gettetrahedron(point pa, point pb, point pc, point pd, 
                               triface *searchtet)
{
  triface spintet;
  int t1ver; 

  if (getedge(pa, pb, searchtet)) {
    spintet = *searchtet;
    while (1) {
      if (apex(spintet) == pc) {
        *searchtet = spintet;
        break;
      }
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    }
    if (apex(*searchtet) == pc) {
      if (oppo(*searchtet) == pd) {
        return 1;
      } else {
        fsymself(*searchtet);
        if (oppo(*searchtet) == pd) {
          return 1;
        }
      }
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// improvequalitybyflips()    Improve the mesh quality by flips.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::improvequalitybyflips()
{
  arraypool *flipqueue, *nextflipqueue, *swapqueue;
  badface *bface, *parybface;
  triface *parytet;
  point *ppt;
  flipconstraints fc;
  REAL *cosdd, ncosdd[6], maxdd;
  long totalremcount, remcount;
  int remflag;
  int n, i, j, k;

  //assert(unflipqueue->objects > 0l);
  flipqueue = new arraypool(sizeof(badface), 10);
  nextflipqueue = new arraypool(sizeof(badface), 10);

  // Backup flip edge options.
  int bakautofliplinklevel = autofliplinklevel;
  int bakfliplinklevel = b->fliplinklevel;
  int bakmaxflipstarsize = b->flipstarsize;

  // Set flip edge options.
  autofliplinklevel = 1; 
  b->fliplinklevel = -1;
  b->flipstarsize = 10; // b->optmaxflipstarsize;

  fc.remove_large_angle = 1;
  fc.unflip = 1;
  fc.collectnewtets = 1;
  fc.checkflipeligibility = 1;

  totalremcount = 0l;

  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  while (flipqueue->objects > 0l) {

    remcount = 0l;

    while (flipqueue->objects > 0l) {
      if (b->verbose > 1) {
        printf("    Improving mesh qualiy by flips [%d]#:  %ld.\n",
               autofliplinklevel, flipqueue->objects);
      }

      for (k = 0; k < flipqueue->objects; k++) {
        bface  = (badface *) fastlookup(flipqueue, k);
        if (gettetrahedron(bface->forg, bface->fdest, bface->fapex,
                           bface->foppo, &bface->tt)) {
          //assert(!ishulltet(bface->tt));
          // There are bad dihedral angles in this tet.
          if (bface->tt.ver != 11) {
            // The dihedral angles are permuted.
            // Here we simply re-compute them. Slow!!.
            ppt = (point *) & (bface->tt.tet[4]);
            tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                           &bface->key, NULL);
            bface->forg = ppt[0];
            bface->fdest = ppt[1];
            bface->fapex = ppt[2];
            bface->foppo = ppt[3];
            bface->tt.ver = 11;
          }
          if (bface->key == 0) {
            // Re-comput the quality values. Due to smoothing operations.
            ppt = (point *) & (bface->tt.tet[4]);
            tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                           &bface->key, NULL);
          }
          cosdd = bface->cent;
          remflag = 0;
          for (i = 0; (i < 6) && !remflag; i++) {
            if (cosdd[i] < cosmaxdihed) {
              // Found a large dihedral angle.
              bface->tt.ver = edge2ver[i]; // Go to the edge.
              fc.cosdihed_in = cosdd[i];
              fc.cosdihed_out = 0.0; // 90 degree.
              n = removeedgebyflips(&(bface->tt), &fc);
              if (n == 2) {
                // Edge is flipped.
                remflag = 1;
                if (fc.cosdihed_out < cosmaxdihed) {
                  // Queue new bad tets for further improvements.
                  for (j = 0; j < cavetetlist->objects; j++) {
                    parytet = (triface *) fastlookup(cavetetlist, j);
                    if (!isdeadtet(*parytet)) {
                      ppt = (point *) & (parytet->tet[4]);
                      // Do not test a hull tet.
                      if (ppt[3] != dummypoint) {
                        tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], ncosdd, 
                                       &maxdd, NULL);
                        if (maxdd < cosmaxdihed) {
                          // There are bad dihedral angles in this tet.
                          nextflipqueue->newindex((void **) &parybface); 
                          parybface->tt.tet = parytet->tet;
                          parybface->tt.ver = 11;
                          parybface->forg = ppt[0];
                          parybface->fdest = ppt[1];
                          parybface->fapex = ppt[2];
                          parybface->foppo = ppt[3];
                          parybface->key = maxdd;
                          for (n = 0; n < 6; n++) {
                            parybface->cent[n] = ncosdd[n];
                          }
                        }
                      } // if (ppt[3] != dummypoint) 
                    }
                  } // j
                } // if (fc.cosdihed_out < cosmaxdihed)
                cavetetlist->restart();
                remcount++;
              }
            }
          } // i          
          if (!remflag) {
            // An unremoved bad tet. Queue it again. 
            unflipqueue->newindex((void **) &parybface);
            *parybface = *bface;
          }
        } // if (gettetrahedron(...))
      } // k

      flipqueue->restart();

      // Swap the two flip queues.
      swapqueue = flipqueue;
      flipqueue = nextflipqueue;
      nextflipqueue = swapqueue;
    } // while (flipqueues->objects > 0)

    if (b->verbose > 1) {
      printf("    Removed %ld bad tets.\n", remcount);
    }
    totalremcount += remcount;

    if (unflipqueue->objects > 0l) {
      //if (autofliplinklevel >= b->optmaxfliplevel) {
      if (autofliplinklevel >= b->optlevel) {
        break;
      }
      autofliplinklevel+=b->fliplinklevelinc;
      //b->flipstarsize = 10 + (1 << (b->optlevel - 1));
    }

    // Swap the two flip queues.
    swapqueue = flipqueue;
    flipqueue = unflipqueue;
    unflipqueue = swapqueue;
  } // while (flipqueues->objects > 0)

  // Restore original flip edge options.
  autofliplinklevel = bakautofliplinklevel;
  b->fliplinklevel = bakfliplinklevel;
  b->flipstarsize = bakmaxflipstarsize;

  delete flipqueue;
  delete nextflipqueue;

  return totalremcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// smoothpoint()    Moving a vertex to improve the mesh quality.             //
//                                                                           //
// 'smtpt' (p) is a point to be smoothed. Generally, it is a Steiner point.  //
// It may be not a vertex of the mesh.                                       //
//                                                                           //
// This routine tries to move 'p' inside its star until a selected objective //
// function over all tetrahedra in the star is improved. The function may be //
// the some quality measures, i.e., aspect ratio, maximum dihedral angel, or //
// simply the volume of the tetrahedra.                                      //
//                                                                           //
// 'linkfacelist' contains the list of link faces of 'p'.  Since a link face //
// has two orientations, ccw or cw, with respect to 'p'.  'ccw' indicates    //
// the orientation is ccw (1) or not (0).                                    //
//                                                                           //
// 'opm' is a structure contains the parameters of the objective function.   //
// It is needed by the evaluation of the function value.                     //
//                                                                           //
// The return value indicates weather the point is smoothed or not.          //
//                                                                           //
// ASSUMPTION: This routine assumes that all link faces are true faces, i.e, //
// no face has 'dummypoint' as its vertex.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::smoothpoint(point smtpt, arraypool *linkfacelist, int ccw,
                            optparameters *opm)
{
  triface *parytet, *parytet1, swaptet;
  point pa, pb, pc;
  REAL fcent[3], startpt[3], nextpt[3], bestpt[3];
  REAL oldval, minval = 0.0, val;
  REAL maxcosd; // oldang, newang;
  REAL ori, diff;
  int numdirs, iter;
  int i, j, k;

  // Decide the number of moving directions.
  numdirs = (int) linkfacelist->objects;
  if (numdirs > opm->numofsearchdirs) {
    numdirs = opm->numofsearchdirs; // Maximum search directions.
  }

  // Set the initial value.
  if (!opm->max_min_volume) {
    assert(opm->initval >= 0.0);
  }
  opm->imprval = opm->initval;
  iter = 0;

  for (i = 0; i < 3; i++) {
    bestpt[i] = startpt[i] = smtpt[i];
  }

  // Iterate until the obj function is not improved.
  while (1) {

    // Find the best next location.
    oldval = opm->imprval;

    for (i = 0; i < numdirs; i++) {
      // Randomly pick a link face (0 <= k <= objects - i - 1).
      k = (int) randomnation(linkfacelist->objects - i);
      parytet = (triface *) fastlookup(linkfacelist, k);
      // Calculate a new position from 'p' to the center of this face.
      pa = org(*parytet);
      pb = dest(*parytet);
      pc = apex(*parytet);
      for (j = 0; j < 3; j++) {
        fcent[j] = (pa[j] + pb[j] + pc[j]) / 3.0;
      }
      for (j = 0; j < 3; j++) {
        nextpt[j] = startpt[j] + opm->searchstep * (fcent[j] - startpt[j]); 
      }
      // Calculate the largest minimum function value for the new location.
      for (j = 0; j < linkfacelist->objects; j++) {
        parytet = (triface *) fastlookup(linkfacelist, j);
        if (ccw) {
          pa = org(*parytet);
          pb = dest(*parytet);
        } else {
          pb = org(*parytet);
          pa = dest(*parytet);
        }
        pc = apex(*parytet);
        ori = orient3d(pa, pb, pc, nextpt);
        if (ori < 0.0) {
          // Calcuate the objective function value. 
          if (opm->max_min_volume) {
            //val = -ori;
            val = - orient3dfast(pa, pb, pc, nextpt);
          } else if (opm->max_min_aspectratio) {
            val = tetaspectratio(pa, pb, pc, nextpt);
          } else if (opm->min_max_dihedangle) {
            tetalldihedral(pa, pb, pc, nextpt, NULL, &maxcosd, NULL);
            if (maxcosd < -1) maxcosd = -1.0; // Rounding.
            val = maxcosd + 1.0; // Make it be positive. 
          } else {
            // Unknown objective function.
            val = 0.0;
          }  
        } else { // ori >= 0.0;
          // An invalid new tet. 
          // This may happen if the mesh contains inverted elements.
          if (opm->max_min_volume) {
            //val = -ori;
            val = - orient3dfast(pa, pb, pc, nextpt);    
          } else {
            // Discard this point.
            break; // j
          }
        } // if (ori >= 0.0)
        // Stop looping when the object value is not improved.
        if (val <= opm->imprval) {
          break; // j
        } else {
          // Remember the smallest improved value.
          if (j == 0) {
            minval = val;
          } else {
            minval = (val < minval) ? val : minval;
          }
        }
      } // j
      if (j == linkfacelist->objects) {
        // The function value has been improved.
        opm->imprval = minval;
        // Save the new location of the point.
        for (j = 0; j < 3; j++) bestpt[j] = nextpt[j];
      }
      // Swap k-th and (object-i-1)-th entries.
      j = linkfacelist->objects - i - 1;
      parytet  = (triface *) fastlookup(linkfacelist, k);
      parytet1 = (triface *) fastlookup(linkfacelist, j);
      swaptet = *parytet1;
      *parytet1 = *parytet;
      *parytet = swaptet;
    } // i

    diff = opm->imprval - oldval;
    if (diff > 0.0) {
      // Is the function value improved effectively?
      if (opm->max_min_volume) {
        //if ((diff / oldval) < b->epsilon) diff = 0.0;  
      } else if (opm->max_min_aspectratio) {
        if ((diff / oldval) < 1e-3) diff = 0.0;
      } else if (opm->min_max_dihedangle) {
        //oldang = acos(oldval - 1.0);
        //newang = acos(opm->imprval - 1.0);
        //if ((oldang - newang) < 0.00174) diff = 0.0; // about 0.1 degree.
      } else {
        // Unknown objective function.
        assert(0); // Not possible.
      }
    }

    if (diff > 0.0) {
      // Yes, move p to the new location and continue.
      for (j = 0; j < 3; j++) startpt[j] = bestpt[j];
      iter++;
      if ((opm->maxiter > 0) && (iter >= opm->maxiter)) {
        // Maximum smoothing iterations reached.
        break;
      }
    } else {
      break;
    }

  } // while (1)

  if (iter > 0) {
    // The point has been smoothed.
    opm->smthiter = iter; // Remember the number of iterations. 
    // The point has been smoothed. Update it to its new position.
    for (i = 0; i < 3; i++) smtpt[i] = startpt[i];
  }

  return iter;
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// improvequalitysmoothing()    Improve mesh quality by smoothing.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::improvequalitybysmoothing(optparameters *opm)
{
  arraypool *flipqueue, *swapqueue;
  triface *parytet;
  badface *bface, *parybface;
  point *ppt;
  long totalsmtcount, smtcount;
  int smtflag;
  int iter, i, j, k;

  //assert(unflipqueue->objects > 0l);
  flipqueue = new arraypool(sizeof(badface), 10);

  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  totalsmtcount = 0l;
  iter = 0;

  while (flipqueue->objects > 0l) {

    smtcount = 0l;

    if (b->verbose > 1) {
      printf("    Improving mesh quality by smoothing [%d]#:  %ld.\n",
             iter, flipqueue->objects);
    }

    for (k = 0; k < flipqueue->objects; k++) {      
      bface  = (badface *) fastlookup(flipqueue, k);
      if (gettetrahedron(bface->forg, bface->fdest, bface->fapex,
                         bface->foppo, &bface->tt)) {
        // Operate on it if it is not in 'unflipqueue'.
        if (!marktested(bface->tt)) {
          // Here we simply re-compute the quality. Since other smoothing
          //   operation may have moved the vertices of this tet.
          ppt = (point *) & (bface->tt.tet[4]);
          tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                         &bface->key, NULL);
          if (bface->key < cossmtdihed) { // if (maxdd < cosslidihed) {
            // It is a sliver. Try to smooth its vertices.
            smtflag = 0;
            opm->initval = bface->key + 1.0; 
            for (i = 0; (i < 4) && !smtflag; i++) {
              if (pointtype(ppt[i]) == FREEVOLVERTEX) {
                getvertexstar(1, ppt[i], cavetetlist, NULL, NULL);
                opm->searchstep = 0.001; // Search step size
                smtflag = smoothpoint(ppt[i], cavetetlist, 1, opm);
                if (smtflag) {
                  while (opm->smthiter == opm->maxiter) {
                    opm->searchstep *= 10.0; // Increase the step size.
                    opm->initval = opm->imprval;
                    opm->smthiter = 0; // reset
                    smoothpoint(ppt[i], cavetetlist, 1, opm);
                  }
                  // This tet is modifed.
                  smtcount++;
                  if ((opm->imprval - 1.0) < cossmtdihed) {
                    // There are slivers in new tets. Queue them.
                    for (j = 0; j < cavetetlist->objects; j++) {
                      parytet = (triface *) fastlookup(cavetetlist, j);
                      assert(!isdeadtet(*parytet));
                      // Operate it if it is not in 'unflipqueue'.
                      if (!marktested(*parytet)) {
                        // Evaluate its quality.
                        // Re-use ppt, bface->key, bface->cent.
                        ppt = (point *) & (parytet->tet[4]);
                        tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], 
                                       bface->cent, &bface->key, NULL);
                        if (bface->key < cossmtdihed) {
                          // A new sliver. Queue it.
                          marktest(*parytet); // It is in unflipqueue.
                          unflipqueue->newindex((void **) &parybface);
                          parybface->tt = *parytet;
                          parybface->forg = ppt[0];
                          parybface->fdest = ppt[1];
                          parybface->fapex = ppt[2];
                          parybface->foppo = ppt[3];
                          parybface->tt.ver = 11; 
                          parybface->key = 0.0;
                        }
                      }
                    } // j
                  } // if ((opm->imprval - 1.0) < cossmtdihed)
                } // if (smtflag)
                cavetetlist->restart();
              } // if (pointtype(ppt[i]) == FREEVOLVERTEX)
            } // i
            if (!smtflag) {
              // Didn't smooth. Queue it again.
              marktest(bface->tt); // It is in unflipqueue.
              unflipqueue->newindex((void **) &parybface);
              parybface->tt = bface->tt;
              parybface->forg = ppt[0];
              parybface->fdest = ppt[1];
              parybface->fapex = ppt[2];
              parybface->foppo = ppt[3];
              parybface->tt.ver = 11;
              parybface->key = 0.0;
            }
	      } // if (maxdd < cosslidihed)
        } // if (!marktested(...))
      } // if (gettetrahedron(...))
    } // k

    flipqueue->restart();

    // Unmark the tets in unflipqueue.
    for (i = 0; i < unflipqueue->objects; i++) {
      bface  = (badface *) fastlookup(unflipqueue, i);
      unmarktest(bface->tt);
    }

    if (b->verbose > 1) {
      printf("    Smooth %ld points.\n", smtcount);
    }
    totalsmtcount += smtcount;

    if (smtcount == 0l) {
      // No point has been smoothed. 
      break;
    } else {
      iter++;
      if (iter == 2) { //if (iter >= b->optpasses) {
        break;
      }
    }

    // Swap the two flip queues.
    swapqueue = flipqueue;
    flipqueue = unflipqueue;
    unflipqueue = swapqueue;
  } // while

  delete flipqueue;

  return totalsmtcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// splitsliver()    Split a sliver.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::splitsliver(triface *slitet, REAL cosd, int chkencflag)
{
  triface *abtets;
  triface searchtet, spintet, *parytet;
  point pa, pb, steinerpt;
  optparameters opm;
  insertvertexflags ivf;
  REAL smtpt[3], midpt[3];
  int success;
  int t1ver;
  int n, i;

  // 'slitet' is [c,d,a,b], where [c,d] has a big dihedral angle. 
  // Go to the opposite edge [a,b].
  edestoppo(*slitet, searchtet); // [a,b,c,d].

  // Do not split a segment.
  if (issubseg(searchtet)) {
    return 0; 
  }

  // Count the number of tets shared at [a,b].
  // Do not split it if it is a hull edge.
  spintet = searchtet;
  n = 0; 
  while (1) {
    if (ishulltet(spintet)) break;
    n++;
    fnextself(spintet);
    if (spintet.tet == searchtet.tet) break;
  }
  if (ishulltet(spintet)) {
    return 0; // It is a hull edge.
  }
  assert(n >= 3);

  // Get all tets at edge [a,b].
  abtets = new triface[n];
  spintet = searchtet;
  for (i = 0; i < n; i++) {
    abtets[i] = spintet;
    fnextself(spintet);
  }

  // Initialize the list of 2n boundary faces.
  for (i = 0; i < n; i++) {    
    eprev(abtets[i], searchtet);
    esymself(searchtet); // [a,p_i,p_i+1].
    cavetetlist->newindex((void **) &parytet);
    *parytet = searchtet;
    enext(abtets[i], searchtet);
    esymself(searchtet); // [p_i,b,p_i+1].
    cavetetlist->newindex((void **) &parytet);
    *parytet = searchtet;
  }

  // Init the Steiner point at the midpoint of edge [a,b].
  pa = org(abtets[0]);
  pb = dest(abtets[0]);
  for (i = 0; i < 3; i++) {
    smtpt[i] = midpt[i] = 0.5 * (pa[i] + pb[i]);
  }

  // Point smooth options.
  opm.min_max_dihedangle = 1;
  opm.initval = cosd + 1.0; // Initial volume is zero.
  opm.numofsearchdirs = 20;
  opm.searchstep = 0.001;  
  opm.maxiter = 100; // Limit the maximum iterations.

  success = smoothpoint(smtpt, cavetetlist, 1, &opm);

  if (success) {
    while (opm.smthiter == opm.maxiter) {
      // It was relocated and the prescribed maximum iteration reached. 
      // Try to increase the search stepsize.
      opm.searchstep *= 10.0;
      //opm.maxiter = 100; // Limit the maximum iterations.
      opm.initval = opm.imprval;
      opm.smthiter = 0; // Init.
      smoothpoint(smtpt, cavetetlist, 1, &opm);  
    }
  } // if (success)

  cavetetlist->restart();

  if (!success) {
    delete [] abtets;
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

  searchtet = abtets[0]; // No need point location.
  if (b->metric) {
    locate(steinerpt, &searchtet); // For size interpolation.
  }

  delete [] abtets;

  ivf.iloc = (int) INSTAR;
  ivf.chkencflag = chkencflag;
  ivf.assignmeshsize = b->metric; 


  if (insertpoint(steinerpt, &searchtet, NULL, NULL, &ivf)) {
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
// removeslivers()    Remove slivers by adding Steiner points.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::removeslivers(int chkencflag)
{
  arraypool *flipqueue, *swapqueue;
  badface *bface, *parybface;
  triface slitet, *parytet;
  point *ppt;
  REAL cosdd[6], maxcosd;
  long totalsptcount, sptcount;
  int iter, i, j, k;

  //assert(unflipqueue->objects > 0l);
  flipqueue = new arraypool(sizeof(badface), 10);

  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  totalsptcount = 0l;
  iter = 0;

  while ((flipqueue->objects > 0l) && (steinerleft != 0)) {

    sptcount = 0l;

    if (b->verbose > 1) {
      printf("    Splitting bad quality tets [%d]#:  %ld.\n",
             iter, flipqueue->objects);
    }

    for (k = 0; (k < flipqueue->objects) && (steinerleft != 0); k++) {      
      bface  = (badface *) fastlookup(flipqueue, k);
      if (gettetrahedron(bface->forg, bface->fdest, bface->fapex,
                         bface->foppo, &bface->tt)) {
        if ((bface->key == 0) || (bface->tt.ver != 11)) {
          // Here we need to re-compute the quality. Since other smoothing
          //   operation may have moved the vertices of this tet.
          ppt = (point *) & (bface->tt.tet[4]);
          tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                         &bface->key, NULL);
        }
        if (bface->key < cosslidihed) { 
          // It is a sliver. Try to split it.
          slitet.tet = bface->tt.tet;
          //cosdd = bface->cent;
          for (j = 0; j < 6; j++) {
            if (bface->cent[j] < cosslidihed) { 
              // Found a large dihedral angle.
              slitet.ver = edge2ver[j]; // Go to the edge.
              if (splitsliver(&slitet, bface->cent[j], chkencflag)) {
                sptcount++;
                break;
              }
            }
          } // j
          if (j < 6) {
            // A sliver is split. Queue new slivers.
            badtetrahedrons->traversalinit();
            parytet = (triface *) badtetrahedrons->traverse();
            while (parytet != NULL) {
              unmarktest2(*parytet);
              ppt = (point *) & (parytet->tet[4]);
              tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], cosdd, 
                             &maxcosd, NULL);
              if (maxcosd < cosslidihed) {
                // A new sliver. Queue it.
                unflipqueue->newindex((void **) &parybface);
                parybface->forg = ppt[0];
                parybface->fdest = ppt[1];
                parybface->fapex = ppt[2];
                parybface->foppo = ppt[3];
                parybface->tt.tet = parytet->tet;
                parybface->tt.ver = 11;
                parybface->key = maxcosd;
                for (i = 0; i < 6; i++) {
                  parybface->cent[i] = cosdd[i];
                }
              }
              parytet = (triface *) badtetrahedrons->traverse();
            }
            badtetrahedrons->restart();
          } else {
            // Didn't split. Queue it again.
            unflipqueue->newindex((void **) &parybface);
            *parybface = *bface;
          } // if (j == 6)
        } // if (bface->key < cosslidihed)
      } // if (gettetrahedron(...))
    } // k

    flipqueue->restart();

    if (b->verbose > 1) {
      printf("    Split %ld tets.\n", sptcount);
    }
    totalsptcount += sptcount;

    if (sptcount == 0l) {
      // No point has been smoothed. 
      break;
    } else {
      iter++;
      if (iter == 2) { //if (iter >= b->optpasses) {
        break;
      }
    }

    // Swap the two flip queues.
    swapqueue = flipqueue;
    flipqueue = unflipqueue;
    unflipqueue = swapqueue;
  } // while

  delete flipqueue;

  return totalsptcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// optimizemesh()    Optimize mesh for specified objective functions.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::optimizemesh()
{
  badface *parybface;
  triface checktet;
  point *ppt;
  int optpasses;
  optparameters opm;
  REAL ncosdd[6], maxdd;
  long totalremcount, remcount;
  long totalsmtcount, smtcount;
  long totalsptcount, sptcount;
  int chkencflag;
  int iter;
  int n;

  if (!b->quiet) {
    printf("Optimizing mesh...\n");
  }

  optpasses = ((1 << b->optlevel) - 1);

  if (b->verbose) {
    printf("  Optimization level  = %d.\n", b->optlevel);
    printf("  Optimization scheme = %d.\n", b->optscheme);
    printf("  Number of iteration = %d.\n", optpasses);
    printf("  Min_Max dihed angle = %g.\n", b->optmaxdihedral);
  }

  totalsmtcount = totalsptcount = totalremcount = 0l;

  cosmaxdihed = cos(b->optmaxdihedral / 180.0 * PI);
  cossmtdihed = cos(b->optminsmtdihed / 180.0 * PI);
  cosslidihed = cos(b->optminslidihed / 180.0 * PI);

  int attrnum = numelemattrib - 1; 

  // Put all bad tetrahedra into array.
  tetrahedrons->traversalinit();
  checktet.tet = tetrahedrontraverse();
  while (checktet.tet != NULL) {
    if (b->convex) { // -c
      // Skip this tet if it lies in the exterior.
      if (elemattribute(checktet.tet, attrnum) == -1.0) {
        checktet.tet = tetrahedrontraverse();
        continue;
      }
    }
    ppt = (point *) & (checktet.tet[4]);
    tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], ncosdd, &maxdd, NULL);
    if (maxdd < cosmaxdihed) {
      // There are bad dihedral angles in this tet.
      unflipqueue->newindex((void **) &parybface); 
      parybface->tt.tet = checktet.tet;
      parybface->tt.ver = 11;
      parybface->forg = ppt[0];
      parybface->fdest = ppt[1];
      parybface->fapex = ppt[2];
      parybface->foppo = ppt[3];
      parybface->key = maxdd;
      for (n = 0; n < 6; n++) {
        parybface->cent[n] = ncosdd[n];
      }
    }
    checktet.tet = tetrahedrontraverse();
  }

  totalremcount = improvequalitybyflips();

  if ((unflipqueue->objects > 0l) && 
      ((b->optscheme & 2) || (b->optscheme & 4))) {
    // The pool is only used by removeslivers().
    badtetrahedrons = new memorypool(sizeof(triface), b->tetrahedraperblock,
                                     sizeof(void *), 0);

    // Smoothing options.
    opm.min_max_dihedangle = 1;
    opm.numofsearchdirs = 10;
    // opm.searchstep = 0.001;  
    opm.maxiter = 30; // Limit the maximum iterations.
    //opm.checkencflag = 4; // Queue affected tets after smoothing.
    chkencflag = 4; // Queue affected tets after splitting a sliver.
    iter = 0;

    while (iter < optpasses) {
      smtcount = sptcount = remcount = 0l;
      if (b->optscheme & 2) {
        smtcount += improvequalitybysmoothing(&opm);
        totalsmtcount += smtcount;
        if (smtcount > 0l) {
          remcount = improvequalitybyflips();
          totalremcount += remcount;
        }
      }
      if (unflipqueue->objects > 0l) {
        if (b->optscheme & 4) {
          sptcount += removeslivers(chkencflag);
          totalsptcount += sptcount;
          if (sptcount > 0l) {
            remcount = improvequalitybyflips();
            totalremcount += remcount;
          }
        }
      }
      if (unflipqueue->objects > 0l) {
        if (remcount > 0l) {
          iter++;
        } else {
          break;
        }
      } else {
        break;
      }
    } // while (iter)

    delete badtetrahedrons;

  }

  if (unflipqueue->objects > 0l) {
    if (b->verbose > 1) {
      printf("    %ld bad tets remained.\n", unflipqueue->objects);
    }
    unflipqueue->restart();
  }

  if (b->verbose) {
    if (totalremcount > 0l) {
      printf("  Removed %ld edges.\n", totalremcount);
    }
    if (totalsmtcount > 0l) {
      printf("  Smoothed %ld points.\n", totalsmtcount);
    }
    if (totalsptcount > 0l) {
      printf("  Split %ld slivers.\n", totalsptcount);
    }
  }
}

////                                                                       ////
////                                                                       ////
//// optimize_cxx /////////////////////////////////////////////////////////////

