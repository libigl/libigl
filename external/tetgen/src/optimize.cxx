#include "../tetgen.h"
//// optimize_cxx /////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// recoverdelaunay()    Recovery the locally Delaunay property.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::recoverdelaunay()
{
  arraypool *flipqueue, *nextflipqueue, *swapqueue;
  badface *bface, *parybface;
  triface tetloop, neightet, *parytet;
  point *ppt;
  flipconstraints fc;
  int i, j;

  if (!b->quiet) {
    printf("Recovering Delaunayness...\n");
  }

  //if (b->verbose) {
  //  printf("  max_flipstarsize = %d.\n", b->optmaxflipstarsize);
  //  printf("  max_fliplinklevel = %d.\n", b->delmaxfliplevel);
  //}

  calc_tetprism_vol = 1;
  tetprism_vol_sum = 0.0; // Initialize it.

  assert(flipstack == NULL);
  assert(unflipqueue->objects == 0l);

  // Put all interior faces of the mesh into 'flipstack'.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != NULL) {
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      // Avoid queue a face twice.
      fsym(tetloop, neightet);
      if (!ishulltet(neightet)) {
        if (!facemarked(neightet)) {
          flippush(flipstack, &tetloop);
        }
      }
    }
    ppt = (point *) &(tetloop.tet[4]);
    tetprism_vol_sum += tetprismvol(ppt[0], ppt[1], ppt[2], ppt[3]);
    tetloop.tet = tetrahedrontraverse();
  }

  if (b->verbose) {
    printf("  Initial obj = %.17g\n", tetprism_vol_sum);
  }

  if (b->verbose > 1) {
    printf("    Recover Delaunay [Lawson] : %ld\n", flippool->items);
  }
  assert(unflipqueue->objects == 0l);

  // First only use the basic Lawson's flip.
  lawsonflip3d(NULL, 4, 0, 0, 1);

  if (b->verbose > 1) {
    printf("    New obj = %.17g\n", tetprism_vol_sum);
  }

  if (unflipqueue->objects == 0l) {
    // The mesh is Delaunay.
    return;
  }

  // Set the common options.
  fc.remove_ndelaunay_edge = 1;
  fc.unflip = 1; // Unflip if the edge is not flipped.
  fc.collectnewtets = 1;

  autofliplinklevel = 1; // Init value.
  b->fliplinklevel = -1;

  // For efficiency reason, we limit the maximium size of the edge star.
  // 'b->optmaxflipstarsize' is set by -OOOOO (5 Os), default is 10.
  int bakmaxflipstarsize = b->flipstarsize;
  b->flipstarsize = 10; //b->optmaxflipstarsize;

  flipqueue = new arraypool(sizeof(badface), 10);
  nextflipqueue = new arraypool(sizeof(badface), 10);


  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  while (flipqueue->objects > 0l) {

    while (flipqueue->objects > 0l) {

      if (b->verbose > 1) {
        printf("    Recover Delaunay [level = %2d] #:  %ld.\n",
               autofliplinklevel, flipqueue->objects);
      }

      for (i = 0; i < flipqueue->objects; i++) {
        bface  = (badface *) fastlookup(flipqueue, i);
        if (getedge(bface->forg, bface->fdest, &bface->tt)) {
          // Remember the the objective value (volume of all tetprisms).
          fc.bak_tetprism_vol = tetprism_vol_sum; 
          if (removeedgebyflips(&(bface->tt), &fc) == 2) {
            if (b->verbose > 2) {
              printf("      Decreased quantity: %.17g.\n", 
                     fc.bak_tetprism_vol - tetprism_vol_sum);
            }
            // Queue new faces for flips.
            for (j = 0; j < cavetetlist->objects; j++) {
              parytet = (triface *) fastlookup(cavetetlist, j);
              // A queued new tet may be dead.
              if (!isdeadtet(*parytet)) {
                for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
                  // Avoid queue a face twice.
                  fsym(*parytet, neightet);
                  if (!facemarked(neightet)) {
                    flippush(flipstack, parytet);
                  }
                } // parytet->ver
              }
            } // j
            cavetetlist->restart();
            // Remove locally non-Delaunay faces. New non-Delaunay edges
            //   may be found. They are saved in 'unflipqueue'.
            lawsonflip3d(NULL, 4, 0, 0, 1);
          } else {
            // Unable to remove this edge. Save it.
            nextflipqueue->newindex((void **) &parybface);
            *parybface = *bface;
          }
        }
      } // i

      flipqueue->restart();

      // Swap the two flip queues.
      swapqueue = flipqueue;
      flipqueue = unflipqueue;
      unflipqueue = swapqueue;
    } // while (flipqueue->objects > 0l)

    if (b->verbose > 1) {
      printf("    New obj = %.17g.\n", tetprism_vol_sum);
    }

    // Swap the two flip queues.
    swapqueue = flipqueue;
    flipqueue = nextflipqueue;
    nextflipqueue = swapqueue;

    if (flipqueue->objects > 0l) {
      // 'b->delmaxfliplevel' is set by -OOOO, default is 1.
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

  b->flipstarsize = bakmaxflipstarsize;

  delete nextflipqueue;
  delete flipqueue;

  calc_tetprism_vol = 0;

  if (b->verbose) {
    printf("  Final  obj  = %.17g\n", tetprism_vol_sum);
  }
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

  if (b->verbose > 2) {
    printf("      Get tet [%d,%d,%d,%d].\n", pointmark(pa), pointmark(pb),
           pointmark(pc), pointmark(pd));
  }

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
                           &maxdd, NULL);
            bface->forg = ppt[0];
            bface->fdest = ppt[1];
            bface->fapex = ppt[2];
            bface->foppo = ppt[3];
            bface->tt.ver = 11;
          }
          cosdd = bface->cent;
          remflag = 0;
          for (i = 0; (i < 6) && !remflag; i++) {
            if (cosdd[i] < cosmaxdihed) {
              // Found a large dihedral angle.
              bface->tt.ver = edge2ver[i]; // Go to the edge.
              if (b->verbose > 2) {
                printf("      Found a large angle [%d,%d,%d,%d] (%g).\n", 
                       pointmark(org(bface->tt)), pointmark(dest(bface->tt)),
                       pointmark(apex(bface->tt)), pointmark(oppo(bface->tt)),
                       acos(cosdd[i]) / PI * 180.0);
              }
              fc.cosdihed_in = cosdd[i];
              fc.cosdihed_out = 0.0; // 90 degree.
              n = removeedgebyflips(&(bface->tt), &fc);
              if (n == 2) {
                // Edge is flipped.
                if (b->verbose > 2) {
                  printf("      Reduced a large angle to %g degree.\n",
                         acos(fc.cosdihed_out) / PI * 180.0);
                }
                remflag = 1;
                if (fc.cosdihed_out < cosmaxdihed) {
                  // Queue new bad tets for further improvements.
                  for (j = 0; j < cavetetlist->objects; j++) {
                    parytet = (triface *) fastlookup(cavetetlist, j);
                    if (!isdeadtet(*parytet)) {
                      ppt = (point *) & (parytet->tet[4]);
                      //if (!marktest2ed(*parytet)) {
                      assert(!marktest2ed(*parytet)); // SELF_CHECK
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
                      } // if (ppt[3] != dummypoint) { 
		      //}
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
// 'of' is a structure contains the parameters of the objective function. It //
// is needed by the evaluation of the function value.                        //
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

  if (b->verbose > 2) {
    printf("      Smooth a point: %ld faces.\n", linkfacelist->objects);
    if (opm->min_max_dihedangle) {
      printf("      Init value = %g (degree).\n", 
             acos(opm->initval - 1.0) / PI * 180.0);
    } else {
      printf("      Init value = %g.\n", opm->initval);
    }
  }

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
            val = -ori;
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
          if (opm->max_min_volume) {
            val = -ori;    
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
        assert(minval > opm->imprval);          
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
    // The point has been smooothed.
    opm->smthiter = iter; // Remember the number of iterations.    
    if (b->verbose > 2) {
      printf("      Smoothed: %d iterations.\n", iter);
      if (opm->min_max_dihedangle) {
        printf("      Fina value = %g (degree).\n", 
               acos(opm->imprval - 1.0) / PI * 180.0);
      } else {
        printf("      Fina value = %g.\n", opm->imprval);
      }
    }
    // The point has been smoothed. Update it to its new position.
    for (i = 0; i < 3; i++) smtpt[i] = startpt[i];

    if (opm->flipflag) {
      // Push all affected faces into 'flipstack'.
      triface starttet, neightet;
      for (i = 0; i < linkfacelist->objects; i++) {
        parytet = (triface *) fastlookup(linkfacelist, i);
        starttet = *parytet;
        for (starttet.ver = 0; starttet.ver < 4; starttet.ver++) {
          fsym(starttet, neightet);
          if (!infected(neightet)) {
            flippush(flipstack, &starttet);
          }
        }
        infect(*parytet);
      }
      for (i = 0; i < linkfacelist->objects; i++) {
        parytet = (triface *) fastlookup(linkfacelist, i);
        uninfect(*parytet);
      }
    } else if (opm->checkencflag) {
      // Push all affected tets into pool.
      badface *bface;
      for (i = 0; i < linkfacelist->objects; i++) {
        parytet = (triface *) fastlookup(linkfacelist, i);
        if (!marktest2ed(*parytet)) {
          marktest2(*parytet); // Only queue it once.
          bface = (badface *) badtetrahedrons->alloc();
          bface->tt = *parytet;
          bface->forg = org(bface->tt);
        }
      }
    }
  } else {
    if (b->verbose > 2) {
      printf("      Not smoothed.\n");
    }
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
  badface *bface, *parybface;
  point *ppt;
  long totalsmtcount, smtcount;
  int smtflag;
  int iter, i, k;

  //assert(unflipqueue->objects > 0l);
  flipqueue = new arraypool(sizeof(badface), 10);

  totalsmtcount = 0l;

  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  iter = 0;

  while (flipqueue->objects > 0l) {

    smtcount = 0l;

    //while (flipqueue->objects > 0l) {
    if (b->verbose > 1) {
      printf("    Improving mesh qualiy by smoothing [%d]#:  %ld.\n",
             iter, flipqueue->objects);
    }

    for (k = 0; k < flipqueue->objects; k++) {      
      bface  = (badface *) fastlookup(flipqueue, k);
      if (gettetrahedron(bface->forg, bface->fdest, bface->fapex,
                         bface->foppo, &bface->tt)) {
        if (!marktested(bface->tt)) {
          // Here we simply re-compute the quality. Since other smoothing
          //   operation may have moved the vertices of this tet.
          ppt = (point *) & (bface->tt.tet[4]);
          tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                         &bface->key, NULL);
          if (bface->key < cossmtdihed) { // if (maxdd < cosslidihed) {
            // It is a sliver. Try to smooth its vertices.
            smtflag = 0;
            //if (opm->min_max_dihedangle) {
            opm->initval = bface->key + 1.0;          
            //opm->checkencflag = 4; // Queue affected tets.            
            //}
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
                  smtcount++;
                }
                cavetetlist->restart();
              }
            } // i
            if (smtflag) {
              // This tet is modifed. 
              smtcount++;
              if ((opm->imprval - 1.0) < cossmtdihed) {
                // Queue new slivers.
                badtetrahedrons->traversalinit();
                bface = badfacetraverse(badtetrahedrons);
                while (bface != NULL) {
                  assert(!isdeadtet(bface->tt));
                  assert(marktest2ed(bface->tt));
                  unmarktest2(bface->tt);
                  if (!marktested(bface->tt)) {
                    ppt = (point *) & (bface->tt.tet[4]);
                    tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                                   &(bface->key), NULL);
                    if (bface->key < cossmtdihed) {
                      // A new sliver. Queue it.
                      marktest(bface->tt); // It is in unflipqueue.
                      bface->forg = ppt[0]; 
                      bface->fdest = ppt[1];
                      bface->fapex = ppt[2];
                      bface->foppo = ppt[3];
                      bface->tt.ver = 11;                    
                      unflipqueue->newindex((void **) &parybface);
                      *parybface = *bface;
                    }
                  }
                  bface = badfacetraverse(badtetrahedrons);
                }
              } else {
                // No new slivers. Only unmark the queued tets.
                badtetrahedrons->traversalinit();
                bface = badfacetraverse(badtetrahedrons);
                while (bface != NULL) {
                  assert(!isdeadtet(bface->tt));
                  assert(marktest2ed(bface->tt));
                  unmarktest2(bface->tt);
                  bface = badfacetraverse(badtetrahedrons);
                }
              }
              badtetrahedrons->restart();
            } else {
              // Didn't smooth. Queue it again.
              // Adjust the vertices for flipping.
              marktest(bface->tt); // It is in unflipqueue.
              bface->forg = ppt[0]; 
              bface->fdest = ppt[1];
              bface->fapex = ppt[2];
              bface->foppo = ppt[3];
              bface->tt.ver = 11;
              unflipqueue->newindex((void **) &parybface);
              *parybface = *bface;
            }
	  } // if (maxdd < cosslidihed)
        } // if (!marktested(...))
      } // gettetrahedron(...)
    } // k

    flipqueue->restart();

    // } // while

    // Unmark the tets in unflipqueue.
    for (i = 0; i < unflipqueue->objects; i++) {
      bface  = (badface *) fastlookup(unflipqueue, i);
      assert(!isdeadtet(bface->tt));
      assert(marktested(bface->tt));
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
  face checkseg;
  point pa, pb, steinerpt;
  optparameters opm;
  insertvertexflags ivf;
  REAL smtpt[3], midpt[3];
  int success;
  int loc;
  int n, i;

  // 'slitet' is [c,d,a,b], where [c,d] has a big dihedral angle. 
  // Go to the opposite edge [a,b].
  eprev(*slitet, searchtet);
  esymself(searchtet);
  enextself(searchtet); // [a,b,c,d].

  // Do not split a segment.
  tsspivot1(searchtet, checkseg);
  if (checkseg.sh != NULL) {
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
    if (b->verbose > 2) {
      printf("      Unable to relocate the initial point.\n");
    }
    delete [] abtets;
    return 0;
  }


  if (steinerleft == 0) {
    // The desired number of Steiner points is reached.
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
  ivf.assignmeshsize = 0; 

  loc = insertvertex(steinerpt, &searchtet, NULL, NULL, &ivf);

  if (loc == (int) INSTAR) {
    // The vertex has been inserted.
    st_volref_count++; //st_inpoly_count++;
    if (steinerleft > 0) steinerleft--;
    return 1;
  } else {
    // The Steiner point is too close to an existing vertex. Reject it.
    pointdealloc(steinerpt);
    return 0;
  }

  delete [] abtets;
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
  point *ppt;
  REAL *cosdd;
  long totalsptcount, sptcount;
  int iter, j, k;

  //assert(unflipqueue->objects > 0l);
  flipqueue = new arraypool(sizeof(badface), 10);

  totalsptcount = 0l;

  // Swap the two flip queues.
  swapqueue = flipqueue;
  flipqueue = unflipqueue;
  unflipqueue = swapqueue;

  iter = 0;

  while (flipqueue->objects > 0l) {

    sptcount = 0l;

    if (b->verbose > 1) {
      printf("    Splitting bad quality tets [%d]#:  %ld.\n",
             iter, flipqueue->objects);
    }

    for (k = 0; k < flipqueue->objects; k++) {      
      bface  = (badface *) fastlookup(flipqueue, k);
      if (gettetrahedron(bface->forg, bface->fdest, bface->fapex,
                         bface->foppo, &bface->tt)) {
        //if (!marktested(bface->tt)) {
          // Here we simply re-compute the quality. Since other smoothing
          //   operation may have moved the vertices of this tet.
          ppt = (point *) & (bface->tt.tet[4]);
          tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                         &bface->key, NULL);
          if (bface->key < cosslidihed) { 
            // It is a sliver. Try to split it.
            cosdd = bface->cent;
            for (j = 0; j < 6; j++) {
              if (cosdd[j] < cosslidihed) { 
                // Found a large dihedral angle.
                bface->tt.ver = edge2ver[j]; // Go to the edge.
                if (b->verbose > 2) {
                  printf("      Found a bad tet [%d,%d,%d,%d] (%g).\n", 
                         pointmark(org(bface->tt)), pointmark(dest(bface->tt)),
                         pointmark(apex(bface->tt)), pointmark(oppo(bface->tt)),
                         acos(cosdd[j]) / PI * 180.0);
                }
                if (splitsliver(&(bface->tt), cosdd[j], chkencflag)) {
                  sptcount++;
                  break;
                }
              }
            } // j
            if (j < 6) {
              // A sliver is split. Queue new slivers.
              badtetrahedrons->traversalinit();
              bface = badfacetraverse(badtetrahedrons);
              while (bface != NULL) {
                assert(!isdeadtet(bface->tt));
                assert(marktest2ed(bface->tt));
                unmarktest2(bface->tt);
                ppt = (point *) & (bface->tt.tet[4]);
                tetalldihedral(ppt[0], ppt[1], ppt[2], ppt[3], bface->cent, 
                               &(bface->key), NULL);
                if (bface->key < cosslidihed) {
                  // A new sliver. Queue it.
                  //marktest(bface->tt); // It is in unflipqueue.
                  bface->forg = ppt[0]; 
                  bface->fdest = ppt[1];
                  bface->fapex = ppt[2];
                  bface->foppo = ppt[3];
                  bface->tt.ver = 11;                    
                  unflipqueue->newindex((void **) &parybface);
                  *parybface = *bface;
                }
                bface = badfacetraverse(badtetrahedrons);
              }
              badtetrahedrons->restart();
            } else {
              // Didn't split. Queue it again.
              // Adjust the vertices for flipping.
              //marktest(bface->tt); // It is in unflipqueue.
              bface->forg = ppt[0]; 
              bface->fdest = ppt[1];
              bface->fapex = ppt[2];
              bface->foppo = ppt[3];
              bface->tt.ver = 11;
              unflipqueue->newindex((void **) &parybface);
              *parybface = *bface;
            } // if (j == 6)
          } // if (bface->key < cosslidihed)
	// } // if (!marktested(bface->tt))
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

  if (b->verbose) {
    printf("  Optimization level  = %d.\n", b->optlevel);
    printf("  Optimization scheme = %d.\n", b->optscheme);
    printf("  Min_Max dihed angle = %g.\n", b->optmaxdihedral);
  }

  optpasses = ((1 << b->optlevel) - 1);

  totalsmtcount = totalsptcount = totalremcount = 0l;

  cosmaxdihed = cos(b->optmaxdihedral / 180.0 * PI);
  cossmtdihed = cos(b->optminsmtdihed / 180.0 * PI);
  cosslidihed = cos(b->optminslidihed / 180.0 * PI);

  // Put all bad tetrahedra into array.
  tetrahedrons->traversalinit();
  checktet.tet = tetrahedrontraverse();
  while (checktet.tet != NULL) {
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

    badtetrahedrons = new memorypool(sizeof(badface), b->tetrahedraperblock,
                                     memorypool::POINTER, 0);

    // Smoothing options.
    opm.min_max_dihedangle = 1;
    opm.numofsearchdirs = 10;
    // opm.searchstep = 0.001;  
    opm.maxiter = 30; // Limit the maximum iterations.
    opm.checkencflag = 4; // Queue affected tets after smoothing.
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
      printf("  Removed %ld bad tets.\n", totalremcount);
    }
    if (totalsmtcount > 0l) {
      printf("  Smoothed %ld points.\n", totalsmtcount);
    }
    if (totalsptcount > 0l) {
      printf("  Split %ld bad tets.\n", totalsptcount);
    }
  }
}

////                                                                       ////
////                                                                       ////
//// optimize_cxx /////////////////////////////////////////////////////////////

