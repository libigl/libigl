#include "../tetgen.h"
//// reconstruct_cxx //////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// carveholes()    Remove tetrahedra not in the mesh domain.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


void tetgenmesh::carveholes()
{
  arraypool *tetarray, *hullarray;
  triface tetloop, neightet, hulltet, *parytet, *parytet1;
  triface openface, casface;
  triface *regiontets;
  face checksh, casingout, casingin, *parysh;
  face checkseg;
  point *ppt, pa, pb, pc, *parypt;
  enum locateresult loc;
  REAL volume;
  long delsegcount, delvertcount, delsteinercount;
  int regioncount;
  int attrnum, attr, maxattr;
  int remflag;
  int i, j, k;

  tetrahedron ptr;
  shellface sptr;

  if (!b->quiet) {
    printf("Removing exterior tetrahedra ...\n");
  }


  // Initialize the pool of exterior tets.
  tetarray = new arraypool(sizeof(triface), 10);
  hullarray = new arraypool(sizeof(triface), 10);

  regiontets = NULL;
  regioncount = 0;
  maxattr = 0; // Choose a small number here.
  //attrnum = in->numberoftetrahedronattributes;
  attrnum = numelemattrib - (b->regionattrib > 0); 
  // Comment: The element region marker is at the end of the list of
  //   the element attributes.

  // Mark as infected any unprotected hull tets.
  tetrahedrons->traversalinit();
  tetloop.ver = 11; // The face opposite to dummypoint.
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    if ((point) tetloop.tet[7] == dummypoint) {
      // Is this side protected by a subface?
      tspivot(tetloop, checksh);
      if (checksh.sh == NULL) {
        infect(tetloop);
        tetarray->newindex((void **) &parytet);
        *parytet = tetloop;
        hullsize--;
        // Add the adjacent tet (not a hull tet) as well.
        // tetloop's face number is 11 & 3 = 3.
        decode(tetloop.tet[3], neightet);
        if (!infected(neightet)) {
          infect(neightet);
          tetarray->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
    }
    tetloop.tet = alltetrahedrontraverse();
  }

  if (in->numberofholes > 0) {
    // Mark as infected any tets inside volume holes.
    for (i = 0; i < 3 * in->numberofholes; i += 3) {
      // Search a tet containing the i-th hole point.
      neightet.tet = NULL;
      randomsample(&(in->holelist[i]), &neightet);
      loc = locate(&(in->holelist[i]), &neightet, 0); 
      if (loc != OUTSIDE) {
        // The tet 'neightet' contain this point.
        if (!infected(neightet)) {
          infect(neightet);
          tetarray->newindex((void **) &parytet);
          *parytet = neightet;
          // Add its adjacent tet if it is not protected.
          tspivot(neightet, checksh);
          if (checksh.sh == NULL) {
            decode(neightet.tet[neightet.ver & 3], tetloop);
            if (!infected(tetloop)) {
              infect(tetloop);
              tetarray->newindex((void **) &parytet);
              *parytet = tetloop;
            }
          } else {
            // It is protected. Check if its adjacent tet is a hull tet.
            decode(neightet.tet[neightet.ver & 3], tetloop);
            if (!infected(tetloop)) {
              if (ishulltet(tetloop)) {
                // It is hull tet, add it into the list. Moreover, the subface
                //   is dead, i.e., both sides are in exterior.
                infect(tetloop);
                tetarray->newindex((void **) &parytet);
                *parytet = tetloop;
                stdissolve(checksh);
                assert(!sinfected(checksh));
                //if (!sinfected(checksh)) {
                  sinfect(checksh); // Only queue it once.
                  subfacstack->newindex((void **) &parysh);
                  *parysh = checksh;
	        //}
                hullsize--;              
              }
            } else {
              // Both sides of this subface are in exterior.
              stdissolve(checksh);
              assert(!sinfected(checksh));
              //if (!sinfected(checksh)) {
                sinfect(checksh); // Only queue it once.
                subfacstack->newindex((void **) &parysh);
                *parysh = checksh;
	      //}
            }
          }
        } // if (!infected(neightet))
      } else {
        // A hole point locates outside of the convex hull.
        if (!b->quiet) {
          printf("Warning:  The %d-th hole point ", i/3 + 1);
          printf("lies outside the convex hull.\n");
        }
      }
    } // i
  }

  if (b->regionattrib && (in->numberofregions > 0)) { // If has -A option.
    // Record the tetrahedra that contains the region points for assigning
    //   region attributes after the holes have been carved.
    regiontets = new triface[in->numberofregions];
    // Mark as marktested any tetrahedra inside volume regions.
    for (i = 0; i < 5 * in->numberofregions; i += 5) {
      // Search a tet containing the i-th region point.
      neightet.tet = NULL;
      randomsample(&(in->regionlist[i]), &neightet);
      loc = locate(&(in->regionlist[i]), &neightet, 0); 
      if (loc != OUTSIDE) {
        regiontets[i/5] = neightet;
        if ((int) in->regionlist[i + 3] > maxattr) {
          maxattr = (int) in->regionlist[i + 3];
        }
      } else {
        if (!b->quiet) {
          printf("Warning:  The %d-th region point ", i/5+1);
          printf("lies outside the convex hull.\n");
        }
        regiontets[i/5].tet = NULL;
      }
    }
  }


  // Find and infect all exterior tets (in concave place and in holes).
  for (i = 0; i < tetarray->objects; i++) {
    parytet = (triface *) fastlookup(tetarray, i);
    // Check its three neighbors if it is not a hull tet.
    if ((point) parytet->tet[7] != dummypoint) {
      j = (parytet->ver & 3); // j is the current face number.
      // Check the neighbors of the other three faces.
      for (k = 0, j++; k < 3; k++, j++) {
        decode(parytet->tet[j % 4], neightet); // neightet may be a hull tet.
        if (!infected(neightet)) {
          // Is neightet protected by a subface.
          tspivot(neightet, checksh);
          if (checksh.sh == NULL) {
            // Not proected. Add it into the list.
            // It should not be a hull tet. Since all unproected hull tets
            //   should have already been added into the list.
            assert(!ishulltet(neightet)); // SELF_CHECK
            infect(neightet);
            tetarray->newindex((void **) &parytet1);
            *parytet1 = neightet;
          } else {
            // It is protected. However, if neightet is a hull tet, it is
            //   also an exterior tet. Moverover, the subface is dead, i.e.,
            //   both sides of it are exterior.
            if ((point) neightet.tet[7] == dummypoint) {
              infect(neightet);
              tetarray->newindex((void **) &parytet1);
              *parytet1 = neightet;
              // Both sides of this subface are exterior.
              stdissolve(checksh);
              // Queue this subface (to be deleted later).
              assert(!sinfected(checksh));
              //if (!sinfected(checksh)) {
                sinfect(checksh); // Only queue it once.
                subfacstack->newindex((void **) &parysh);
                *parysh = checksh;
	      //}
              hullsize--;
            }
          }
        } else {
          // Both sides of this face are in exterior.
          // Check if there is a subface.
          tspivot(neightet, checksh);
          if (checksh.sh != NULL) {
            if (!sinfected(checksh)) {
              sinfect(checksh); // Only queue it once.
              subfacstack->newindex((void **) &parysh);
              *parysh = checksh;
	    }
          }
        }
      } // j, k
    }
  } // i

  if (b->regionattrib && (in->numberofregions > 0)) {
    // Re-check saved region tets to see if they lie outside.
    for (i = 0; i < in->numberofregions; i++) {
      if (infected(regiontets[i])) {
        if (b->verbose) {
          printf("Warning:  The %d-th region point ", i+1);
          printf("lies in the exterior of the domain.\n");
        }
        regiontets[i].tet = NULL;
      }
    }
  }


if (!b->convex) {

  // Create new hull tets. 
  // Update point-to-tet map, segment-to-tet map, and subface-to-tet map.
  for (i = 0; i < tetarray->objects; i++) {
    parytet = (triface *) fastlookup(tetarray, i);
    if ((point) parytet->tet[7] != dummypoint) {
      // We must check all four adjacent tets.
      for (j = 0; j < 4; j++) {
        decode(parytet->tet[j], tetloop);
        if (!infected(tetloop)) {
          // This face becomes a hull face.
          tspivot(tetloop, checksh);
          assert(checksh.sh != NULL); // SELF_CHECK
          maketetrahedron(&hulltet);
          pa = org(tetloop);
          pb = dest(tetloop);
          pc = apex(tetloop);
          setvertices(hulltet, pb, pa, pc, dummypoint);
          bond(tetloop, hulltet);
          // Update the subface-to-tet map.
          sesymself(checksh);
          tsbond(hulltet, checksh);
          // Update the segment-to-tet map.
          for (k = 0; k < 3; k++) {
            tsspivot1(tetloop, checkseg);
            if (checkseg.sh != NULL) {
              tssbond1(hulltet, checkseg);
              sstbond1(checkseg, hulltet);
            }
            enextself(tetloop);
            eprevself(hulltet);
          }
          // Update the point-to-tet map.
          ptr = encode(tetloop);
          setpoint2tet(pa, ptr);
          setpoint2tet(pb, ptr);
          setpoint2tet(pc, ptr);
          // Save this hull tet in list.
          hullarray->newindex((void **) &parytet1);
          *parytet1 = hulltet;
        }
      } // j
    } else {
      // It is a hull tet. Clear the adjacent hull tets' connections to it.
      // Our data structure ensures that the 3rd face opposites dummypoint.
      for (j = 0; j < 3; j++) {
        decode(parytet->tet[j], neightet);
        if (neightet.tet != NULL) {
          assert(ishulltet(neightet));
          if (!infected(neightet)) {
            neightet.tet[neightet.ver & 3] = NULL;
          }
        }
      } // j
    }
  } // i

  // Update the hull size.
  hullsize += hullarray->objects;

  // Remove all exterior tetrahedra (including infected hull tets).
  for (i = 0; i < tetarray->objects; i++) {
    parytet = (triface *) fastlookup(tetarray, i);
    tetrahedrondealloc(parytet->tet);
  } // i

  tetarray->restart(); 


  if (subfacstack->objects > 0) {
    // Remove all subfaces which do not attach to any tetrahedron.
    //   Segments which are not attached to any subfaces and tets
    //   are deleted too.
    delsegcount = 0;
    for (i = 0; i < subfacstack->objects; i++) {
      parysh = (face *) fastlookup(subfacstack, i);
      if (i == 0) {
        if (b->verbose) {
          printf("Warning:  Removing an open face (%d, %d, %d)\n",
                 pointmark(sorg(*parysh)), pointmark(sdest(*parysh)),
                 pointmark(sapex(*parysh)));
        }
      }
      // Dissolve this subface from face links.
      for (j = 0; j < 3; j++) {         
        spivot(*parysh, casingout);
        sspivot(*parysh, checkseg);
        if (casingout.sh != NULL) {
          casingin = casingout;
          while (1) {
            spivot(casingin, checksh);
            if (checksh.sh == parysh->sh) break;
            casingin = checksh;
          }
          if (casingin.sh != casingout.sh) {
            // Update the link: ... -> casingin -> casingout ->...
            sbond1(casingin, casingout);
          } else {
            // Only one subface at this edge is left.
            sdissolve(casingout);
          }
          if (checkseg.sh != NULL) {
            // Make sure the segment does not connect to a dead one.
            ssbond(casingout, checkseg);
          }
        } else {
          if (checkseg.sh != NULL) {
            // The segment is also dead.
            if (delsegcount == 0) {
              if (b->verbose) {
                printf("Warning:  Removing a dangling segment (%d, %d)\n",
                       pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
              }
            }
            shellfacedealloc(subsegs, checkseg.sh);
            delsegcount++;
          }
        }
        senextself(*parysh);
      } // j
      // Delete this subface.
      shellfacedealloc(subfaces, parysh->sh);
    } // i
    if (b->verbose) {
      printf("  Deleted %ld subfaces.\n", subfacstack->objects);
      if (delsegcount > 0) {
        printf("  Deleted %ld segments.\n", delsegcount);
      }
    }
    subfacstack->restart();
  }


  // Some vertices may be not belong to any tet. Mark them.
  delvertcount = unuverts;
  delsteinercount = 0l;
  points->traversalinit();
  pa = pointtraverse();
  while (pa != NULL) {
    if (pointtype(pa) != UNUSEDVERTEX) {
      remflag = 0;
      decode(point2tet(pa), neightet);
      if ((neightet.tet == NULL) || (neightet.tet[4] == NULL)) {
        remflag = 1; // It's a dead tet.
      } else {
        // Check if this tet contains pa.
        ppt = (point *) &(neightet.tet[4]);
        if (!((ppt[0] == pa) || (ppt[1] == pa) || 
              (ppt[2] == pa) || (ppt[3] == pa))) {
          remflag = 1; // It's a wrong pointer.
        }
      }
      if (remflag) {
        // Found an exterior vertex.
        if (pointmark(pa) > 
              (in->numberofpoints - (in->firstnumber ? 0 : 1))) {
          if (pointtype(pa) == FREESEGVERTEX) {
            st_segref_count--;
          } else if (pointtype(pa) == FREEFACETVERTEX) {
            st_facref_count--;
          } else {
            assert(pointtype(pa) == FREEVOLVERTEX);
            st_volref_count--; //st_inpoly_count--;
          }
          delsteinercount++; // A Steiner point.
          if (steinerleft > 0) steinerleft++;
        }
        setpointtype(pa, UNUSEDVERTEX);
        unuverts++;
      } else {
        // This vertex survived. 
        if (b->nobisect && (b->nobisect_param > 1)) { // -Y2
          // Queue it if it is a Steiner point.
          if ((pointtype(pa) == FREESEGVERTEX) ||
              (pointtype(pa) == FREEFACETVERTEX) ||
              (pointtype(pa) == FREEVOLVERTEX)) {
            subvertstack->newindex((void **) &parypt);
            *parypt = pa;
          }
        }
      }
    }
    pa = pointtraverse();
  }

  if (b->verbose) {
    if (unuverts > delvertcount) {
      if (delsteinercount > 0l) {
        if (unuverts > (delvertcount + delsteinercount)) {
          printf("  Removed %ld exterior input vertices.\n", 
                 unuverts - delvertcount - delsteinercount);
        }
        printf("  Removed %ld exterior Steiner vertices.\n", delsteinercount);
      } else {
        printf("  Removed %ld exterior input vertices.\n", 
               unuverts - delvertcount);
      }
    }
  }


  // Connect new hull tets.
  for (i = 0; i < hullarray->objects; i++) {
    parytet = (triface *) fastlookup(hullarray, i);
    hulltet = *parytet;
    for (j = 0; j < 3; j++) {
      esym(hulltet, neightet);
      if (neightet.tet[neightet.ver & 3] == NULL) {
        tspivot(hulltet, checksh);
        assert(checksh.sh != NULL);
        // Get the next subface in the same face ring of checksh. It must
        //   exist, otherwise, checksh is either a dangling subface (which
        //   should be removed already), or it is not a hull face.
        sfnext(checksh, casingout);
        assert(casingout.sh != NULL);
        // Go to the hull side.
        sesymself(casingout);
        stpivot(casingout, casface);
        assert(ishulltet(casface));
        esymself(casface);
        assert(casface.tet[casface.ver & 3] == NULL);
        // Bond the two hull tets together.
        bond(neightet, casface);
      }
      enextself(hulltet);
    }
  }

} else {  // '-c' option is set.


  long bak_subface_count = subfaces->items;
  long bak_segment_count = subsegs->items;

  // In this case, we regard every hull face/edge is a subface/segment.
  for (i = 0; i < tetarray->objects; i++) {
    parytet = (triface *) fastlookup(tetarray, i);
    // Only need the hull tet to find convex hull faces.
    if ((point) parytet->tet[7] == dummypoint) {
      hulltet.tet = parytet->tet;
      hulltet.ver = 3; // The hull face.
      tspivot(hulltet, checksh); // SELF_CHECK
      if (checksh.sh == NULL) {
        // Create a subface.
        makeshellface(subfaces, &checksh);
        pa = org(hulltet);
        pb = dest(hulltet);
        pc = apex(hulltet);
        setsorg(checksh, pa);
        setsdest(checksh, pb);
        setsapex(checksh, pc);
        // Create the point-to-subface map.
        sptr = sencode(checksh);
        setpoint2sh(pa, sptr);
        setpoint2sh(pb, sptr);
        setpoint2sh(pc, sptr);
      }
      // Insert this subface. 
      // Note: Even the subface is already exist, it may have been 
      //   disconnected from its adjacent tets.
      tsbond(hulltet, checksh);
      fsym(hulltet, neightet);
      assert(infected(neightet));
      sesymself(checksh);
      tsbond(neightet, checksh);
      sesymself(checksh);
      // Create three segments.
      for (j = 0; j < 3; j++) {
        tsspivot1(hulltet, checkseg);
        if (checkseg.sh == NULL) {
          // Create a segment.
          makeshellface(subsegs, &checkseg);
          pa = org(hulltet);
          pb = dest(hulltet);
          setshvertices(checkseg, pa, pb, NULL);
          // Insert the segment into the mesh.
          tetloop = hulltet;
          pc = apex(hulltet);
          checksh.sh = NULL;
          while (1) {
            tssbond1(tetloop, checkseg);
            tspivot(tetloop, checksh);
            if (checksh.sh != NULL) {
              ssbond1(checksh, checkseg);
              sbond1(checkseg, checksh);
            }
            fnextself(tetloop);
            if (apex(tetloop) == pc) break;
          }
          sstbond1(checkseg, tetloop);
        }
        enextself(hulltet);
      }
      // Save this hull tet in list.
      hullarray->newindex((void **) &parytet1);
      *parytet1 = hulltet;
    }
  } // i

  hullsize += hullarray->objects;

  if (subfacstack->objects > 0) {
    // Uninfect the collected exterior subfaces.
    for (i = 0; i < subfacstack->objects; i++) {
      parysh = (face *) fastlookup(subfacstack, i);
      suninfect(*parysh);
    }
  }

  if (b->regionattrib) {
    // Only the hull tets need to be uninfected.
    for (i = 0; i < hullarray->objects; i++) {
      parytet = (triface *) fastlookup(hullarray, i);
      uninfect(*parytet);
    }
  } else {
    // Uninfect all collected tets.
    for (i = 0; i < tetarray->objects; i++) {
      parytet = (triface *) fastlookup(tetarray, i);
      uninfect(*parytet);
    }
  }

  tetarray->restart();

  if (b->verbose) {
    printf("  Created %ld convex hull boundary faces.\n", 
           subfaces->items - bak_subface_count);
    printf("  Created %ld convex hull boundary edges.\n", 
           subsegs->items - bak_segment_count);
  }

} // if (b->convex)


  // Set region attributes (the -A option).
  if (b->regionattrib) {
    if (!b->quiet) {
      printf("Spreading region attributes.\n");
    }

    // If has user-defined region attributes.
    if (in->numberofregions > 0) {
      // Spread region attributes.
      for (i = 0; i < 5 * in->numberofregions; i += 5) {
        if (regiontets[i/5].tet != NULL) {
          attr = (int) in->regionlist[i + 3];
          volume = in->regionlist[i + 4];
          tetarray->restart(); // Re-use this array.
          infect(regiontets[i/5]);
          tetarray->newindex((void **) &parytet);
          *parytet = regiontets[i/5];
          // Collect and set attrs for all tets of this region.
          for (j = 0; j < tetarray->objects; j++) {
            parytet = (triface *) fastlookup(tetarray, j);
            tetloop = *parytet;
            setelemattribute(tetloop.tet, attrnum, attr);
            if (b->varvolume) { // If has -a option.
              setvolumebound(tetloop.tet, volume);
            }
            for (k = 0; k < 4; k++) {
              decode(tetloop.tet[k], neightet);
              // Is the adjacent already checked?
              if (!infected(neightet)) {
                // Is this side protected by a subface?
                tspivot(neightet, checksh);
                if (checksh.sh == NULL) {
                  infect(neightet);
                  tetarray->newindex((void **) &parytet);
                  *parytet = neightet;
                }
              }
            } // k
          } // j
          regioncount++;
        } // if (regiontets[i/5].tet != NULL)
      } // i
    }

    // Set attributes for all tetrahedra.
    attr = maxattr + 1;
    tetrahedrons->traversalinit();
    tetloop.tet = tetrahedrontraverse();
    while (tetloop.tet != (tetrahedron *) NULL) {
      if (!infected(tetloop)) {
        // An unmarked region.
        tetarray->restart(); // Re-use this array.
        infect(tetloop);
        tetarray->newindex((void **) &parytet);
        *parytet = tetloop;
        // Find and mark all tets.
        for (j = 0; j < tetarray->objects; j++) {
          parytet = (triface *) fastlookup(tetarray, j);
          tetloop = *parytet;
          setelemattribute(tetloop.tet, attrnum, attr);
          for (k = 0; k < 4; k++) {
            decode(tetloop.tet[k], neightet);
            // Is the adjacent tet already checked?
            if (!infected(neightet)) {
              // Is this side protected by a subface?
              tspivot(neightet, checksh);
              if (checksh.sh == NULL) {
                infect(neightet);
                tetarray->newindex((void **) &parytet);
                *parytet = neightet;
              }
            }
          } // k
        } // j
        attr++; // Increase the attribute.
        regioncount++;
      }
      tetloop.tet = tetrahedrontraverse();
    }
    // Until here, every tet has a region attribute.

    // Uninfect processed tets.
    tetrahedrons->traversalinit();
    tetloop.tet = tetrahedrontraverse();
    while (tetloop.tet != (tetrahedron *) NULL) {
      uninfect(tetloop);
      tetloop.tet = tetrahedrontraverse();
    }


    if (b->verbose) {
      assert(regioncount > 0);
      if (regioncount > 1) {
        printf("  Found %d subdomains.\n", regioncount);
      } else {
        printf("  Found 1 domain.\n");
      }
    }
  } // if (b->regionattrib)

  if (b->regionattrib && (in->numberofregions > 0)) { // If has -A option.
    delete [] regiontets;
  }
  delete tetarray;
  delete hullarray;

if (!b->convex) {

  // The mesh is non-convex now.
  nonconvex = 1;


  // Push all hull tets into 'flipstack'.
  tetrahedrons->traversalinit();
  tetloop.ver = 11; // The face opposite to dummypoint.
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    if ((point) tetloop.tet[7] == dummypoint) {
      flippush(flipstack, &tetloop);
    }
    tetloop.tet = alltetrahedrontraverse();
  }

  // Peel "slivers" off the hull.
  lawsonflip3d(NULL, 4, 1, 0, 0);

  if (b->verbose && (opt_sliver_peels > 0l)) {
    printf("  Peeled %ld hull slivers.\n", opt_sliver_peels);
  }

}

}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// reconstructmesh()    Reconstruct a tetrahedral mesh.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::reconstructmesh()
{
  tetrahedron *ver2tetarray;
  point *idx2verlist;
  triface tetloop, checktet, prevchktet;
  triface hulltet, face1, face2;
  tetrahedron tptr;
  face subloop, neighsh, nextsh;
  face segloop;
  shellface sptr;
  point p[4], q[3];
  REAL ori, attrib, volume;
  REAL angtol, ang;
  int eextras, marker = 0;
  int bondflag;
  int idx, i, j, k;

  if (!b->quiet) {
    printf("Reconstructing mesh ...\n");
  }
  // Default assume the mesh is non-convex.
  nonconvex = 1;

  // Create a map from indices to vertices.
  makeindex2pointmap(idx2verlist);

  // Allocate an array that maps each vertex to its adjacent tets.
  ver2tetarray = new tetrahedron[in->numberofpoints + 1];
  for (i = 0; i < in->numberofpoints; i++) {
    ver2tetarray[i] = NULL;
  }

  // Create the tetrahedra and connect those that share a common face.
  for (i = 0; i < in->numberoftetrahedra; i++) {
    // Get the four vertices.
    idx = i * in->numberofcorners;
    for (j = 0; j < 4; j++) {
      p[j] = idx2verlist[in->tetrahedronlist[idx++]];
      setpointtype(p[j], VOLVERTEX); // initial type.
    }
    // Check the orientation.
    ori = orient3d(p[0], p[1], p[2], p[3]);
    if (ori > 0.0) {
      // Swap the first two vertices.
      q[0] = p[0]; p[0] = p[1]; p[1] = q[0];
    } else if (ori == 0.0) {
      if (!b->quiet) {
        printf("Warning:  Tet #%d is degenerate.\n", i + in->firstnumber);
      }
    }
    // Create a new tetrahedron.
    maketetrahedron(&tetloop); // tetloop.ver = 11.
    setvertices(tetloop, p[0], p[1], p[2], p[3]);
    // Set element attributes if they exist.
    for (j = 0; j < in->numberoftetrahedronattributes; j++) {
      idx = i * in->numberoftetrahedronattributes;
      attrib = in->tetrahedronattributelist[idx + j];
      setelemattribute(tetloop.tet, j, attrib);
    }
    // If -a switch is used (with no number follows) Set a volume
    //   constraint if it exists.
    if (b->varvolume) {
      if (in->tetrahedronvolumelist != (REAL *) NULL) {
        volume = in->tetrahedronvolumelist[i];
      } else {
        volume = -1.0;
      }
      setvolumebound(tetloop.tet, volume);
    }
    // Try connecting this tet to others that share the common faces.
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      p[3] = oppo(tetloop);
      // Look for other tets having this vertex.
      idx = pointmark(p[3]);
      tptr = ver2tetarray[idx];
      // Link the current tet to the next one in the stack.
      tetloop.tet[8 + tetloop.ver] = tptr;
      // Push the current tet onto the stack.
      ver2tetarray[idx] = encode(tetloop);
      decode(tptr, checktet);
      if (checktet.tet != NULL) {
        p[0] =  org(tetloop); // a
        p[1] = dest(tetloop); // b
        p[2] = apex(tetloop); // c
        prevchktet = tetloop;
        do {
          assert(checktet.ver < 4); // SELF_CHECK
          q[0] =  org(checktet); // a'
          q[1] = dest(checktet); // b'
          q[2] = apex(checktet); // c'
          // Check the three faces at 'd' in 'checktet'.
          bondflag = 0;
          for (j = 0; j < 3; j++) {
            // Go to the face [b',a',d], or [c',b',d], or [a',c',d].
            esym(checktet, face2);
            if (face2.tet[face2.ver & 3] == NULL) {
              k = ((j + 1) % 3);
              if (q[k] == p[0]) {   // b', c', a' = a
                if (q[j] == p[1]) { // a', b', c' = b
                  // [#,#,d] is matched to [b,a,d].
                  esym(tetloop, face1);
                  bond(face1, face2);
                  bondflag++;
                }
              }
              if (q[k] == p[1]) {   // b',c',a' = b
                if (q[j] == p[2]) { // a',b',c' = c
                  // [#,#,d] is matched to [c,b,d].
                  enext(tetloop, face1);
                  esymself(face1);
                  bond(face1, face2);
                  bondflag++;
                }
              }
              if (q[k] == p[2]) {   // b',c',a' = c
                if (q[j] == p[0]) { // a',b',c' = a
                  // [#,#,d] is matched to [a,c,d].
                  eprev(tetloop, face1);
                  esymself(face1);
                  bond(face1, face2);
                  bondflag++;
                }
              }
            } else {
              bondflag++;
            }
            enextself(checktet);
          } // j
          // Go to the next tet in the link.
          tptr = checktet.tet[8 + checktet.ver];
          if (bondflag == 3) {
            // All three faces at d in 'checktet' have been connected.
            // It can be removed from the link.            
            prevchktet.tet[8 + prevchktet.ver] = tptr;
          } else {
            // Bakup the previous tet in the link.
            prevchktet = checktet;
          }
          decode(tptr, checktet);
        } while (checktet.tet != NULL);
      } // if (checktet.tet != NULL)
    } // for (tetloop.ver = 0; ...
  } // i

  // Remember a tet of the mesh.
  recenttet = tetloop;

  // Create hull tets, create the point-to-tet map, and clean up the
  //   temporary spaces used in each tet. 
  hullsize = tetrahedrons->items;

  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    tptr = encode(tetloop);
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      if (tetloop.tet[tetloop.ver] == NULL) {
        // Create a hull tet.
        maketetrahedron(&hulltet);
        p[0] =  org(tetloop);
        p[1] = dest(tetloop);
        p[2] = apex(tetloop);
        setvertices(hulltet, p[1], p[0], p[2], dummypoint);
        bond(tetloop, hulltet);
        // Try connecting this to others that share common hull edges.
        for (j = 0; j < 3; j++) {
          fsym(hulltet, face2);
          while (1) {
            if (face2.tet == NULL) break;
            esymself(face2);
            if (apex(face2) == dummypoint) break;
            fsymself(face2);
          }
          if (face2.tet != NULL) {
            // Found an adjacent hull tet.
            assert(face2.tet[face2.ver & 3] == NULL);
            esym(hulltet, face1);
            bond(face1, face2);
          }
          enextself(hulltet);
        }
        //hullsize++;
      }
      // Create the point-to-tet map.
      setpoint2tet((point) (tetloop.tet[4 + tetloop.ver]), tptr);
      // Clean the temporary used space.
      tetloop.tet[8 + tetloop.ver] = NULL;
    }
    tetloop.tet = tetrahedrontraverse();
  }

  hullsize = tetrahedrons->items - hullsize;

  // Subfaces will be inserted into the mesh. 
  if (in->trifacelist != NULL) {
    // A .face file is given. It may contain boundary faces. Insert them.
    for (i = 0; i < in->numberoftrifaces; i++) {
      // Is it a subface?
      if (in->trifacemarkerlist != NULL) {
        marker = in->trifacemarkerlist[i];
      } else {
        // Face markers are not available. Assume all of them are subfaces.
        marker = 1;
      }
      if (marker > 0) {
        idx = i * 3;
        for (j = 0; j < 3; j++) {
          p[j] = idx2verlist[in->trifacelist[idx++]];
        }
        // Search the subface.
        bondflag = 0;
        // Make sure all vertices are in the mesh. Avoid crash.
        for (j = 0; j < 3; j++) {
          decode(point2tet(p[j]), checktet);
          if (checktet.tet == NULL) break;
        }
        if ((j == 3) && getedge(p[0], p[1], &checktet)) {
          tetloop = checktet;
          q[2] = apex(checktet);
          while (1) {
            if (apex(tetloop) == p[2]) {
              // Found the face.
              // Check if there exist a subface already?
              tspivot(tetloop, neighsh); 
              if (neighsh.sh != NULL) {
                // Found a duplicated subface. 
                // This happens when the mesh was generated by other mesher.
                bondflag = 0;
              } else {
                bondflag = 1;
              }
              break;
            }
            fnextself(tetloop);
            if (apex(tetloop) == q[2]) break;
          }
        }
        if (bondflag) {
          // Create a new subface.
          makeshellface(subfaces, &subloop);
          setshvertices(subloop, p[0], p[1], p[2]);
          // Create the point-to-subface map.
          sptr = sencode(subloop);
          for (j = 0; j < 3; j++) {
            setpointtype(p[j], FACETVERTEX); // initial type.
            setpoint2sh(p[j], sptr);
          }
          if (in->trifacemarkerlist != NULL) {
            setshellmark(subloop, in->trifacemarkerlist[i]);
          }
          // Insert the subface into the mesh.
          tsbond(tetloop, subloop);
          fsymself(tetloop);
          sesymself(subloop);
          tsbond(tetloop, subloop);
        } else {
          if (!b->quiet) {
            if (neighsh.sh == NULL) {
              printf("Warning:  Subface #%d [%d,%d,%d] is missing.\n", 
                     i + in->firstnumber, pointmark(p[0]), pointmark(p[1]),
                     pointmark(p[2]));
            } else {
              printf("Warning: Ignore a dunplicated subface #%d [%d,%d,%d].\n", 
                     i + in->firstnumber, pointmark(p[0]), pointmark(p[1]),
                     pointmark(p[2]));
            }
          }
        } // if (bondflag)
      } // if (marker > 0)
    } // i
  } // if (in->trifacelist)

    // Indentify subfaces from the mesh.
    // Create subfaces for hull faces (if they're not subface yet) and
    //   interior faces which separate two different materials.
    eextras = in->numberoftetrahedronattributes;
    tetrahedrons->traversalinit();
    tetloop.tet = tetrahedrontraverse();
    while (tetloop.tet != (tetrahedron *) NULL) {
      for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
        tspivot(tetloop, neighsh);
        if (neighsh.sh == NULL) {
          bondflag = 0;
          fsym(tetloop, checktet);
          if (ishulltet(checktet)) {
            bondflag = 1;  // A hull face.
          } else {
            if (eextras > 0) {
              if (elemattribute(tetloop.tet, eextras - 1) !=
                  elemattribute(checktet.tet, eextras - 1)) {
                bondflag = 1; // An interior interface.
              }
            }
          }
          if (bondflag) {
            // Create a new subface.
            makeshellface(subfaces, &subloop);
            p[0] = org(tetloop);
            p[1] = dest(tetloop);
            p[2] = apex(tetloop);
            setshvertices(subloop, p[0], p[1], p[2]);
            // Create the point-to-subface map.
            sptr = sencode(subloop);
            for (j = 0; j < 3; j++) {
              setpointtype(p[j], FACETVERTEX); // initial type.
              setpoint2sh(p[j], sptr);
            }
            setshellmark(subloop, 0); // Default marker.
            // Insert the subface into the mesh.
            tsbond(tetloop, subloop);
            sesymself(subloop);
            tsbond(checktet, subloop);
          } // if (bondflag)
        } // if (neighsh.sh == NULL)
      }
      tetloop.tet = tetrahedrontraverse();
    }

  // Connect subfaces together. 
  subfaces->traversalinit();
  subloop.shver = 0;
  subloop.sh = shellfacetraverse(subfaces);
  while (subloop.sh != (shellface *) NULL) {
    for (i = 0; i < 3; i++) {
      spivot(subloop, neighsh);
      if (neighsh.sh == NULL) {
        // Form a subface ring by linking all subfaces at this edge.
        // Traversing all faces of the tets at this edge.
        stpivot(subloop, tetloop);
        q[2] = apex(tetloop);
        neighsh = subloop;
        while (1) {
          fnextself(tetloop);
          tspivot(tetloop, nextsh);
          if (nextsh.sh != NULL) {
            // Link neighsh <= nextsh.
            sbond1(neighsh, nextsh);
            neighsh = nextsh;
          }
          if (apex(tetloop) == q[2]) {
            assert(nextsh.sh == subloop.sh); // It's a ring.
            break;
          }
        } // while (1)
      } // if (neighsh.sh == NULL)
      senextself(subloop);
    }
    subloop.sh = shellfacetraverse(subfaces);
  }

  //if (b->verbose) {
  //  printf("  Created %ld subfaces.\n", subfaces->items);
  //}

  // Segments will be introudced. 
  if (in->edgelist != NULL) {
    // A .edge file is given. It may contain boundary edges. Insert them.
    for (i = 0; i < in->numberofedges; i++) {
      // Is it a segment?
      if (in->edgemarkerlist != NULL) {
        marker = in->edgemarkerlist[i];
      } else {
        // Edge markers are not available. Assume all of them are segments.
        marker = 1;
      }
      if (marker != 0) { 
        // Insert a segment.
        idx = i * 2;
        for (j = 0; j < 2; j++) {
          p[j] = idx2verlist[in->edgelist[idx++]];
        }
        // Make sure all vertices are in the mesh. Avoid crash.
        for (j = 0; j < 2; j++) {
          decode(point2tet(p[j]), checktet);
          if (checktet.tet == NULL) break;
        }
        // Search the segment.
        if ((j == 2) && getedge(p[0], p[1], &checktet)) {
          // Create a new subface.
          makeshellface(subsegs, &segloop);
          setshvertices(segloop, p[0], p[1], NULL);
          // Create the point-to-segment map.
          sptr = sencode(segloop);
          for (j = 0; j < 2; j++) {
            setpointtype(p[j], RIDGEVERTEX); // initial type.
            setpoint2sh(p[j], sptr);
          }
          if (in->edgemarkerlist != NULL) {
            setshellmark(segloop, marker);
          }
          // Insert the segment into the mesh.
          tetloop = checktet;
          q[2] = apex(checktet);
          subloop.sh = NULL;
          while (1) {
            tssbond1(tetloop, segloop);
            tspivot(tetloop, subloop);
            if (subloop.sh != NULL) {
              ssbond1(subloop, segloop);
              sbond1(segloop, subloop);
            }
            fnextself(tetloop);
            if (apex(tetloop) == q[2]) break;
          } // while (1)
          // Remember an adjacent tet for this segment.
          sstbond1(segloop, tetloop);
        } else {
          if (!b->quiet) {
            printf("Warning:  Segment #%d [%d,%d] is missing.\n", 
                   i + in->firstnumber, pointmark(p[0]), pointmark(p[1]));
          }
        }
      } // if (marker != 0)
    } // i
  } // if (in->edgelist)

    // Identify segments from the mesh. 
    // Create segments for non-manifold edges (which are shared by more 
    //   than two subfaces), and for non-coplanar edges, i.e., two subfaces
    //   form an dihedral angle > 'b->facet_ang_tol' (degree).
    angtol = b->facet_ang_tol / 180.0 * PI;
    subfaces->traversalinit();
    subloop.shver = 0;
    subloop.sh = shellfacetraverse(subfaces);
    while (subloop.sh != (shellface *) NULL) {
      for (i = 0; i < 3; i++) {
        sspivot(subloop, segloop);
        if (segloop.sh == NULL) {
          // Check if this edge is a segment.
          bondflag = 0;
          // Counter the number of subfaces at this edge.
          idx = 0;
          nextsh = subloop;
          while (1) {
            idx++;
            spivotself(nextsh);
            if (nextsh.sh == subloop.sh) break;
          }
          if (idx != 2) {
            // It's a non-manifold edge. Insert a segment.
            p[0] = sorg(subloop);
            p[1] = sdest(subloop);
            bondflag = 1;
          } else {
            // Check the dihedral angle formed by the two subfaces.
            spivot(subloop, neighsh);
            p[0] = sorg(subloop);
            p[1] = sdest(subloop);
            p[2] = sapex(subloop);
            p[3] = sapex(neighsh);
            ang = facedihedral(p[0], p[1], p[2], p[3]);
            if (ang > PI) ang = 2 * PI - ang;
            if (ang < angtol) {
              bondflag = 1;
            }
          }
          if (bondflag) {
            // Create a new subface.
            makeshellface(subsegs, &segloop);
            setshvertices(segloop, p[0], p[1], NULL);
            // Create the point-to-segment map.
            sptr = sencode(segloop);
            for (j = 0; j < 2; j++) {
              setpointtype(p[j], RIDGEVERTEX); // initial type.
              setpoint2sh(p[j], sptr);
            }
            setshellmark(segloop, marker);
            // Insert the subface into the mesh.
            stpivot(subloop, tetloop);
            q[2] = apex(tetloop);
            while (1) {
              tssbond1(tetloop, segloop);
              tspivot(tetloop, neighsh);
              if (neighsh.sh != NULL) {
                ssbond1(neighsh, segloop);
              }
              fnextself(tetloop);
              if (apex(tetloop) == q[2]) break;
            } // while (1)
            // Remember an adjacent tet for this segment.
            sstbond1(segloop, tetloop);
            sbond1(segloop, subloop);
          } // if (bondflag)
        } // if (neighsh.sh == NULL)
        senextself(subloop);
      }
      subloop.sh = shellfacetraverse(subfaces);
    }

  // Remember the number of input segments.
  insegments = subsegs->items;

  //if (b->verbose) {
  //  printf("  Created %ld segments.\n", subsegs->items);
  //}

  // Set global flags.
  checksubsegflag = 1;
  checksubfaceflag = 1;
  //nonconvex = 1; 

  delete [] idx2verlist;
  delete [] ver2tetarray;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// scoutpoint()    Search a point in mesh.                                   //
//                                                                           //
// This function searches the point in a mesh whose domain may be not convex.//
// In case of a convex domain, the locate() function is sufficient.          //
//                                                                           //
// If 'randflag' is used, randomly select a start searching tet.  Otherwise, //
// start searching directly from 'searchtet'.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::scoutpoint(point searchpt, triface *searchtet, int randflag)
{
  point pa, pb, pc, pd;
  enum locateresult loc = OUTSIDE;
  REAL vol, ori1, ori2, ori3, ori4;
  int iter;

  if (searchtet->tet == NULL) {
    *searchtet = recenttet;
  }

  iter = 0;
  while (1) {
    // Randonmly select a good starting tet.
    if (randflag) {
      randomsample(searchpt, searchtet);
    }
    loc = locate(searchpt, searchtet, 0);
    if (loc == OUTSIDE) {
      // Not found. This happens when the mesh is not convex.
      if (!randflag) break;
      iter++;
      if (iter > 3) {
        searchtet->tet = NULL;
        break;
      }
    } else {
      // Found the point.
      break;
    }
  } // while (1)

  if (loc != OUTSIDE) {
    // Round the result of location.
    pa = org(*searchtet);
    pb = dest(*searchtet);
    pc = apex(*searchtet);
    pd = oppo(*searchtet);
    vol = orient3d(pa, pb, pc, pd);
    ori1 = orient3d(pa, pb, pc, searchpt);
    ori2 = orient3d(pb, pa, pd, searchpt);
    ori3 = orient3d(pc, pb, pd, searchpt);
    ori4 = orient3d(pa, pc, pd, searchpt);
    if (fabs(ori1 / vol) < b->epsilon) ori1 = 0;
    if (fabs(ori2 / vol) < b->epsilon) ori2 = 0;
    if (fabs(ori3 / vol) < b->epsilon) ori3 = 0;
    if (fabs(ori4 / vol) < b->epsilon) ori4 = 0;
  } else { // if (loc == OUTSIDE) {
    // Do a brute force search for the point (with rounding).
    tetrahedrons->traversalinit();
    searchtet->tet = tetrahedrontraverse();
    while (searchtet->tet != NULL) {
      pa = org(*searchtet);
      pb = dest(*searchtet);
      pc = apex(*searchtet);
      pd = oppo(*searchtet);

      vol = orient3d(pa, pb, pc, pd); 
      assert(vol < 0); // vol != 0

      ori1 = orient3d(pa, pb, pc, searchpt);
      if (fabs(ori1 / vol) < b->epsilon) ori1 = 0; // Rounding.
      if (ori1 <= 0) {
        ori2 = orient3d(pb, pa, pd, searchpt);
        if (fabs(ori2 / vol) < b->epsilon) ori2 = 0;
        if (ori2 <= 0) {
          ori3 = orient3d(pc, pb, pd, searchpt);
          if (fabs(ori3 / vol) < b->epsilon) ori3 = 0;
          if (ori3 <= 0) {
            ori4 = orient3d(pa, pc, pd, searchpt);
            if (fabs(ori4 / vol) < b->epsilon) ori4 = 0;
            if (ori4 <= 0) {
              // Found the tet. Return its location. 
              break;
            } // ori4
          } // ori3
        } // ori2
      } // ori1

      searchtet->tet = bgm->tetrahedrontraverse();
    } // while (searchtet->tet != NULL)
  }

  if (searchtet->tet != NULL) {
    // Return the point location.
    if (ori1 == 0) { // on face [a,b,c]
      if (ori2 == 0) { // on edge [a,b].
        if (ori3 == 0) { // on vertex [b].
          assert(ori4 != 0);
          enextself(*searchtet); // [b,c,a,d]
          loc = ONVERTEX;
        } else {
          if (ori4 == 0) { // on vertex [a]
            loc =  ONVERTEX; // [a,b,c,d]
          } else {    
            loc =  ONEDGE; // [a,b,c,d]
          }
        }
      } else { // ori2 != 0
        if (ori3 == 0) { // on edge [b,c]
          if (ori4 == 0) { // on vertex [c]
            eprevself(*searchtet); // [c,a,b,d]
            loc =  ONVERTEX;
          } else {
            enextself(*searchtet); // [b,c,a,d]
            loc =  ONEDGE;
          }
        } else { // ori3 != 0
          if (ori4 == 0) { // on edge [c,a]
            eprevself(*searchtet); // [c,a,b,d]
            loc =  ONEDGE;
          } else {
            loc =  ONFACE;
          }
        }
      }
    } else { // ori1 != 0
      if (ori2 == 0) { // on face [b,a,d]
        esymself(*searchtet); // [b,a,d,c]
        if (ori3 == 0) { // on edge [b,d]
          eprevself(*searchtet); // [d,b,a,c]
          if (ori4 == 0) { // on vertex [d]                      
            loc =  ONVERTEX;
          } else {
            loc =  ONEDGE;
          }
        } else { // ori3 != 0
          if (ori4 == 0) { // on edge [a,d]
            enextself(*searchtet); // [a,d,b,c]
            loc =  ONEDGE;
          } else {
            loc =  ONFACE;
          }
        }
      } else { // ori2 != 0
        if (ori3 == 0) { // on face [c,b,d]
          enextself(*searchtet);
          esymself(*searchtet);
          if (ori4 == 0) { // on edge [c,d]
            eprevself(*searchtet);
            loc =  ONEDGE;
          } else {
            loc =  ONFACE;
          }
        } else {
          if (ori4 == 0) { // on face [a,c,d]
            eprevself(*searchtet);
            esymself(*searchtet);
            loc =  ONFACE;
          } else { // inside tet [a,b,c,d]
            loc =  INTETRAHEDRON;
          } // ori4
        } // ori3
      } // ori2
    } // ori1
  } else {
    loc = OUTSIDE;
  }

  return (int) loc;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// getpointmeshsize()    Interpolate the mesh size at given point.           //
//                                                                           //
// 'iloc' indicates the location of the point w.r.t. 'searchtet'.  The size  //
// is obtained by linear interpolation on the vertices of the tet.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL tetgenmesh::getpointmeshsize(point searchpt, triface *searchtet, int iloc)
{
  point *pts, pa, pb, pc;
  REAL volume, vol[4], wei[4];
  REAL size;
  int i;

  size = 0;

  if (iloc == (int) INTETRAHEDRON) {
    pts = (point *) &(searchtet->tet[4]);
    assert(pts[3] != dummypoint);
    // Only do interpolation if all vertices have non-zero sizes.
    if ((pts[0][pointmtrindex] > 0) && (pts[1][pointmtrindex] > 0) &&
        (pts[2][pointmtrindex] > 0) && (pts[3][pointmtrindex] > 0)) {
      // P1 interpolation.
      volume = orient3d(pts[0], pts[1], pts[2], pts[3]);
      vol[0] = orient3d(searchpt, pts[1], pts[2], pts[3]);
      vol[1] = orient3d(pts[0], searchpt, pts[2], pts[3]);
      vol[2] = orient3d(pts[0], pts[1], searchpt, pts[3]);
      vol[3] = orient3d(pts[0], pts[1], pts[2], searchpt);
      for (i = 0; i < 4; i++) {
        wei[i] = fabs(vol[i] / volume);
        size += (wei[i] * pts[i][pointmtrindex]);
      }
    }
  } else if (iloc == (int) ONFACE) {
    pa = org(*searchtet);
    pb = dest(*searchtet);
    pc = apex(*searchtet);
    if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0) &&
        (pc[pointmtrindex] > 0)) {
      volume = triarea(pa, pb, pc);
      vol[0] = triarea(searchpt, pb, pc);
      vol[1] = triarea(pa, searchpt, pc);
      vol[2] = triarea(pa, pb, searchpt);
      size = (vol[0] / volume) * pa[pointmtrindex]
           + (vol[1] / volume) * pb[pointmtrindex]
           + (vol[2] / volume) * pc[pointmtrindex];
    }
  } else if (iloc == (int) ONEDGE) {
    pa = org(*searchtet);
    pb = dest(*searchtet);
    if ((pa[pointmtrindex] > 0) && (pb[pointmtrindex] > 0)) {
      volume = distance(pa, pb);
      vol[0] = distance(searchpt, pb);
      vol[1] = distance(pa, searchpt);
      size = (vol[0] / volume) * pa[pointmtrindex]
           + (vol[1] / volume) * pb[pointmtrindex];
    }
  } else if (iloc == (int) ONVERTEX) {
    pa = org(*searchtet);
    if (pa[pointmtrindex] > 0) {
      size = pa[pointmtrindex];
    }
  }

  return size;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// interpolatemeshsize()    Interpolate the mesh size from a background mesh //
//                          (source) to the current mesh (destination).      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::interpolatemeshsize()
{
  triface searchtet;
  point ploop;
  REAL minval = 0.0, maxval = 0.0;
  int iloc;
  int count;

  if (!b->quiet) {
    printf("Interpolating mesh size ...\n");
  }
  count = 0; // Count the number of interpolated points.

  // Interpolate sizes for all points in the current mesh.
  points->traversalinit();
  ploop = pointtraverse();
  while (ploop != NULL) {
    // Search a tet in bgm which containing this point.
    searchtet.tet = NULL;
    iloc = bgm->scoutpoint(ploop, &searchtet, 1); // randflag = 1
    if (iloc != (int) OUTSIDE) {
      // Interpolate the mesh size.
      ploop[pointmtrindex] = bgm->getpointmeshsize(ploop, &searchtet, iloc);
      setpoint2bgmtet(ploop, bgm->encode(searchtet));
      if (count == 0) {
        // This is the first interpolated point.
        minval = maxval = ploop[pointmtrindex];
      } else {
        if (ploop[pointmtrindex] < minval) {
          minval = ploop[pointmtrindex];
        }
        if (ploop[pointmtrindex] > maxval) {
          maxval = ploop[pointmtrindex];
        }
      }
      count++;
    } else {
      if (!b->quiet) {
        printf("Warnning:  Failed to locate point %d in source mesh.\n",
               pointmark(ploop));
      }
    }
    ploop = pointtraverse();
  }

  if (b->verbose) {
    printf("  Interoplated %d points.\n", count);
    printf("  Size rangle [%.17g, %.17g].\n", minval, maxval);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertconstrainedpoints()    Insert a list of points into the mesh.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::insertconstrainedpoints(tetgenio *addio)
{
  triface searchtet, spintet;
  face checksh, *splitsh;
  face checkseg, *splitseg;
  point newpt;
  insertvertexflags ivf;
  REAL x, y, z, w;
  int attribindex, mtrindex;
  int randflag;
  int count, index;
  int loc;
  int i, j;

  if (!b->quiet) {
    printf("Inserting constrained points ...\n");
  }

  randflag = 1; // Randomly select start tet for point location. 
  count = 0;
  index = 0;
  attribindex = 0;
  mtrindex = 0;

  for (i = 0; i < addio->numberofpoints; i++) {
    makepoint(&newpt, UNUSEDVERTEX);
    x = newpt[0] = addio->pointlist[index++];
    y = newpt[1] = addio->pointlist[index++];
    z = newpt[2] = addio->pointlist[index++];
    // Read the point attributes. (Including point weights.)
    for (j = 0; j < addio->numberofpointattributes; j++) {
      newpt[3 + j] = addio->pointattributelist[attribindex++];
    }
    // Read the point metric tensor.
    for (j = 0; j < addio->numberofpointmtrs; j++) {
      newpt[pointmtrindex + j] = addio->pointmtrlist[mtrindex++];
    }
    if (b->weighted) { // -w option
      if (addio->numberofpointattributes > 0) {
        // The first point attribute is its weight.
        w = newpt[3];
      } else {
        // No given weight available. Default choose the maximum
        //   absolute value among its coordinates.        
        w = fabs(x);
        if (w < fabs(y)) w = fabs(y);
        if (w < fabs(z)) w = fabs(z);
      }
      if (b->weighted_param == 0) {
        newpt[3] = x * x + y * y + z * z - w; // Weighted DT.
      } else { // -w1 option
        newpt[3] = w;  // Regular tetrahedralization.
      }
    }

    // Find the location of the inserted point.
    searchtet.tet = NULL;
    ivf.iloc = scoutpoint(newpt, &searchtet, randflag);
    if (ivf.iloc != (int) OUTSIDE) {
      // Found the point. 
      // Initialize the insertion parameters. 
      if (b->psc) {
        ivf.bowywat = 0;   // Do not enlarge the initial cavity.
        ivf.validflag = 0; // Do not validate the initial cavity.
      } else {
        ivf.bowywat = 3;   // Use the "Bowyer-Watson" algorithm to form cavity.
        ivf.validflag = 1; // Validate the B-W cavity.
      }
      ivf.lawson = 3;    // ???
      ivf.rejflag = 0;   // ???
      ivf.chkencflag = 0;
      ivf.sloc = ivf.iloc;
      ivf.sbowywat = ivf.bowywat;  // Surface mesh options.
      ivf.splitbdflag = 1;
      ivf.respectbdflag = 1;
      ivf.assignmeshsize = b->metric;

      splitsh = NULL;
      splitseg = NULL;

      // Set the right point type.
      if (ivf.iloc == (int) ONEDGE) {
        tsspivot1(searchtet, checkseg);
        if (checkseg.sh != NULL) {
          setpointtype(newpt, RIDGEVERTEX);
          spivot(checkseg, checksh);
          splitsh = &checksh;
          splitseg = &checkseg;        
        } else {
          // Check if it is a subface edge.
          spintet = searchtet;
          while (1) {
            tspivot(spintet, checksh);
            if (checksh.sh != NULL) {
              setpointtype(newpt, FACETVERTEX);
              splitsh = &checksh;
              break;
            }
            fnextself(spintet);
            if (spintet.tet == searchtet.tet) break;
          }
        }
      } else if (ivf.iloc == (int) ONFACE) {
        tspivot(searchtet, checksh);
        if (checksh.sh != NULL) {
          setpointtype(newpt, FACETVERTEX);
          splitsh = &checksh;
        }
      } else {
        setpointtype(newpt, VOLVERTEX);
      }

      // Insert the vertex.
      loc = insertvertex(newpt, &searchtet, splitsh, splitseg, &ivf);

      if (loc == ivf.iloc) {
        // The point has been inserted.
        lawsonflip3d(newpt, 4, 0, ivf.chkencflag, 0);
        count++;
      } else {
        if (!b->quiet) {
          printf("Warning:  Failed to insert point #%d. Ignored.\n", i);
        }
        pointdealloc(newpt);
      }
    } else {
      if (!b->quiet) {
        printf("Warning:  Can't locate add point #%d. Ignored.\n", i);
      }
      pointdealloc(newpt);
    }
  } // i

  if (b->verbose) {
    printf("  Inserted %d of %d vertices.\n", count, addio->numberofpoints);
  }
}

////                                                                       ////
////                                                                       ////
//// reconstruct_cxx //////////////////////////////////////////////////////////

