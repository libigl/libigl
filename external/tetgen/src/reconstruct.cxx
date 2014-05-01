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
  triface tetloop, neightet, *parytet, *parytet1;
  triface *regiontets = NULL;
  face checksh, *parysh;
  face checkseg;
  point ptloop, *parypt;
  int t1ver;
  int i, j, k;

  if (!b->quiet) {
    if (b->convex) {
      printf("Marking exterior tetrahedra ...\n");
    } else {
      printf("Removing exterior tetrahedra ...\n");
    }
  }

  // Initialize the pool of exterior tets.
  tetarray = new arraypool(sizeof(triface), 10);
  hullarray = new arraypool(sizeof(triface), 10);

  // Collect unprotected tets and hull tets.
  tetrahedrons->traversalinit();
  tetloop.ver = 11; // The face opposite to dummypoint.
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    if (ishulltet(tetloop)) {
      // Is this side protected by a subface?
      if (!issubface(tetloop)) {
        // Collect an unprotected hull tet and tet.
        infect(tetloop);
        hullarray->newindex((void **) &parytet);
        *parytet = tetloop;
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
      if (locate(&(in->holelist[i]), &neightet) != OUTSIDE) {
        // The tet 'neightet' contain this point.
        if (!infected(neightet)) {
          infect(neightet);
          tetarray->newindex((void **) &parytet);
          *parytet = neightet;
          // Add its adjacent tet if it is not protected.
          if (!issubface(neightet)) {
            decode(neightet.tet[neightet.ver & 3], tetloop);
            if (!infected(tetloop)) {
              infect(tetloop);
              if (ishulltet(tetloop)) {
                hullarray->newindex((void **) &parytet);
              } else {
                tetarray->newindex((void **) &parytet);
              }
              *parytet = tetloop;
            }
          }
          else {
            // It is protected. Check if its adjacent tet is a hull tet.
            decode(neightet.tet[neightet.ver & 3], tetloop);
            if (ishulltet(tetloop)) {
              // It is hull tet, add it into the list. Moreover, the subface
              //   is dead, i.e., both sides are in exterior.
              if (!infected(tetloop)) {
                infect(tetloop);
                hullarray->newindex((void **) &parytet);
                *parytet = tetloop;
              }
            }
            if (infected(tetloop)) {
              // Both sides of this subface are in exterior.
              tspivot(neightet, checksh);
              sinfect(checksh); // Only queue it once.
              subfacstack->newindex((void **) &parysh);
              *parysh = checksh;
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
  } // if (in->numberofholes > 0)

  if (b->regionattrib && (in->numberofregions > 0)) { // -A option.
    // Record the tetrahedra that contains the region points for assigning
    //   region attributes after the holes have been carved.
    regiontets = new triface[in->numberofregions];
    // Mark as marktested any tetrahedra inside volume regions.
    for (i = 0; i < 5 * in->numberofregions; i += 5) {
      // Search a tet containing the i-th region point.
      neightet.tet = NULL;
      randomsample(&(in->regionlist[i]), &neightet);
      if (locate(&(in->regionlist[i]), &neightet) != OUTSIDE) {
        regiontets[i/5] = neightet;
      } else {
        if (!b->quiet) {
          printf("Warning:  The %d-th region point ", i/5+1);
          printf("lies outside the convex hull.\n");
        }
        regiontets[i/5].tet = NULL;
      }
    }
  }

  // Collect all exterior tets (in concave place and in holes).
  for (i = 0; i < tetarray->objects; i++) {
    parytet = (triface *) fastlookup(tetarray, i);
    j = (parytet->ver & 3); // j is the current face number.
    // Check the other three adjacent tets.
    for (k = 1; k < 4; k++) {
      decode(parytet->tet[(j + k) % 4], neightet); 
      // neightet may be a hull tet.
      if (!infected(neightet)) {
        // Is neightet protected by a subface.
        if (!issubface(neightet)) {
          // Not proected. Collect it. (It must not be a hull tet).
          infect(neightet);
          tetarray->newindex((void **) &parytet1);
          *parytet1 = neightet;
        } else {
          // Protected. Check if it is a hull tet.
          if (ishulltet(neightet)) {
            // A hull tet. Collect it.
            infect(neightet);
            hullarray->newindex((void **) &parytet1);
            *parytet1 = neightet;
            // Both sides of this subface are exterior.
            tspivot(neightet, checksh);
            // Queue this subface (to be deleted later).
            assert(!sinfected(checksh));
            sinfect(checksh); // Only queue it once.
            subfacstack->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
      } else {
        // Both sides of this face are in exterior.
        // If there is a subface. It should be collected.
        if (issubface(neightet)) {
          tspivot(neightet, checksh);
          if (!sinfected(checksh)) {
            sinfect(checksh);
            subfacstack->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
      }
    } // j, k
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

  // Collect vertices which point to infected tets. These vertices
  //   may get deleted after the removal of exterior tets.
  //   If -Y1 option is used, collect all Steiner points for removal.
  //   The lists 'cavetetvertlist' and 'subvertstack' are re-used.
  points->traversalinit();
  ptloop = pointtraverse();
  while (ptloop != NULL) {
    if ((pointtype(ptloop) != UNUSEDVERTEX) &&
        (pointtype(ptloop) != DUPLICATEDVERTEX)) {
      decode(point2tet(ptloop), neightet);
      if (infected(neightet)) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = ptloop;
      }
      if (b->nobisect && (b->nobisect_param > 0)) { // -Y1
        // Queue it if it is a Steiner point.
        if (pointmark(ptloop) > 
              (in->numberofpoints - (in->firstnumber ? 0 : 1))) {
          subvertstack->newindex((void **) &parypt);
          *parypt = ptloop;
        }
      }
    }
    ptloop = pointtraverse();
  }

  if (!b->convex && (tetarray->objects > 0l)) { // No -c option.
    // Remove exterior tets. Hull tets are updated.
    arraypool *newhullfacearray;
    triface hulltet, casface;
    point pa, pb, pc;

    newhullfacearray = new arraypool(sizeof(triface), 10);

    // Create and save new hull tets.
    for (i = 0; i < tetarray->objects; i++) {
      parytet = (triface *) fastlookup(tetarray, i);
      for (j = 0; j < 4; j++) {
        decode(parytet->tet[j], tetloop);
        if (!infected(tetloop)) {
          // Found a new hull face (must be a subface).
          tspivot(tetloop, checksh);
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
            if (issubseg(tetloop)) {
              tsspivot1(tetloop, checkseg);
              tssbond1(hulltet, checkseg);
              sstbond1(checkseg, hulltet);
            }
            enextself(tetloop);
            eprevself(hulltet);
          }
          // Update the point-to-tet map.
          setpoint2tet(pa, (tetrahedron) tetloop.tet);
          setpoint2tet(pb, (tetrahedron) tetloop.tet);
          setpoint2tet(pc, (tetrahedron) tetloop.tet);
          // Save the exterior tet at this hull face. It still holds pointer
          //   to the adjacent interior tet. Use it to connect new hull tets. 
          newhullfacearray->newindex((void **) &parytet1);
          parytet1->tet = parytet->tet;
          parytet1->ver = j;
        } // if (!infected(tetloop))
      } // j
    } // i

    // Connect new hull tets.
    for (i = 0; i < newhullfacearray->objects; i++) {
      parytet = (triface *) fastlookup(newhullfacearray, i);
      fsym(*parytet, neightet);
      // Get the new hull tet.
      fsym(neightet, hulltet);
      for (j = 0; j < 3; j++) {
        esym(hulltet, casface);
        if (casface.tet[casface.ver & 3] == NULL) {
          // Since the boundary of the domain may not be a manifold, we
          //   find the adjacent hull face by traversing the tets in the
          //   exterior (which are all infected tets).
          neightet = *parytet;
          while (1) {
            fnextself(neightet);
            if (!infected(neightet)) break;
          }
          if (!ishulltet(neightet)) {
            // An interior tet. Get the new hull tet.
            fsymself(neightet);
            esymself(neightet);
          } 
          // Bond them together.
          bond(casface, neightet);
        }
        enextself(hulltet);
        enextself(*parytet);
      } // j
    } // i

    if (subfacstack->objects > 0l) {
      // Remove all subfaces which do not attach to any tetrahedron.
      //   Segments which are not attached to any subfaces and tets
      //   are deleted too.
      face casingout, casingin;
      long delsegcount = 0l;

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
    } // if (subfacstack->objects > 0l)

    if (cavetetvertlist->objects > 0l) {
      // Some vertices may lie in exterior. Marke them as UNUSEDVERTEX.
      long delvertcount = unuverts;
      long delsteinercount = 0l;

      for (i = 0; i < cavetetvertlist->objects; i++) {
        parypt = (point *) fastlookup(cavetetvertlist, i);
        decode(point2tet(*parypt), neightet);
        if (infected(neightet)) {
          // Found an exterior vertex.
          if (pointmark(*parypt) > 
                (in->numberofpoints - (in->firstnumber ? 0 : 1))) {
            // A Steiner point.
            if (pointtype(*parypt) == FREESEGVERTEX) {
              st_segref_count--;
            } else if (pointtype(*parypt) == FREEFACETVERTEX) {
              st_facref_count--;
            } else {
              assert(pointtype(*parypt) == FREEVOLVERTEX);
              st_volref_count--;
            }
            delsteinercount++;
            if (steinerleft > 0) steinerleft++;
          }
          setpointtype(*parypt, UNUSEDVERTEX);
          unuverts++;
        }
      }

      if (b->verbose) {
        if (unuverts > delvertcount) {
          if (delsteinercount > 0l) {
            if (unuverts > (delvertcount + delsteinercount)) {
              printf("  Removed %ld exterior input vertices.\n", 
                     unuverts - delvertcount - delsteinercount);
            }
            printf("  Removed %ld exterior Steiner vertices.\n", 
                   delsteinercount);
          } else {
            printf("  Removed %ld exterior input vertices.\n", 
                   unuverts - delvertcount);
          }
        }
      }
      cavetetvertlist->restart();
      // Comment: 'subvertstack' will be cleaned in routine
      //   suppresssteinerpoints().
    } // if (cavetetvertlist->objects > 0l)

    // Update the hull size.
    hullsize += (newhullfacearray->objects - hullarray->objects);

    // Delete all exterior tets and old hull tets.
    for (i = 0; i < tetarray->objects; i++) {
      parytet = (triface *) fastlookup(tetarray, i);
      tetrahedrondealloc(parytet->tet);
    }
    tetarray->restart();

    for (i = 0; i < hullarray->objects; i++) {
      parytet = (triface *) fastlookup(hullarray, i);
      tetrahedrondealloc(parytet->tet);
    }
    hullarray->restart();

    delete newhullfacearray;
  } // if (!b->convex && (tetarray->objects > 0l))

  if (b->convex && (tetarray->objects > 0l)) { // With -c option
    // In this case, all exterior tets get a region marker '-1'.
    assert(b->regionattrib > 0); // -A option must be enabled.
    int attrnum = numelemattrib - 1;

    for (i = 0; i < tetarray->objects; i++) {
      parytet = (triface *) fastlookup(tetarray, i);
      setelemattribute(parytet->tet, attrnum, -1);
    }
    tetarray->restart();

    for (i = 0; i < hullarray->objects; i++) {
      parytet = (triface *) fastlookup(hullarray, i);
      uninfect(*parytet);
    }
    hullarray->restart();

    if (subfacstack->objects > 0l) {
      for (i = 0; i < subfacstack->objects; i++) {
        parysh = (face *) fastlookup(subfacstack, i);
        suninfect(*parysh);
      }
      subfacstack->restart();
    }

    if (cavetetvertlist->objects > 0l) {
      cavetetvertlist->restart();
    }
  } // if (b->convex && (tetarray->objects > 0l))

  if (b->regionattrib) { // With -A option.
    if (!b->quiet) {
      printf("Spreading region attributes.\n");
    }
    REAL volume;
    int attr, maxattr = 0; // Choose a small number here.
    int attrnum = numelemattrib - 1; 
    // Comment: The element region marker is at the end of the list of
    //   the element attributes.
    int regioncount = 0;

    // If has user-defined region attributes.
    if (in->numberofregions > 0) {
      // Spread region attributes.
      for (i = 0; i < 5 * in->numberofregions; i += 5) {
        if (regiontets[i/5].tet != NULL) {
          attr = (int) in->regionlist[i + 3];
          if (attr > maxattr) {
            maxattr = attr;
          }
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
                if (!issubface(neightet)) {
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
              if (!issubface(neightet)) {
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
      //assert(regioncount > 0);
      if (regioncount > 1) {
        printf("  Found %d subdomains.\n", regioncount);
      } else {
        printf("  Found %d domain.\n", regioncount);
      }
    }
  } // if (b->regionattrib)

  if (regiontets != NULL) {
    delete [] regiontets;
  }
  delete tetarray;
  delete hullarray;

  if (!b->convex) { // No -c option
    // The mesh is non-convex now.
    nonconvex = 1;

    // Push all hull tets into 'flipstack'.
    tetrahedrons->traversalinit();
    tetloop.ver = 11; // The face opposite to dummypoint.
    tetloop.tet = alltetrahedrontraverse();
    while (tetloop.tet != (tetrahedron *) NULL) {
      if ((point) tetloop.tet[7] == dummypoint) {
        fsym(tetloop, neightet);
        flippush(flipstack, &neightet);
      }
      tetloop.tet = alltetrahedrontraverse();
    }

    flipconstraints fc;
    fc.enqflag = 2;
    long sliver_peel_count = lawsonflip3d(&fc);

    if (sliver_peel_count > 0l) {
      if (b->verbose) {
        printf("  Removed %ld hull slivers.\n", sliver_peel_count);
      }
    }
    unflipqueue->restart();
  } // if (!b->convex)
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
  int t1ver;
  int idx, i, j, k;

  if (!b->quiet) {
    printf("Reconstructing mesh ...\n");
  }

  if (b->convex) { // -c option.
    // Assume the mesh is convex. Exterior tets have region attribute -1.
    assert(in->numberoftetrahedronattributes > 0);
  } else {
    // Assume the mesh is non-convex.
    nonconvex = 1;
  }

  // Create a map from indices to vertices. 
  makeindex2pointmap(idx2verlist);
  // 'idx2verlist' has length 'in->numberofpoints + 1'.
  if (in->firstnumber == 1) {
    idx2verlist[0] = dummypoint; // Let 0th-entry be dummypoint.
  }

  // Allocate an array that maps each vertex to its adjacent tets.
  ver2tetarray = new tetrahedron[in->numberofpoints + 1];
  //for (i = 0; i < in->numberofpoints + 1; i++) {
  for (i = in->firstnumber; i < in->numberofpoints + in->firstnumber; i++) {
    setpointtype(idx2verlist[i], VOLVERTEX); // initial type.
    ver2tetarray[i] = NULL;
  }

  // Create the tetrahedra and connect those that share a common face.
  for (i = 0; i < in->numberoftetrahedra; i++) {
    // Get the four vertices.
    idx = i * in->numberofcorners;
    for (j = 0; j < 4; j++) {
      p[j] = idx2verlist[in->tetrahedronlist[idx++]];
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
          // A hull face.
          if (!b->convex) {
            bondflag = 1;  // Insert a hull subface.
          }
        } else {
          if (eextras > 0) {
            if (elemattribute(tetloop.tet, eextras - 1) !=
                elemattribute(checktet.tet, eextras - 1)) {
              bondflag = 1; // Insert an interior interface.
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


  // Segments will be introduced. 
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
          spivot(subloop, neighsh);
          if (shellmark(subloop) != shellmark(neighsh)) {
            // It's an interior interface. Insert a segment.
            p[0] = sorg(subloop);
            p[1] = sdest(subloop);
            bondflag = 1;
          } else {
            if (!b->convex) {
              // Check the dihedral angle formed by the two subfaces.
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
          }
        }
        if (bondflag) {
          // Create a new segment.
          makeshellface(subsegs, &segloop);
          setshvertices(segloop, p[0], p[1], NULL);
          // Create the point-to-segment map.
          sptr = sencode(segloop);
          for (j = 0; j < 2; j++) {
            setpointtype(p[j], RIDGEVERTEX); // initial type.
            setpoint2sh(p[j], sptr);
          }
          setshellmark(segloop, 0); // Initially has no marker.
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
    } // i
    subloop.sh = shellfacetraverse(subfaces);
  }

  // Remember the number of input segments.
  insegments = subsegs->items;

  if (!b->nobisect || checkconstraints) {
    // Mark Steiner points on segments and facets.
    //   - all vertices which remaining type FEACTVERTEX become
    //     Steiner points in facets (= FREEFACERVERTEX).
    //   - vertices on segment need to be checked.
    face* segperverlist;
    int* idx2seglist;
    face parentseg, nextseg;
    verttype vt;
    REAL area, len, l1, l2;
    int fmarker;

    makepoint2submap(subsegs, idx2seglist, segperverlist);

    points->traversalinit();
    point ptloop = pointtraverse();
    while (ptloop != NULL) {
      vt = pointtype(ptloop);
      if (vt == VOLVERTEX) {
        setpointtype(ptloop, FREEVOLVERTEX);
        st_volref_count++;
      } else if (vt == FACETVERTEX) {
        setpointtype(ptloop, FREEFACETVERTEX);
        st_facref_count++;
      } else if (vt == RIDGEVERTEX) {
        idx = pointmark(ptloop) - in->firstnumber;
        if ((idx2seglist[idx + 1] - idx2seglist[idx]) == 2) {
          i = idx2seglist[idx];
          parentseg = segperverlist[i];
          nextseg = segperverlist[i + 1];
          sesymself(nextseg);
          p[0] = sorg(nextseg);
          p[1] = sdest(parentseg);
          // Check if three points p[0], ptloop, p[2] are (nearly) collinear.
          len = distance(p[0], p[1]);
          l1 = distance(p[0], ptloop);
          l2 = distance(ptloop, p[1]);
          if (((l1 + l2 - len) / len) < b->epsilon) {
            // They are (nearly) collinear.
            setpointtype(ptloop, FREESEGVERTEX);
            // Connect nextseg and parentseg together at ptloop.
            senextself(nextseg);
            senext2self(parentseg);
            sbond(nextseg, parentseg);
            st_segref_count++;
          }
        }
      }
      ptloop = pointtraverse();
    }

    // Are there area constraints?
    if (b->quality && (in->facetconstraintlist != (REAL *) NULL)) {
      // Set maximum area constraints on facets.
      for (i = 0; i < in->numberoffacetconstraints; i++) {
        fmarker = (int) in->facetconstraintlist[i * 2];
        area = in->facetconstraintlist[i * 2 + 1];
        subfaces->traversalinit();
        subloop.sh = shellfacetraverse(subfaces);
        while (subloop.sh != NULL) {
          if (shellmark(subloop) == fmarker) {
            setareabound(subloop, area);
          }
          subloop.sh = shellfacetraverse(subfaces);
        }
      }
    }

    // Are there length constraints?
    if (b->quality && (in->segmentconstraintlist != (REAL *) NULL)) {
      // Set maximum length constraints on segments.
      int e1, e2;
      for (i = 0; i < in->numberofsegmentconstraints; i++) {
        e1 = (int) in->segmentconstraintlist[i * 3];
        e2 = (int) in->segmentconstraintlist[i * 3 + 1];
        len = in->segmentconstraintlist[i * 3 + 2];
        // Search for edge [e1, e2].
        idx = e1 - in->firstnumber;
        for (j = idx2seglist[idx]; j <  idx2seglist[idx + 1]; j++) {
          parentseg = segperverlist[j];
          if (pointmark(sdest(parentseg)) == e2) {
            setareabound(parentseg, len);
            break;
          }
        }
      }
    }

    delete [] idx2seglist;
    delete [] segperverlist;
  }


  // Set global flags.
  checksubsegflag = 1;
  checksubfaceflag = 1;

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
  REAL vol, ori1, ori2 = 0, ori3 = 0, ori4 = 0;
  int t1ver;


  // Randomly select a good starting tet.
  if (randflag) {
    randomsample(searchpt, searchtet);
  } else {
    if (searchtet->tet == NULL) {
      *searchtet = recenttet;
    }
  }
  loc = locate(searchpt, searchtet);

  if (loc == OUTSIDE) {
    if (b->convex) { // -c option
      // The point lies outside of the convex hull.
      return (int) loc;
    }
    // Test if it lies nearly on the hull face.
    // Reuse vol, ori1.
    pa = org(*searchtet);
    pb = dest(*searchtet);
    pc = apex(*searchtet);
    vol = triarea(pa, pb, pc);
    ori1 = orient3dfast(pa, pb, pc, searchpt);
    if (fabs(ori1 / vol) < b->epsilon) {
      loc = ONFACE; // On face (or on edge, or on vertex).
      fsymself(*searchtet);
    }
  }

  if (loc != OUTSIDE) {
    // Round the result of location.
    pa = org(*searchtet);
    pb = dest(*searchtet);
    pc = apex(*searchtet);
    pd = oppo(*searchtet);
    vol = orient3dfast(pa, pb, pc, pd);
    ori1 = orient3dfast(pa, pb, pc, searchpt);
    ori2 = orient3dfast(pb, pa, pd, searchpt);
    ori3 = orient3dfast(pc, pb, pd, searchpt);
    ori4 = orient3dfast(pa, pc, pd, searchpt);
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

      vol = orient3dfast(pa, pb, pc, pd); 
      if (vol < 0) {
        ori1 = orient3dfast(pa, pb, pc, searchpt);
        if (fabs(ori1 / vol) < b->epsilon) ori1 = 0; // Rounding.
        if (ori1 <= 0) {
          ori2 = orient3dfast(pb, pa, pd, searchpt);
          if (fabs(ori2 / vol) < b->epsilon) ori2 = 0;
          if (ori2 <= 0) {
            ori3 = orient3dfast(pc, pb, pd, searchpt);
            if (fabs(ori3 / vol) < b->epsilon) ori3 = 0;
            if (ori3 <= 0) {
              ori4 = orient3dfast(pa, pc, pd, searchpt);
              if (fabs(ori4 / vol) < b->epsilon) ori4 = 0;
              if (ori4 <= 0) {
                // Found the tet. Return its location. 
                break;
              } // ori4
            } // ori3
          } // ori2
        } // ori1
      }

      searchtet->tet = tetrahedrontraverse();
    } // while (searchtet->tet != NULL)
    nonregularcount++;  // Re-use this counter.
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
      volume = orient3dfast(pts[0], pts[1], pts[2], pts[3]);
      vol[0] = orient3dfast(searchpt, pts[1], pts[2], pts[3]);
      vol[1] = orient3dfast(pts[0], searchpt, pts[2], pts[3]);
      vol[2] = orient3dfast(pts[0], pts[1], searchpt, pts[3]);
      vol[3] = orient3dfast(pts[0], pts[1], pts[2], searchpt);
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

  long bak_nonregularcount = nonregularcount;
  nonregularcount = 0l; // Count the number of (slow) global searches.
  long baksmaples = bgm->samples;
  bgm->samples = 3l;
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
    if (nonregularcount > 0l) {
      printf("  Performed %ld brute-force searches.\n", nonregularcount);
    }
    printf("  Size rangle [%.17g, %.17g].\n", minval, maxval);
  }

  bgm->samples = baksmaples;
  nonregularcount = bak_nonregularcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertconstrainedpoints()    Insert a list of points into the mesh.       //
//                                                                           //
// Assumption:  The bounding box of the insert point set should be no larger //
// than the bounding box of the mesh.  (Required by point sorting).          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::insertconstrainedpoints(point *insertarray, int arylen,
                                         int rejflag)
{
  triface searchtet, spintet;
  face splitsh;
  face splitseg;
  insertvertexflags ivf;
  flipconstraints fc;
  int randflag = 0;
  int t1ver;
  int i;

  if (b->verbose) {
    printf("  Inserting %d constrained points\n", arylen);
  }

  if (b->no_sort) { // -b/1 option.
    if (b->verbose) {
      printf("  Using the input order.\n"); 
    }
  } else {
    if (b->verbose) {
      printf("  Permuting vertices.\n"); 
    }
    point swappoint;
    int randindex;
    srand(arylen);
    for (i = 0; i < arylen; i++) {
      randindex = rand() % (i + 1); 
      swappoint = insertarray[i];
      insertarray[i] = insertarray[randindex];
      insertarray[randindex] = swappoint;
    }
    if (b->brio_hilbert) { // -b1 option
      if (b->verbose) {
        printf("  Sorting vertices.\n"); 
      }
      hilbert_init(in->mesh_dim);
      int ngroup = 0; 
      brio_multiscale_sort(insertarray, arylen, b->brio_threshold, 
                           b->brio_ratio, &ngroup);
    } else { // -b0 option.
      randflag = 1;
    } // if (!b->brio_hilbert)
  } // if (!b->no_sort)

  long bak_nonregularcount = nonregularcount;
  nonregularcount = 0l;
  long baksmaples = samples;
  samples = 3l; // Use at least 3 samples. Updated in randomsample().

  long bak_seg_count = st_segref_count;
  long bak_fac_count = st_facref_count;
  long bak_vol_count = st_volref_count;

  // Initialize the insertion parameters. 
  if (b->incrflip) { // -l option
    // Use incremental flip algorithm.
    ivf.bowywat = 0; 
    ivf.lawson = 1;
    ivf.validflag = 0; // No need to validate the cavity.
    fc.enqflag = 2;
  } else {
    // Use Bowyer-Watson algorithm.
    ivf.bowywat = 1; 
    ivf.lawson = 0;
    ivf.validflag = 1; // Validate the B-W cavity.
  }
  ivf.rejflag = rejflag;
  ivf.chkencflag = 0; 
  ivf.sloc = (int) INSTAR;
  ivf.sbowywat = 3; 
  ivf.splitbdflag = 1;
  ivf.respectbdflag = 1;
  ivf.assignmeshsize = b->metric;

  encseglist = new arraypool(sizeof(face), 8);
  encshlist = new arraypool(sizeof(badface), 8);

  // Insert the points.
  for (i = 0; i < arylen; i++) {
    // Find the location of the inserted point.
    // Do not use 'recenttet', since the mesh may be non-convex.
    searchtet.tet = NULL; 
    ivf.iloc = scoutpoint(insertarray[i], &searchtet, randflag);

    // Decide the right type for this point.
    setpointtype(insertarray[i], FREEVOLVERTEX); // Default.
    splitsh.sh = NULL;
    splitseg.sh = NULL;
    if (ivf.iloc == (int) ONEDGE) {
      if (issubseg(searchtet)) {
        tsspivot1(searchtet, splitseg);
        setpointtype(insertarray[i], FREESEGVERTEX);
        //ivf.rejflag = 0;
      } else {
        // Check if it is a subface edge.
        spintet = searchtet;
        while (1) {
          if (issubface(spintet)) {
            tspivot(spintet, splitsh);
            setpointtype(insertarray[i], FREEFACETVERTEX);
            //ivf.rejflag |= 1;
            break;
          }
          fnextself(spintet);
          if (spintet.tet == searchtet.tet) break;
        }
      }
    } else if (ivf.iloc == (int) ONFACE) {
      if (issubface(searchtet)) {
        tspivot(searchtet, splitsh);
        setpointtype(insertarray[i], FREEFACETVERTEX);
        //ivf.rejflag |= 1;
      }
    }

    // Now insert the point.
    if (insertpoint(insertarray[i], &searchtet, &splitsh, &splitseg, &ivf)) {
      if (flipstack != NULL) {
        // There are queued faces. Use flips to recover Delaunayness.
        lawsonflip3d(&fc);
        // There may be unflippable edges. Ignore them.
        unflipqueue->restart();
      }
      // Update the Steiner counters.
      if (pointtype(insertarray[i]) == FREESEGVERTEX) {
        st_segref_count++;
      } else if (pointtype(insertarray[i]) == FREEFACETVERTEX) {
        st_facref_count++;
      } else {
        st_volref_count++;
      }
    } else {
      // Point is not inserted.
      //pointdealloc(insertarray[i]);
      setpointtype(insertarray[i], UNUSEDVERTEX);
      unuverts++;
      encseglist->restart();
      encshlist->restart();
    }
  } // i

  delete encseglist;
  delete encshlist;

  if (b->verbose) {
    printf("  Inserted %ld (%ld, %ld, %ld) vertices.\n", 
           st_segref_count + st_facref_count + st_volref_count - 
           (bak_seg_count + bak_fac_count + bak_vol_count),
           st_segref_count - bak_seg_count, st_facref_count - bak_fac_count,
           st_volref_count - bak_vol_count);
    if (nonregularcount > 0l) {
      printf("  Performed %ld brute-force searches.\n", nonregularcount);
    }
  }

  nonregularcount = bak_nonregularcount;
  samples = baksmaples; 
}

void tetgenmesh::insertconstrainedpoints(tetgenio *addio)
{
  point *insertarray, newpt;
  REAL x, y, z, w;
  int index, attribindex, mtrindex;
  int arylen, i, j;

  if (!b->quiet) {
    printf("Inserting constrained points ...\n");
  }

  insertarray = new point[addio->numberofpoints];
  arylen = 0;
  index = 0;
  attribindex = 0;
  mtrindex = 0;

  for (i = 0; i < addio->numberofpoints; i++) {
    x = addio->pointlist[index++];
    y = addio->pointlist[index++];
    z = addio->pointlist[index++];
    // Test if this point lies inside the bounding box.
    if ((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax) ||
        (z < zmin) || (z > zmax)) {
      if (b->verbose) {
        printf("Warning:  Point #%d lies outside the bounding box. Ignored\n",
               i + in->firstnumber);
      }
      continue;
    }
    makepoint(&newpt, UNUSEDVERTEX);
    newpt[0] = x;
    newpt[1] = y;
    newpt[2] = z;
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
    insertarray[arylen] = newpt;
    arylen++;
  } // i

  // Insert the points.
  int rejflag = 0;  // Do not check encroachment.
  if (b->metric) { // -m option.
    rejflag |= 4; // Reject it if it lies in some protecting balls.
  }

  insertconstrainedpoints(insertarray, arylen, rejflag);

  delete [] insertarray;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// meshcoarsening()    Deleting (selected) vertices.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::collectremovepoints(arraypool *remptlist)
{
  point ptloop, *parypt;
  verttype vt;

  // If a mesh sizing function is given. Collect vertices whose mesh size
  //   is greater than its smallest edge length.
  if (b->metric) { // -m option
    REAL len, smlen;
    int i;
    points->traversalinit();
    ptloop = pointtraverse();
    while (ptloop != NULL) {
      if (ptloop[pointmtrindex] > 0) {
        // Get the smallest edge length at this vertex.
        getvertexstar(1, ptloop, cavetetlist, cavetetvertlist, NULL);
        parypt = (point *) fastlookup(cavetetvertlist, 0);
        smlen = distance(ptloop, *parypt);
        for (i = 1; i < cavetetvertlist->objects; i++) {
          parypt = (point *) fastlookup(cavetetvertlist, i);
          len = distance(ptloop, *parypt);
          if (len < smlen) {
            smlen = len;
          }
        }
        cavetetvertlist->restart();
        cavetetlist->restart();
        if (smlen < ptloop[pointmtrindex]) {
          pinfect(ptloop);
          remptlist->newindex((void **) &parypt);
          *parypt = ptloop;
        }
      }
      ptloop = pointtraverse();
    }
    if (b->verbose > 1) {
      printf("    Coarsen %ld oversized points.\n", remptlist->objects); 
    }
  }

  // If 'in->pointmarkerlist' exists, Collect vertices with markers '-1'.
  if (in->pointmarkerlist != NULL) {
    long bak_count = remptlist->objects;
    points->traversalinit();
    ptloop = pointtraverse();
    int index = 0;
    while (ptloop != NULL) {
      if (index < in->numberofpoints) {
        if (in->pointmarkerlist[index] == -1) {
          pinfect(ptloop);
          remptlist->newindex((void **) &parypt);
          *parypt = ptloop;
        }
      } else {
        // Remaining are not input points. Stop here.
        break; 
      }
      index++;
      ptloop = pointtraverse();
    }
    if (b->verbose > 1) {
      printf("    Coarsen %ld marked points.\n", remptlist->objects - bak_count); 
    }
  } // if (in->pointmarkerlist != NULL)

  if (b->coarsen_param > 0) { // -R1/#
    // Remove a coarsen_percent number of interior points.
    assert((b->coarsen_percent > 0) && (b->coarsen_percent <= 1.0));
    if (b->verbose > 1) {
      printf("    Coarsen %g percent of interior points.\n", 
             b->coarsen_percent * 100.0);
    }
    arraypool *intptlist = new arraypool(sizeof(point *), 10);
    // Count the total number of interior points.
    points->traversalinit();
    ptloop = pointtraverse();
    while (ptloop != NULL) {
      vt = pointtype(ptloop);
      if ((vt == VOLVERTEX) || (vt == FREEVOLVERTEX) || 
          (vt == FREEFACETVERTEX) || (vt == FREESEGVERTEX)) {
        intptlist->newindex((void **) &parypt);
        *parypt = ptloop;
      }
      ptloop = pointtraverse();
    }
    if (intptlist->objects > 0l) {
      // Sort the list of points randomly.
      point *parypt_i, swappt;
      int randindex, i;
      srand(intptlist->objects);
      for (i = 0; i < intptlist->objects; i++) {
        randindex = rand() % (i + 1); // randomnation(i + 1);
        parypt_i = (point *) fastlookup(intptlist, i); 
        parypt = (point *) fastlookup(intptlist, randindex);
        // Swap this two points.
        swappt = *parypt_i;
        *parypt_i = *parypt;
        *parypt = swappt;
      }
      int remcount = (int) ((REAL) intptlist->objects * b->coarsen_percent);
      // Return the first remcount points.
      for (i = 0; i < remcount; i++) {
        parypt_i = (point *) fastlookup(intptlist, i);
        if (!pinfected(*parypt_i)) {
          pinfected(*parypt_i);
          remptlist->newindex((void **) &parypt);
          *parypt = *parypt_i;
        }
      }
    }
    delete intptlist;
  }

  // Unmark all collected vertices.
  for (int i = 0; i < remptlist->objects; i++) {
    parypt = (point *) fastlookup(remptlist, i);
    puninfect(*parypt);
  }
}

void tetgenmesh::meshcoarsening()
{
  arraypool *remptlist;

  if (!b->quiet) {
    printf("Mesh coarsening ...\n");
  }

  // Collect the set of points to be removed
  remptlist = new arraypool(sizeof(point *), 10);
  collectremovepoints(remptlist);

  if (remptlist->objects == 0l) {
    delete remptlist;
    return;
  }

  if (b->verbose) {
    if (remptlist->objects > 0l) {
      printf("  Removing %ld points...\n", remptlist->objects);
    }
  }

  point *parypt, *plastpt;
  long ms = remptlist->objects;
  int nit = 0; 
  int bak_fliplinklevel = b->fliplinklevel;
  b->fliplinklevel = -1;
  autofliplinklevel = 1; // Init value.
  int i;

  while (1) {
  
    if (b->verbose > 1) {
      printf("    Removing points [%s level = %2d] #:  %ld.\n", 
             (b->fliplinklevel > 0) ? "fixed" : "auto",
             (b->fliplinklevel > 0) ? b->fliplinklevel : autofliplinklevel,
             remptlist->objects);
    }

    // Remove the list of points.
    for (i = 0; i < remptlist->objects; i++) {
      parypt = (point *) fastlookup(remptlist, i);
      assert(pointtype(*parypt) != UNUSEDVERTEX);
      if (removevertexbyflips(*parypt)) {
        // Move the last entry to the current place.
        plastpt = (point *) fastlookup(remptlist, remptlist->objects - 1);
        *parypt = *plastpt;
        remptlist->objects--;
        i--;
      }
    }

    if (remptlist->objects > 0l) {
      if (b->fliplinklevel >= 0) {
        break; // We have tried all levels.
      }
      if (remptlist->objects == ms) {
        nit++;
        if (nit >= 3) {
          // Do the last round with unbounded flip link level.
          b->fliplinklevel = 100000;
        }
      } else {
        ms = remptlist->objects;
        if (nit > 0) {
          nit--;
        }
      }
      autofliplinklevel+=b->fliplinklevelinc;
    } else {
      // All points are removed.
      break;
    }
  } // while (1)

  if (remptlist->objects > 0l) {
    if (b->verbose) {
      printf("  %ld points are not removed !\n", remptlist->objects);
    }
  }

  b->fliplinklevel = bak_fliplinklevel;
  delete remptlist;
}

////                                                                       ////
////                                                                       ////
//// reconstruct_cxx //////////////////////////////////////////////////////////

