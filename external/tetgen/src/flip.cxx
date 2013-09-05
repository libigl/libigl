#include "../tetgen.h"
//// flip_cxx /////////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

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
  badface *newflipface;

  if (!facemarked(*flipface)) {
    newflipface = (badface *) flippool->alloc();
    newflipface->tt = *flipface;
    markface(newflipface->tt);
    // Push this face into stack.
    newflipface->nextitem = fstack;
    fstack = newflipface;
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip23()    Perform a 2-to-3 flip (face-to-edge flip).                    //
//                                                                           //
// 'fliptets' is an array of tetrahedra.  On input it contains two tets      //
// [a,b,c,d] and [b,a,c,e]. It returns three new tets: [e,d,a,b], [e,d,b,c], //
// [e,d,c,a]. The face [a,b,c] is removed, and the edge [d,e] is created.    //
//                                                                           //
// If 'hullflag' > 0, hull tets may be involved in this flip, i.e., one of   //
// the five vertices may be 'dummypoint'. There are two canonical cases:     //
//   (1) d is 'dummypoint', then all three new tets are hull tets.  If e is  //
//       'dummypoint', we reconfigure e to d, i.e., turn it up-side down.    //
//   (2) c is 'dummypoint', then two new tets: [e,d,b,c] and [e,d,c,a], are  //
//       hull tets. If a or b is 'dummypoint', we reconfigure it to c, i.e., //
//       rotate the three input tets counterclockwisely (right-hand rule)    //
//       until a or b is in c's position.                                    //
//                                                                           //
// If 'flipflag > 0', faces on the convex hull of the five vertices might    //
// need to be flipped, e.g., for incremental DT construction or mesh quality //
// improvement. They will be queued in 'flipstack'.                          //
//                                                                           //
// If 'flipflag = 1', it is in the process of incrmental flip DT algorithm,  //
// and we assume that 'd' must be the newly inserted vertex.  In such case,  //
// only the link faces at 'd', i.e., three faces [a,b,e], [b,c,e], and [c,a, //
// e] needs to be queued ([Edelsbrunner & Shah'1996] and [M\"ucke'1998]).    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip23(triface* fliptets, int hullflag, int flipflag, 
                        int chkencflag)
{
  triface topcastets[3], botcastets[3];
  triface newface, casface;
  face checksh;
  face checkseg;
  badface *bface; // used by chkencflag
  point pa, pb, pc, pd, pe;
  REAL volneg[2], volpos[3], vol_diff; // volumes of involved tet-prisms.
  REAL attrib, volume;
  int dummyflag = 0;  // range = {-1, 0, 1, 2}.
  int i;

  if (hullflag > 0) {
    // Check if e is dummypoint.
    if (oppo(fliptets[1]) == dummypoint) {
      // Swap the two old tets.
      newface = fliptets[0];
      fliptets[0] = fliptets[1];
      fliptets[1] = newface;
      dummyflag = -1;  // d is dummypoint.
    } else {
      // Check if either a or b is dummypoint.
      if (org(fliptets[0]) == dummypoint) {
        dummyflag = 1; // a is dummypoint.
        enextself(fliptets[0]);
        eprevself(fliptets[1]);
      } else if (dest(fliptets[0]) == dummypoint) {
        dummyflag = 2; // b is dummypoint.
        eprevself(fliptets[0]);
        enextself(fliptets[1]);
      } else {
        dummyflag = 0; // either c or d may be dummypoint.
      }
    }
  }

  pa = org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

  if (b->verbose > 3) {
    printf("        flip 2-to-3: (%d, %d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip23count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    fnext(fliptets[0], topcastets[i]);
    enextself(fliptets[0]);
  }
  for (i = 0; i < 3; i++) {
    fnext(fliptets[1], botcastets[i]);
    eprevself(fliptets[1]);
  }

  // Re-use fliptets[0] and fliptets[1].
  fliptets[0].ver = 11;
  fliptets[1].ver = 11;
  setelemmarker(fliptets[0].tet, 0); // Clear all flags.
  setelemmarker(fliptets[1].tet, 0);
  // NOTE: the element attributes and volume constraint remain unchanged.
  if (checksubsegflag) {
    // Dealloc the space to subsegments.
    if (fliptets[0].tet[8] != NULL) {
      tet2segpool->dealloc((shellface *) fliptets[0].tet[8]);
      fliptets[0].tet[8] = NULL;
    }
    if (fliptets[1].tet[8] != NULL) {
      tet2segpool->dealloc((shellface *) fliptets[1].tet[8]);
      fliptets[1].tet[8] = NULL;
    }
  }
  if (checksubfaceflag) {
    // Dealloc the space to subfaces.
    if (fliptets[0].tet[9] != NULL) {
      tet2subpool->dealloc((shellface *) fliptets[0].tet[9]);
      fliptets[0].tet[9] = NULL;
    }
    if (fliptets[1].tet[9] != NULL) {
      tet2subpool->dealloc((shellface *) fliptets[1].tet[9]);
      fliptets[1].tet[9] = NULL;
    }
  }
  // Create a new tet.
  maketetrahedron(&(fliptets[2]));
  // The new tet have the same attributes from the old tet.
  for (i = 0; i < numelemattrib; i++) {
    attrib = elemattribute(fliptets[0].tet, i);
    setelemattribute(fliptets[2].tet, i, attrib);
  }
  if (b->varvolume) {
    volume = volumebound(fliptets[0].tet);
    setvolumebound(fliptets[2].tet, volume);
  }

  if (hullflag > 0) {
    // Check if d is dummytet.
    if (pd != dummypoint) {
      setvertices(fliptets[0], pe, pd, pa, pb); // [e,d,a,b] *
      setvertices(fliptets[1], pe, pd, pb, pc); // [e,d,b,c] *
      // Check if c is dummypoint.
      if (pc != dummypoint) {
        setvertices(fliptets[2], pe, pd, pc, pa);  // [e,d,c,a] *
      } else {
        setvertices(fliptets[2], pd, pe, pa, pc); // [d,e,a,c]
        esymself(fliptets[2]);                    // [e,d,c,a] *
      }
      // The hullsize does not change.
    } else {
      // d is dummypoint.
      setvertices(fliptets[0], pa, pb, pe, pd); // [a,b,e,d]
      setvertices(fliptets[1], pb, pc, pe, pd); // [b,c,e,d]
      setvertices(fliptets[2], pc, pa, pe, pd); // [c,a,e,d]
      // Adjust the faces to [e,d,a,b], [e,d,b,c], [e,d,c,a] *
      for (i = 0; i < 3; i++) {
        eprevesymself(fliptets[i]);
        enextself(fliptets[i]);
      }
      // We deleted one hull tet, and created three hull tets.
      hullsize += 2;
    }
  } else {
    setvertices(fliptets[0], pe, pd, pa, pb); // [e,d,a,b] *
    setvertices(fliptets[1], pe, pd, pb, pc); // [e,d,b,c] *
    setvertices(fliptets[2], pe, pd, pc, pa); // [e,d,c,a] *
  }

  if (calc_tetprism_vol) {
    if (pd != dummypoint) {
      if (pc != dummypoint) {
        volpos[0] = tetprismvol(pe, pd, pa, pb);
        volpos[1] = tetprismvol(pe, pd, pb, pc);
        volpos[2] = tetprismvol(pe, pd, pc, pa);
        volneg[0] = tetprismvol(pa, pb, pc, pd);
        volneg[1] = tetprismvol(pb, pa, pc, pe);
      } else { // pc == dummypoint
        volpos[0] = tetprismvol(pe, pd, pa, pb);
        volpos[1] = 0.;
        volpos[2] = 0.;
        volneg[0] = 0.;
        volneg[1] = 0.;
      }
    } else { // pd == dummypoint.
      volpos[0] = 0.;
      volpos[1] = 0.;
      volpos[2] = 0.;
      volneg[0] = 0.;
      volneg[1] = tetprismvol(pb, pa, pc, pe);
    }
    vol_diff = volpos[0] + volpos[1] + volpos[2] - volneg[0] - volneg[1];
    tetprism_vol_sum  += vol_diff; // Update the total sum.
  } // if (check_tetprism_vol_diff)

  // Bond three new tets together.
  for (i = 0; i < 3; i++) {
    esym(fliptets[i], newface);
    bond(newface, fliptets[(i + 1) % 3]);
  }
  // Bond to top outer boundary faces (at [a,b,c,d]).
  for (i = 0; i < 3; i++) {
    enextesym(fliptets[i], newface);
    eprevself(newface); // At edges [b,a], [c,b], [a,c].
    bond(newface, topcastets[i]);
  }
  // Bond bottom outer boundary faces (at [b,a,c,e]).
  for (i = 0; i < 3; i++) {
    eprevesym(fliptets[i], newface);
    enextself(newface); // At edges [a,b], [b,c], [c,a].
    bond(newface, botcastets[i]);
  }

  // Bond 15 subsegments if there are.
  if (checksubsegflag) {
    // The middle three: [a,b], [b,c], [c,a].
    for (i = 0; i < 3; i++) {
      tsspivot1(topcastets[i], checkseg);
      if (checkseg.sh != NULL) {
        enextesym(fliptets[i], newface);
        eprevself(newface); // At edges [b,a], [c,b], [a,c].
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
    }
    // The top three: [d,a], [d,b], [d,c]. Two tets per edge.
    for (i = 0; i < 3; i++) {
      eprev(topcastets[i], casface);
      tsspivot1(casface, checkseg);
      if (checkseg.sh != NULL) {
        enext(fliptets[i], newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        esym(fliptets[(i + 2) % 3], newface);
        eprevself(newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
    }
    // The bot three: [a,e], [b,e], [c,e]. Two tets per edge.
    for (i = 0; i < 3; i++) {
      enext(botcastets[i], casface);
      tsspivot1(casface, checkseg);
      if (checkseg.sh != NULL) {
        eprev(fliptets[i], newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        esym(fliptets[(i + 2) % 3], newface);
        enextself(newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
    }
  }

  // Bond 6 subfaces if there are.
  if (checksubfaceflag) {
    for (i = 0; i < 3; i++) {
      tspivot(topcastets[i], checksh);
      if (checksh.sh != NULL) {
        enextesym(fliptets[i], newface);
        eprevself(newface); // At edge [b,a], [c,b], [a,c].
        sesymself(checksh);
        tsbond(newface, checksh);
        if (chkencflag & 2) {
          if (!smarktest2ed(checksh)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = checksh;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checksh); // An alive badface
          }
        }
      }
    }
    for (i = 0; i < 3; i++) {
      tspivot(botcastets[i], checksh);
      if (checksh.sh != NULL) {
        eprevesym(fliptets[i], newface);
        enextself(newface); // At edge [a,b], [b,c], [c,a]
        sesymself(checksh);
        tsbond(newface, checksh);
        if (chkencflag & 2) {
          if (!smarktest2ed(checksh)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = checksh;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checksh); // An alive badface
          }
        }
      }
    }
  }

  if (chkencflag & 4) {
    // Put three new tets into check list.
    for (i = 0; i < 3; i++) {
      if (!marktest2ed(fliptets[i])) {
        bface = (badface *) badtetrahedrons->alloc();
        bface->tt = fliptets[i];
        marktest2(bface->tt);
        bface->forg = org(fliptets[i]);
      }
    }
  }

  // Update the point-to-tet map.
  setpoint2tet(pa, encode(fliptets[0]));
  setpoint2tet(pb, encode(fliptets[0]));
  setpoint2tet(pc, encode(fliptets[1]));
  setpoint2tet(pd, encode(fliptets[0]));
  setpoint2tet(pe, encode(fliptets[0]));

  if (hullflag > 0) {
    if (dummyflag != 0) {
      // Restore the original position of the points (for flipnm()).
      if (dummyflag == -1) { 
        // Reverse the edge.
        for (i = 0; i < 3; i++) {
          esymself(fliptets[i]);
        }
        // Swap the last two new tets.
        newface = fliptets[1];
        fliptets[1] = fliptets[2];
        fliptets[2] = newface;
      } else {
        // either a or b were swapped.
        if (dummyflag == 1) {
          // a is dummypoint.
          newface = fliptets[0];
          fliptets[0] = fliptets[2];
          fliptets[2] = fliptets[1];
          fliptets[1] = newface;
        } else { // dummyflag == 2
          // b is dummypoint.
          newface = fliptets[0];
          fliptets[0] = fliptets[1];
          fliptets[1] = fliptets[2];
          fliptets[2] = newface;
        }
      }
    }
  }

  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    for (i = 0; i < 3; i++) {
      eprevesym(fliptets[i], newface);
      flippush(flipstack, &newface);
    }
    if (flipflag > 1) {
      for (i = 0; i < 3; i++) {
        enextesym(fliptets[i], newface);
        flippush(flipstack, &newface);
      }
    }
  }

  recenttet = fliptets[0];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip32()    Perform a 3-to-2 flip (edge-to-face flip).                    //
//                                                                           //
// 'fliptets' is an array of three tetrahedra. On input, it contains three   //
// tets: [e,d,a,b], [e,d,b,c], and [e,d,c,a]. It returns tw tets: [a,b,c,d], //
// and [b,a,c,e]. The edge [e,d] is replaced by the face [a,b,c].            //
//                                                                           //
// If 'hullflag' > 0, hull tets may be involved in this flip, i.e., one of   //
// the five vertices may be 'dummypoint'. There are two canonical cases:     //
//   (1) d is 'dummypoint', then [a,b,c,d] is hull tet. If e is 'dummypoint',//
//       we reconfigure e to d, i.e., turnover it.                           //
//   (2) c is 'dummypoint' then both [a,b,c,d] and [b,a,c,e] are hull tets.  //
//       If a or b is 'dummypoint', we reconfigure it to c, i.e., rotate the //
//       three old tets counterclockwisely (right-hand rule) until a or b    //
//       is in c's position.                                                 //
//                                                                           //
// If 'flipflag > 0', faces on the convex hull of the five vertices might    //
// need to be flipped, e.g., for incremental DT construction or mesh quality //
// improvement. They will be queued in 'flipstack'.                          //
//                                                                           //
// If 'flipflag = 1', it is in the process of incrmental flip DT algorithm,  //
// and we assume that 'a' must be the newly inserted vertex.  In such case,  //
// only the link faces at 'a', i.e., two faces [c,b,d] and [b,c,e] needs to  //
// be queued ( [Edelsbrunner & Shah'1996] and [M\"ucke'1998]).               //
//                                                                           //
// If 'checksubfaceflag' is on (global variable), and assume [e,d] is not a  //
// segment. There may be two (interior) subfaces sharing at [e,d], which are //
// [e,d,p] and [e,d,q], where the pair (p,q) may be either (a,b), or (b,c),  //
// or (c,a)  In such case, a 2-to-2 flip is performed on these two subfaces  //
// and two new subfaces [p,q,e] and [p,q,d] are created. They are inserted   //
// back into the tetrahedralization.  However, it is possible that the new   //
// subface ([p,q,e] or [p,q,d] already exists. In such case, we just delete  //
// the conflict subface. As a result, either 'd' or 'e' is removed from the  //
// surface mesh.  A better solution would be to detect and perform a 3-to-1  //
// flip to remove 'd' or 'e' (see also 2011-11-15).                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip32(triface* fliptets, int hullflag, int flipflag,
                        int chkencflag)
{
  triface topcastets[3], botcastets[3];
  triface newface, casface;
  face checksh; 
  face checkseg; 
  badface *bface; // used by chkencflag
  point pa, pb, pc, pd, pe;
  REAL volneg[3], volpos[2], vol_diff; // volumes of involved tet-prisms.
  REAL attrib, volume;
  int dummyflag = 0;  // Rangle = {-1, 0, 1, 2}
  int i, j;

  // For 2-to-2 flip (subfaces).
  face flipshs[3], flipfaces[2];
  point rempt;
  int spivot = -1, scount = 0;

  if (hullflag > 0) {
    // Check if e is 'dummypoint'.
    if (org(fliptets[0]) == dummypoint) {
      // Reverse the edge.
      for (i = 0; i < 3; i++) {
        esymself(fliptets[i]);
      }
      // Swap the last two tets.
      newface = fliptets[1];
      fliptets[1] = fliptets[2];
      fliptets[2] = newface;
      dummyflag = -1; // e is dummypoint.
    } else {
      // Check if a or b is the 'dummypoint'.
      if (apex(fliptets[0]) == dummypoint) { 
        dummyflag = 1;  // a is dummypoint.
        newface = fliptets[0];
        fliptets[0] = fliptets[1];
        fliptets[1] = fliptets[2];
        fliptets[2] = newface;
      } else if (apex(fliptets[1]) == dummypoint) {
        dummyflag = 2;  // b is dummypoint.
        newface = fliptets[0];
        fliptets[0] = fliptets[2];
        fliptets[2] = fliptets[1];
        fliptets[1] = newface;
      } else {
        dummyflag = 0;  // either c or d may be dummypoint.
      }
    }
  }

  pa = apex(fliptets[0]);
  pb = apex(fliptets[1]);
  pc = apex(fliptets[2]);
  pd = dest(fliptets[0]);
  pe = org(fliptets[0]);

  if (b->verbose > 3) {
    printf("        flip 3-to-2: (%d, %d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pe));
  }
  flip32count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    enextesym(fliptets[i], casface);
    eprevself(casface);
    fsym(casface, topcastets[i]);
  }
  for (i = 0; i < 3; i++) {
    eprevesym(fliptets[i], casface);
    enextself(casface);
    fsym(casface, botcastets[i]);
  }

  if (checksubfaceflag) {
    // Check if there are interior subfaces at the edge [e,d].
    spivot = -1;
    scount = 0;
    for (i = 0; i < 3; i++) {
      tspivot(fliptets[i], flipshs[i]);
      if (flipshs[i].sh != NULL) {
        if (b->verbose > 3) {
          printf("        Found an interior subface (%d, %d, %d).\n",
                 pointmark(sorg(flipshs[i])), pointmark(sdest(flipshs[i])),
                 pointmark(sapex(flipshs[i])));
        }
        stdissolve(flipshs[i]); // Disconnect the sub-tet bond.
        scount++;
      } else {
        spivot = i;
      }
    }
  }

  // Re-use fliptets[0] and fliptets[1].
  fliptets[0].ver = 11;
  fliptets[1].ver = 11;
  setelemmarker(fliptets[0].tet, 0); // Clear all flags.
  setelemmarker(fliptets[1].tet, 0);
  // NOTE: the element attributes and volume constraint must be set correctly.
  if (checksubfaceflag) {
    if (scount > 0) {
      // There are two subfaces involved in this flip. The three tets are
      //   separated into two different regions, one may be exterior. The
      //   first region has two tets, and the second region has only one.
      //   The two created tets must be in the same region as the first region. 
      //   The element attributes and volume constraint must be set correctly.
      //assert(spivot != -1);
      // The tet fliptets[spivot] is in the first region.
      for (j = 0; j < 2; j++) {
        for (i = 0; i < numelemattrib; i++) {
          attrib = elemattribute(fliptets[spivot].tet, i);
          setelemattribute(fliptets[j].tet, i, attrib);
        }
        if (b->varvolume) {
          volume = volumebound(fliptets[spivot].tet);
          setvolumebound(fliptets[j].tet, volume);
        }
      }
    }
  }
  if (checksubsegflag) {
    // Dealloc the space to subsegments.
    if (fliptets[0].tet[8] != NULL) {
      tet2segpool->dealloc((shellface *) fliptets[0].tet[8]);
      fliptets[0].tet[8] = NULL;
    }
    if (fliptets[1].tet[8] != NULL) {
      tet2segpool->dealloc((shellface *) fliptets[1].tet[8]);
      fliptets[1].tet[8] = NULL;
    }
  }
  if (checksubfaceflag) {
    // Dealloc the space to subfaces.
    if (fliptets[0].tet[9] != NULL) {
      tet2subpool->dealloc((shellface *) fliptets[0].tet[9]);
      fliptets[0].tet[9] = NULL;
    }
    if (fliptets[1].tet[9] != NULL) {
      tet2subpool->dealloc((shellface *) fliptets[1].tet[9]);
      fliptets[1].tet[9] = NULL;
    }
  }

  // Delete an old tet.
  tetrahedrondealloc(fliptets[2].tet);

  if (hullflag > 0) {
    // Check if c is dummypointc.
    if (pc != dummypoint) {
      // Check if d is dummypoint.
      if (pd != dummypoint) {
        // No hull tet is involved.
      } else {
        // We deleted three hull tets, and created one hull tet.
        hullsize -= 2;
      }
      setvertices(fliptets[0], pa, pb, pc, pd);
      setvertices(fliptets[1], pb, pa, pc, pe);
    } else {
      // c is dummypoint. The two new tets are hull tets.
      setvertices(fliptets[0], pb, pa, pd, pc);
      setvertices(fliptets[1], pa, pb, pe, pc);
      // Adjust badc -> abcd.
      esymself(fliptets[0]);
      // Adjust abec -> bace.
      esymself(fliptets[1]);
      // The hullsize does not changle.
    }
  } else {
    setvertices(fliptets[0], pa, pb, pc, pd);
    setvertices(fliptets[1], pb, pa, pc, pe);
  }

  if (calc_tetprism_vol) {
    if (pc != dummypoint) {
      if (pd != dummypoint) {
        volneg[0] = tetprismvol(pe, pd, pa, pb);
        volneg[1] = tetprismvol(pe, pd, pb, pc);
        volneg[2] = tetprismvol(pe, pd, pc, pa);
        volpos[0] = tetprismvol(pa, pb, pc, pd);
        volpos[1] = tetprismvol(pb, pa, pc, pe);
      } else { // pd == dummypoint
        volneg[0] = 0.;
        volneg[1] = 0.;
        volneg[2] = 0.;
        volpos[0] = 0.;
        volpos[1] = tetprismvol(pb, pa, pc, pe);
      }
    } else { // pc == dummypoint.
      volneg[0] = tetprismvol(pe, pd, pa, pb);
      volneg[1] = 0.;
      volneg[2] = 0.;
      volpos[0] = 0.;
      volpos[1] = 0.;
    }
    vol_diff = volpos[0] + volpos[1] - volneg[0] - volneg[1] - volneg[2];
    tetprism_vol_sum  += vol_diff; // Update the total sum.
  }

  // Bond abcd <==> bace.
  bond(fliptets[0], fliptets[1]);
  // Bond new faces to top outer boundary faces (at abcd).
  for (i = 0; i < 3; i++) {
    esym(fliptets[0], newface);
    bond(newface, topcastets[i]);
    enextself(fliptets[0]);
  }
  // Bond new faces to bottom outer boundary faces (at bace).
  for (i = 0; i < 3; i++) {
    esym(fliptets[1], newface);
    bond(newface, botcastets[i]);
    eprevself(fliptets[1]);
  }

  if (checksubsegflag) {
    // Bond segments to new (flipped) tets.
    for (i = 0; i < 3; i++) {
      tsspivot1(topcastets[i], checkseg);
      if (checkseg.sh != NULL) {
        tssbond1(fliptets[0], checkseg);
        sstbond1(checkseg, fliptets[0]);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      enextself(fliptets[0]);
    }
    // The three top edges.
    for (i = 0; i < 3; i++) {
      esym(fliptets[0], newface);
      eprevself(newface); // edge b->d, c->d, a->d.
      enext(topcastets[i], casface);
      tsspivot1(casface, checkseg);
      if (checkseg.sh != NULL) {
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      enextself(fliptets[0]);
    }
    // Process the bottom tet bace.
    for (i = 0; i < 3; i++) {
      tsspivot1(botcastets[i], checkseg);
      if (checkseg.sh != NULL) {
        tssbond1(fliptets[1], checkseg);
        sstbond1(checkseg, fliptets[1]);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      eprevself(fliptets[1]);
    }
    // The three bot edges.
    for (i = 0; i < 3; i++) {
      esym(fliptets[1], newface);
      enextself(newface); // edge b<-e, c<-e, a<-e.
      eprev(botcastets[i], casface);
      tsspivot1(casface, checkseg);
      if (checkseg.sh != NULL) {
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      eprevself(fliptets[1]);
    }
  }

  if (checksubfaceflag) {
    // Bond the top three casing subfaces.
    for (i = 0; i < 3; i++) {
      tspivot(topcastets[i], checksh);
      if (checksh.sh != NULL) {
        esym(fliptets[0], newface); // At edge [b,a], [c,b], [a,c]
        sesymself(checksh);
        tsbond(newface, checksh);
        if (chkencflag & 2) {
          if (!smarktest2ed(checksh)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = checksh;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checksh); // An alive badface
          }
        }
      }
      enextself(fliptets[0]);
    }
    // Bond the bottom three casing subfaces.
    for (i = 0; i < 3; i++) {
      tspivot(botcastets[i], checksh);
      if (checksh.sh != NULL) {
        esym(fliptets[1], newface); // // At edge [a,b], [b,c], [c,a]
        sesymself(checksh);
        tsbond(newface, checksh);
        if (chkencflag & 2) {
          if (!smarktest2ed(checksh)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = checksh;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checksh); // An alive badface
          }
        }
      }
      eprevself(fliptets[1]);
    }
  }

  if (checksubfaceflag) {
    if (scount > 0) {
      assert(spivot != -1); // spivot = i, in {0,1,2}
      // Perform a 2-to-2 flip in subfaces.
      flipfaces[0] = flipshs[(spivot + 1) % 3];
      flipfaces[1] = flipshs[(spivot + 2) % 3];
      sesymself(flipfaces[1]);
      flip22(flipfaces, 0, chkencflag);
      // Connect the flipped subfaces to flipped tets.
      // First go to the corresponding flipping edge.
      //   Re-use top- and botcastets[0].
      topcastets[0] = fliptets[0];
      botcastets[0] = fliptets[1];
      for (i = 0; i < ((spivot + 1) % 3); i++) {
        enextself(topcastets[0]);
        eprevself(botcastets[0]);
      }
      // Connect the top subface to the top tets.
      esymself(topcastets[0]);
      sesymself(flipfaces[0]);
      // Check if there already exists a subface.
      tspivot(topcastets[0], checksh);
      if (checksh.sh == NULL) {
        tsbond(topcastets[0], flipfaces[0]);
        fsymself(topcastets[0]);
        sesymself(flipfaces[0]);
        tsbond(topcastets[0], flipfaces[0]);
      } else {
        // Found two subfaces are duplicated at the same tet face. 
        //   Due to the same reason explained below.
        assert(sapex(checksh) == sapex(flipfaces[0]));
        sspivot(checksh, checkseg);
        assert(checkseg.sh == NULL);
        // Delete the two duplicated subfaces.
        rempt = sapex(checksh);
        if (b->verbose > 2) {
          printf("      Remove vertex %d from surface.\n", pointmark(rempt));
        }
        // Make sure we do not delete a Steiner points in segment.
        assert(pointtype(rempt) == FREEFACETVERTEX);
        setpointtype(rempt, FREEVOLVERTEX);
        // Re-use flipshs.
        //spivot(checksh, flipshs[0]);
        flipshs[0] = checksh; 
        spivotself(flipshs[0]);
        if (flipshs[0].sh == flipfaces[0].sh) {
          sesym(checksh, flipshs[0]);
          spivotself(flipshs[0]);
        }
        assert(flipshs[0].sh != flipfaces[0].sh);
        //spivot(flipfaces[0], flipshs[1]);
        flipshs[1] = flipfaces[0];
        spivotself(flipshs[1]);
        if (flipshs[1].sh == checksh.sh) {
          sesym(flipfaces[0], flipshs[1]);
          spivotself(flipshs[1]);
        }
        assert(flipshs[1].sh != checksh.sh);
        // Bond the two subfaces together.
        sbond(flipshs[0], flipshs[1]);
        // Detach 'checksh' from the adjacent tets.
        tsdissolve(topcastets[0]);
        fsymself(topcastets[0]);
        tsdissolve(topcastets[0]);
        // Delete the two duplicated subfaces.
        shellfacedealloc(subfaces, checksh.sh);
        shellfacedealloc(subfaces, flipfaces[0].sh);
      }
      // // Push topcastets[0] into queue for checking new sliver.
      // assert(oppo(topcastets[0]) != dummypoint);
      // flippush(&(topcastets[0]), oppo(topcastets[0]));
      // Connect the bot subface to the bottom tets.
      esymself(botcastets[0]);
      sesymself(flipfaces[1]);
      // Check if there already exists a subface.
      tspivot(botcastets[0], checksh);
      if (checksh.sh == NULL) {
        tsbond(botcastets[0], flipfaces[1]);
        fsymself(botcastets[0]);
        sesymself(flipfaces[1]);
        tsbond(botcastets[0], flipfaces[1]);
      } else {
        // Found two subfaces are duplicated at the same tet face. 
        assert(sapex(checksh) == sapex(flipfaces[1]));
        // This happens in case when a Steiner point is not exactly coplanar
        //   or collinear with the subface or subedge where it was added.
        //   See figs illustrated in 2011-11-09.
        sspivot(checksh, checkseg);
        assert(checkseg.sh == NULL); 
        // Since the edge [p,q] is not a segment, both subfaces must be 
        //   removed. The effect is that the Steiner point is removed from
        //   the surface triangulation.
        // Delete the two duplicated subfaces.
        rempt = sapex(checksh);
        if (b->verbose > 2) {
          printf("      Remove vertex %d from surface.\n", pointmark(rempt));
        }
        // Make sure we do not delete a Steiner points in segment.
        assert(pointtype(rempt) == FREEFACETVERTEX);
        setpointtype(rempt, FREEVOLVERTEX);
        // Re-use flipshs.
        //spivot(checksh, flipshs[0]);
        flipshs[0] = checksh;
        spivotself(flipshs[0]);
        if (flipshs[0].sh == flipfaces[1].sh) {
          sesym(checksh, flipshs[0]);
          spivotself(flipshs[0]);
        }
        assert(flipshs[0].sh != flipfaces[1].sh);
        //spivot(flipfaces[1], flipshs[1]);
        flipshs[1] = flipfaces[1];
        spivotself(flipshs[1]);
        if (flipshs[1].sh == checksh.sh) {
          sesym(flipfaces[1], flipshs[1]);
          spivotself(flipshs[1]);
        }
        assert(flipshs[1].sh != checksh.sh);
        // Bond the two subfaces together.
        sbond(flipshs[0], flipshs[1]);
        // Detach 'checksh' from the adjacent tets.
        tsdissolve(botcastets[0]);
        fsymself(botcastets[0]);
        tsdissolve(botcastets[0]);
        // Delete the two duplicated subfaces.
        shellfacedealloc(subfaces, checksh.sh);
        shellfacedealloc(subfaces, flipfaces[1].sh);
      }
      // // Push botcastets[0] into queue for checking new sliver.
      // assert(oppo(botcastets[0]) != dummypoint);
      // flippush(&(botcastets[0]), oppo(botcastets[0]));
    }
  }

  if (chkencflag & 4) {
    // Put two new tets into check list.
    for (i = 0; i < 2; i++) {
      if (!marktest2ed(fliptets[i])) {
        bface = (badface *) badtetrahedrons->alloc();
        bface->tt = fliptets[i];
        marktest2(bface->tt);
        bface->forg = org(fliptets[i]);
      }
    }
  }

  setpoint2tet(pa, encode(fliptets[0]));
  setpoint2tet(pb, encode(fliptets[0]));
  setpoint2tet(pc, encode(fliptets[0]));
  setpoint2tet(pd, encode(fliptets[0]));
  setpoint2tet(pe, encode(fliptets[1]));

  if (hullflag > 0) {
    if (dummyflag != 0) {
      // Restore the original position of the points (for flipnm()).
      if (dummyflag == -1) {
        // e were dummypoint. Swap the two new tets.
        newface = fliptets[0];
        fliptets[0] = fliptets[1];
        fliptets[1] = newface;
      } else {
        // a or b was dummypoint.
        if (dummyflag == 1) {
          eprevself(fliptets[0]);
          enextself(fliptets[1]);
        } else { // dummyflag == 2
          enextself(fliptets[0]);
          eprevself(fliptets[1]);
        }
      }
    }
  }
  
  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.
    // pa = org(fliptets[0]); // 'a' may be a new vertex.
    enextesym(fliptets[0], newface);
    flippush(flipstack, &newface);
    eprevesym(fliptets[1], newface);
    flippush(flipstack, &newface);
    if (flipflag > 1) {
      //pb = dest(fliptets[0]);
      eprevesym(fliptets[0], newface);
      flippush(flipstack, &newface);
      enextesym(fliptets[1], newface);
      flippush(flipstack, &newface);
      //pc = apex(fliptets[0]);
      esym(fliptets[0], newface);
      flippush(flipstack, &newface);
      esym(fliptets[1], newface);
      flippush(flipstack, &newface);
    }
  }

  recenttet = fliptets[0];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip41()    Perform a 4-to-1 flip (Remove a vertex).                      //
//                                                                           //
// 'fliptets' is an array of four tetrahedra in the star of the removing     //
// vertex 'p'. Let the four vertices in the star of p be a, b, c, and d. The //
// four tets in 'fliptets' are: [p,d,a,b], [p,d,b,c], [p,d,c,a], and [a,b,c, //
// p].  On return, 'fliptets[0]' is the new tet [a,b,c,d].                   //
//                                                                           //
// If 'hullflag' is set (> 0), one of the four vertices may be 'duumypoint'. //
// The 'hullsize' may be changed.                                            //
//                                                                           //
// If 'checksubface' flag is set (>0), it is possible that there are three   //
// interior subfaces connecting at p.  If so, a 3-to-1 flip is performed to  //
// to remove p from the surface triangulation.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip41(triface* fliptets, int hullflag, int flipflag,
                        int chkencflag)
{
  triface topcastets[3], botcastet;
  triface newface, neightet;
  face flipshs[4];
  face checksh;
  face checkseg;
  point pa, pb, pc, pd, pp;
  badface *bface; // used by chkencflag
  REAL volneg[4], volpos[1], vol_diff; // volumes of involved tet-prisms.
  int dummyflag = 0; // in {0, 1, 2, 3, 4}
  int spivot = -1, scount = 0;
  int i;

  pa =  org(fliptets[3]);
  pb = dest(fliptets[3]);
  pc = apex(fliptets[3]);
  pd = dest(fliptets[0]);
  pp =  org(fliptets[0]); // The removing vertex.

  if (b->verbose > 3) {
    printf("        flip 4-to-1: (%d, %d, %d, %d, %d)\n", pointmark(pa),
           pointmark(pb), pointmark(pc), pointmark(pd), pointmark(pp));
  }
  // flip41count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    enext(fliptets[i], topcastets[i]);
    fnextself(topcastets[i]); // [d,a,b,#], [d,b,c,#], [d,c,a,#]
    enextself(topcastets[i]); // [a,b,d,#], [b,c,d,#], [c,a,d,#]
  }
  fsym(fliptets[3], botcastet); // [b,a,c,#]

  if (checksubfaceflag) {
    // Check if there are three subfaces at 'p'.
    //   Re-use 'newface'.
    spivot = -1;
    scount = 0;
    for (i = 0; i < 3; i++) {
      fnext(fliptets[3], newface); // [a,b,p,d],[b,c,p,d],[c,a,p,d].
      tspivot(newface, flipshs[i]);
      if (flipshs[i].sh != NULL) {
        spivot = i; // Remember this subface.
        scount++;
      }
      enextself(fliptets[3]);
    }
    if (scount > 0) {
      // There are three subfaces connecting at p.
      if (scount < 3) {
        // The new subface is one of {[a,b,d], [b,c,d], [c,a,d]}.
        assert(scount == 1); // spivot >= 0
        // Go to the tet containing the three subfaces.
        fsym(topcastets[spivot], neightet);
        // Get the three subfaces connecting at p.
        for (i = 0; i < 3; i++) {
          esym(neightet, newface);
          tspivot(newface, flipshs[i]);
          assert(flipshs[i].sh != NULL);
          eprevself(neightet);
        }
      } else {
        spivot = 3; // The new subface is [a,b,c].
      }
    }
  } // if (checksubfaceflag)

  // Re-use fliptets[0] for [a,b,c,d].
  fliptets[0].ver = 11;
  setelemmarker(fliptets[0].tet, 0); // Clean all flags.
  // NOTE: the element attributes and volume constraint remain unchanged.
  if (checksubsegflag) {
    // Dealloc the space to subsegments.
    if (fliptets[0].tet[8] != NULL) {
      tet2segpool->dealloc((shellface *) fliptets[0].tet[8]);
      fliptets[0].tet[8] = NULL;
    }
  }
  if (checksubfaceflag) {
    // Dealloc the space to subfaces.
    if (fliptets[0].tet[9] != NULL) {
      tet2subpool->dealloc((shellface *) fliptets[0].tet[9]);
      fliptets[0].tet[9] = NULL;
    }
  }

  // Delete the other three tets.
  for (i = 1; i < 4; i++) {
    tetrahedrondealloc(fliptets[i].tet);
  }

  // Mark the point pp as unused.
  setpointtype(pp, UNUSEDVERTEX);
  unuverts++;

  // Create the new tet [a,b,c,d].
  if (hullflag > 0) {
    // One of the four vertices may be 'dummypoint'.
    if (pa == dummypoint) {
      // pa is dummypoint.
      setvertices(fliptets[0], pc, pb, pd, pa);
      esymself(fliptets[0]);  // [b,c,a,d]
      eprevself(fliptets[0]); // [a,b,c,d]
      dummyflag = 1;
    } else if (pb == dummypoint) {
      setvertices(fliptets[0], pa, pc, pd, pb);
      esymself(fliptets[0]);  // [c,a,b,d]
      enextself(fliptets[0]); // [a,b,c,d]
      dummyflag = 2;
    } else if (pc == dummypoint) {
      setvertices(fliptets[0], pb, pa, pd, pc);
      esymself(fliptets[0]);  // [a,b,c,d]
      dummyflag = 3;
    } else if (pd == dummypoint) {
      setvertices(fliptets[0], pa, pb, pc, pd);
      dummyflag = 4;
    } else {
      setvertices(fliptets[0], pa, pb, pc, pd);
      dummyflag = 0;
    }
    if (dummyflag > 0) {
      // We delete 3 hull tets, and create 1 hull tet.
      hullsize -= 2;
    }
  } else {
    setvertices(fliptets[0], pa, pb, pc, pd);
  }

  if (calc_tetprism_vol) {
    if (dummyflag > 0) {
      if (pa == dummypoint) {
        volneg[0] = 0.;
        volneg[1] = tetprismvol(pp, pd, pb, pc);
        volneg[2] = 0.;
        volneg[3] = 0.;
      } else if (pb == dummypoint) {
        volneg[0] = 0.;
        volneg[1] = 0.;
        volneg[2] = tetprismvol(pp, pd, pc, pa);
        volneg[3] = 0.;
      } else if (pc == dummypoint) {
        volneg[0] = tetprismvol(pp, pd, pa, pb);
        volneg[1] = 0.;
        volneg[2] = 0.;
        volneg[3] = 0.;
      } else { // pd == dummypoint
        volneg[0] = 0.;
        volneg[1] = 0.;
        volneg[2] = 0.;
        volneg[3] = tetprismvol(pa, pb, pc, pp);
      }
      volpos[0] = 0.;
    } else {
      volneg[0] = tetprismvol(pp, pd, pa, pb);
      volneg[1] = tetprismvol(pp, pd, pb, pc);
      volneg[2] = tetprismvol(pp, pd, pc, pa);
      volneg[3] = tetprismvol(pa, pb, pc, pp);
      volpos[0] = tetprismvol(pa, pb, pc, pd);
    }
    vol_diff = volpos[0] - volneg[0] - volneg[1] - volneg[2] - volneg[3];
    tetprism_vol_sum  += vol_diff; // Update the total sum.
  }

  // Bond the new tet to adjacent tets.
  for (i = 0; i < 3; i++) {
    esym(fliptets[0], newface); // At faces [b,a,d], [c,b,d], [a,c,d].
    bond(newface, topcastets[i]);
    enextself(fliptets[0]);
  }
  bond(fliptets[0], botcastet);

  if (checksubsegflag) {
    // Bond 6 segments (at edges of [a,b,c,d]) if there there are.
    for (i = 0; i < 3; i++) {
      eprev(topcastets[i], newface); // At edges [d,a],[d,b],[d,c].
      tsspivot1(newface, checkseg);
      if (checkseg.sh != NULL) {
        esym(fliptets[0], newface);
        enextself(newface); // At edges [a,d], [b,d], [c,d].
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      enextself(fliptets[0]);
    }
    for (i = 0; i < 3; i++) {
      tsspivot1(topcastets[i], checkseg); // At edges [a,b],[b,c],[c,a].
      if (checkseg.sh != NULL) {
        tssbond1(fliptets[0], checkseg);
        sstbond1(checkseg, fliptets[0]);
        if (chkencflag & 1) {
          // Skip it if it has already queued.
          if (!smarktest2ed(checkseg)) {
            bface = (badface *) badsubsegs->alloc();
            bface->ss = checkseg;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checkseg); // An alive badface.
          }
        }
      }
      enextself(fliptets[0]);
    }
  }

  if (checksubfaceflag) {
    // Bond 4 subfaces (at faces of [a,b,c,d]) if there are.
    for (i = 0; i < 3; i++) {
      tspivot(topcastets[i], checksh); // At faces [a,b,d],[b,c,d],[c,a,d]
      if (checksh.sh != NULL) {
        esym(fliptets[0], newface); // At faces [b,a,d],[c,b,d],[a,c,d]
        sesymself(checksh);
        tsbond(newface, checksh);
        if (chkencflag & 2) {
          if (!smarktest2ed(checksh)) {
            bface = (badface *) badsubfacs->alloc();
            bface->ss = checksh;
            smarktest2(bface->ss); // Only queue it once.
            bface->forg = sorg(checksh); // An alive badface
          }
        }
      }
      enextself(fliptets[0]);
    }
    tspivot(botcastet, checksh); // At face [b,a,c]
    if (checksh.sh != NULL) {
      sesymself(checksh);
      tsbond(fliptets[0], checksh);
      if (chkencflag & 2) {
        if (!smarktest2ed(checksh)) {
          bface = (badface *) badsubfacs->alloc();
          bface->ss = checksh;
          smarktest2(bface->ss); // Only queue it once.
          bface->forg = sorg(checksh); // An alive badface
        }
      }
    }
  }

  if (checksubfaceflag) {
    if (spivot >= 0) {
      // Perform a 3-to-1 flip in surface triangulation.
      // Depending on the value of 'spivot', the three subfaces are:
      //   - 0: [a,b,p], [b,d,p], [d,a,p]
      //   - 1: [b,c,p], [c,d,p], [d,b,p] 
      //   - 2: [c,a,p], [a,d,p], [d,c,p] 
      //   - 3: [a,b,p], [b,c,p], [c,a,p]
      // Adjust the three subfaces such that their origins are p, i.e., 
      //   - 3: [p,a,b], [p,b,c], [p,c,a]. (Required by the flip31()).
      for (i = 0; i < 3; i++) {
        senext2self(flipshs[i]);
      }
      flip31(flipshs, 0);
      // Delete the three old subfaces.
      for (i = 0; i < 3; i++) {
        shellfacedealloc(subfaces, flipshs[i].sh);
      }
      if (spivot < 3) {
        // // Bond the new subface to the new tet [a,b,c,d].
        tsbond(topcastets[spivot], flipshs[3]);
        fsym(topcastets[spivot], newface);
        sesym(flipshs[3], checksh);
        tsbond(newface, checksh);
      } else {
        // Bound the new subface [a,b,c] to the new tet [a,b,c,d].
        tsbond(fliptets[0], flipshs[3]);
        fsym(fliptets[0], newface);
        sesym(flipshs[3], checksh);
        tsbond(newface, checksh);
      }
    } // if (spivot > 0)
  } // if (checksubfaceflag)

  if (chkencflag & 4) {
    // Put the new tet into check list.
    if (!marktest2ed(fliptets[0])) {
      bface = (badface *) badtetrahedrons->alloc();
      bface->tt = fliptets[0];
      marktest2(bface->tt);
      bface->forg = org(fliptets[0]);
    }
  }

  // Update the point-to-tet map.
  setpoint2tet(pa, encode(fliptets[0]));
  setpoint2tet(pb, encode(fliptets[0]));
  setpoint2tet(pc, encode(fliptets[0]));
  setpoint2tet(pd, encode(fliptets[0]));

  if (flipflag > 0) {
    // Queue faces which may be locally non-Delaunay.
    for (i = 0; i < 3; i++) {
      esym(fliptets[0], newface);
      flippush(flipstack, &newface);
      enextself(fliptets[0]);
    }
    flippush(flipstack, &(fliptets[0]));
  }

  recenttet = fliptets[0];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipnm()    Try to flip an edge through a sequence of elementary flips.   //
//                                                                           //
// 'abtets' is an array of 'n' tets in the star of edge [a,b].These tets are //
// ordered in a counterclockwise cycle with respect to the vector a->b, i.e.,//
// use the right-hand rule.                                                  //
//                                                                           //
// 'level' (>= 0) indicates the current link level. If 'level > 0', we are   //
// flipping a link edge of an edge [a',b'],  and 'abedgepivot' indicates     //
// which link edge, i.e., [c',b'] or [a',c'], is [a,b]  These two parameters //
// allow us to determine the new tets after a 3-to-2 flip, i.e., tets that   //
// do not inside the reduced star of edge [a',b'].                           //
//                                                                           //
// If the flag 'fc->unflip' is set, this routine un-does the flips performed //
// in flipnm([a,b]) so that the mesh is returned to its original state       //
// before doing the flipnm([a,b]) operation.                                 //
//                                                                           //
// The return value is an integer nn, where nn <= n.  If nn is 2, then the   //
// edge is flipped.  The first and the second tets in 'abtets' are new tets. //
// Otherwise, nn > 2, the edge is not flipped, and nn is the number of tets  //
// in the current star of [a,b].                                             //
//                                                                           //
// ASSUMPTIONS:                                                              //
//  - Neither a nor b is 'dummypoint'.                                       //
//  - [a,b] must not be a segment.                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::flipnm(triface* abtets, int n, int level, int abedgepivot,
                       flipconstraints* fc)
{
  triface fliptets[3], spintet, flipedge;
  triface *tmpabtets, *parytet;
  face checksh;
  face checkseg, *paryseg;
  point pa, pb, pc, pd, pe, pf;
  point tmppts[3];
  REAL abovept[3];
  REAL ori, ori1, ori2;
  int reducflag, rejflag;
  int hullflag;
  int reflexlinkedgecount;
  int edgepivot;
  int n1, nn;
  int i, j;

  pa = org(abtets[0]);
  pb = dest(abtets[0]);

  if (b->verbose > 2) {
    printf("      flipnm(%d): (%d, %d) - n(%d), e(%d).\n", level, pointmark(pa),
           pointmark(pb), n, abedgepivot);
  }

  if (n > 3) {
    // Try to reduce the size of the Star(ab) by flipping a face in it. 
    reflexlinkedgecount = 0;

    for (i = 0; i < n; i++) {
      // Let the face of 'abtets[i]' be [a,b,c].
      if (checksubfaceflag) {
        // Do not flip this face if it is a constraining face.
        tspivot(abtets[i], checksh);
        if (checksh.sh != NULL) {
          continue; // Skip a subface.
        }
      }
      // Do not flip this face if it is involved in two Stars.
      if ((elemcounter(abtets[i]) > 1) ||
          (elemcounter(abtets[(i - 1 + n) % n]) > 1)) {
        continue;
      }
      pc = apex(abtets[i]); 
      pd = apex(abtets[(i + 1) % n]);
      pe = apex(abtets[(i - 1 + n) % n]);
      if ((pd == dummypoint) || (pe == dummypoint)) {
        // [a,b,c] is a hull face, it is not flipable.
        continue;
      }
      if (checkinverttetflag) {
        // The mesh contains inverted (or degenerated) elements.
        // Only do check if both elements are valid.
        if (pc != dummypoint) {
          ori = orient3d(pa, pb, pc, pd);
          if (ori < 0) {
            ori = orient3d(pb, pa, pc, pe);
          }
          if (ori >= 0) {
            continue; // An invalid tet.
          }
        } else {
          continue;
        }
      } // if (checkinverttetflag)

      reducflag = 0; // Not reducible.

      hullflag = (pc == dummypoint); // pc may be dummypoint.
      if (hullflag == 0) {
        ori = orient3d(pb, pc, pd, pe); // Is [b,c] locally convex?
        if (ori > 0) {
          ori = orient3d(pc, pa, pd, pe); // Is [c,a] locally convex?
          if (ori > 0) {
            // Test if [a,b] is locally convex OR flat.
            ori = orient3d(pa, pb, pd, pe);
            if (ori > 0) {
              // Found a 2-to-3 flip: [a,b,c] => [e,d]
              reducflag = 1;
            } else if (ori == 0) {
              // [a,b] is flat.
              if (n == 4) {
                // The "flat" tet can be removed immedately by a 3-to-2 flip.
                reducflag = 1;
              }
            }
          }
        }
        if (!reducflag) {
          reflexlinkedgecount++;
        }
      } else {
        // 'c' is dummypoint.
        if (n == 4) {
          // Let the vertex opposite to 'c' is 'f'.
          // A 4-to-4 flip is possible if the two tets [d,e,f,a] and [e,d,f,b]
          //   are valid tets. 
          // Note: When the mesh is not convex, it is possible that [a,b] is
          //   locally non-convex (at hull faces [a,b,e] and [b,a,d]).
          //   In this case, an edge flip [a,b] to [e,d] is still possible.
          pf = apex(abtets[(i + 2) % n]);
          assert(pf != dummypoint);
          ori = orient3d(pd, pe, pf, pa);
          if (ori < 0) {
            ori = orient3d(pe, pd, pf, pb);
            if (ori < 0) {
              // Found a 4-to-4 flip: [a,b] => [e,d]
              reducflag = 1;
              ori = 0; // Signal as a 4-to-4 flip (like a co-palanar case).
            }
          }
        }
      } // if (hullflag)

      if (reducflag) {
        // [a,b,c] could be removed by a 2-to-3 flip.
        rejflag = 0;
        if (fc != NULL) {
          // Check if the flip can be performed.
          rejflag = checkflipeligibility(1, pa, pb, pc, pd, pe, level,
                                         abedgepivot, fc);
        }
        if (!rejflag) {
          // Do flip: [a,b,c] => [e,d].
          fliptets[0] = abtets[i];
          fsym(fliptets[0], fliptets[1]); // abtets[i-1].
          flip23(fliptets, hullflag, 0, 0);

          // Shrink the array 'abtets', maintain the original order.
          //   Two tets 'abtets[i-1] ([a,b,e,c])' and 'abtets[i] ([a,b,c,d])'
          //   are flipped, i.e., they do not in Star(ab) anymore. 
          //   'fliptets[0]' ([e,d,a,b]) is in Star(ab), it is saved in
          //   'abtets[i-1]' (adjust it to be [a,b,e,d]), see below: 
          // 
          //            before                   after
          //     [0] |___________|        [0] |___________| 
          //     ... |___________|        ... |___________|
          //   [i-1] |_[a,b,e,c]_|      [i-1] |_[a,b,e,d]_|
          //     [i] |_[a,b,c,d]_| -->    [i] |_[a,b,d,#]_|
          //   [i+1] |_[a,b,d,#]_|      [i+1] |_[a,b,#,*]_|
          //     ... |___________|        ... |___________|
          //   [n-2] |___________|      [n-2] |___________| 
          //   [n-1] |___________|      [n-1] |_[i]_2-t-3_|
          //
          eprevself(fliptets[0]);
          esymself(fliptets[0]);
          enextself(fliptets[0]); // [a,b,e,d]
          // Increase the counter of this new tet (it is in Star(ab)).
          increaseelemcounter(fliptets[0]); //marktest(fliptets[0]);
          abtets[(i - 1 + n) % n] = fliptets[0];
          for (j = i; j < n - 1; j++) {
            abtets[j] = abtets[j + 1];  // Upshift
          }
          // The last entry 'abtets[n-1]' is empty. It is used in two ways:
          //   (i) it remebers the vertex 'c' (in 'abtets[n-1].tet'), and
          //  (ii) it remebers the position [i] where this flip took place.
          // These informations let us to either undo this flip or recover
          //   the original edge link (for collecting new created tets).
          //abtets[n - 1] = fliptets[1]; // [e,d,b,c] is remebered.
          abtets[n - 1].tet = (tetrahedron *) pc;
          abtets[n - 1].ver = 0; // Clear it.
          // 'abtets[n - 1].ver' is in range [0,11] -- only uses 4 bits.
          // Use the 5th bit in 'abtets[n - 1].ver' to signal a 2-to-3 flip.
          abtets[n - 1].ver |= (1 << 4);
          // The poisition [i] of this flip is saved above the 7th bit.
          abtets[n - 1].ver |= (i << 6);

          if (fc->collectnewtets) {
            // Push the two new tets [e,d,b,c] and [e,d,c,a] into a stack.
            //   Re-use the global array 'cavetetlist'.
            for (j = 1; j < 3; j++) {
              cavetetlist->newindex((void **) &parytet);
              *parytet = fliptets[j]; // fliptets[1], fliptets[2].
            }
          }

          // Star(ab) is reduced. Try to flip the edge [a,b].
          nn = flipnm(abtets, n - 1, level, abedgepivot, fc);

          if (nn > 2) {
            // The edge is not flipped.
            if (fc->unflip || (ori == 0)) {
              // Undo the previous 2-to-3 flip, i.e., do a 3-to-2 flip to 
              //   transform [e,d] => [a,b,c].
              // 'ori == 0' means that the previous flip created a degenrated
              //   tet. It must be removed. 
              // Remeber that 'abtets[i-1]' is [a,b,e,d]. We can use it to
              //   find another two tets [e,d,b,c] and [e,d,c,a].
              fliptets[0] = abtets[(i-1 + (n-1)) % (n-1)]; // [a,b,e,d]
              eprevself(fliptets[0]);
              esymself(fliptets[0]);
              enextself(fliptets[0]); // [e,d,a,b]
              fnext(fliptets[0], fliptets[1]); // [1] is [e,d,b,c]
              fnext(fliptets[1], fliptets[2]); // [2] is [e,d,c,a]
              assert(apex(fliptets[0]) == oppo(fliptets[2])); // SELF_CHECK
              // Restore the two original tets in Star(ab). 
              flip32(fliptets, hullflag, 0, 0);
              // Marktest the two restored tets in Star(ab).
              for (j = 0; j < 2; j++) {
                increaseelemcounter(fliptets[j]); //marktest(fliptets[j]);
              }
              // Expand the array 'abtets', maintain the original order.
              for (j = n - 2; j>= i; j--) {
                abtets[j + 1] = abtets[j];  // Downshift
              }
              // Insert the two new tets 'fliptets[0]' [a,b,c,d] and 
              //  'fliptets[1]' [b,a,c,e] into the (i-1)-th and i-th entries, 
              //  respectively.
              esym(fliptets[1], abtets[(i - 1 + n) % n]); // [a,b,e,c]
              abtets[i] = fliptets[0]; // [a,b,c,d]
              nn++;
              if (fc->collectnewtets) {
                // Pop two (flipped) tets from the stack.
                cavetetlist->objects -= 2;
              }
            } // if (upflip || (ori == 0))
          } // if (nn > 2)

          if (nn == 2) { //if ((nn == 2) || !fullsearch) {
            // The edge has been flipped.
            return nn;
          }
          if (!fc->unflip) {
            // The flips are not reversed. The current Star(ab) can not be
            //   further reduced. Return its size (# of tets).
            return nn; 
          }
          // unflip is set. 
          // Continue the search for flips.
        } else {
          if (b->verbose > 2) {
            printf("      -- Reject a 2-to-3 flip at star face (%d, %d, %d)",
                   pointmark(pa), pointmark(pb), pointmark(pc));
            printf(", link (%d)\n", level);
          }
          if (fc != NULL) {
            fc->rejf23count++;
          }
        } // if (rejflag)
      } // if (reducflag)
    } // i

    // The Star(ab) is not reduced. 
    if (reflexlinkedgecount > 0) {
      // There are reflex edges in the Link(ab).
      if (((b->fliplinklevel < 0) && (level < autofliplinklevel)) || 
          ((b->fliplinklevel >= 0) && (level < b->fliplinklevel))) {
        // Record the largest level.
        if ((level + 1) > maxfliplinklevel) {
          maxfliplinklevel = level + 1;
        }
        if (fc != NULL) {
          // Increase the link level counter.
          if ((level + 1) > fc->maxflippedlinklevelcount) {
            fc->maxflippedlinklevelcount = level + 1;
          }
        }
        // Try to reduce the Star(ab) by flipping a reflex edge in Link(ab).
        for (i = 0; i < n; i++) {
          // Do not flip this face [a,b,c] if there are two Stars involved.
          if ((elemcounter(abtets[i]) > 1) ||
              (elemcounter(abtets[(i - 1 + n) % n]) > 1)) {
            continue;
          }
          pc = apex(abtets[i]);
          if (pc == dummypoint) {
            continue; // [a,b,dummypoint] is a hull edge.
          }
          pd = apex(abtets[(i + 1) % n]);
          pe = apex(abtets[(i - 1 + n) % n]);
          if ((pd == dummypoint) || (pe == dummypoint)) {
            continue; // [a,b,c] is a hull face.
          }
          if (checkinverttetflag) {
            // The mesh contains inverted (or degenerated) elements.
            // Only do check if both elements are valid.
            // assert(pc != dummypoint);
            ori = orient3d(pa, pb, pc, pd);
            if (ori < 0) {
              ori = orient3d(pb, pa, pc, pe);
            }
            if (ori >= 0) {
              continue; // An invalid tet.
            }
          } // if (checkinverttetflag)

          edgepivot = 0; // No edge is selected yet.

          // Test if [b,c] is locally convex or flat.
          ori = orient3d(pb, pc, pd, pe);
          if (ori <= 0) {
            // Select the edge [c,b].
            enext(abtets[i], flipedge); // [b,c,a,d]
            edgepivot = 1;
          }
          if (!edgepivot) {
            // Test if [c,a] is locally convex or flat.
            ori = orient3d(pc, pa, pd, pe);
            if (ori <= 0) {
              // Select the edge [a,c].
              eprev(abtets[i], flipedge); // [c,a,b,d].
              edgepivot = 2;
            }
          }

          if (!edgepivot) continue;

          // An edge is selected.
          if (checksubsegflag) {
            // Do not flip it if it is a segment.
            tsspivot1(flipedge, checkseg);
            if (checkseg.sh != NULL) {
              if (b->verbose > 2) {
                printf("      -- Can't flip a link(%d) segment (%d, %d).\n",
                  level, pointmark(org(flipedge)), pointmark(dest(flipedge)));
              }
              if (fc != NULL) {
                fc->encsegcount++;
                if (fc->collectencsegflag) {
                  if (!sinfected(checkseg)) {
                    // Queue this segment in list.
                    sinfect(checkseg);                
                    caveencseglist->newindex((void **) &paryseg);
                    *paryseg = checkseg;
                  }
                }
              }
              continue;
            }
          }

          // Try to flip the selected edge ([c,b] or [a,c]).
          esymself(flipedge); 
          // Count the number of tets at the edge.
          n1 = 0;
          j = 0; // Sum of the star counters.
          spintet = flipedge;
          while (1) {
            n1++;
            j += (elemcounter(spintet)); //if (marktested(spintet)) j++;
            fnextself(spintet);
            if (spintet.tet == flipedge.tet) break;
          }
          assert(n1 >= 3);
          if (j > 2) {
            // The Star(flipedge) overlaps other Stars.
            continue; // Do not flip this edge.
          }
          // Only two tets can be marktested.
          assert(j == 2); 

          flipstarcount++;
          // Record the maximum star size.
          if (n1 > maxflipstarsize) {
            maxflipstarsize = n1;
          }
          if ((b->flipstarsize > 0) && (n1 > b->flipstarsize)) {
            // The star size exceeds the given limit (-LL__).
            skpflipstarcount++;
            continue; // Do not flip it.
          }

          // Allocate spaces for Star(flipedge).
          tmpabtets = new triface[n1];
          // Form the Star(flipedge).
          j = 0;
          spintet = flipedge;
          while (1) {
            tmpabtets[j] = spintet;
            // Increase the star counter of this tet.
            increaseelemcounter(tmpabtets[j]); 
            j++;
            fnextself(spintet);
            if (spintet.tet == flipedge.tet) break;
          }
          // SELF_CHECK BEGIN
          // These two tets are inside both of the Stars.
          assert(elemcounter(tmpabtets[0]) == 2);
          assert(elemcounter(tmpabtets[1]) == 2);
          // Marktest the tets in Star(flipedge) but not in Star(ab).
          for (j = 2; j < n1; j++) {
            assert(elemcounter(tmpabtets[j]) == 1);
            //marktest(tmpabtets[j]);
          }

          // Try to flip the selected edge away.
          nn = flipnm(tmpabtets, n1, level + 1, edgepivot, fc);

          if (nn == 2) {
            // The edge is flipped. Star(ab) is reduced.
            // Shrink the array 'abtets', maintain the original order.
            if (edgepivot == 1) {
              // 'tmpabtets[0]' is [d,a,e,b] => contains [a,b].
              spintet = tmpabtets[0]; // [d,a,e,b]
              enextself(spintet);
              esymself(spintet);
              enextself(spintet); // [a,b,e,d]
            } else {
              // 'tmpabtets[1]' is [b,d,e,a] => contains [a,b].
              spintet = tmpabtets[1]; // [b,d,e,a]
              eprevself(spintet);
              esymself(spintet);
              eprevself(spintet); // [a,b,e,d]
            } // edgepivot == 2
            //assert(!marktested(spintet)); // It's a new tet.
            assert(elemcounter(spintet) == 0);
            //marktest(spintet); // It is in Star(ab).
            increaseelemcounter(spintet);
            // Put the new tet at [i-1]-th entry.
            abtets[(i - 1 + n) % n] = spintet;
            for (j = i; j < n - 1; j++) {
              abtets[j] = abtets[j + 1];  // Upshift
            }
            // Remember the flips in the last entry of the array 'abtets'.
            // They can be used to recover the flipped edge.
            abtets[n - 1].tet = (tetrahedron *) tmpabtets; // The star(fedge).
            abtets[n - 1].ver = 0; // Clear it.
            // Use the 1st and 2nd bit to save 'edgepivot' (1 or 2).
            abtets[n - 1].ver |= edgepivot;
            // Use the 6th bit to signal this n1-to-m1 flip.
            abtets[n - 1].ver |= (1 << 5); 
            // The poisition [i] of this flip is saved from 7th to 19th bit.
            abtets[n - 1].ver |= (i << 6);
            // The size of the star 'n1' is saved from 20th bit.
            abtets[n - 1].ver |= (n1 << 19);

            // Remember the flipped link vertex 'c'. It can be used to recover
            //   the original edge link of [a,b], and to collect new tets.
            tmpabtets[0].tet = (tetrahedron *) pc;
            tmpabtets[0].ver = (1 << 5); // Flag it as a vertex handle.

            // Continue to flip the edge [a,b].
            nn = flipnm(abtets, n - 1, level, abedgepivot, fc);

            if (nn > 2) {
              // The edge is not flipped.
              if (fc->unflip) {
                // Recover the flipped edge ([c,b] or [a,c]).
                assert(nn == (n - 1));
                // The sequence of flips are saved in 'tmpabtets'. 
                // abtets[(i-1) % (n-1)] is [a,b,e,d], i.e., the tet created by
                //   the flipping of edge [c,b] or [a,c].It must still exist in
                //   Star(ab). It is the start tet to recover the flipped edge.
                if (edgepivot == 1) { 
                  // The flip edge is [c,b].
                  tmpabtets[0] = abtets[((i-1)+(n-1))%(n-1)]; // [a,b,e,d]
                  eprevself(tmpabtets[0]);
                  esymself(tmpabtets[0]);
                  eprevself(tmpabtets[0]); // [d,a,e,b]
                  fsym(tmpabtets[0], tmpabtets[1]); // [a,d,e,c]
                } else {
                  // The flip edge is [a,c].
                  tmpabtets[1] = abtets[((i-1)+(n-1))%(n-1)]; // [a,b,e,d]
                  enextself(tmpabtets[1]);
                  esymself(tmpabtets[1]);
                  enextself(tmpabtets[1]); // [b,d,e,a]
                  fsym(tmpabtets[1], tmpabtets[0]); // [d,b,e,c]
                } // if (edgepivot == 2)

                // Recover the flipped edge ([c,b] or [a,c]).
                flipnm_post(tmpabtets, n1, 2, edgepivot, fc);

                // Insert the two recovered tets into Star(ab).
                for (j = n - 2; j >= i; j--) {
                  abtets[j + 1] = abtets[j];  // Downshift
                }
                if (edgepivot == 1) {
                  // tmpabtets[0] is [c,b,d,a] ==> contains [a,b]
                  // tmpabtets[1] is [c,b,a,e] ==> contains [a,b]
                  // tmpabtets[2] is [c,b,e,d]
                  fliptets[0] = tmpabtets[1];
                  enextself(fliptets[0]);
                  esymself(fliptets[0]); // [a,b,e,c]
                  fliptets[1] = tmpabtets[0];
                  esymself(fliptets[1]);
                  eprevself(fliptets[1]); // [a,b,c,d]
                } else {
                  // tmpabtets[0] is [a,c,d,b] ==> contains [a,b]
                  // tmpabtets[1] is [a,c,b,e] ==> contains [a,b]
                  // tmpabtets[2] is [a,c,e,d]
                  fliptets[0] = tmpabtets[1];
                  eprevself(fliptets[0]);
                  esymself(fliptets[0]); // [a,b,e,c]
                  fliptets[1] = tmpabtets[0];
                  esymself(fliptets[1]);
                  enextself(fliptets[1]); // [a,b,c,d]
                } // edgepivot == 2
                for (j = 0; j < 2; j++) {
                  assert(elemcounter(fliptets[j]) == 0); // SELF_CHECK
                  increaseelemcounter(fliptets[j]);
                }
                // Insert the two recovered tets into Star(ab).
                abtets[(i - 1 + n) % n] = fliptets[0];
                abtets[i] = fliptets[1];
                nn++;
                // Release the allocated spaces.
                delete [] tmpabtets;
              } // if (unflip)
            } // if (nn > 2)

            if (nn == 2) { //if ((nn == 2) || !fullsearch) {
              // The edge has been flipped.
              return nn;
            }
            if (!fc->unflip) {
              // The flips are not reversed. The current Star(ab) can not be
              //   further reduced. Return its size (# of tets).
              return nn; 
            }
            // unflip is set. 
            // Continue the search for flips.
          } else {
            // The seclected edge is not flipped.
            if (fc->unflip) {
              // The memory should already be freed.
              assert(nn == n1);
            } else {
              // Release the memory used in this attempted flip.
              flipnm_post(tmpabtets, n1, nn, edgepivot, fc);
            }
            // Decrease the star counters of tets in Star(flipedge).
            for (j = 0; j < nn; j++) {
              assert(elemcounter(tmpabtets[j]) > 0); // SELF_CHECK
              decreaseelemcounter(tmpabtets[j]);
            }
            // Release the allocated spaces.
            delete [] tmpabtets;
          }
        } // i
      } else {
        if (b->verbose > 2) {
          printf("      -- Maximal link level (%d) reached at edge (%d, %d).\n",
                 level, pointmark(org(abtets[0])), pointmark(dest(abtets[0])));
        }
        if (fc != NULL) {
          fc->misfliplinklevelcount++;
        }
      } // if (level...)
    } // if (reflexlinkedgecount > 0)
  } else {
    // Check if a 3-to-2 flip is possible.
    pc = apex(abtets[0]);
    pd = apex(abtets[1]);
    pe = apex(abtets[2]);

    // Check if one of them is dummypoint. If so, we rearrange the vertices
    //   c, d, and e into p0, p1, and p2, such that p2 is the dummypoint.
    hullflag = 0;
    if (pc == dummypoint) {
      hullflag = 1;
      tmppts[0] = pd;
      tmppts[1] = pe;
      tmppts[2] = pc;
    } else if (pd == dummypoint) {
      hullflag = 1;
      tmppts[0] = pe;
      tmppts[1] = pc;
      tmppts[2] = pd;
    } else if (pe == dummypoint) {
      hullflag = 1;
      tmppts[0] = pc;
      tmppts[1] = pd;
      tmppts[2] = pe;
    } else {
      tmppts[0] = pc;
      tmppts[1] = pd;
      tmppts[2] = pe;
    }

    reducflag = 0;
    rejflag = 0;

    if (checkinverttetflag) {
      // Only do flip if no tet is inverted (or degenerated).
      if (hullflag == 0) {
        ori = orient3d(pa, pb, pc, pd);
        if (ori < 0) {
          ori = orient3d(pa, pb, pd, pe);
          if (ori < 0) {
            ori = orient3d(pa, pb, pe, pc);
          }
        }
      } else {
        ori = orient3d(pa, pb, tmppts[0], tmppts[1]);
      }
      if (ori >= 0) {
        if (b->verbose > 2) {
          printf("      -- Hit a non-valid tet (%d, %d) - (%d, %d, %d)",
                 pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd),
                 pointmark(pe));
          printf(" at link(%d)\n", level);
        }
        return 3;
      }
    } // if (checkinverttetflag)

    if (hullflag == 0) {
      // Make sure that no inverted tet will be created, i.e. the new tets
      //   [d,c,e,a] and [c,d,e,b] must be valid tets. 
      ori = orient3d(pd, pc, pe, pa);
      if (ori < 0) {
        ori = orient3d(pc, pd, pe, pb);
        if (ori < 0) {
          reducflag = 1;
        }
      } else {
        if (b->verbose > 2) {
          printf("      -- Hit a chrismastree (%d, %d) - (%d, %d, %d)",
                 pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd),
                 pointmark(pe));
          printf(" at link(%d)\n", level);
        }
        if (fc != NULL) {
          fc->chrismastreecount++;
        }
      }
    } else {
      // [a,b] is a hull edge. Moreover, the tet [a,b,p0,p1] is a hull tet
      //   ([a,b,p0] and [a,b,p1] are two hull faces). 
      //   This can happen when it is in the middle of a 4-to-4 flip.
      //   Note that [a,b] may even be a non-convex hull edge.
      if (!nonconvex) {
        // [a,b], [a,b,p0] and [a,b,p1] are on the convex hull.
        ori = orient3d(pa, pb, tmppts[0], tmppts[1]);
        if (ori == 0) {
          // They four vertices are coplanar. A 2-to-2 flip is possible if
          //   [a,b] and [p0,p1] are intersecting each other. 
          // NOTE: The following test is not robust, should be replaced in 
          //   the future. 2011-12-01.
          calculateabovepoint4(pa, pb, tmppts[0], tmppts[1]);
          for (j = 0; j < 3; j++) {
            abovept[j] = dummypoint[j];
          }
          // Make sure that no inverted face will be created, i.e., [p1,p0,
          //   abvpt,pa] and [p0,p1,abvpt,pb] must be valid tets.
          ori1 = orient3d(tmppts[0], tmppts[1], abovept, pa);
          ori2 = orient3d(tmppts[0], tmppts[1], abovept, pb);
          if (ori1 * ori2 < 0) {
            reducflag = 1; // Flipable.
          }
          if (!reducflag) {
            if (b->verbose > 2) {
              printf("      -- Hit a degenerate chrismastree (%d, %d)",
                     pointmark(pa), pointmark(pb));
              printf(" - (%d, %d, -1) at link(%d)\n", 
                     pointmark(tmppts[0]), pointmark(tmppts[1]), level);
            }
            if (fc != NULL) {
              fc->chrismastreecount++;
            }
          }
        } else {
          if (b->verbose > 2) {
            printf("      -- Hit a convex hull edge (%d, %d) at link(%d).\n",
                   pointmark(pa), pointmark(pb), level);
          }
          if (fc != NULL) {
            fc->convexhulledgecount++;
          }
        }
      } else { // if (nonconvex)
        // [a,b,p0] and [a,b,p1] must be two subfaces.
        // Since [a,b] is not a segment. A 3-to-2 flip (including a 2-to-2
        //   flip) is possible.
        // Here we only do flip if there are exactly three tets containing
        //   the edge [p0,p1]. In this case, the other two tets at [p0,p1]
        //   (not [a,b,p0,p1]) must be valid. Since they already exist.
        for (j = 0; j < 3; j++) {
          if (apex(abtets[j]) == dummypoint) {
            flipedge = abtets[(j + 1) % 3]; // [a,b,p0,p1].
            break;
          }
        }
        // assert(j < 3);
        eprevself(flipedge);
        esymself(flipedge);
        enextself(flipedge); // [p0,p1,a,b].
        assert(apex(flipedge) == pa);
        spintet = flipedge;
        j = 0;
        while (1) {
          j++;
          fnextself(spintet);
          if (spintet.tet == flipedge.tet) break;
        }
        if (j == 3) {
          reducflag = 1;
        } else {
          if (b->verbose > 2) {
            printf("      -- Hit a hull edge (%d, %d) at link(%d).\n",
                   pointmark(pa), pointmark(pb), level);
          }
          //if (fc != NULL) {
          //  fc->convexhulledgecount++;
          //}
        }
      }
    } // if (hullflag == 1)

    if (reducflag) {
      // A 3-to-2 flip is possible.
      if (checksubfaceflag) {
        // This edge (must not be a segment) can be flipped ONLY IF it belongs
        //   to either 0 or 2 subfaces.  In the latter case, a 2-to-2 flip in 
        //   the surface mesh will be automatically performed within the 
        //   3-to-2 flip.
        nn = 0;
        for (j = 0; j < 3; j++) {
          tspivot(abtets[j], checksh);
          if (checksh.sh != NULL) {
            nn++; // Found a subface.
          }
        }
        assert(nn < 3);
        if (nn == 1) {
          // Found only 1 subface containing this edge. This can happen in 
          //   the boundary recovery phase. The neighbor subface is not yet 
          //   recovered. This edge should not be flipped at this moment.
          rejflag = 1; 
        }
      }
      if (!rejflag && (fc != NULL)) {
        // Here we must exchange 'a' and 'b'. Since in the check... function,
        //   we assume the following point sequence, 'a,b,c,d,e', where
        //   the face [a,b,c] will be flipped and the edge [e,d] will be
        //   created. The two new tets are [a,b,c,d] and [b,a,c,e]. 
        rejflag = checkflipeligibility(2, tmppts[0], tmppts[1], tmppts[2], 
                                       pb, pa, level, abedgepivot, fc);
      }
      if (!rejflag) {
        // Do flip: [a,b] => [c,d,e]
        flip32(abtets, hullflag, 0, 0);
        sucflipstarcount++;
        if (fc->remove_ndelaunay_edge) {
          if (level == 0) {
            // It is the desired removing edge.
            if (tetprism_vol_sum >= fc->bak_tetprism_vol) {
              if (b->verbose > 2) {
                printf("      -- Reject to flip (%d, %d) at link(%d)\n",
                       pointmark(pa), pointmark(pb), level);
                printf("         due to an increased volume (%.17g).\n",
                       tetprism_vol_sum - fc->bak_tetprism_vol);
              }
              // flip back: [c,d,e] => [a,b].
              flip23(abtets, hullflag, 0, 0);
              // Increase the element counter -- They are in cavity.
              for (j = 0; j < 3; j++) {
                increaseelemcounter(abtets[j]); 
              }
              return 3;
            }
          } // if (level == 0)
        }
        if (fc->collectnewtets) {
          // Collect new tets.
          if (level == 0) {
            // Push the two new tets into stack.
            for (j = 0; j < 2; j++) {
              cavetetlist->newindex((void **) &parytet);
              *parytet = abtets[j];
            }
          } else {
            // Only one of the new tets is collected. The other one is inside
            //   the reduced edge star. 'abedgepivot' is either '1' or '2'.
            cavetetlist->newindex((void **) &parytet);
            if (abedgepivot == 1) { // [c,b]
              *parytet = abtets[1];
            } else { 
              assert(abedgepivot == 2); // [a,c]
              *parytet = abtets[0];
            }
          }
        } // if (fc->collectnewtets)
        return 2;
      } else {
        if (b->verbose > 2) {
          printf("      -- Reject a 3-to-2 flip (%d, %d) at link(%d).\n",
                 pointmark(pa), pointmark(pb), level);
        }
        if (fc != NULL) {
          fc->rejf32count++;
        }
      } // if (rejflag)
    } // if (reducflag)
  } // if (n == 3)

  // The current (reduced) Star size.
  return n;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipnm_post()    Post process a n-to-m flip.                              //
//                                                                           //
// IMPORTANT: This routine only works when there is no other flip operation  //
// is done after flipnm([a,b]) which attempts to remove an edge [a,b].       //
//                                                                           //
// 'abtets' is an array of 'n' (>= 3) tets which are in the original star of //
// [a,b] before flipnm([a,b]).  'nn' (< n) is the value returned by flipnm.  //
// If 'nn == 2', the edge [a,b] has been flipped. 'abtets[0]' and 'abtets[1]'//
// are [c,d,e,b] and [d,c,e,a], i.e., a 2-to-3 flip can recover the edge [a, //
// b] and its initial Star([a,b]).  If 'nn >= 3' edge [a,b] still exists in  //
// current mesh and 'nn' is the current number of tets in Star([a,b]).       //
//                                                                           //
// Each 'abtets[i]', where nn <= i < n, saves either a 2-to-3 flip or a      //
// flipnm([p1,p2]) operation ([p1,p2] != [a,b]) which created the tet        //
// 'abtets[t-1]', where '0 <= t <= i'.  These information can be used to     //
// undo the flips performed in flipnm([a,b]) or to collect new tets created  //
// by the flipnm([a,b]) operation.                                           //
//                                                                           //
// Default, this routine only walks through the flips and frees the spaces   //
// allocated during the flipnm([a,b]) operation.                             //
//                                                                           //
// If the flag 'fc->unflip' is set, this routine un-does the flips performed //
// in flipnm([a,b]) so that the mesh is returned to its original state       //
// before doing the flipnm([a,b]) operation.                                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::flipnm_post(triface* abtets, int n, int nn, int abedgepivot,
                            flipconstraints* fc)
{
  triface fliptets[3], flipface;
  triface *tmpabtets;
  int fliptype;
  int edgepivot;
  int t, n1;
  int i, j;


  if (nn == 2) {
    // The edge [a,b] has been flipped.
    // 'abtets[0]' is [c,d,e,b] or [#,#,#,b].
    // 'abtets[1]' is [d,c,e,a] or [#,#,#,a].
    if (fc->unflip) {
      // Do a 2-to-3 flip to recover the edge [a,b]. There may be hull tets.
      flip23(abtets, 1, 0, 0);
      if (fc->collectnewtets) {
        // Pop up new (flipped) tets from the stack.
        if (abedgepivot == 0) {
          // Two new tets were collected.
          cavetetlist->objects -= 2;
        } else {
          // Only one of the two new tets was collected.
          cavetetlist->objects -= 1;
        }
      }
    } 
    // The initial size of Star(ab) is 3.
    nn++;
  } 

  // Walk through the performed flips.
  for (i = nn; i < n; i++) {
    // At the beginning of each step 'i', the size of the Star([a,b]) is 'i'.
    // At the end of this step, the size of the Star([a,b]) is 'i+1'.
    // The sizes of the Link([a,b]) are the same.
    fliptype = ((abtets[i].ver >> 4) & 3); // 0, 1, or 2.
    if (fliptype == 1) {
      // It was a 2-to-3 flip: [a,b,c]->[e,d].
      t = (abtets[i].ver >> 6);
      assert(t <= i);
      if (fc->unflip) {
        if (b->verbose > 2) {
          printf("      Recover a 2-to-3 flip at f[%d].\n", t);
        }
        // 'abtets[(t-1)%i]' is the tet [a,b,e,d] in current Star(ab), i.e.,
        //   it is created by a 2-to-3 flip [a,b,c] => [e,d].
        fliptets[0] = abtets[((t - 1) + i) % i]; // [a,b,e,d]
        eprevself(fliptets[0]);
        esymself(fliptets[0]);
        enextself(fliptets[0]); // [e,d,a,b]
        fnext(fliptets[0], fliptets[1]); // [e,d,b,c]
        fnext(fliptets[1], fliptets[2]); // [e,d,c,a]
        // Do a 3-to-2 flip: [e,d] => [a,b,c].
        // NOTE: hull tets may be invloved.
        flip32(fliptets, 1, 0, 0);
        // Expand the array 'abtets', maintain the original order.
        // The new array length is (i+1).
        for (j = i - 1; j >= t; j--) {
          abtets[j + 1] = abtets[j];  // Downshift
        }
        // The tet abtets[(t-1)%i] is deleted. Insert the two new tets 
        //   'fliptets[0]' [a,b,c,d] and 'fliptets[1]' [b,a,c,e] into
        //   the (t-1)-th and t-th entries, respectively.
        esym(fliptets[1], abtets[((t-1) + (i+1)) % (i+1)]); // [a,b,e,c]
        abtets[t] = fliptets[0]; // [a,b,c,d]
        if (fc->collectnewtets) {
          // Pop up two (flipped) tets from the stack.
          cavetetlist->objects -= 2;
        }
      } 
    } else if (fliptype == 2) {
      tmpabtets = (triface *) (abtets[i].tet);
      n1 = ((abtets[i].ver >> 19) & 8191); // \sum_{i=0^12}{2^i} = 8191
      edgepivot = (abtets[i].ver & 3); 
      t = ((abtets[i].ver >> 6) & 8191);
      assert(t <= i);
      if (fc->unflip) {        
        if (b->verbose > 2) {
          printf("      Recover a %d-to-m flip at e[%d] of f[%d].\n", n1, 
                 edgepivot, t);
        }
        // Recover the flipped edge ([c,b] or [a,c]).
        // abtets[(t - 1 + i) % i] is [a,b,e,d], i.e., the tet created by
        //   the flipping of edge [c,b] or [a,c]. It must still exist in
        //   Star(ab). Use it to recover the flipped edge.
        if (edgepivot == 1) { 
          // The flip edge is [c,b].
          tmpabtets[0] = abtets[(t - 1 + i) % i]; // [a,b,e,d]
          eprevself(tmpabtets[0]);
          esymself(tmpabtets[0]);
          eprevself(tmpabtets[0]); // [d,a,e,b]
          fsym(tmpabtets[0], tmpabtets[1]); // [a,d,e,c]
        } else {
          // The flip edge is [a,c].
          tmpabtets[1] = abtets[(t - 1 + i) % i]; // [a,b,e,d]
          enextself(tmpabtets[1]);
          esymself(tmpabtets[1]);
          enextself(tmpabtets[1]); // [b,d,e,a]
          fsym(tmpabtets[1], tmpabtets[0]); // [d,b,e,c]
        } // if (edgepivot == 2)

        // Do a n1-to-m1 flip to recover the flipped edge ([c,b] or [a,c]).
        flipnm_post(tmpabtets, n1, 2, edgepivot, fc);

        // Insert the two recovered tets into the original Star(ab).
        for (j = i - 1; j >= t; j--) {
          abtets[j + 1] = abtets[j];  // Downshift
        }
        if (edgepivot == 1) {
          // tmpabtets[0] is [c,b,d,a] ==> contains [a,b]
          // tmpabtets[1] is [c,b,a,e] ==> contains [a,b]
          // tmpabtets[2] is [c,b,e,d]
          fliptets[0] = tmpabtets[1];
          enextself(fliptets[0]);
          esymself(fliptets[0]); // [a,b,e,c]
          fliptets[1] = tmpabtets[0];
          esymself(fliptets[1]);
          eprevself(fliptets[1]); // [a,b,c,d]
        } else {
          // tmpabtets[0] is [a,c,d,b] ==> contains [a,b]
          // tmpabtets[1] is [a,c,b,e] ==> contains [a,b]
          // tmpabtets[2] is [a,c,e,d]
          fliptets[0] = tmpabtets[1];
          eprevself(fliptets[0]);
          esymself(fliptets[0]); // [a,b,e,c]
          fliptets[1] = tmpabtets[0];
          esymself(fliptets[1]);
          enextself(fliptets[1]); // [a,b,c,d]
        } // edgepivot == 2
        // Insert the two recovered tets into Star(ab).
        abtets[((t-1) + (i+1)) % (i+1)] = fliptets[0];
        abtets[t] = fliptets[1];
      } 
      else {
        // Only free the spaces.
        flipnm_post(tmpabtets, n1, 2, edgepivot, fc);
      } // if (!unflip)
      if (b->verbose > 2) {
        printf("      Release %d spaces at f[%d].\n", n1, i);
      }
      delete [] tmpabtets;
    } else {
      assert(fliptype == 0); // Not a saved flip.
      assert(0); // Should be not possible.
    }
  } // i

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// lawsonflip3d()    A three-dimensional Lawson's flip algorithm.            //
//                                                                           //
// The basic idea of Lawson's algorithm is to flip every face of the triang- //
// ulation which is not locally Delaunay until no such face exists, then the //
// triangulation is a DT. However, in 3D, it is common that a face which is  //
// not locally Delaunay and is not flippable. Hence, Laowson's algorithm may //
// get stuck. It is still an open problem, whether there exists a flip algo- //
// rithm which has a guarantee to create a DT in 3D.                         //
//                                                                           //
// If only one vertex is added into a DT, then Lawson's flip algorithm is    //
// guaranteed to transform it into a new DT [Joe'91]. Moreover, an arbitrary //
// order of flips is sufficient [Edelsbrunner & Shah'96].                    //
//                                                                           //
// In practice, it is desired to remove not locally Delaunay faces by flips  //
// as many as possible. For this purpose, a second queue is used to store    //
// the not locally Delaunay faces which are not flippable, and try them at a //
// later time.                                                               //
//                                                                           //
// If 'newpt' (p) is not NULL, it is a new vertex just inserted into the     //
// tetrahedralization T.                                                     //
//                                                                           //
// 'flipflag' indicates the property of the tetrahedralization 'T' which     //
// does not include 'p' yet.                                                 //
//                                                                           //
// If 'peelsliverflag' is set, the purpose of calling Lawson's flip is to    //
// remove "hull slivers". This flag only works with a non-convex mesh, i.e., //
// the mesh must contains boundaries (segments and subfaces).                //
//                                                                           //
// 'chkencflag' indicates whether segments, subfaces, and tets should be     //
//  checked (for encroaching and quality) after flips.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

long tetgenmesh::lawsonflip3d(point newpt, int flipflag, int peelsliverflag, 
                              int chkencflag, int flipedgeflag)
{
  badface *popface, *bface;
  triface fliptets[5], baktets[2];
  triface fliptet, neightet, *parytet;
  face checksh, *parysh;
  face checkseg, *paryseg;
  point *ppt, pd, pe, pf;
  long flipcount;
  REAL sign, ori;
  int convflag;
  int n, i;

  // For removing hull slivers.
  face neighsh; 
  point p1, p2;
  point pa, pb, pc, rempt;
  REAL ang; 
  long tetpeelcount; 
  int remflag;

  flipconstraints fc;

  if (b->verbose > 2) {
    printf("      Lawson flip %ld faces.\n", flippool->items);
  }

  flipcount = flip23count + flip32count + flip44count;
  tetpeelcount = opt_sliver_peels;

  if (flipedgeflag) {
    fc.remove_ndelaunay_edge = 1;
    fc.unflip = 1; // Unflip if the edge is not flipped.
    fc.collectnewtets = 1;
    assert(cavetetlist->objects == 0l);
    assert(calc_tetprism_vol == 1); // Swith on.
  } else {
    assert(unflipqueue->objects == 0); // The second queue must be empty.
  }

  while (1) {

    while (flipstack != (badface *) NULL) {

      // Pop a face from the stack.
      popface = flipstack;
      flipstack = flipstack->nextitem; // The next top item in stack.
      fliptet = popface->tt;
      flippool->dealloc((void *) popface);

      // Skip it if it is a dead tet (destroyed by previous flips).
      if (isdeadtet(fliptet)) continue;
      // Skip it if it is not the same tet as we saved.
      if (!facemarked(fliptet)) continue;

      unmarkface(fliptet);

      if (ishulltet(fliptet)) {
        // It is a hull tet.
        if (((flipflag == 4) || peelsliverflag) && !b->convex) {
          fliptet.ver = epivot[fliptet.ver & 3];
          if (oppo(fliptet) == dummypoint) {
            // It's a hull face (oppo(fliptet) == dummypoint).
            // Check if there exists a "hull sliver".
            fsymself(fliptet);
            tspivot(fliptet, checksh);
            assert(checksh.sh != NULL);
            for (i = 0; i < 3; i++) {
              sspivot(checksh, checkseg);
              if (checkseg.sh == NULL) {
                spivot(checksh, neighsh);
                assert(neighsh.sh != NULL);
                if (sorg(checksh) != sdest(neighsh)) {
                  sesymself(neighsh);
                }
                stpivot(neighsh, neightet);
                if (neightet.tet == fliptet.tet) {
                  // Found a hull sliver 'neightet' [d,e,a,b], where [d,e,a] 
                  //   and [e,d,b] are two hull faces. Normally, a 3-to-2 flip
                  //   (including a 2-to-2 flip on hull subfaces) can remove 
                  //   this hull sliver.
                  // A special case is the existence of a hull tet [b,a,d,-1]
                  //   or [a,b,e,-1]. It was creared by a previous hull tet
                  //   removal. Moreover, 'd' or 'e' might be Steiner points
                  //   on segments [a,b]. In this case, eithe [a,d],[b,d] or 
                  //   [a,e],[b,e] are subsegments. If so, a 4-to-1 flip
                  //   (including a 3-to-1, and maybe a 2-to-1 flip) should be
                  //   applied to remove an exterior vertex.
                  // First check if the face [b,a,d] is a hull face.
                  eprev(neightet, fliptets[0]);
                  esymself(fliptets[0]);  // [d,a,b,e]
                  enextself(fliptets[0]); // [a,b,d,e]
                  fsymself(fliptets[0]);  // [b,a,d,#]
                  if (oppo(fliptets[0]) != dummypoint) {
                    // Second check if the face [a,b,e] is a hull face.
                    enext(neightet, fliptets[0]);
                    esymself(fliptets[0]);  // [a,e,b,d]
                    eprevself(fliptets[0]); // [b,a,e,d]
                    fsymself(fliptets[0]);  // [b,a,e,#]
                  }

                  if (oppo(fliptets[0]) != dummypoint) {
                    // Make sure we do not create an "inverted triangle" in the
                    //   boundary, i.e., in exactly planar case, d and e must
                    //   lie in the different sides of the edge [a,b]. 
                    // If the dihedral angle formed by [a,b,e] and [a,b,d] is
                    //   larger than 90 degree, we can remove [a,b,e,d].  
                    fliptets[0] = neightet; // [e,d,a,b]
                    eprevself(fliptets[0]);
                    esymself(fliptets[0]);
                    enextself(fliptets[0]); // [a,b,e,d].
                    pa = org(fliptets[0]);
                    pb = dest(fliptets[0]);
                    p1 = apex(fliptets[0]); // pe
                    p2 = oppo(fliptets[0]); // pd
                    ang = facedihedral(pa, pb, p1, p2);
                    ang *= 2.0;
                    if (ang > PI) {
                      if (b->verbose > 2) {
                        printf("      Remove a hull sliver (%d, %d, %d, %d).\n",
                          pointmark(org(fliptet)), pointmark(dest(fliptet)),
                          pointmark(apex(fliptet)), pointmark(oppo(fliptet)));
                      }
                      // Remove the ill tet from bounday.
                      fliptets[0] = neightet;          // [e,d,a,b]
                      fnext(fliptets[0], fliptets[1]); // [e,d,b,c]
                      fnext(fliptets[1], fliptets[2]); // [e,d,c,a]
                      // FOR DEBUG
                      fnext(fliptets[2], fliptets[3]);
                      assert(fliptets[3].tet == neightet.tet);
                      assert(oppo(fliptets[1]) == dummypoint);
                      // Do a 3-to-2 flip to remove the ill tet. Two hull tets
                      // are removed toether. Two hull subfaces are flipped.
                      flip32(fliptets, 1, flipflag, 0);
                      // Update counters.
                      flip32count--;
                      flip22count--;
                      opt_sliver_peels++;
                    }
                  } else {
                    // There exists a thrid hull tet at vertex.
                    rempt = apex(fliptets[0]);
                    if (pointmark(rempt) > 
                      (in->numberofpoints - (in->firstnumber ? 0 : 1))) {
                      if (pointtype(rempt) == FREESEGVERTEX) {
                        st_segref_count--;
                      } else if (pointtype(rempt) == FREEFACETVERTEX) {
                        st_facref_count--;
                      } else {
                        assert(0); // Impossible.
                      }
                      if (b->verbose > 2) {
                        printf("      Remove an exterior Steiner vertex %d.\n", 
                               pointmark(rempt));
                      }
                      if (removevertexbyflips(rempt)) {
                        // exsteinercount++;
                      } else {
                        assert(0); // Not possible.
                      }
                    } else {
                      //if (b->verbose > 2) {
                      //  printf("      Remove an exterior input vertex %d.\n", 
                      //         pointmark(rempt));
                      //}
                      // Comment: We do not remove an input point.
                    }
                  }
                  break;
                }
              } // if (checkseg.sh == NULL)
              senextself(checksh);
            } // i
          } else {
            // It's a hull edge.
            assert(apex(fliptet) == dummypoint);
            if (!peelsliverflag) { 
              // The hull edge may be not locally Delaunay. Put interior
              //   faces at this edge into 'flipstack' for flipping.
              neightet = fliptet;  // [a,b,c,d] ('c' is dummypoint).
              fnextself(neightet); // [a,b,d,#1] ([a,b,d] is a hull face).
              while (1) {
                fnextself(neightet); // [a,b,#1,#2]
                if (oppo(neightet) != dummypoint) {
                  // It is an interior face.
                  flippush(flipstack, &neightet);
                } else {
                  // We assume the interior of the domain is connected.
                  // Hence we can hit hull faces only twice.
                  break;
                }
              } // while (1)
            } // if (!peelsliverflag)
          }
        } // if ((flipflag == 4) || peelsliverflag)

        // Do not flip a hull face/edge UNLESS it is in the process of
        //   incrementally creating a DT in which the convex hull may be
        //   enlarged by the flips (when p lies outside of it).
        if (flipflag != 1) {
          continue;
        }
      } // if (ishulltet(fliptet))

      if (peelsliverflag) {
        continue; // Only check hull tets.
      }

      // Let 'fliptet' be [a,b,c,d], the face [a,b,c] is the flip face.
      // Get its opposite tet [b,a,c,e].
      fsym(fliptet, neightet);

      if (ishulltet(neightet)) {
        // It is a hull tet.
        if (flipflag == 1) {
          // Check if the new point is visible by the hull face.
          ppt = (point *) neightet.tet;
          ori = orient3d(ppt[4], ppt[5], ppt[6], newpt); 
          if (ori < 0) {
            // Visible. Perform a 2-to-3 flip on the flip face.
            fliptets[0] = fliptet;  // [a,b,c,d], d = newpt.
            fliptets[1] = neightet; // [b,a,c,e], c = dummypoint.
            flip23(fliptets, 1, flipflag, chkencflag); // flip a hull tet.
            //recenttet = fliptets[0];
          } else if (ori == 0) {
            // Handle degenerate case ori == 0.
            if (oppo(neightet) == newpt) {
              // Two hull tets have the same base face.
              if (b->verbose > 2) {
                printf("      Close an open face (%d, %d, %d)\n", 
                       pointmark(org(fliptet)), pointmark(dest(fliptet)),
                       pointmark(apex(fliptet)));
              }
              // The following code connect adjacent tets at corresponding
              //   sides of the two hull tets. It is hard to understand.
              //   See an example in 2011-11-11.
              // First infect the two hull tets (they will be deleted).
              infect(fliptet);
              infect(neightet);
              // Connect the actual adjacent tets.
              for (i = 0; i < 3; i++) {
                fnext(fliptet, fliptets[0]);
                fnext(neightet, fliptets[1]);
                if (!infected(fliptets[0])) {
                  assert(!infected(fliptets[1]));
                  bond(fliptets[0], fliptets[1]);
                  // Update the point-to-tet map.
                  pa = org(fliptet);
                  pb = dest(fliptet);
                  setpoint2tet(pa, encode(fliptets[0]));
                  setpoint2tet(pb, encode(fliptets[0]));
                  // Remeber a recent tet for point location.
                  recenttet = fliptets[0];
                  // apex(fliptets[0]) is the new point. The opposite face may
                  // be not locally Delaunay. Put it in flip stack.
                  assert(apex(fliptets[0]) == newpt); // SELF_CHECK
                  esymself(fliptets[0]);
                  flippush(flipstack, &(fliptets[0]));
                  assert(apex(fliptets[1]) == newpt); // SELF_CHECK
                  esymself(fliptets[1]);
                  flippush(flipstack, &(fliptets[1]));                  
                }
                enextself(fliptet);
                eprevself(neightet);
              }
              // Delete the two tets.
              tetrahedrondealloc(fliptet.tet);
              tetrahedrondealloc(neightet.tet);
              // Update the hull size.
              hullsize -= 2;
            }
          }
        } // if (flipflag == 1)
 
        continue; // Do not flip a hull face.
      } // if (ishulltet(neightet))

      if (ishulltet(fliptet)) {
        continue; // Do not flip a hull tet.
      }
      
      if ((flipflag == 3) || (flipflag == 4)) {
        if (checksubfaceflag) {
          // Do not flip a subface.
          tspivot(fliptet, checksh);
          if (checksh.sh != NULL) {
            if (chkencflag & 2) {
              // Mesh refinement.
              // Put this subface into list.
              if (!smarktest2ed(checksh)) {
                bface = (badface *) badsubfacs->alloc();
                bface->ss = checksh;
                smarktest2(checksh); // Only queue it once.
                bface->forg = sorg(checksh); // An alive badface.
              }
            }
            continue;
          }
        }
      } // if ((flipflag == 3) || (flipflag == 4))

      ppt = (point *) fliptet.tet;
      pe = oppo(neightet);

      sign = insphere_s(ppt[4], ppt[5], ppt[6], ppt[7], pe);

      if (sign < 0) {
        if (b->verbose > 3) {
          printf("        A non-Delaunay face (%d, %d, %d) - %d, %d\n",
                 pointmark(org(fliptet)), pointmark(dest(fliptet)),
                 pointmark(apex(fliptet)), pointmark(oppo(fliptet)),
                 pointmark(pe));
        }

        // Try to flip this face.
        pd = oppo(fliptet);
        // Check the convexity of its three edges.
        convflag = 1;
        for (i = 0; i < 3; i++) {
          p1 = org(fliptet);
          p2 = dest(fliptet);
          ori = orient3d(p1, p2, pd, pe); 
          if (ori < 0) {
            // A locally non-convex edge.
            convflag = -1;
            break;  
          } else if (ori == 0) {
            // A locally flat edge.
            convflag = 0;
            break;
          }
          enextself(fliptet);
        }

        if (convflag > 0) {
          // A 2-to-3 flip is found.
          fliptets[0] = fliptet; // abcd, d may be the new vertex.
          fliptets[1] = neightet; // bace.
          if ((flipflag == 1) || (flipflag == 2)) { // CDT boundary recovery.
            if (checksubfaceflag) {
              // Check if a subface will be flipped.
              tspivot(fliptets[0], checksh);
              if (checksh.sh != NULL) {
                assert(flipflag < 3); // 1 or 2.
                // It is updateing a conforming DT or a CDT.
                if (b->verbose > 3) {
                  printf("        Queue a flipped subface (%d, %d, %d).\n",
                         pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                         pointmark(sapex(checksh)));
                }
                for (i = 0; i < 2; i++) {
                  tsdissolve(fliptets[i]); // Disconnect the tet->sub bond.
                }
                stdissolve(checksh); // Disconnect the sub->tet bond.
                // Add the missing subface into list.
                subfacstack->newindex((void **) &parysh);
                *parysh = checksh;
              } // if (checksh.sh != NULL)
            }
          } // if ((flipflag == 1) || (flipflag == 2))
          flip23(fliptets, 0, flipflag, chkencflag);
          //recenttet = fliptets[0]; // for point location.
        } else {
          // The edge ('fliptet') is non-convex or flat.
          if ((flipflag == 3) || (flipflag == 4)) {
            // Do not flip a subsegment.
            tsspivot1(fliptet, checkseg);
            if (checkseg.sh != NULL) {
              if (b->verbose > 3) {
                printf("        Found a non-Delaunay segment (%d, %d).\n",
                       pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
              }
              // Comment: this should be only possible when a new Steiner
              //   point is inserted on a segment nearby.
              if (chkencflag & 1) {
                // Put this segment into list.
                if (!smarktest2ed(checkseg)) {
                  bface = (badface *) badsubsegs->alloc();
                  bface->ss = checkseg;
                  smarktest2(checkseg); // Only queue it once.
                  bface->forg = sorg(checkseg); // An alive badface.
                }
              }
              continue;
            }
          }

          // A 3-to-2 or 4-to-4 may be possible.
          esym(fliptet, fliptets[0]); // [b,a,d,c]
          // assert(apex(fliptets[0]) == pd);
          n = 0;
          do {
            fnext(fliptets[n], fliptets[n + 1]);
            n++;
          } while ((fliptets[n].tet != fliptet.tet) && (n < 5));

          if (n == 3) {
            // Found a 3-to-2 flip.
            if ((flipflag == 1) || (flipflag == 2)) { // CDT boundary recovery.
              if (checksubsegflag) {
                // Check if the flip edge is subsegment.
                tsspivot1(fliptets[0], checkseg);
                if (checkseg.sh != NULL) {
                  if (!sinfected(checkseg)) {
                    // This subsegment will be flipped. Queue it.
                    if (b->verbose > 3) {
                      printf("        Queue a flipped segment (%d, %d).\n",
                        pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
                    }
                    sinfect(checkseg);  // Only save it once.
                    subsegstack->newindex((void **) &paryseg);
                    *paryseg = checkseg;
                  }
                  // Clean tet-to-seg pointers.
                  for (i = 0; i < 3; i++) {
                    tssdissolve1(fliptets[i]);
                  }
                  // Clean the seg-to-tet pointer.
                  sstdissolve1(checkseg);
                }
              }
              if (checksubfaceflag) {
                // Check if there are subfaces to be flipped.
                for (i = 0; i < 3; i++) {
                  tspivot(fliptets[i], checksh);
                  if (checksh.sh != NULL) {//if (flipshs[i].sh != NULL) {
                    if (b->verbose > 2) {
                      printf("        Queue a flipped subface (%d, %d, %d).\n",
                        pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                        pointmark(sapex(checksh)));
                    }
                    tsdissolve(fliptets[i]); // Disconnect the tet->sub bond.
                    stdissolve(checksh); // Disconnect the sub->tet bond.
                    // Add the missing subface into list.
                    subfacstack->newindex((void **) &parysh);
                    *parysh = checksh;
                  }
                }
              }
            } // if ((flipflag == 1) || (flipflag == 2))

            // Now flip the edge.
            flip32(fliptets, 0, flipflag, chkencflag);
            //recenttet = fliptets[0]; // for point location.
          } else {
            // There are more than 3 tets shared at this edge.
            if ((n == 4) && (convflag < 0)) {
              // Check if a 4-to-4 flip is possible.
              pf = apex(fliptets[3]);
              if (pf == dummypoint) {
                // It is a non-convex hull edge shared by four tets (two hull
                //   tets and two interior tets). 
                // Let the two interior tets be [a,b,c,d] and [b,a,c,e] where
                //   [a,b] be the hull edge, [a,b,c] be the interior face.
                //   [a,b,d] and [a,b,e] are two hull faces.
                //   A 4-to-4 flip is possible if the two new tets [e,d,b,c]
                //   and [e,d,c,a] are valid tets.
                // Current status:
                //   'fliptets[0]' is [a,b,e,c]
                //   'fliptets[1]' is [a,b,c,d]
                //   'fliptets[2]' is [a,b,d,f] (hull tet)
                //   'fliptets[3]' is [a,b,f,e] (hull tet)
                pa =  org(fliptets[1]);
                pb = dest(fliptets[1]);
                pc = apex(fliptets[1]);
                p1 = oppo(fliptets[1]); // pd
                p2 = apex(fliptets[0]); // pe
                ori = orient3d(p2, p1, pb, pc);
                if (ori < 0) {
                  ori = orient3d(p2, p1, pc, pa);
                  if (ori < 0) {
                    convflag = -2; // A 4-to-4 flip is possible.
                  }
                }
              }
	    } // if ((n == 4) && (convflag < 0))
            if ((n == 4) && ((convflag == 0) || (convflag == -2))) {
              // Found a 4-to-4 flip.
              if (b->verbose > 3) {
                printf("        A 4-to-4 flip (%d, %d) - (%d, %d).\n",
                       pointmark(org(fliptet)), pointmark(dest(fliptet)),
                       pointmark(pd), pointmark(pe));
              }
              if ((flipflag == 1) || (flipflag == 2)) { // CDT boundary recovery
                if (checksubsegflag) {
                  // Check if the flip edge is subsegment.
                  tsspivot1(fliptets[0], checkseg);
                  if (checkseg.sh != NULL) {
                    if (!sinfected(checkseg)) {
                      // This subsegment will be flipped. Queue it.
                      if (b->verbose > 3) {
                        printf("        Queue a flipped segment (%d, %d).\n",
                         pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
                      }
                      sinfect(checkseg);  // Only save it once.
                      subsegstack->newindex((void **) &paryseg);
                      *paryseg = checkseg;
                    }
                    // Clean the tet-to-seg pointers.
                    for (i = 0; i < 4; i++) {
                      tssdissolve1(fliptets[i]);
                    }
                    // Clean the seg-to-tet pointer.
                    sstdissolve1(checkseg);
                  }
                }
                if (checksubfaceflag) {
                  // Check if there are subfaces to be flipped.
                  for (i = 0; i < 4; i++) {
                    tspivot(fliptets[i], checksh);
                    if (checksh.sh != NULL) {
                      if (b->verbose > 3) {
                        printf("        Queue a flipped subface (%d,%d,%d).\n",
                          pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                          pointmark(sapex(checksh)));
                      }
                      tsdissolve(fliptets[i]); // Disconnect the tet->sub bond.
                      stdissolve(checksh); // Disconnect the sub->tet bond.
                      // Add the missing subface into list.
                      subfacstack->newindex((void **) &parysh);
                      *parysh = checksh;
                    }
                  }
                }
              } // if ((flipflag == 1) || (flipflag == 2))

              // First do a 2-to-3 flip.
              // Comment: This flip temporarily creates either a degenerated
              //   tet (convflag == 0) or an inverted tet (convflag < 0).
              //   It is removed by the followed 3-to-2 flip.
              fliptets[0] = fliptet; // tet abcd, d is the new vertex.
              baktets[0] = fliptets[2];
              baktets[1] = fliptets[3];
              // The flip may involve hull tets.
              flip23(fliptets, 1, flipflag, chkencflag);
              // Then do a 3-to-2 flip.
              enextesymself(fliptets[0]);  // fliptets[0] is edab.
              eprevself(fliptets[0]); // tet badc, d is the new vertex.
              fliptets[1] = baktets[0];
              fliptets[2] = baktets[1];
              flip32(fliptets, 1, flipflag, chkencflag);
              flip23count--;
              flip32count--;
              flip44count++;
              //recenttet = fliptets[0]; // for point location.
            } else {
              // This edge is shared by more than 4 tets.
              if (b->verbose > 2) {
                printf("        An unflippable non-Delaunay edge (%d,%d).\n",
                       pointmark(org(fliptet)), pointmark(dest(fliptet)));
              }
              remflag = 0;
              if (flipedgeflag == 2) {
                // Try to flip this edge by my edge flip algorithm.
                // Remember the the objective value (volume of all tetprisms).
                fc.bak_tetprism_vol = tetprism_vol_sum; 
                if (removeedgebyflips(&fliptet, &fc) == 2) {
                  if (b->verbose > 2) {
                    printf("      Decreased quantity: %.17g.\n", 
                           fc.bak_tetprism_vol - tetprism_vol_sum);
                  }
                  // Queue new faces in flipstack.
                  for (i = 0; i < cavetetlist->objects; i++) {
                    parytet = (triface *) fastlookup(cavetetlist, i);
                    if (!isdeadtet(*parytet)) { // Skip a dead tet.
                      for (parytet->ver = 0; parytet->ver < 4; parytet->ver++) {
                        // Avoid queue a face twice.
                        fsym(*parytet, neightet);
                        if (!facemarked(neightet)) {
                          //flippush(flipstack, parytet);
                          bface = (badface *) flippool->alloc();
                          bface->tt = *parytet;
                          markface(bface->tt);
                          bface->forg = org(bface->tt); // An alive badface.
                          bface->fdest = dest(bface->tt); 
                          bface->fapex = apex(bface->tt); 
                          // bface->foppo = oppo(bface->tt);
                          // Push this face into stack.
                          bface->nextitem = flipstack;
                          flipstack = bface;
                        }
                      } // parytet->ver
                    }
                  } // i
                  cavetetlist->restart();
                  remflag = 1;
                }
              }
              if (!remflag) {
                // Found an unflippable non-Delaunay edge.
                if (flipedgeflag > 0) { // if (flipflag > 1) {
                  // Save this face (of the edge) in a second queue.
                  unflipqueue->newindex((void **) &bface);
                  bface->tt = fliptet;
                  bface->forg = org(fliptet);
                  bface->fdest = dest(fliptet);
                  bface->fapex = apex(fliptet); // FOR DEBUG.
                }
              }
            }
          } // if (n > 3)
        } // if (convflag <= 0)
      } // if (sign < 0)

    } // while (flipstack != NULL)


    break;

  } // while (1)


  if (b->verbose > 2) {
    printf("      Total %ld flips", flip23count + flip32count + flip44count
           - flipcount);
    if ((flipflag == 4) || peelsliverflag) {
      printf(", %ld sliver peels", opt_sliver_peels - tetpeelcount);
    }
    printf("\n");
  }


  return flip23count + flip32count + flip44count - flipcount;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertvertex()    Insert a point into current tetrahedralization.         //
//                                                                           //
// The Bowyer-Watson (B-W) algorithm is used to add a new point p into the   //
// tetrahedralization T. It first finds a "cavity", denoted as C, in T,  C   //
// consists of tetrahedra in T that "conflict" with p.  If T is a Delaunay   //
// tetrahedralization, then all boundary faces (triangles) of C are visible  //
// by p, i.e.,C is star-shaped. We can insert p into T by first deleting all //
// tetrahedra in C, then creating new tetrahedra formed by boundary faces of //
// C and p.  If T is not a DT, then C may be not star-shaped.  It must be    //
// modified so that it becomes star-shaped.                                  //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::insertvertex(point insertpt, triface *searchtet, face *splitsh,
                             face *splitseg, insertvertexflags *ivf)
{
  arraypool *swaplist; // for updating cavity.
  triface *cavetet, spintet, neightet, neineitet, *parytet;
  triface oldtet, newtet, newneitet;
  face checksh, *parysh, neighsh, spinsh;
  face checkseg, *paryseg;
  point *pts, pa, pb, pc, *parypt;
  badface *bface;
  enum locateresult loc;
  REAL sign, ori;
  REAL rd, cent[3];
  REAL attrib, volume;
  long cutcount, cutshcount, tetcount = 0;
  long bakhullsize;
  bool enqflag;
  int i, j, k, s;

  int rejptflag, encptflag; // for protecting balls.
  int bgmloc;


  if (b->verbose > 2) {
    printf("      Insert point %d\n", pointmark(insertpt));
  }


  // Locate the point.
  loc = OUTSIDE; // Set a default value.

  if (searchtet->tet != NULL) {
    loc = (enum locateresult) ivf->iloc;
  }

  if (loc == OUTSIDE) {
    tetcount = ptloc_count; // Count the number of visited tets.
    if (searchtet->tet == NULL) {
      if (!b->weighted) {
        if (b->brio_hilbert) { // -b
          *searchtet = recenttet;
        } else { // -b0
          randomsample(insertpt, searchtet);
        }
      } else {
        // There may exist dangling vertex. 
        *searchtet = recenttet;
      }
    }
    // Locate the point.
    loc = locate(insertpt, searchtet, ivf->chkencflag); 
    if (b->verbose > 3) {
        printf("        Walk distance (# tets): %ld\n", ptloc_count-tetcount);
    }
    if (ptloc_max_count < (ptloc_count - tetcount)) {
      ptloc_max_count = (ptloc_count - tetcount);
    }
  }

  if (b->verbose > 3) {
    printf("        Located tet (%d, %d, %d, %d).\n",
           pointmark(org(*searchtet)), pointmark(dest(*searchtet)), 
           pointmark(apex(*searchtet)), pointmark(oppo(*searchtet)));
  }


  if (b->weighted) {
    if (loc != OUTSIDE) {
      // Check if this vertex is regular.
      pts = (point *) searchtet->tet;
      assert(pts[7] != dummypoint);
      sign = orient4d_s(pts[4], pts[5], pts[6], pts[7], insertpt,
                        pts[4][3], pts[5][3], pts[6][3], pts[7][3],
                        insertpt[3]);
      if (sign > 0) {
        // This new vertex does not lie below the lower hull. Skip it.
        if (b->verbose > 1) {
	  printf("    Point #%d is non-regular, skipped.\n",
                 pointmark(insertpt));
	}
        setpointtype(insertpt, NREGULARVERTEX);
        nonregularcount++;
        return NONREGULAR;
      }
    }
  }

  // Create the initial cavity C(p) which contains all tetrahedra directly
  //   intersect with p.

  // Remember the current hullsize. It is used to restore the hullsize
  //   if the new point is rejected for insertion.
  bakhullsize = hullsize;

  if (loc == OUTSIDE) {
    if (b->verbose > 3) {
      printf("        Outside hull.\n");
    }
    // The current hull will be enlarged.
    // Add four adjacent boundary tets into list.
    for (i = 0; i < 4; i++) {
      decode(searchtet->tet[i], neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    if ((point) searchtet->tet[7] == dummypoint) hullsize--;
    // tetrahedrondealloc(searchtet->tet);
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;
    flip14count++;
  } else if (loc == INTETRAHEDRON) {
    if (b->verbose > 3) {
      printf("        Inside tet.\n");
    }
    // Add four adjacent boundary tets into list.
    for (i = 0; i < 4; i++) {
      decode(searchtet->tet[i], neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    // tetrahedrondealloc(searchtet->tet);
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;
    flip14count++;
  } else if (loc == ONFACE) {
    if (b->verbose > 3) {
      printf("        On face.\n");
    }
    // Add six adjacent boundary tets into list.
    j = (searchtet->ver & 3); // The current face number.
    for (i = 1; i < 4; i++) { 
      decode(searchtet->tet[(j + i) % 4], neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    decode(searchtet->tet[j], spintet);
    j = (spintet.ver & 3); // The current face number.
    for (i = 1; i < 4; i++) {
      decode(spintet.tet[(j + i) % 4], neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    if ((point) spintet.tet[7] == dummypoint) hullsize--;
    if ((point) searchtet->tet[7] == dummypoint) hullsize--;
    // tetrahedrondealloc(spintet.tet);
    infect(spintet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = spintet;
    // tetrahedrondealloc(searchtet->tet);
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;
    flip26count++;

    if (ivf->splitbdflag) { //if (bowywat > 2) {
      if (splitsh != NULL) {
        // Create the initial sub-cavity sC(p).
        smarktest(*splitsh);
        caveshlist->newindex((void **) &parysh);
        *parysh = *splitsh;
      }
    } // if (splitbdflag)
  } else if (loc == ONEDGE) {
    if (b->verbose > 3) {
      printf("        On edge.\n");
    }
    // Add all adjacent boundary tets into list.
    spintet = *searchtet;
    while (1) {
      enextesym(spintet, neightet);
      fsymself(neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
      eprevesym(spintet, neightet);
      fsymself(neightet);
      neightet.ver = epivot[neightet.ver & 3];
      cavetetlist->newindex((void **) &parytet);
      *parytet = neightet;
      if ((point) spintet.tet[7] == dummypoint) hullsize--;
      // tetrahedrondealloc(spintet.tet);
      infect(spintet);
      caveoldtetlist->newindex((void **) &parytet);
      *parytet = spintet;
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)
    flipn2ncount++;

    if (ivf->splitbdflag) { //if (bowywat > 2) {
      // Create the initial sub-cavity sC(p).
      if (splitseg != NULL) {
        smarktest(*splitseg);
        splitseg->shver = 0;
        spivot(*splitseg, *splitsh);
      }
      if (splitsh != NULL) {
        if (splitsh->sh != NULL) {
          // Collect all subfaces share at this edge.
          pa = sorg(*splitsh);
          neighsh = *splitsh;
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
            if (neighsh.sh == splitsh->sh) break;
            if (neighsh.sh == NULL) break;
          } // while (1)
        } // if (not a dangling segment)
      }
    } // if (splitbdflag)
  } else if (loc == INSTAR) {
    if (b->verbose > 3) {
      printf("        Inside star.\n");
    }
    // We assume that all tets in the star are given in 'caveoldtetlist',
    //   and they are all infected.
    assert(caveoldtetlist->objects > 0);
    // Collect the boundary faces of the star.
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      // Check its 4 neighbor tets.
      for (j = 0; j < 4; j++) {
        decode(cavetet->tet[j], neightet);
        if (!infected(neightet)) {
          // It's a boundary face.
          neightet.ver = epivot[neightet.ver & 3];
          cavetetlist->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
    }
  } else if (loc == ONVERTEX) {
    pa = org(*searchtet);
    if (b->verbose > 3) {
      printf("        On vertex %d.\n", pointmark(pa));
    }
    if (insertpt != pa) {
      // Remember it is a duplicated point.
      setpointtype(insertpt, DUPLICATEDVERTEX);
      // Set a pointer to the point it duplicates.
      setpoint2ppt(insertpt, pa);
    }
    // The point already exist. Do nothing and return.
    return (int) loc;
  } else if (loc == ENCSUBFACE) {
    if (b->verbose > 3) {
      printf("        Beyond boundary.\n");
    }
    // The vertex lies outside of the region boundary.
    if (ivf->rejflag & 2) {
      // Check if this vertex lies very close to the boundary face.
      // This case needs to be handled due to the rounding off error.
      tspivot(*searchtet, checksh);
      assert(checksh.sh != NULL);
      pa = sorg(checksh);
      pb = sdest(checksh);
      pc = sapex(checksh);
      ori = orient3d(pa, pb, pc, insertpt);
      // Re-use cent[3].
      cent[0] = distance(pa, insertpt);
      cent[1] = distance(pb, insertpt);
      cent[2] = distance(pb, insertpt);
      // Choose the largest distance.
      if (cent[0] < cent[1]) cent[0] = cent[1];
      if (cent[0] < cent[2]) cent[0] = cent[2];
      if (fabs(ori) / (cent[0] * cent[0] * cent[0]) < b->epsilon) {
        // A nearly co-planar subface. We treat this case as coplanar, so
        //   the insertion point does not lie outside of the domain. 
        // Queue an encroached subface.
        // Calculate the circumcenter of this subface (for refinement).
        circumsphere(pa, pb, pc, NULL, cent, &rd);
        encshlist->newindex((void **) &bface);
        bface->ss = checksh;
        bface->forg = pa; // Not a dad one.
        for (j = 0; j < 3; j++) bface->cent[j] = cent[j];
        bface->key = rd;
        return (int) ENCSUBFACE;
      }
    }
    // Treated it as outside
    loc = OUTSIDE;
    return (int) loc;
  } else {
    assert(0); // Unknown type.
  }


  if (ivf->assignmeshsize) {
    // Assign mesh size for the new point.
    if (bgm != NULL) {
      // Interpolate the mesh size from the background mesh. 
      pa = org(*searchtet);
      bgm->decode(point2bgmtet(pa), neightet); // neightet is in 'bgm'!
      bgmloc = bgm->scoutpoint(insertpt, &neightet, 0); // randflag = 0
      if (bgmloc != (int) OUTSIDE) {
        insertpt[pointmtrindex] =  
          bgm->getpointmeshsize(insertpt, &neightet, bgmloc);
        setpoint2bgmtet(insertpt, bgm->encode(neightet));
      }
    } else {
      insertpt[pointmtrindex] = getpointmeshsize(insertpt,searchtet,(int)loc);
    }
  } // if (assignmeshsize)

  if (ivf->validflag) { //if (bowywat > 2) { 
    // Validate the initial C(p). Enlarge it at a face which is not visible
    //   by p. This removes (interior) slivers. Re-use 'cavebdrylist'.
    tetcount = 0l;

    for (i = 0; i < cavetetlist->objects; i++) {
      cavetet = (triface *) fastlookup(cavetetlist, i);
      // Other expansions may make this face inside C(p).
      if (!infected(*cavetet)) {
        pc = apex(*cavetet);
        // Do valid if it is a face (not a hull edge).
        if (pc != dummypoint) {
          pa = org(*cavetet);
          pb = dest(*cavetet);
          ori = orient3d(pa, pb, pc, insertpt); 
          if (ori <= 0) {
            // An invalid face. Enlarge the cavity.
              if (b->verbose > 3) {
                printf("        Enlarge cavity at (%d, %d, %d)\n",
                       pointmark(pa), pointmark(pb), pointmark(pc));
              }
              // Add the other three faces into list.
              j = (cavetet->ver & 3); // The current face number.
              for (k = 1; k < 4; k++) { 
                decode(cavetet->tet[(j + k) % 4], neightet);
                neightet.ver = epivot[neightet.ver & 3];
                cavetetlist->newindex((void **) &parytet);
                *parytet = neightet;
              }
              if ((point) cavetet->tet[7] == dummypoint) hullsize--;
              infect(*cavetet);
              caveoldtetlist->newindex((void **) &parytet);
              *parytet = *cavetet;
              tetcount++;
          } else {
            // A valid face.
            cavebdrylist->newindex((void **) &parytet);
            *parytet = *cavetet;
          }
        } else {
          // A hull edge is valid.
          cavebdrylist->newindex((void **) &parytet);
          *parytet = *cavetet;
        }
      } // if (!infected(*cavetet))
    } // i

    if (tetcount > 0l) {
      // The cavity has been enlarged. Update it.
      cavetetlist->restart();
      for (i = 0; i < cavebdrylist->objects; i++) {
        cavetet = (triface *) fastlookup(cavebdrylist, i);
        if (!infected(*cavetet)) {
          cavetetlist->newindex((void **) &parytet);
          *parytet = *cavetet;
        }
      } // i
    } // if (tetcount)

    cavebdrylist->restart();
    tetcount = 0l;
  } // if (bowywat > 2)

  // Update the cavity C(p) using the Bowyer-Watson approach (bowywat > 0). 

  for (i = 0; i < cavetetlist->objects; i++) {
    // 'cavetet' is an adjacent tet at outside of the cavity.
    cavetet = (triface *) fastlookup(cavetetlist, i);
    // The tet may be tested and included in the (enlarged) cavity.
    if (!infected(*cavetet)) {
      // Check for two possible cases for this tet: 
      //   (1) It is a cavity tet, or
      //   (2) it is a cavity boundary face.
      // In case (1), this tet is grabbed in the cavity and three adjacent 
      //   tets on other faces of this tet are added into 'cavetetlist'.
      enqflag = false;
      if (!marktested(*cavetet)) {
        if (ivf->bowywat) {
          // Do Delaunay (in-sphere) test.
          pts = (point *) cavetet->tet;
          if (pts[7] != dummypoint) {
            // A volume tet. Operate on it.
            if (b->weighted) {
              sign = orient4d_s(pts[4], pts[5], pts[6], pts[7], insertpt,
                                pts[4][3], pts[5][3], pts[6][3], pts[7][3],
                                insertpt[3]);
            } else {
              sign = insphere_s(pts[4], pts[5], pts[6], pts[7], insertpt);
            }
            enqflag = (sign < 0.0);
          } else {
            if (!nonconvex) {
              // Test if this hull face is visible by the new point. 
              ori = orient3d(pts[4], pts[5], pts[6], insertpt); 
              if (ori < 0) {
                // A visible hull face. 
                //if (!nonconvex) { 
                // Include it in the cavity. The convex hull will be enlarged.
                enqflag = true; // (ori < 0.0);
		//}
              } else if (ori == 0.0) {
                // A coplanar hull face. We need to test if this hull face is
                //   Delaunay or not. We test if the adjacent tet (not faked)
                //   of this hull face is Delaunay or not.
                neightet = *cavetet;
                neightet.ver = 3; // The face opposite to dummypoint.
                fsym(neightet, neineitet);
                if (!infected(neineitet)) {
                  if (!marktested(neineitet)) {
                    // Do Delaunay test on this tet.
                    pts = (point *) neineitet.tet;
                    assert(pts[7] != dummypoint);
                    if (b->weighted) {
                      sign = orient4d_s(pts[4],pts[5],pts[6],pts[7], insertpt,
                                        pts[4][3], pts[5][3], pts[6][3], 
                                        pts[7][3], insertpt[3]);
                    } else {
                      sign = insphere_s(pts[4],pts[5],pts[6],pts[7], insertpt);
                    }
                    enqflag = (sign < 0.0);
                  } else {
                    // The adjacent tet has been tested (marktested), and it
                    //   is Delaunay (not get infected). Hence the the hull
                    //   face is Delaunay as well.
                    // enqflag = false;
                  }
                } else {
                  // The adjacent tet is non-Delaunay. The hull face is non-
                  //   Delaunay as well. Include it in the cavity.
                  enqflag = true;
                } // if (!infected(neineitet))
              } // if (ori == 0.0)
            } else {
              // A hull face (must be a subface).
              assert(checksubfaceflag);
              assert(ivf->validflag);
              // We FIRST include it in the initial cavity if the adjacent tet
              //   (not faked) of this hull face is not Delaunay wrt p.
              //   Whether it belongs to the final cavity will be determined
              //   during the validation process. 'validflag'.
              neightet = *cavetet;
              neightet.ver = 3; // The face opposite to dummypoint.
              fsym(neightet, neineitet);
              if (!infected(neineitet)) {
                if (!marktested(neineitet)) {
                  // Do Delaunay test on this tet.
                  pts = (point *) neineitet.tet;
                  assert(pts[7] != dummypoint);
                  if (b->weighted) {
                    sign = orient4d_s(pts[4],pts[5],pts[6],pts[7], insertpt,
                                      pts[4][3], pts[5][3], pts[6][3], 
                                      pts[7][3], insertpt[3]);
                  } else {
                    sign = insphere_s(pts[4],pts[5],pts[6],pts[7], insertpt);
                  }
                  enqflag = (sign < 0.0);
                } else {
                  // The adjacent tet has been tested (marktested), and it
                  //   is Delaunay (not get infected). Hence the the hull
                  //   face is Delaunay as well.
                  // enqflag = false;
                } // if (marktested(neineitet))
              } else {
                // The adjacent tet is non-Delaunay. The hull face is non-
                //   Delaunay as well. Include it in the cavity.
                enqflag = true;
              } // if (infected(neineitet))
            } // if (nonconvex)
          } // if (pts[7] != dummypoint)
        } // if (bowywat)
        marktest(*cavetet); // Only test it once.
      } // if (!marktested(*cavetet))

      if (enqflag) {
        // Found a tet in the cavity. Put other three faces in check list.
        k = (cavetet->ver & 3); // The current face number
        for (j = 1; j < 4; j++) {
          decode(cavetet->tet[(j + k) % 4], neightet);
          neightet.ver = epivot[neightet.ver & 3];
          cavetetlist->newindex((void **) &parytet);
          *parytet = neightet;
        }
        if ((point) cavetet->tet[7] == dummypoint) hullsize--;
        // tetrahedrondealloc(cavetet->tet);
        infect(*cavetet);
        caveoldtetlist->newindex((void **) &parytet);
        *parytet = *cavetet;
      } else {
        // Found a boundary face of the cavity. It may be a face of a hull
        //   tet which contains 'dummypoint'. Choose the edge in the face 
        //   such that its endpoints are not 'dummypoint', while its apex
        //   may be 'dummypoint'.
        //j = (cavetet->ver & 3); // j is the face number.
        //cavetet->ver = epivot[j]; // [4,5,2,11]
        cavebdrylist->newindex((void **) &parytet);
        *parytet = *cavetet;
      }
    } // if (!infected(*cavetet))
  } // i

  if (b->verbose > 3) {
    printf("        Initial cavity size: %ld tets, %ld faces.\n", 
           caveoldtetlist->objects, cavebdrylist->objects);
  }


  if (checksubsegflag) {
    // Collect all segments of C(p).
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      for (j = 0; j < 6; j++) {
        cavetet->ver = edge2ver[j];
        tsspivot1(*cavetet, checkseg);
        if (checkseg.sh != NULL) {
          if (!sinfected(checkseg)) {
            sinfect(checkseg);
            cavetetseglist->newindex((void **) &paryseg);
            *paryseg = checkseg;
          }
        }
      }
    }
    // Uninfect collected segments.
    for (i = 0; i < cavetetseglist->objects; i++) {
      checkseg = * (face *) fastlookup(cavetetseglist, i);
      suninfect(checkseg);
    }
  } // if (checksubsegflag)

  if (checksubfaceflag) {
    // Collect all subfaces of C(p).
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      oldtet = *cavetet;
      for (oldtet.ver = 0; oldtet.ver < 4; oldtet.ver++) {
        tspivot(oldtet, checksh);
        if (checksh.sh != NULL) {
          if (!sinfected(checksh)) {
            sinfect(checksh);
            cavetetshlist->newindex((void **) &parysh);
            *parysh = checksh;
          }
        }
      }
    }
    // Uninfect collected subfaces.
    for (i = 0; i < cavetetshlist->objects; i++) {
      checksh = * (face *) fastlookup(cavetetshlist, i);
      suninfect(checksh);
    }
  } // if (checksubfaceflag)

  if (ivf->rejflag & 1) {
    // Reject insertion of this point if it encroaches upon any segment.
    for (i = 0; i < cavetetseglist->objects; i++) {
      checkseg = * (face *) fastlookup(cavetetseglist, i);
      pa = sorg(checkseg);
      pb = sdest(checkseg);
      if (checkseg4encroach(pa, pb, insertpt)) {
        if (b->verbose > 3) {
          printf("        Found an encroached seg (%d, %d).\n",
                 pointmark(pa), pointmark(pb));
        }
        encseglist->newindex((void **) &paryseg);
        *paryseg = checkseg;
      }
    } // i
    if (encseglist->objects > 0) {
      if (b->verbose > 3) {
        printf("        Found %ld encroached segments. Reject it.\n",
               encseglist->objects);
      }
      for (i = 0; i < caveoldtetlist->objects; i++) {
        cavetet = (triface *) fastlookup(caveoldtetlist, i);
        uninfect(*cavetet);
        unmarktest(*cavetet);
      }
      for (i = 0; i < cavebdrylist->objects; i++) {
        cavetet = (triface *) fastlookup(cavebdrylist, i);
        unmarktest(*cavetet); // Unmark it.
      }
      // Clear working lists.
      cavetetlist->restart();
      cavebdrylist->restart();
      caveoldtetlist->restart();
      cavetetseglist->restart();
      cavetetshlist->restart();
      if (ivf->splitbdflag) { //if (bowywat > 2) {
        if (splitseg != NULL) {
          sunmarktest(*splitseg);
        }
        for (i = 0; i < caveshlist->objects; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          assert(smarktested(*parysh));
          sunmarktest(*parysh);
        }
        caveshlist->restart();
        cavesegshlist->restart();
      }
      // Restore the hullsize.
      hullsize = bakhullsize;
      return (int) ENCSEGMENT;
    }
  } // if (reject & 1)

  if (ivf->rejflag & 2) {
    // Reject insertion of this point if it encroaches upon any subface.
    for (i = 0; i < cavetetshlist->objects; i++) {
      checksh = * (face *) fastlookup(cavetetshlist, i);
      pa = sorg(checksh);
      pb = sdest(checksh);
      pc = sapex(checksh);
      if (checkfac4encroach(pa, pb, pc, insertpt, cent, &rd)) {
        if (b->verbose > 3) {
          printf("        Found an encroached subface (%d, %d, %d).\n",
                 pointmark(pa), pointmark(pb), pointmark(pc));
        }
        encshlist->newindex((void **) &bface);
        bface->ss = checksh;
        bface->forg = pa; // Not a dad one.
        for (j = 0; j < 3; j++) bface->cent[j] = cent[j];
        bface->key = rd;
      }
    } // i
    if (encshlist->objects > 0) {
      if (b->verbose > 3) {
        printf("        Found %ld encroached subfaces. Reject it.\n",
               encshlist->objects);
      }
      for (i = 0; i < caveoldtetlist->objects; i++) {
        cavetet = (triface *) fastlookup(caveoldtetlist, i);
        uninfect(*cavetet);
        unmarktest(*cavetet);
      }
      for (i = 0; i < cavebdrylist->objects; i++) {
        cavetet = (triface *) fastlookup(cavebdrylist, i);
        unmarktest(*cavetet); // Unmark it.
      }
      cavetetlist->restart();
      cavebdrylist->restart();
      caveoldtetlist->restart();
      cavetetseglist->restart();
      cavetetshlist->restart();
      if (ivf->splitbdflag) { //if (bowywat > 2) {
        if (splitseg != NULL) {
          sunmarktest(*splitseg);
        }
        for (i = 0; i < caveshlist->objects; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          assert(smarktested(*parysh));
          sunmarktest(*parysh);
        }
        caveshlist->restart();
        cavesegshlist->restart();
      }
      // Restore the hullsize.
      hullsize = bakhullsize;
      return (int) ENCSUBFACE;
    }
  } // if (reject & 2)

  if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
    // Update the sC(p). 
    // We have already 'smarktested' the subfaces which directly intersect
    //   with p in 'caveshlist'. From them, we 'smarktest' their neighboring
    //   subfaces which are included in C(p). Do not across a segment.
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      assert(smarktested(*parysh));
      checksh = *parysh;
      for (j = 0; j < 3; j++) {
        sspivot(checksh, checkseg);
        if (checkseg.sh == NULL) {
          spivot(checksh, neighsh);
          assert(neighsh.sh != NULL);
          if (!smarktested(neighsh)) {
            stpivot(neighsh, neightet);
            if (infected(neightet)) {
              fsymself(neightet);
              if (infected(neightet)) {
                // This subface is inside C(p). 
                // Check if its diametrical circumsphere encloses 'p'.
                pa = sorg(neighsh);
                pb = sdest(neighsh);
                pc = sapex(neighsh);
                sign = incircle3d(pa, pb, pc, insertpt);
                if (sign < 0) {
                  smarktest(neighsh);
                  caveshlist->newindex((void **) &parysh);
                  *parysh = neighsh;
                }
              }
            }
          }
        }
        senextself(checksh);
      } // j
    } // i
    if (b->verbose > 3) {
      printf("        Initial subcavity size: %ld subfacess.\n", 
             caveshlist->objects);
    }
  }

  cutcount = 0l;

  if (ivf->validflag) { 
    //if (bowywat > 1) { // if (bowywat == 2 || bowywat == 3) {
    // T is a CT. Validation is needed (fig/dump-cavity-case8).
    cavetetlist->restart(); // Re-use it.

    //if (splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
    if (ivf->respectbdflag) {
      // The initial cavity may include subfaces which are not on the facets
      //   being splitting. Find them and make them as boundary of C(p).      
      // Comment: We have already 'smarktested' the subfaces in sC(p).
      //   It is needed by 'splitbdflag'.
      for (i = 0; i < cavetetshlist->objects; i++) {
        parysh = (face *) fastlookup(cavetetshlist, i);
        stpivot(*parysh, neightet);
        if (infected(neightet)) {
          fsymself(neightet);
          if (infected(neightet)) {
            if (!smarktested(*parysh)) {
              if (b->verbose > 3) {
                printf("        Found a subface (%d, %d, %d) inside cavity.\n", 
                       pointmark(sorg(*parysh)), pointmark(sdest(*parysh)),
                       pointmark(sapex(*parysh)));
              }
              // It is possible that this face is a boundary subface.
              // Check if it is a hull face.
              assert(apex(neightet) != dummypoint);
              if (oppo(neightet) != dummypoint) {
                fsymself(neightet);
              }
              if (oppo(neightet) != dummypoint) {
                pa = org(neightet);
                pb = dest(neightet);
                pc = apex(neightet);
                ori = orient3d(pa, pb, pc, insertpt);
                if (ori < 0) {
                  // A visible face, get its neighbor face.
                  fsymself(neightet);
                  ori = -ori; // It must be invisible by p.
                }
              } else {
                // A hull tet. It needs to be cut.
                ori = 1;
              }
              // Cut this tet if it is either invisible by or coplanar with p.
              if (ori >= 0) {
                if (b->verbose > 3) {
                  printf("        Cut tet (%d, %d, %d, %d)\n", 
                         pointmark(org(neightet)), pointmark(dest(neightet)),
                         pointmark(apex(neightet)), pointmark(oppo(neightet)));
                }
                uninfect(neightet);
                unmarktest(neightet);
                cutcount++;
                neightet.ver = epivot[neightet.ver & 3];
                cavebdrylist->newindex((void **) &parytet);
                *parytet = neightet;
                // Add three new faces to find new boundaries.
                for (j = 0; j < 3; j++) {
                  esym(neightet, neineitet);
                  neineitet.ver = epivot[neineitet.ver & 3];
                  cavebdrylist->newindex((void **) &parytet);
                  *parytet = neineitet;
                  enextself(neightet);
                }
                // Update hullsize.
                if (oppo(neightet) == dummypoint) hullsize++;
              } // if (ori >= 0) 
            }
          }
        }
      } // i

      // The initial cavity may include segments in its interior. We need to
      //   Update the cavity so that these segments are on the boundary of
      //   the cavity.
      for (i = 0; i < cavetetseglist->objects; i++) {
        paryseg = (face *) fastlookup(cavetetseglist, i);
        // Check this segment if it is not a splitting segment.
        if (!smarktested(*paryseg)) {
          sstpivot1(*paryseg, neightet);
          spintet = neightet;
          while (1) {
            if (!infected(spintet)) break;
            fnextself(spintet);
            if (spintet.tet == neightet.tet) break;
          }
          if (infected(spintet)) {
            if (b->verbose > 3) {
              printf("        Found an interior segment (%d, %d).\n", 
                     pointmark(sorg(*paryseg)), pointmark(sdest(*paryseg)));
            }
            // Find an adjacent tet at this segment such that both faces
            //   at this segment are not visible by p.
            pa = org(neightet);
            pb = dest(neightet);
            spintet = neightet;
            j = 0;
            while (1) {
              // Check if this face is visible by p.
              pc = apex(spintet);
              if (pc != dummypoint) {
                ori = orient3d(pa, pb, pc, insertpt);
                if (ori >= 0) {
                  // Not visible. Check another face in this tet.
                  esym(spintet, neineitet);
                  pc = apex(neineitet);
                  if (pc != dummypoint) {
                    ori = orient3d(pb, pa, pc, insertpt);
                    if (ori >= 0) {
                      // Not visible. Found this face.
                      j = 1; // Flag that it is found.
                      break;
                    }
                  }
                }
              } else {
              }
              fnextself(spintet);
              if (spintet.tet == neightet.tet) break;
            }
            if (j == 0) {
              // Not found such a face.
              assert(0); // debug this case.
            }
            neightet = spintet;
            if (b->verbose > 3) {
               printf("        Cut tet (%d, %d, %d, %d)\n", 
                      pointmark(org(neightet)), pointmark(dest(neightet)),
                      pointmark(apex(neightet)), pointmark(oppo(neightet)));
            }
            uninfect(neightet);
            unmarktest(neightet);
            cutcount++;
            neightet.ver = epivot[neightet.ver & 3];
            cavebdrylist->newindex((void **) &parytet);
            *parytet = neightet;
            // Add three new faces to find new boundaries.
            for (j = 0; j < 3; j++) {
              esym(neightet, neineitet);
              neineitet.ver = epivot[neineitet.ver & 3];
              cavebdrylist->newindex((void **) &parytet);
              *parytet = neineitet;
              enextself(neightet);
            }
            // Update hullsize.
            //if (oppo(neightet) == dummypoint) hullsize++;
            if ((point) (neightet.tet[7]) == dummypoint) hullsize++;
          }
        }
      } // i
    } // if (bowywat > 2)

    // Update the cavity by removing invisible faces until it is star-shaped.
    for (i = 0; i < cavebdrylist->objects; i++) {
      cavetet = (triface *) fastlookup(cavebdrylist, i);
      // 'cavetet' is an exterior tet adjacent to the cavity.      
      assert(cavetet->ver == epivot[cavetet->ver & 3]); // SELF_CHECK
      // It must be not inside the cavity (since we only cut tets).
      assert(!infected(*cavetet));
      // Check if its neighbor is inside C(p).
      fsym(*cavetet, neightet);
      if (infected(neightet)) {        
        if (apex(*cavetet) != dummypoint) {
          // It is a cavity boundary face. Check its visibility.
          if (oppo(neightet) != dummypoint) {
            pa = org(*cavetet);
            pb = dest(*cavetet);
            pc = apex(*cavetet);
            ori = orient3d(pa, pb, pc, insertpt); 
            enqflag = (ori > 0);
            // Comment: if ori == 0 (coplanar case), we also cut the tet.
          } else {
            // It is a hull face. And its adjacent tet (at inside of the 
            //   domain) has been cut from the cavity. Cut it as well.
            //assert(nonconvex);
            enqflag = false;
          }
        } else {
          enqflag = true; // A hull edge.
        }
        if (enqflag) {
          // This face is valid, save it.
          cavetetlist->newindex((void **) &parytet);
          *parytet = *cavetet; 
        } else {
          if (b->verbose > 3) {
            printf("        Cut tet (%d, %d, %d, %d)\n", 
                   pointmark(org(neightet)), pointmark(dest(neightet)),
                   pointmark(apex(neightet)), pointmark(oppo(neightet)));
          }
          uninfect(neightet);
          unmarktest(neightet);
          cutcount++;
          // Add three new faces to find new boundaries.
          for (j = 0; j < 3; j++) {
            esym(neightet, neineitet);
            neineitet.ver = epivot[neineitet.ver & 3];
            cavebdrylist->newindex((void **) &parytet);
            *parytet = neineitet;
            enextself(neightet);
          }
          // Update the hullsize.
          if (oppo(neightet) == dummypoint) hullsize++;
          // 'cavetet' is not on the cavity boundary anymore.
          unmarktest(*cavetet);
        }
      } else {
        // 'cavetet' is not on the cavity boundary anymore.
        unmarktest(*cavetet);
      }
    } // i

    if (cutcount > 0) {
      // The cavity has been updated.

      // Update the cavity boundary faces.
      cavebdrylist->restart();
      for (i = 0; i < cavetetlist->objects; i++) {
        cavetet = (triface *) fastlookup(cavetetlist, i);
        // 'cavetet' was an exterior tet adjacent to the cavity.
        assert(cavetet->ver == epivot[cavetet->ver & 3]); // SELF_CHECK
        assert(!infected(*cavetet));
        fsym(*cavetet, neightet);
        if (infected(neightet)) {
          // It is a cavity boundary face.
          cavebdrylist->newindex((void **) &parytet);
          *parytet = *cavetet;
        } else {
          // Not a cavity boundary face.
          unmarktest(*cavetet);
        }
      }

      // Update the list of old tets.
      cavetetlist->restart();
      for (i = 0; i < caveoldtetlist->objects; i++) {
        cavetet = (triface *) fastlookup(caveoldtetlist, i);
        if (infected(*cavetet)) {
          cavetetlist->newindex((void **) &parytet);
          *parytet = *cavetet;
        }
      }
      // Swap 'cavetetlist' and 'caveoldtetlist'.
      swaplist = caveoldtetlist;
      caveoldtetlist = cavetetlist;
      cavetetlist = swaplist;

      // The cavity should contain at least one tet.
      if (caveoldtetlist->objects == 0l) {
        assert(cavebdrylist->objects == 0l);
        cavetetlist->restart();
        cavebdrylist->restart();
        caveoldtetlist->restart();
        cavetetseglist->restart();
        cavetetshlist->restart();
        if (ivf->splitbdflag) { //if (bowywat > 2) {
          if (splitseg != NULL) {
            sunmarktest(*splitseg);
          }
          for (i = 0; i < caveshlist->objects; i++) {
            parysh = (face *) fastlookup(caveshlist, i);
            assert(smarktested(*parysh));
            sunmarktest(*parysh);
          }
          caveshlist->restart();
          cavesegshlist->restart();
        }
        // Restore the hullsize.
        hullsize = bakhullsize;
        return (int) BADELEMENT;
      }

      if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
        cutshcount = 0;
        // Update the sub-cavity sC(p).
        for (i = 0; i < caveshlist->objects; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          if (smarktested(*parysh)) {
            enqflag = false;
            stpivot(*parysh, neightet);
            if (infected(neightet)) {
              fsymself(neightet);
              if (infected(neightet)) {
                enqflag = true;
              }
            }
            if (!enqflag) {
              if (b->verbose > 3) {
                printf("        Cut subface (%d, %d, %d).\n",
                       pointmark(sorg(*parysh)), pointmark(sdest(*parysh)),
                       pointmark(sapex(*parysh)));
              }
              sunmarktest(*parysh);
              // Use the last entry of this array to fill this entry.
              j = caveshlist->objects - 1;
              checksh = * (face *) fastlookup(caveshlist, j);
              *parysh = checksh;
              cutshcount++;
              caveshlist->objects--; // The list is shrinked.
              i--;
            }
          }
        }

        if (cutshcount > 0) {
          i = 0; // Count the number of invalid subfaces/segments.
          // Valid the updated sub-cavity sC(p).
          if (loc == ONFACE) {
            if (splitsh != NULL) {
              // The to-be split subface should be in sC(p).
              if (!smarktested(*splitsh)) i++;
            }
          } else if (loc == ONEDGE) {
            if (splitseg != NULL) {
              // The to-be split segment should be in sC(p).
              if (!smarktested(*splitseg)) i++;
            }
            if (splitsh != NULL) {
              // All subfaces at this edge should be in sC(p).
              pa = sorg(*splitsh);
              neighsh = *splitsh;
              while (1) {
                // Adjust the origin of its edge to be 'pa'.
                if (sorg(neighsh) != pa) {
                  sesymself(neighsh);
                }
                // Add this face into list (in B-W cavity).
                if (!smarktested(neighsh)) i++;
                // Go to the next face at the edge.
                spivotself(neighsh);
                // Stop if all faces at the edge have been visited.
                if (neighsh.sh == splitsh->sh) break;
                if (neighsh.sh == NULL) break;
              } // while (1)
            }
          }

          if (i > 0) {
            // The updated sC(p) is invalid. Do not insert this vertex.
            if (b->verbose > 3) {
              printf("        Found %d invalid items. Reject it.\n", i);
            }
            for (i = 0; i < caveoldtetlist->objects; i++) {
              cavetet = (triface *) fastlookup(caveoldtetlist, i);
              uninfect(*cavetet);
              unmarktest(*cavetet);
            }
            for (i = 0; i < cavebdrylist->objects; i++) {
              cavetet = (triface *) fastlookup(cavebdrylist, i);
              unmarktest(*cavetet); // Unmark it.
            }
            cavetetlist->restart();
            cavebdrylist->restart();
            caveoldtetlist->restart();
            cavetetseglist->restart();
            cavetetshlist->restart();
            if (ivf->splitbdflag) { //if (bowywat > 2) {
              if (splitseg != NULL) {
                sunmarktest(*splitseg);
              }
              for (i = 0; i < caveshlist->objects; i++) {
                parysh = (face *) fastlookup(caveshlist, i);
                assert(smarktested(*parysh));
                sunmarktest(*parysh);
              }
              caveshlist->restart();
              cavesegshlist->restart();
            }
            // Restore the hullsize.
            hullsize = bakhullsize;
            return (int) BADELEMENT;
          }
        } // if (cutshcount > 0)
      } // if (bowywat > 2)

    } // if (cutcount > 0)

  } // if (validflag) // if (bowywat > 1)

  if (b->verbose > 3) {
    printf("        Final cavity: %ld tets, %ld faces.", 
           caveoldtetlist->objects, cavebdrylist->objects);
    if (cutcount > 0l) {
      printf(" Updated %ld times.", cutcount);
    }
    printf("\n");
  }


  if (ivf->refineflag) {
    // The new point is inserted by Delaunay refinement, i.e., it is the 
    //   circumcenter of a tetrahedron, or a subface, or a segment.
    //   Do not insert this point if the tetrahedron, or subface, or segment
    //   is not inside the final cavity.
    rejptflag = 0;
    if (ivf->refineflag == 1) {
      // The new point is the circumcenter of a tetrahedron.
      assert(!isdeadtet(ivf->refinetet));
      if (!infected(ivf->refinetet)) {
        rejrefinetetcount++;
        rejptflag = 1;
      }
    } else if (ivf->refineflag == 2) {
      // The new point is the circumcenter of a subface.
      assert(ivf->refinesh.sh != NULL);
      if (!smarktested(ivf->refinesh)) {
        rejrefineshcount++;
        rejptflag = 1;
      }
    }
    if (rejptflag) {
      if (b->verbose > 2) {
        printf("      Point %d does not refine its element. Rejected.\n",
               pointmark(insertpt));
      }
      // Restore the original status.
      for (i = 0; i < caveoldtetlist->objects; i++) {
        cavetet = (triface *) fastlookup(caveoldtetlist, i);
        uninfect(*cavetet);
        unmarktest(*cavetet);
      }
      for (i = 0; i < cavebdrylist->objects; i++) {
        cavetet = (triface *) fastlookup(cavebdrylist, i);
        unmarktest(*cavetet); // Unmark it.
      }
      // Clear working lists.
      cavetetlist->restart();
      cavebdrylist->restart();
      caveoldtetlist->restart();
      cavetetshlist->restart();
      cavetetseglist->restart();
      cavetetvertlist->restart();
      if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
        if (splitseg != NULL) {
          sunmarktest(*splitseg);
        }
        for (i = 0; i < caveshlist->objects; i++) {
          parysh = (face *) fastlookup(caveshlist, i);
          assert(smarktested(*parysh));
          sunmarktest(*parysh);
        }
        caveshlist->restart();
        cavesegshlist->restart();
      }

      // Restore the hullsize.
      hullsize = bakhullsize;
      loc = BADELEMENT;
      return (int) loc;
    } // if (rejptflag)
  } // if (ivf->refineflag)

  rejptflag = (ivf->rejflag & 4);
  encptflag = 0;

  if (b->weighted || b->plc || rejptflag) {
    // Collect all vertices of C(p).
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      assert(infected(*cavetet));
      pts = (point *) &(cavetet->tet[4]);
      for (j = 0; j < 4; j++) {
        if (pts[j] != dummypoint) {
          if (!pinfected(pts[j])) {
            pinfect(pts[j]);
            cavetetvertlist->newindex((void **) &parypt);
            *parypt = pts[j];
          }
        }
      } // j
    } // i
    if (b->verbose > 3) {
      printf("        %ld cavity vertices.\n", cavetetvertlist->objects);
    }
    // Uninfect all collected (cavity) vertices.
    for (i = 0; i < cavetetvertlist->objects; i++) {
      parypt = (point *) fastlookup(cavetetvertlist, i);
      puninfect(*parypt);
    }
    if (b->plc || rejptflag) {
      // Check if p is too close to an existing vertex.
      pts = NULL;
      for (i = 0; i < cavetetvertlist->objects; i++) {
        parypt = (point *) fastlookup(cavetetvertlist, i);
        rd = distance(*parypt, insertpt);
        // Is the point very close to an existing point?
        if (rd < b->minedgelength) {
          pts = parypt; 
          break;
        }
        if (rejptflag) {
          // Is the point encroaches upon an existing point?
          if (rd < (*parypt)[pointmtrindex]) {
            // The point lies inside the protection ball.
            if (b->verbose > 2) {
              printf("      Point %d lies in protball of %d. Rejected.\n",
                     pointmark(insertpt), pointmark(*parypt));
            }
            pts = parypt;
            encptflag = 1; 
            break;
          }
        }
      } // i
      if (pts != NULL) {
        // p is too close to *pts.
        if (ivf->iloc != (int) INSTAR) {
          if (pointmark(insertpt) <= in->numberofpoints) {
            // It's an input point.
            if (!b->quiet) {
              printf("Warning:  Point %d is replaced by point %d.\n",
                     pointmark(insertpt), pointmark(*pts));
            }
            // Count the number of duplicated points.
            dupverts++;
          } else { // It's a Steiner point.
            if (b->verbose) {
              if (!rejptflag) {
                printf("Warning:  Reject a Steiner point %d (close to %d).\n",
                       pointmark(insertpt), pointmark(*pts));
              }
            }
          }
          // Remember it is a duplicated point.
          setpointtype(insertpt, DUPLICATEDVERTEX);
          // Set a pointer to the point it duplicates.
          setpoint2ppt(insertpt, *pts);

          // Restore the original status.
          for (i = 0; i < caveoldtetlist->objects; i++) {
            cavetet = (triface *) fastlookup(caveoldtetlist, i);
            uninfect(*cavetet);
            unmarktest(*cavetet);
          }
          for (i = 0; i < cavebdrylist->objects; i++) {
            cavetet = (triface *) fastlookup(cavebdrylist, i);
            unmarktest(*cavetet); // Unmark it.
          }
          // Clear working lists.
          cavetetlist->restart();
          cavebdrylist->restart();
          caveoldtetlist->restart();
          cavetetshlist->restart();
          cavetetseglist->restart();
          cavetetvertlist->restart();
          if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
            if (splitseg != NULL) {
              sunmarktest(*splitseg);
            }
            for (i = 0; i < caveshlist->objects; i++) {
              parysh = (face *) fastlookup(caveshlist, i);
              assert(smarktested(*parysh));
              sunmarktest(*parysh);
            }
            caveshlist->restart();
            cavesegshlist->restart();
          }

          // Restore the hullsize.
          hullsize = bakhullsize;
          if (!encptflag) {
            loc = NEARVERTEX;
          } else {
            loc = ENCVERTEX; 
          }
          return (int) loc;
        } else {  // (iloc == (int) INSTAR)
          // The cavity is guaranteed to be valid by the caller of this 
          //   function. We still insert this vertex.
          if (b->verbose) {
            printf("Warning:  The Steiner point %d is very close to %d.\n",
                   pointmark(insertpt), pointmark(*pts));
          }
        }
      } // if (pts != NULL)
    } 
  } 

  // The new point will be inserted.
  totaldeadtets += caveoldtetlist->objects;
  totalbowatcavsize += cavebdrylist->objects;
  if (maxbowatcavsize < cavebdrylist->objects) {
    maxbowatcavsize = cavebdrylist->objects;
  }

  // Before re-mesh C(p). Process the segments and subfaces which are on the
  //   boundary of C(p). Make sure that each such segment or subface is
  //   connecting to a tet outside C(p). So we can re-connect them to the
  //   new tets inside the C(p) later.

  if (checksubsegflag) {
    for (i = 0; i < cavetetseglist->objects; i++) {
      paryseg = (face *) fastlookup(cavetetseglist, i);
      // Operate on it if it is not the splitting segment, i.e., in sC(p).
      if (!smarktested(*paryseg)) {
        // Check if the segment is inside the cavity.
        //   'j' counts the num of adjacent tets of this seg.
        //   'k' counts the num of adjacent tets which are 'sinfected'.
        j = k = 0;
        sstpivot1(*paryseg, neightet);
        spintet = neightet;
        while (1) {
          j++;
          if (!infected(spintet)) {
            neineitet =  spintet; // An outer tet. Remember it.
          } else {
            k++; // An in tet.
          }
          fnextself(spintet);
          if (spintet.tet == neightet.tet) break;
        }
        // assert(j > 0);
        if (k == 0) {
          // The segment is not connect to C(p) anymore. Remove it by
          //   Replacing it by the last entry of this list.
          s = cavetetseglist->objects - 1;
          checkseg = * (face *) fastlookup(cavetetseglist, s);
          *paryseg = checkseg;
          cavetetseglist->objects--;
          i--;
        } else if (k < j) {
          // The segment is on the boundary of C(p).
          sstbond1(*paryseg, neineitet);
        } else { // k == j
          // The segment is inside C(p).
          if (!ivf->splitbdflag) {//if (bowywat < 3) { // if (bowywat == 2) {
            checkseg = *paryseg;
            if (b->verbose > 3) {
              printf("        Queueing a missing seg (%d, %d)\n", 
	             pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
            }
            sinfect(checkseg); // Flag it as an interior segment.
            caveencseglist->newindex((void **) &paryseg);
            *paryseg = checkseg;
          } else {
            assert(0); // Not possible.
          }
        }
      } else { 
        // assert(smarktested(*paryseg));
        // Flag it as an interior segment. Do not queue it, since it will
        //   be deleted after the segment splitting.
        sinfect(*paryseg);
      }
    } // i
    if (b->verbose > 3) {
      printf("        %ld (%ld) cavity (interior) segments.\n", 
             cavetetseglist->objects, caveencseglist->objects);
    }
  } // if (checksubsegflag)

  if (checksubfaceflag) {
    for (i = 0; i < cavetetshlist->objects; i++) {
      parysh = (face *) fastlookup(cavetetshlist, i);
      // Operate on it if it is not inside the sub-cavity sC(p).
      if (!smarktested(*parysh)) {
        // Check if this subface is inside the cavity.
        k = 0;
        for (j = 0; j < 2; j++) {
          stpivot(*parysh, neightet);
          if (!infected(neightet)) {
            checksh = *parysh; // Remeber this side.
          } else {
            k++;
          }
          sesymself(*parysh);
        }
        if (k == 0) {
          // The subface is not connected to C(p). Remove it.
          s = cavetetshlist->objects - 1;
          checksh = * (face *) fastlookup(cavetetshlist, s);
          *parysh = checksh;
          cavetetshlist->objects--;
          i--;
        } else if (k == 1) {
          // This side is the outer boundary of C(p).
          *parysh = checksh;
        } else { // k == 2
          if (!ivf->splitbdflag) { //if (bowywat < 3) { // if (bowywat == 2) {
            checksh = *parysh;
            if (b->verbose > 3) {
              printf("        Queueing a missing subface (%d, %d, %d)\n", 
                     pointmark(sorg(checksh)), pointmark(sdest(checksh)),
                     pointmark(sapex(checksh)));
            }
            sinfect(checksh); // Flag it.
            caveencshlist->newindex((void **) &parysh);
            *parysh = checksh;
          } else {
            assert(0); // Not possible.
          }
        }
      } else {
        // assert(smarktested(*parysh));
        // Flag it as an interior subface. Do not queue it. It will be
        //   deleted after the facet point insertion.
        sinfect(*parysh);
      }
    } // i
    if (b->verbose > 3) {
      printf("        %ld (%ld) cavity (interior) subfaces.\n", 
             cavetetshlist->objects, caveencshlist->objects);
    }
  } // if (checksubfaceflag) {

  // Create new tetrahedra to fill the cavity.

  for (i = 0; i < cavebdrylist->objects; i++) {
    cavetet = (triface *) fastlookup(cavebdrylist, i);
    neightet = *cavetet;
    assert(!infected(neightet));
    unmarktest(neightet); // Unmark it.
    // Get the oldtet (inside the cavity).
    fsym(neightet, oldtet);
    if (apex(neightet) != dummypoint) {
      // Create a new tet in the cavity (see Fig. bowyerwatson 1 or 3).
      maketetrahedron(&newtet);
      setorg(newtet, dest(neightet));
      setdest(newtet, org(neightet));
      setapex(newtet, apex(neightet));
      setoppo(newtet, insertpt);
    } else {
      // Create a new hull tet (see Fig. bowyerwatson 2).
      hullsize++;
      maketetrahedron(&newtet);
      setorg(newtet, org(neightet));
      setdest(newtet, dest(neightet));
      setapex(newtet, insertpt);
      setoppo(newtet, dummypoint); // It must opposite to face 3.
      // Adjust back to the cavity bounday face.
      esymself(newtet);
    }
    // The new tet inherits attribtes from the old tet.
    for (j = 0; j < numelemattrib; j++) {
      attrib = elemattribute(oldtet.tet, j);
      setelemattribute(newtet.tet, j, attrib);
    }
    if (b->varvolume) {
      volume = volumebound(oldtet.tet);
      setvolumebound(newtet.tet, volume);
    }
    // Connect newtet <==> neightet, this also disconnect the old bond.
    bond(newtet, neightet);
    // oldtet still connects to neightet.
    *cavetet = oldtet; // *cavetet = newtet;
  } // i

  // Set a handle for speeding point location.
  recenttet = newtet;
  //setpoint2tet(insertpt, encode(newtet));
  setpoint2tet(insertpt, (tetrahedron) (newtet.tet));

  if (ivf->lawson > 1) { // if (lawson == 2 || lawson == 3) {
    // Re-use this list to save new interior cavity faces.
    cavetetlist->restart();
  }

  // Connect adjacent new tetrahedra together.
  for (i = 0; i < cavebdrylist->objects; i++) {
    cavetet = (triface *) fastlookup(cavebdrylist, i);
    // cavtet is an oldtet, get the newtet at this face.
    oldtet = *cavetet;
    fsym(oldtet, neightet);
    fsym(neightet, newtet);
    // Comment: oldtet and newtet must be at the same directed edge.
    // Connect the three other faces of this newtet.
    for (j = 0; j < 3; j++) {
      esym(newtet, neightet); // Go to the face.
      if (neightet.tet[neightet.ver & 3] == NULL) {
        // Find the adjacent face of this newtet.
        spintet = oldtet;
        while (1) {
          fnextself(spintet);
          if (!infected(spintet)) break;
        }
        fsym(spintet, newneitet);
        esymself(newneitet);
        assert(newneitet.tet[newneitet.ver & 3] == NULL); // FOR DEBUG
        bond(neightet, newneitet);
        if (ivf->lawson > 1) {
          // We are updateing a CDT. Queue the internal face.
          //   See also fig/dump-cavity-case13, -case21.
          cavetetlist->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
      //setpoint2tet(org(newtet), encode(newtet));
      setpoint2tet(org(newtet), (tetrahedron) (newtet.tet));
      enextself(newtet);
      enextself(oldtet);
    }
    *cavetet = newtet; // Save the new tet.
  } // i

  if (checksubfaceflag) {
    // Connect subfaces on the boundary of the cavity to the new tets.
    for (i = 0; i < cavetetshlist->objects; i++) {
      parysh = (face *) fastlookup(cavetetshlist, i);
      // Connect it if it is not a missing subface.
      if (!sinfected(*parysh)) {
        stpivot(*parysh, neightet);
        fsym(neightet, spintet);
        sesymself(*parysh);
        tsbond(spintet, *parysh);
      }
    }
  }

  if (checksubsegflag) {
    // Connect segments on the boundary of the cavity to the new tets.
    for (i = 0; i < cavetetseglist->objects; i++) {
      paryseg = (face *) fastlookup(cavetetseglist, i);
      // Connect it if it is not a missing segment.
      if (!sinfected(*paryseg)) {
        sstpivot1(*paryseg, neightet);
        spintet = neightet;
        while (1) {
          tssbond1(spintet, *paryseg);
          fnextself(spintet);
          if (spintet.tet == neightet.tet) break;
        }
      }
    }
  }

  if (splitsh != NULL) {
    // Split a subface or a segment.
    sinsertvertex(insertpt, splitsh, splitseg, ivf->sloc, ivf->sbowywat);
  }

  if (checksubfaceflag) {
    if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
      // Recover new subfaces in C(p).
      for (i = 0; i < caveshbdlist->objects; i++) {
        // Get an old subface at edge [a, b].
        parysh = (face *) fastlookup(caveshbdlist, i);
        spivot(*parysh, checksh); // The new subface [a, b, p].
        // Do not recover a deleted new face (degenerated).
        if (checksh.sh[3] != NULL) {
          // Note that the old subface still connects to adjacent old tets 
          //   of C(p), which still connect to the tets outside C(p).
          stpivot(*parysh, neightet);
          assert(infected(neightet));
          // Find the adjacent tet containing the edge [a,b] outside C(p).
          spintet = neightet;
          while (1) {
            fnextself(spintet);
            if (!infected(spintet)) break;
            assert(spintet.tet != neightet.tet);
          }
          // The adjacent tet connects to a new tet in C(p).
          fsym(spintet, neightet);
          assert(!infected(neightet));
          // Find the tet containing the face [a, b, p].
          spintet = neightet;
          while (1) {
            fnextself(spintet);
            if (apex(spintet) == insertpt) break;
            assert(spintet.tet != neightet.tet);
          }
          // Adjust the edge direction in spintet and checksh.
          if (sorg(checksh) != org(spintet)) {
            sesymself(checksh);
            assert(sorg(checksh) == org(spintet));
          }
          assert(sdest(checksh) == dest(spintet));
          // Connect the subface to two adjacent tets.
          tsbond(spintet, checksh);
          fsymself(spintet);
          sesymself(checksh);
          tsbond(spintet, checksh);
        } // if (checksh.sh[3] != NULL)
      }
      // There should be no missing interior subfaces in C(p).
      assert(caveencshlist->objects == 0l);
    } else { 
      // bowywat = 1 or bowywat = 2.
      // The Boundary reocvery phase.
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
      // Put all interior subfaces into stack for recovery.
      for (i = 0; i < caveencshlist->objects; i++) {
        parysh = (face *) fastlookup(caveencshlist, i);
        assert(sinfected(*parysh));
        // Some subfaces inside C(p) might be split in sinsertvertex().
        //   Only queue those faces which are not split.
        if (!smarktested(*parysh)) {
          if (b->verbose > 3) {
            printf("        Queue a missing subface (%d, %d, %d) x%lx.\n",
                   pointmark(sorg(*parysh)), pointmark(sdest(*parysh)),
                   pointmark(sapex(*parysh)), (uintptr_t) parysh->sh);
          }
          checksh = *parysh;
          suninfect(checksh);
          stdissolve(checksh); // Detach connections to old tets.
          subfacstack->newindex((void **) &parysh);
          *parysh = checksh;
        }
      }
    }
  } // if (checksubfaceflag)

  if (checksubsegflag) {
    if (ivf->splitbdflag) { //if (bowywat > 2) { // if (bowywat == 3) {
      if (splitseg != NULL) {
        // Recover the two new subsegments in C(p).
        for (i = 0; i < cavesegshlist->objects; i++) {
          paryseg = (face *) fastlookup(cavesegshlist, i);
          // Insert this subsegment into C(p).
          checkseg = *paryseg;
          // Get the adjacent new subface.
          checkseg.shver = 0;
          spivot(checkseg, checksh);
          if (checksh.sh != NULL) {
            // Get the adjacent new tetrahedron.
            stpivot(checksh, neightet);
          } else {
            // It's a dangling segment.
            pa = sorg(checkseg);
            pb = sdest(checkseg);
            point2tetorg(pa, neightet);
            finddirection(&neightet, pb);
            assert(dest(neightet) == pb);
          }
          assert(!infected(neightet));
          sstbond1(checkseg, neightet);
          spintet = neightet;
          while (1) {
            tssbond1(spintet, checkseg);
            fnextself(spintet);
            if (spintet.tet == neightet.tet) break;
          }
        }
      } // if (splitseg != NULL)
      // There should be no interior segment in C(p).
      assert(caveencseglist->objects == 0l);
    } else {
      // bowywat == 1 or bowywat == 2;
      // The Boundary Recovery Phase.  
      // Queue missing segments in C(p) for recovery.
      if (splitseg != NULL) {
        // Queue two new subsegments in C(p) for recovery.
        for (i = 0; i < cavesegshlist->objects; i++) {
          paryseg = (face *) fastlookup(cavesegshlist, i);
          if (b->verbose > 3) {
            printf("        Queue new subseg (%d, %d)\n", 
                   pointmark(sorg(*paryseg)), pointmark(sdest(*paryseg)));
          }
          checkseg = *paryseg;
          //sstdissolve1(checkseg); // It has not been connected yet.
          s = randomnation(subsegstack->objects + 1);
          subsegstack->newindex((void **) &paryseg);
          *paryseg = * (face *) fastlookup(subsegstack, s); 
          paryseg = (face *) fastlookup(subsegstack, s);
          *paryseg = checkseg;
        }
      } // if (splitseg != NULL)
      for (i = 0; i < caveencseglist->objects; i++) {
        paryseg = (face *) fastlookup(caveencseglist, i);
        assert(sinfected(*paryseg));
        if (!smarktested(*paryseg)) { // It may be split.
          if (b->verbose > 3) {
            printf("        Queue a missing segment (%d, %d).\n",
                   pointmark(sorg(*paryseg)), pointmark(sdest(*paryseg)));
          }
          checkseg = *paryseg;
          suninfect(checkseg);
          sstdissolve1(checkseg); // Detach connections to old tets.
          s = randomnation(subsegstack->objects + 1);
          subsegstack->newindex((void **) &paryseg);
          *paryseg = * (face *) fastlookup(subsegstack, s); 
          paryseg = (face *) fastlookup(subsegstack, s);
          *paryseg = checkseg;
        }
      }
    }
  } // if (checksubsegflag)

  if (b->plc || b->weighted) {
    // Some vertices may be completed inside the cavity. They must be
    //   detected and added to recovering list.
    if (b->plc) {
      tetcount = subvertstack->objects; // Re-use tetcount;
    }
    // Since every "live" vertex must contain a pointer to a non-dead
    //   tetrahedron, we can check for each vertex this pointer.
    for (i = 0; i < cavetetvertlist->objects; i++) {
      pts = (point *) fastlookup(cavetetvertlist, i);
      decode(point2tet(*pts), *searchtet);
      assert(searchtet->tet != NULL); // No tet has been deleted yet.
      if (infected(*searchtet)) {
        if (b->weighted) {
          if (b->verbose > 1) {
            printf("    Point #%d is non-regular after the insertion of #%d.\n",
                   pointmark(*pts), pointmark(insertpt));
	  }
          setpointtype(*pts, NREGULARVERTEX);
          nonregularcount++;
        } else {
          if (b->verbose > 3) {
            printf("        Queue a dangling vertex %d.\n", pointmark(*pts));
          }
          subvertstack->newindex((void **) &parypt);
          *parypt = *pts;
        }
      }
    }
    if (b->plc) {
      if (subvertstack->objects > tetcount) {
        // There are missing vertices after inserting the new point.
        printf("DBG: Insert %d. Found %ld interior vertices.\n",
               pointmark(insertpt), subvertstack->objects);
        assert(0); // NEED TO DEBUG.
      }
    }
  }

  if (ivf->chkencflag & 1) {
    // Queue all segment outside C(p).
    for (i = 0; i < cavetetseglist->objects; i++) {
      paryseg = (face *) fastlookup(cavetetseglist, i);
      // Skip if it is the split segment.
      if (!sinfected(*paryseg)) {
        // Skip it if it has already queued.
        if (!smarktest2ed(*paryseg)) {
          bface = (badface *) badsubsegs->alloc();
          bface->ss = *paryseg;
          smarktest2(bface->ss); // Only queue it once.
          bface->forg = sorg(*paryseg); // An alive badface.
        }
      }
    }
    if (splitseg != NULL) {
      // Queue the two new subsegments inside C(p).
      for (i = 0; i < cavesegshlist->objects; i++) {
        paryseg = (face *) fastlookup(cavesegshlist, i);
        bface = (badface *) badsubsegs->alloc();
        bface->ss = *paryseg;
        smarktest2(bface->ss); // Only queue it once.
        bface->forg = sorg(*paryseg); // An alive badface.
      }
    }
  } // if (chkencflag & 1)

  if (ivf->chkencflag & 2) {
    // Queue all subfaces outside C(p).
    for (i = 0; i < cavetetshlist->objects; i++) {
      parysh = (face *) fastlookup(cavetetshlist, i);
      // Skip if it is a split subface.
      if (!sinfected(*parysh)) {
        // Skip it if it has already queued.
        if (!smarktest2ed(*parysh)) {
          bface = (badface *) badsubfacs->alloc();
          bface->ss = *parysh;
          smarktest2(bface->ss); // Only queue it once.
          bface->forg = sorg(*parysh); // An alive badface.
          //bface->fdest = sdest(*parysh);
          //bface->fapex = sapex(*parysh);
        }
      }
    }
    // Queue all new subfaces inside C(p).
    for (i = 0; i < caveshbdlist->objects; i++) {
      // Get an old subface at edge [a, b].
      parysh = (face *) fastlookup(caveshbdlist, i);
      spivot(*parysh, checksh); // checksh is a new subface [a, b, p].
      // Do not recover a deleted new face (degenerated).
      if (checksh.sh[3] != NULL) {
        //assert(!smarktest2ed(checksh));
        bface = (badface *) badsubfacs->alloc();
        bface->ss = checksh;
        smarktest2(bface->ss); // Only queue it once.
        bface->forg = sorg(checksh); // An alive badface.
      }
    }
  } // if (chkencflag & 2)

  if (ivf->chkencflag & 4) {
    // Queue all new tetrahedra in C(p).
    for (i = 0; i < cavebdrylist->objects; i++) {
      cavetet = (triface *) fastlookup(cavebdrylist, i);
      //assert(!marktest2ed(*cavetet));
      bface = (badface *) badtetrahedrons->alloc();
      bface->tt = *cavetet;
      marktest2(bface->tt);
      bface->forg = org(*cavetet);
    }
  }

  // C(p) is re-meshed successfully. 

  // Deleted the old tets in C(p).
  for (i = 0; i < caveoldtetlist->objects; i++) {
    searchtet = (triface *) fastlookup(caveoldtetlist, i);
    tetrahedrondealloc(searchtet->tet);
  }

  if (splitsh != NULL) {
    // Delete the old subfaces in sC(p).
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      if (checksubfaceflag) {//if (bowywat == 2) {
        // It is possible that this subface still connects to adjacent
        //   tets which are not in C(p). If so, clear connections in the
        //   adjacent tets at this subface.
        stpivot(*parysh, neightet);
        if (neightet.tet != NULL) {
          if (neightet.tet[4] != NULL) {
            // Found an adjacent tet. It must be not in C(p).
            assert(!infected(neightet));
            tsdissolve(neightet);
            fsymself(neightet);
            assert(!infected(neightet));
            tsdissolve(neightet);
          }
        }
      }
      shellfacedealloc(subfaces, parysh->sh);
    }
    if (splitseg != NULL) {
      // Delete the old segment in sC(p).
      shellfacedealloc(subsegs, splitseg->sh);
    }
  }

  if (ivf->lawson) {
    for (i = 0; i < cavebdrylist->objects; i++) {
      searchtet = (triface *) fastlookup(cavebdrylist, i);
      //flippush(flipstack, searchtet, insertpt);
      flippush(flipstack, searchtet);
    }
    if (ivf->lawson > 1) {
      for (i = 0; i < cavetetlist->objects; i++) {
        searchtet = (triface *) fastlookup(cavetetlist, i);
        //flippush(flipstack, searchtet, oppo(*searchtet));
        flippush(flipstack, searchtet);
      }
    }
  }

  // The vertex should already have a type.
  assert(pointtype(insertpt) != UNUSEDVERTEX);



  // Clean the working lists.

  caveoldtetlist->restart();
  cavebdrylist->restart();
  cavetetlist->restart();

  if (checksubsegflag) {
    cavetetseglist->restart();
    caveencseglist->restart();
  }

  if (checksubfaceflag) {
    cavetetshlist->restart();
    caveencshlist->restart();
  }
  
  if (b->plc || b->weighted) {
    cavetetvertlist->restart();
  }
  
  if (splitsh != NULL) {
    caveshlist->restart();
    caveshbdlist->restart();
    cavesegshlist->restart();
  }

  return (int) loc;
}

////                                                                       ////
////                                                                       ////
//// flip_cxx /////////////////////////////////////////////////////////////////

