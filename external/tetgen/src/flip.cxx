#include "../tetgen.h"
//// flip_cxx /////////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flip23()    Perform a 2-to-3 flip (face-to-edge flip).                    //
//                                                                           //
// 'fliptets' is an array of three tets (handles), where the [0] and [1] are  //
// [a,b,c,d] and [b,a,c,e].  The three new tets: [e,d,a,b], [e,d,b,c], and   //
// [e,d,c,a] are returned in [0], [1], and [2] of 'fliptets'.  As a result,  //
// The face [a,b,c] is removed, and the edge [d,e] is created.               //
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
// If 'fc->enqflag' is set, convex hull faces will be queued for flipping.   //
// In particular, if 'fc->enqflag' is 1, it is called by incrementalflip()   //
// after the insertion of a new point.  It is assumed that 'd' is the new    //
// point. IN this case, only link faces of 'd' are queued.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip23(triface* fliptets, int hullflag, flipconstraints *fc)
{
  triface topcastets[3], botcastets[3];
  triface newface, casface;
  point pa, pb, pc, pd, pe;
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

  pa =  org(fliptets[0]);
  pb = dest(fliptets[0]);
  pc = apex(fliptets[0]);
  pd = oppo(fliptets[0]);
  pe = oppo(fliptets[1]);

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

  if (fc->remove_ndelaunay_edge) { // calc_tetprism_vol
    REAL volneg[2], volpos[3], vol_diff;
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
    fc->tetprism_vol_sum  += vol_diff; // Update the total sum.
  }

  // Bond three new tets together.
  for (i = 0; i < 3; i++) {
    esym(fliptets[i], newface);
    bond(newface, fliptets[(i + 1) % 3]);
  }
  // Bond to top outer boundary faces (at [a,b,c,d]).
  for (i = 0; i < 3; i++) {
    eorgoppo(fliptets[i], newface); // At edges [b,a], [c,b], [a,c].
    bond(newface, topcastets[i]);
  }
  // Bond bottom outer boundary faces (at [b,a,c,e]).
  for (i = 0; i < 3; i++) {
    edestoppo(fliptets[i], newface); // At edges [a,b], [b,c], [c,a].
    bond(newface, botcastets[i]);
  }

  if (checksubsegflag) {
    // Bond subsegments if there are.
    // Each new tet has 5 edges to be checked (except the edge [e,d]). 
    face checkseg;
    // The middle three: [a,b], [b,c], [c,a].
    for (i = 0; i < 3; i++) {      
      if (issubseg(topcastets[i])) {
        tsspivot1(topcastets[i], checkseg);
        eorgoppo(fliptets[i], newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
    }
    // The top three: [d,a], [d,b], [d,c]. Two tets per edge.
    for (i = 0; i < 3; i++) {
      eprev(topcastets[i], casface);      
      if (issubseg(casface)) {
        tsspivot1(casface, checkseg);
        enext(fliptets[i], newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        esym(fliptets[(i + 2) % 3], newface);
        eprevself(newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
    }
    // The bot three: [a,e], [b,e], [c,e]. Two tets per edge.
    for (i = 0; i < 3; i++) {
      enext(botcastets[i], casface);
      if (issubseg(casface)) {
        tsspivot1(casface, checkseg);
        eprev(fliptets[i], newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        esym(fliptets[(i + 2) % 3], newface);
        enextself(newface);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
    }
  } // if (checksubsegflag)

  if (checksubfaceflag) {
    // Bond 6 subfaces if there are.
    face checksh;
    for (i = 0; i < 3; i++) {      
      if (issubface(topcastets[i])) {
        tspivot(topcastets[i], checksh);
        eorgoppo(fliptets[i], newface);
        sesymself(checksh);
        tsbond(newface, checksh);
        if (fc->chkencflag & 2) {
          enqueuesubface(badsubfacs, &checksh);
        }
      }
    }
    for (i = 0; i < 3; i++) {
      if (issubface(botcastets[i])) {
        tspivot(botcastets[i], checksh);
        edestoppo(fliptets[i], newface);
        sesymself(checksh);
        tsbond(newface, checksh);
        if (fc->chkencflag & 2) {
          enqueuesubface(badsubfacs, &checksh);
        }
      }
    }
  } // if (checksubfaceflag)

  if (fc->chkencflag & 4) {
    // Put three new tets into check list.
    for (i = 0; i < 3; i++) {
      enqueuetetrahedron(&(fliptets[i]));
    }
  }

  // Update the point-to-tet map.
  setpoint2tet(pa, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pb, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pc, (tetrahedron) fliptets[1].tet);
  setpoint2tet(pd, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pe, (tetrahedron) fliptets[0].tet);

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

  if (fc->enqflag > 0) {
    // Queue faces which may be locally non-Delaunay.  
    for (i = 0; i < 3; i++) {
      eprevesym(fliptets[i], newface);
      flippush(flipstack, &newface);
    }
    if (fc->enqflag > 1) {
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
// 'fliptets' is an array of three tets (handles),  which are [e,d,a,b],     //
// [e,d,b,c], and [e,d,c,a].  The two new tets: [a,b,c,d] and [b,a,c,e] are  //
// returned in [0] and [1] of 'fliptets'.  As a result, the edge [e,d] is    //
// replaced by the face [a,b,c].                                             //
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
// If 'fc->enqflag' is set, convex hull faces will be queued for flipping.   //
// In particular, if 'fc->enqflag' is 1, it is called by incrementalflip()   //
// after the insertion of a new point.  It is assumed that 'a' is the new    //
// point. In this case, only link faces of 'a' are queued.                   //
//                                                                           //
// If 'checksubfaceflag' is on (global variable), and assume [e,d] is not a  //
// segment. There may be two (interior) subfaces sharing at [e,d], which are //
// [e,d,p] and [e,d,q], where the pair (p,q) may be either (a,b), or (b,c),  //
// or (c,a)  In such case, a 2-to-2 flip is performed on these two subfaces  //
// and two new subfaces [p,q,e] and [p,q,d] are created. They are inserted   //
// back into the tetrahedralization.                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip32(triface* fliptets, int hullflag, flipconstraints *fc)
{
  triface topcastets[3], botcastets[3];
  triface newface, casface;
  face flipshs[3]; 
  face checkseg; 
  point pa, pb, pc, pd, pe;
  REAL attrib, volume;
  int dummyflag = 0;  // Rangle = {-1, 0, 1, 2}
  int spivot = -1, scount = 0; // for flip22()
  int t1ver; 
  int i, j;

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

  flip32count++;

  // Get the outer boundary faces.
  for (i = 0; i < 3; i++) {
    eorgoppo(fliptets[i], casface);
    fsym(casface, topcastets[i]);
  }
  for (i = 0; i < 3; i++) {
    edestoppo(fliptets[i], casface);
    fsym(casface, botcastets[i]);
  }

  if (checksubfaceflag) {
    // Check if there are interior subfaces at the edge [e,d].
    for (i = 0; i < 3; i++) {
      tspivot(fliptets[i], flipshs[i]);
      if (flipshs[i].sh != NULL) {
        // Found an interior subface.
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
  if (checksubfaceflag) {
    if (scount > 0) {
      // The element attributes and volume constraint must be set correctly.
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
      // The hullsize does not change.
    }
  } else {
    setvertices(fliptets[0], pa, pb, pc, pd);
    setvertices(fliptets[1], pb, pa, pc, pe);
  }

  if (fc->remove_ndelaunay_edge) { // calc_tetprism_vol
    REAL volneg[3], volpos[2], vol_diff;
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
    fc->tetprism_vol_sum  += vol_diff; // Update the total sum.
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
    // Bond 9 segments to new (flipped) tets.
    for (i = 0; i < 3; i++) { // edges a->b, b->c, c->a.      
      if (issubseg(topcastets[i])) {
        tsspivot1(topcastets[i], checkseg);
        tssbond1(fliptets[0], checkseg);
        sstbond1(checkseg, fliptets[0]);
        tssbond1(fliptets[1], checkseg);
        sstbond1(checkseg, fliptets[1]);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
      enextself(fliptets[0]);
      eprevself(fliptets[1]);
    }
    // The three top edges.
    for (i = 0; i < 3; i++) { // edges b->d, c->d, a->d.
      esym(fliptets[0], newface);
      eprevself(newface); 
      enext(topcastets[i], casface);      
      if (issubseg(casface)) {
        tsspivot1(casface, checkseg);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
      enextself(fliptets[0]);
    }
    // The three bot edges.
    for (i = 0; i < 3; i++) { // edges b<-e, c<-e, a<-e.
      esym(fliptets[1], newface);
      enextself(newface); 
      eprev(botcastets[i], casface);
      if (issubseg(casface)) {
        tsspivot1(casface, checkseg);
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
      eprevself(fliptets[1]);
    }
  } // if (checksubsegflag)

  if (checksubfaceflag) {
    face checksh;
    // Bond the top three casing subfaces.
    for (i = 0; i < 3; i++) { // At edges [b,a], [c,b], [a,c]
      if (issubface(topcastets[i])) {
        tspivot(topcastets[i], checksh);
        esym(fliptets[0], newface); 
        sesymself(checksh);
        tsbond(newface, checksh);
        if (fc->chkencflag & 2) {
          enqueuesubface(badsubfacs, &checksh);
        }
      }
      enextself(fliptets[0]);
    }
    // Bond the bottom three casing subfaces.
    for (i = 0; i < 3; i++) { // At edges [a,b], [b,c], [c,a]
      if (issubface(botcastets[i])) {
        tspivot(botcastets[i], checksh);
        esym(fliptets[1], newface); 
        sesymself(checksh);
        tsbond(newface, checksh);
        if (fc->chkencflag & 2) {
          enqueuesubface(badsubfacs, &checksh);
        }
      }
      eprevself(fliptets[1]);
    }

    if (scount > 0) {
      face flipfaces[2];
      // Perform a 2-to-2 flip in subfaces.
      flipfaces[0] = flipshs[(spivot + 1) % 3];
      flipfaces[1] = flipshs[(spivot + 2) % 3];
      sesymself(flipfaces[1]);
      flip22(flipfaces, 0, fc->chkencflag);
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
        // An invalid 2-to-2 flip. Report a bug.
        terminatetetgen(this, 2); 
      }
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
        // An invalid 2-to-2 flip. Report a bug.
        terminatetetgen(this, 2);  
      }
    } // if (scount > 0)
  } // if (checksubfaceflag)

  if (fc->chkencflag & 4) {
    // Put two new tets into check list.
    for (i = 0; i < 2; i++) {
      enqueuetetrahedron(&(fliptets[i]));
    }
  }

  setpoint2tet(pa, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pb, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pc, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pd, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pe, (tetrahedron) fliptets[1].tet);

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
  
  if (fc->enqflag > 0) {
    // Queue faces which may be locally non-Delaunay.
    // pa = org(fliptets[0]); // 'a' may be a new vertex.
    enextesym(fliptets[0], newface);
    flippush(flipstack, &newface);
    eprevesym(fliptets[1], newface);
    flippush(flipstack, &newface);
    if (fc->enqflag > 1) {
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
// If 'hullflag' is set (> 0), one of the five vertices may be 'dummypoint'. //
// The 'hullsize' may be changed.  Note that p may be dummypoint.  In this   //
// case, four hull tets are replaced by one real tet.                        //
//                                                                           //
// If 'checksubface' flag is set (>0), it is possible that there are three   //
// interior subfaces connecting at p.  If so, a 3-to-1 flip is performed to  //
// to remove p from the surface triangulation.                               //
//                                                                           //
// If it is called by the routine incrementalflip(), we assume that d is the //
// newly inserted vertex.                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::flip41(triface* fliptets, int hullflag, flipconstraints *fc)
{
  triface topcastets[3], botcastet;
  triface newface, neightet;
  face flipshs[4];
  point pa, pb, pc, pd, pp;
  int dummyflag = 0; // in {0, 1, 2, 3, 4}
  int spivot = -1, scount = 0;
  int t1ver; 
  int i;

  pa =  org(fliptets[3]);
  pb = dest(fliptets[3]);
  pc = apex(fliptets[3]);
  pd = dest(fliptets[0]);
  pp =  org(fliptets[0]); // The removing vertex.

  flip41count++;

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

  if (pp != dummypoint) {
    // Mark the point pp as unused.
    setpointtype(pp, UNUSEDVERTEX);
    unuverts++;
  }

  // Create the new tet [a,b,c,d].
  if (hullflag > 0) {
    // One of the five vertices may be 'dummypoint'.
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
      if (pp == dummypoint) {
        dummyflag = -1;
      } else {
        dummyflag = 0;
      }
    }
    if (dummyflag > 0) {
      // We deleted 3 hull tets, and create 1 hull tet.
      hullsize -= 2;
    } else if (dummyflag < 0) {
      // We deleted 4 hull tets.
      hullsize -= 4;
      // meshedges does not change.
    }
  } else {
    setvertices(fliptets[0], pa, pb, pc, pd);
  }

  if (fc->remove_ndelaunay_edge) { // calc_tetprism_vol
    REAL volneg[4], volpos[1], vol_diff;
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
    } else if (dummyflag < 0) {
      volneg[0] = 0.;
      volneg[1] = 0.;
      volneg[2] = 0.;
      volneg[3] = 0.;
      volpos[0] = tetprismvol(pa, pb, pc, pd);
    } else {
      volneg[0] = tetprismvol(pp, pd, pa, pb);
      volneg[1] = tetprismvol(pp, pd, pb, pc);
      volneg[2] = tetprismvol(pp, pd, pc, pa);
      volneg[3] = tetprismvol(pa, pb, pc, pp);
      volpos[0] = tetprismvol(pa, pb, pc, pd);
    }
    vol_diff = volpos[0] - volneg[0] - volneg[1] - volneg[2] - volneg[3];
    fc->tetprism_vol_sum  += vol_diff; // Update the total sum.
  }

  // Bond the new tet to adjacent tets.
  for (i = 0; i < 3; i++) {
    esym(fliptets[0], newface); // At faces [b,a,d], [c,b,d], [a,c,d].
    bond(newface, topcastets[i]);
    enextself(fliptets[0]);
  }
  bond(fliptets[0], botcastet);

  if (checksubsegflag) {
    face checkseg;
    // Bond 6 segments (at edges of [a,b,c,d]) if there there are.
    for (i = 0; i < 3; i++) {
      eprev(topcastets[i], newface); // At edges [d,a],[d,b],[d,c].
      if (issubseg(newface)) {
        tsspivot1(newface, checkseg);
        esym(fliptets[0], newface);
        enextself(newface); // At edges [a,d], [b,d], [c,d].
        tssbond1(newface, checkseg);
        sstbond1(checkseg, newface);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
      enextself(fliptets[0]);
    }
    for (i = 0; i < 3; i++) {
      if (issubseg(topcastets[i])) {
        tsspivot1(topcastets[i], checkseg); // At edges [a,b],[b,c],[c,a].
        tssbond1(fliptets[0], checkseg);
        sstbond1(checkseg, fliptets[0]);
        if (fc->chkencflag & 1) {
          enqueuesubface(badsubsegs, &checkseg);
        }
      }
      enextself(fliptets[0]);
    }
  }

  if (checksubfaceflag) {
    face checksh;
    // Bond 4 subfaces (at faces of [a,b,c,d]) if there are.
    for (i = 0; i < 3; i++) {
      if (issubface(topcastets[i])) {
        tspivot(topcastets[i], checksh); // At faces [a,b,d],[b,c,d],[c,a,d]
        esym(fliptets[0], newface); // At faces [b,a,d],[c,b,d],[a,c,d]
        sesymself(checksh);
        tsbond(newface, checksh);
        if (fc->chkencflag & 2) {
          enqueuesubface(badsubfacs, &checksh);
        }
      }
      enextself(fliptets[0]);
    }
    if (issubface(botcastet)) {
      tspivot(botcastet, checksh); // At face [b,a,c]
      sesymself(checksh);
      tsbond(fliptets[0], checksh);
      if (fc->chkencflag & 2) {
        enqueuesubface(badsubfacs, &checksh);
      }
    }

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

  if (fc->chkencflag & 4) {
    enqueuetetrahedron(&(fliptets[0]));
  }

  // Update the point-to-tet map.
  setpoint2tet(pa, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pb, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pc, (tetrahedron) fliptets[0].tet);
  setpoint2tet(pd, (tetrahedron) fliptets[0].tet);

  if (fc->enqflag > 0) {
    // Queue faces which may be locally non-Delaunay.
    flippush(flipstack, &(fliptets[0])); // [a,b,c] (opposite to new point).
    if (fc->enqflag > 1) {
      for (i = 0; i < 3; i++) {
        esym(fliptets[0], newface);
        flippush(flipstack, &newface);
        enextself(fliptets[0]);
      }
    }
  }

  recenttet = fliptets[0];
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipnm()    Flip an edge through a sequence of elementary flips.          //
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
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::flipnm(triface* abtets, int n, int level, int abedgepivot,
                       flipconstraints* fc)
{
  triface fliptets[3], spintet, flipedge;
  triface *tmpabtets, *parytet;
  point pa, pb, pc, pd, pe, pf;
  REAL ori;
  int hullflag, hulledgeflag;
  int reducflag, rejflag;
  int reflexlinkedgecount;
  int edgepivot;
  int n1, nn;
  int t1ver;
  int i, j;

  pa = org(abtets[0]);
  pb = dest(abtets[0]);

  if (n > 3) {
    // Try to reduce the size of the Star(ab) by flipping a face in it. 
    reflexlinkedgecount = 0;

    for (i = 0; i < n; i++) {
      // Let the face of 'abtets[i]' be [a,b,c].
      if (checksubfaceflag) {
        if (issubface(abtets[i])) {
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
        continue; // [a,b,c] is a hull face.
      }


      // Decide whether [a,b,c] is flippable or not.
      reducflag = 0; 

      hullflag = (pc == dummypoint); // pc may be dummypoint.
      hulledgeflag = 0;
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
                // The "flat" tet can be removed immediately by a 3-to-2 flip.
                reducflag = 1;
                // Check if [e,d] is a hull edge.
                pf = apex(abtets[(i + 2) % n]);
                hulledgeflag = (pf == dummypoint);
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
              ori = 0; // Signal as a 4-to-4 flip (like a co-planar case).
              hulledgeflag = 1; // [e,d] is a hull edge.
            }
          }
        }
      } // if (hullflag)

      if (reducflag) {
        if (nonconvex && hulledgeflag) {
          // We will create a hull edge [e,d]. Make sure it does not exist.
          if (getedge(pe, pd, &spintet)) {
            // The 2-to-3 flip is not a topological valid flip. 
            reducflag = 0;
          }
        }
      }

      if (reducflag) {
        // [a,b,c] could be removed by a 2-to-3 flip.
        rejflag = 0;
        if (fc->checkflipeligibility) {
          // Check if the flip can be performed.
          rejflag = checkflipeligibility(1, pa, pb, pc, pd, pe, level,
                                         abedgepivot, fc);
        }
        if (!rejflag) {
          // Do flip: [a,b,c] => [e,d].
          fliptets[0] = abtets[i];
          fsym(fliptets[0], fliptets[1]); // abtets[i-1].
          flip23(fliptets, hullflag, fc);

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
          edestoppoself(fliptets[0]); // [a,b,e,d]
          // Increase the counter of this new tet (it is in Star(ab)).
          increaseelemcounter(fliptets[0]); 
          abtets[(i - 1 + n) % n] = fliptets[0];
          for (j = i; j < n - 1; j++) {
            abtets[j] = abtets[j + 1];  // Upshift
          }
          // The last entry 'abtets[n-1]' is empty. It is used in two ways:
          //   (i) it remembers the vertex 'c' (in 'abtets[n-1].tet'), and
          //  (ii) it remembers the position [i] where this flip took place.
          // These informations let us to either undo this flip or recover
          //   the original edge link (for collecting new created tets).
          //abtets[n - 1] = fliptets[1]; // [e,d,b,c] is remembered.
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

          if (nn == 2) {
            // The edge has been flipped.
            return nn;
          } else { // if (nn > 2)
            // The edge is not flipped.
            if (fc->unflip || (ori == 0)) {
              // Undo the previous 2-to-3 flip, i.e., do a 3-to-2 flip to 
              //   transform [e,d] => [a,b,c].
              // 'ori == 0' means that the previous flip created a degenerated
              //   tet. It must be removed. 
              // Remember that 'abtets[i-1]' is [a,b,e,d]. We can use it to
              //   find another two tets [e,d,b,c] and [e,d,c,a].
              fliptets[0] = abtets[(i-1 + (n-1)) % (n-1)]; // [a,b,e,d]
              edestoppoself(fliptets[0]); // [e,d,a,b]
              fnext(fliptets[0], fliptets[1]); // [1] is [e,d,b,c]
              fnext(fliptets[1], fliptets[2]); // [2] is [e,d,c,a]
              assert(apex(fliptets[0]) == oppo(fliptets[2])); // SELF_CHECK
              // Restore the two original tets in Star(ab). 
              flip32(fliptets, hullflag, fc);
              // Marktest the two restored tets in Star(ab).
              for (j = 0; j < 2; j++) {
                increaseelemcounter(fliptets[j]);
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
            } // if (unflip || (ori == 0))
          } // if (nn > 2)

          if (!fc->unflip) {
            // The flips are not reversed. The current Star(ab) can not be
            //   further reduced. Return its current size (# of tets).
            return nn; 
          }
          // unflip is set. 
          // Continue the search for flips.
        }
      } // if (reducflag)
    } // i

    // The Star(ab) is not reduced. 
    if (reflexlinkedgecount > 0) {
      // There are reflex edges in the Link(ab).
      if (((b->fliplinklevel < 0) && (level < autofliplinklevel)) || 
          ((b->fliplinklevel >= 0) && (level < b->fliplinklevel))) {
        // Try to reduce the Star(ab) by flipping a reflex edge in Link(ab).
        for (i = 0; i < n; i++) {
          // Do not flip this face [a,b,c] if there are two Stars involved.
          if ((elemcounter(abtets[i]) > 1) ||
              (elemcounter(abtets[(i - 1 + n) % n]) > 1)) {
            continue;
          }
          pc = apex(abtets[i]);
          if (pc == dummypoint) {
            continue; // [a,b] is a hull edge.
          }
          pd = apex(abtets[(i + 1) % n]);
          pe = apex(abtets[(i - 1 + n) % n]);
          if ((pd == dummypoint) || (pe == dummypoint)) {
            continue; // [a,b,c] is a hull face.
          }


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
            if (issubseg(flipedge)) {
              if (fc->collectencsegflag) {
                face checkseg, *paryseg;
                tsspivot1(flipedge, checkseg);
                if (!sinfected(checkseg)) {
                  // Queue this segment in list.
                  sinfect(checkseg);                
                  caveencseglist->newindex((void **) &paryseg);
                  *paryseg = checkseg;
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
            j += (elemcounter(spintet)); 
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

          if ((b->flipstarsize > 0) && (n1 > b->flipstarsize)) {
            // The star size exceeds the given limit.
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
            assert(elemcounter(spintet) == 0); // It's a new tet.
            increaseelemcounter(spintet); // It is in Star(ab).
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

            if (nn == 2) { 
              // The edge has been flipped.
              return nn;
            } else { // if (nn > 2) {
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

            if (!fc->unflip) {
              // The flips are not reversed. The current Star(ab) can not be
              //   further reduced. Return its size (# of tets).
              return nn; 
            }
            // unflip is set. 
            // Continue the search for flips.
          } else {
            // The selected edge is not flipped.
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
      } // if (level...)
    } // if (reflexlinkedgecount > 0)
  } else {
    // Check if a 3-to-2 flip is possible.
    // Let the three apexes be c, d,and e. Hull tets may be involved. If so, 
    //   we rearrange them such that the vertex e is dummypoint. 
    hullflag = 0;

    if (apex(abtets[0]) == dummypoint) {
      pc = apex(abtets[1]);
      pd = apex(abtets[2]);
      pe = apex(abtets[0]);
      hullflag = 1;
    } else if (apex(abtets[1]) == dummypoint) {
      pc = apex(abtets[2]);
      pd = apex(abtets[0]);
      pe = apex(abtets[1]);
      hullflag = 2;
    } else {
      pc = apex(abtets[0]);
      pd = apex(abtets[1]);
      pe = apex(abtets[2]);
      hullflag = (pe == dummypoint) ? 3 : 0;
    }

    reducflag = 0;
    rejflag = 0;


    if (hullflag == 0) {
      // Make sure that no inverted tet will be created, i.e. the new tets
      //   [d,c,e,a] and [c,d,e,b] must be valid tets. 
      ori = orient3d(pd, pc, pe, pa);
      if (ori < 0) {
        ori = orient3d(pc, pd, pe, pb);
        if (ori < 0) {
          reducflag = 1;
        }
      }
    } else {
      // [a,b] is a hull edge.
      //   Note: This can happen when it is in the middle of a 4-to-4 flip.
      //   Note: [a,b] may even be a non-convex hull edge.
      if (!nonconvex) {
        //  The mesh is convex, only do flip if it is a coplanar hull edge.
        ori = orient3d(pa, pb, pc, pd); 
        if (ori == 0) {
          reducflag = 1;
        }
      } else { // nonconvex
        reducflag = 1;
      }
      if (reducflag == 1) {
        // [a,b], [a,b,c] and [a,b,d] are on the convex hull.
        // Make sure that no inverted tet will be created.
        point searchpt = NULL, chkpt;
        REAL bigvol = 0.0, ori1, ori2;
        // Search an interior vertex which is an apex of edge [c,d].
        //   In principle, it can be arbitrary interior vertex.  To avoid
        //   numerical issue, we choose the vertex which belongs to a tet
        //   't' at edge [c,d] and 't' has the biggest volume.  
        fliptets[0] = abtets[hullflag % 3]; // [a,b,c,d].
        eorgoppoself(fliptets[0]);  // [d,c,b,a]
        spintet = fliptets[0];
        while (1) {
          fnextself(spintet);
          chkpt = oppo(spintet);
          if (chkpt == pb) break;
          if ((chkpt != dummypoint) && (apex(spintet) != dummypoint)) {
            ori = -orient3d(pd, pc, apex(spintet), chkpt);
            assert(ori > 0);
            if (ori > bigvol) {
              bigvol = ori;
              searchpt = chkpt;
            }
          }
        }
        if (searchpt != NULL) { 
          // Now valid the configuration.
          ori1 = orient3d(pd, pc, searchpt, pa);
          ori2 = orient3d(pd, pc, searchpt, pb);
          if (ori1 * ori2 >= 0.0) {
            reducflag = 0; // Not valid. 
          } else {
            ori1 = orient3d(pa, pb, searchpt, pc);
            ori2 = orient3d(pa, pb, searchpt, pd);
            if (ori1 * ori2 >= 0.0) {
              reducflag = 0; // Not valid.
            }
          }
        } else {
          // No valid searchpt is found.
          reducflag = 0; // Do not flip it.
        }
      } // if (reducflag == 1)
    } // if (hullflag == 1)

    if (reducflag) {
      // A 3-to-2 flip is possible.
      if (checksubfaceflag) {
        // This edge (must not be a segment) can be flipped ONLY IF it belongs
        //   to either 0 or 2 subfaces.  In the latter case, a 2-to-2 flip in 
        //   the surface mesh will be automatically performed within the 
        //   3-to-2 flip.
        nn = 0;
        edgepivot = -1; // Re-use it.
        for (j = 0; j < 3; j++) {
          if (issubface(abtets[j])) {
            nn++; // Found a subface.
          } else {
            edgepivot = j;
          }
        }
        assert(nn < 3);
        if (nn == 1) {
          // Found only 1 subface containing this edge. This can happen in 
          //   the boundary recovery phase. The neighbor subface is not yet 
          //   recovered. This edge should not be flipped at this moment.
          rejflag = 1; 
        } else if (nn == 2) {
          // Found two subfaces. A 2-to-2 flip is possible. Validate it.
          // Below we check if the two faces [p,q,a] and [p,q,b] are subfaces.
          eorgoppo(abtets[(edgepivot + 1) % 3], spintet); // [q,p,b,a]
          if (issubface(spintet)) {
            rejflag = 1; // Conflict to a 2-to-2 flip.
          } else {
            esymself(spintet);
            if (issubface(spintet)) {
              rejflag = 1; // Conflict to a 2-to-2 flip.
            }
          }
        }
      }
      if (!rejflag && fc->checkflipeligibility) {
        // Here we must exchange 'a' and 'b'. Since in the check... function,
        //   we assume the following point sequence, 'a,b,c,d,e', where
        //   the face [a,b,c] will be flipped and the edge [e,d] will be
        //   created. The two new tets are [a,b,c,d] and [b,a,c,e]. 
        rejflag = checkflipeligibility(2, pc, pd, pe, pb, pa, level, 
                                       abedgepivot, fc);
      }
      if (!rejflag) {
        // Do flip: [a,b] => [c,d,e]
        flip32(abtets, hullflag, fc);
        if (fc->remove_ndelaunay_edge) {
          if (level == 0) {
            // It is the desired removing edge. Check if we have improved
            //   the objective function.
            if ((fc->tetprism_vol_sum >= 0.0) ||
                (fabs(fc->tetprism_vol_sum) < fc->bak_tetprism_vol)) {
              // No improvement! flip back: [c,d,e] => [a,b].
              flip23(abtets, hullflag, fc);
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
      }
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
      flip23(abtets, 1, fc);
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
        flip32(fliptets, 1, fc);
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
    }
  } // i

  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertpoint()    Insert a point into current tetrahedralization.          //
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
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::insertpoint(point insertpt, triface *searchtet, face *splitsh,
                            face *splitseg, insertvertexflags *ivf)
{
  arraypool *swaplist;
  triface *cavetet, spintet, neightet, neineitet, *parytet;
  triface oldtet, newtet, newneitet;
  face checksh, neighsh, *parysh;
  face checkseg, *paryseg;
  point *pts, pa, pb, pc, *parypt;
  enum locateresult loc = OUTSIDE;
  REAL sign, ori;
  REAL attrib, volume;
  bool enqflag;
  int t1ver;
  int i, j, k, s;

  if (b->verbose > 2) {
    printf("      Insert point %d\n", pointmark(insertpt));
  }

  // Locate the point.
  if (searchtet->tet != NULL) {
    loc = (enum locateresult) ivf->iloc;
  }

  if (loc == OUTSIDE) {
    if (searchtet->tet == NULL) {
      if (!b->weighted) {
        randomsample(insertpt, searchtet);
      } else {
        // Weighted DT. There may exist dangling vertex. 
        *searchtet = recenttet;
      }
    }
    // Locate the point.
    loc = locate(insertpt, searchtet); 
  }

  ivf->iloc = (int) loc; // The return value.

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
        setpointtype(insertpt, NREGULARVERTEX);
        nonregularcount++;
        ivf->iloc = (int) NONREGULAR;
        return 0;
      }
    }
  }

  // Create the initial cavity C(p) which contains all tetrahedra that
  //   intersect p. It may include 1, 2, or n tetrahedra.
  // If p lies on a segment or subface, also create the initial sub-cavity
  //   sC(p) which contains all subfaces (and segment) which intersect p.

  if (loc == OUTSIDE) {
    flip14count++;
    // The current hull will be enlarged.
    // Add four adjacent boundary tets into list.
    for (i = 0; i < 4; i++) {
      decode(searchtet->tet[i], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;
  } else if (loc == INTETRAHEDRON) {
    flip14count++;
    // Add four adjacent boundary tets into list.
    for (i = 0; i < 4; i++) {
      decode(searchtet->tet[i], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;
  } else if (loc == ONFACE) {
    flip26count++;
    // Add six adjacent boundary tets into list.
    j = (searchtet->ver & 3); // The current face number.
    for (i = 1; i < 4; i++) { 
      decode(searchtet->tet[(j + i) % 4], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    decode(searchtet->tet[j], spintet);
    j = (spintet.ver & 3); // The current face number.
    for (i = 1; i < 4; i++) {
      decode(spintet.tet[(j + i) % 4], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
    }
    infect(spintet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = spintet;
    infect(*searchtet);
    caveoldtetlist->newindex((void **) &parytet);
    *parytet = *searchtet;

    if (ivf->splitbdflag) { 
      if ((splitsh != NULL) && (splitsh->sh != NULL)) {
        // Create the initial sub-cavity sC(p).
        smarktest(*splitsh);
        caveshlist->newindex((void **) &parysh);
        *parysh = *splitsh;
      }
    } // if (splitbdflag)
  } else if (loc == ONEDGE) {
    flipn2ncount++;
    // Add all adjacent boundary tets into list.
    spintet = *searchtet;
    while (1) {
      eorgoppo(spintet, neightet);
      decode(neightet.tet[neightet.ver & 3], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
      edestoppo(spintet, neightet);
      decode(neightet.tet[neightet.ver & 3], neightet);
      neightet.ver = epivot[neightet.ver];
      cavebdrylist->newindex((void **) &parytet);
      *parytet = neightet;
      infect(spintet);
      caveoldtetlist->newindex((void **) &parytet);
      *parytet = spintet;
      fnextself(spintet);
      if (spintet.tet == searchtet->tet) break;
    } // while (1)

    if (ivf->splitbdflag) {
      // Create the initial sub-cavity sC(p).
      if ((splitseg != NULL) && (splitseg->sh != NULL)) {
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
          neightet.ver = epivot[neightet.ver];
          cavebdrylist->newindex((void **) &parytet);
          *parytet = neightet;
        }
      }
    }
  } else if (loc == ONVERTEX) {
    // The point already exist. Do nothing and return.
    return 0;
  } 


  if (ivf->assignmeshsize) {
    // Assign mesh size for the new point.
    if (bgm != NULL) {
      // Interpolate the mesh size from the background mesh. 
      bgm->decode(point2bgmtet(org(*searchtet)), neightet);
      int bgmloc = (int) bgm->scoutpoint(insertpt, &neightet, 0);
      if (bgmloc != (int) OUTSIDE) {
        insertpt[pointmtrindex] =  
          bgm->getpointmeshsize(insertpt, &neightet, bgmloc);
        setpoint2bgmtet(insertpt, bgm->encode(neightet));
      }
    } else {
      insertpt[pointmtrindex] = getpointmeshsize(insertpt,searchtet,(int)loc);
    }
  } // if (assignmeshsize)

  if (ivf->bowywat) {
    // Update the cavity C(p) using the Bowyer-Watson algorithm.
    swaplist = cavetetlist;
    cavetetlist = cavebdrylist;
    cavebdrylist = swaplist;
    for (i = 0; i < cavetetlist->objects; i++) {
      // 'cavetet' is an adjacent tet at outside of the cavity.
      cavetet = (triface *) fastlookup(cavetetlist, i);
      // The tet may be tested and included in the (enlarged) cavity.
      if (!infected(*cavetet)) {
        // Check for two possible cases for this tet: 
        //   (1) It is a cavity tet, or
        //   (2) it is a cavity boundary face.
        enqflag = false;
        if (!marktested(*cavetet)) {
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
                decode(cavetet->tet[3], neineitet);
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
                  } 
                } else {
                  // The adjacent tet is non-Delaunay. The hull face is non-
                  //   Delaunay as well. Include it in the cavity.
                  enqflag = true;
                } // if (!infected(neineitet))
              } // if (ori == 0.0)
            } else {
              // A hull face (must be a subface).
              // We FIRST include it in the initial cavity if the adjacent tet
              //   (not faked) of this hull face is not Delaunay wrt p.
              //   Whether it belongs to the final cavity will be determined
              //   during the validation process. 'validflag'.
              decode(cavetet->tet[3], neineitet);
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
                } 
              } else {
                // The adjacent tet is non-Delaunay. The hull face is non-
                //   Delaunay as well. Include it in the cavity.
                enqflag = true;
              } // if (infected(neineitet))
            } // if (nonconvex)
          } // if (pts[7] != dummypoint)
          marktest(*cavetet); // Only test it once.
        } // if (!marktested(*cavetet))

        if (enqflag) {
          // Found a tet in the cavity. Put other three faces in check list.
          k = (cavetet->ver & 3); // The current face number
          for (j = 1; j < 4; j++) {
            decode(cavetet->tet[(j + k) % 4], neightet);
            cavetetlist->newindex((void **) &parytet);
            *parytet = neightet;
          }
          infect(*cavetet);
          caveoldtetlist->newindex((void **) &parytet);
          *parytet = *cavetet;
        } else {
          // Found a boundary face of the cavity. 
          cavetet->ver = epivot[cavetet->ver];
          cavebdrylist->newindex((void **) &parytet);
          *parytet = *cavetet;
        }
      } // if (!infected(*cavetet))
    } // i

    cavetetlist->restart(); // Clear the working list.
  } // if (ivf->bowywat)

  if (checksubsegflag) {
    // Collect all segments of C(p).
    shellface *ssptr;
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      if ((ssptr = (shellface*) cavetet->tet[8]) != NULL) {
        for (j = 0; j < 6; j++) {
          if (ssptr[j]) {
            sdecode(ssptr[j], checkseg);
            if (!sinfected(checkseg)) {
              sinfect(checkseg);
              cavetetseglist->newindex((void **) &paryseg);
              *paryseg = checkseg;
            }
          }
        } // j
      }
    } // i
    // Uninfect collected segments.
    for (i = 0; i < cavetetseglist->objects; i++) {
      paryseg = (face *) fastlookup(cavetetseglist, i);
      suninfect(*paryseg);
    }

    if (ivf->rejflag & 1) {
      // Reject this point if it encroaches upon any segment.
      face *paryseg1;
      for (i = 0; i < cavetetseglist->objects; i++) {
        paryseg1 = (face *) fastlookup(cavetetseglist, i);
        if (checkseg4encroach((point) paryseg1->sh[3], (point) paryseg1->sh[4], 
                              insertpt)) {
          encseglist->newindex((void **) &paryseg);
          *paryseg = *paryseg1;
        }
      } // i
      if (encseglist->objects > 0) {
        insertpoint_abort(splitseg, ivf);
        ivf->iloc = (int) ENCSEGMENT;
        return 0;
      }
    }
  } // if (checksubsegflag)

  if (checksubfaceflag) {
    // Collect all subfaces of C(p).
    shellface *sptr;
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      if ((sptr = (shellface*) cavetet->tet[9]) != NULL) {
        for (j = 0; j < 4; j++) {
          if (sptr[j]) {
            sdecode(sptr[j], checksh);
            if (!sinfected(checksh)) {
              sinfect(checksh);
              cavetetshlist->newindex((void **) &parysh);
              *parysh = checksh;
            }
          }
        } // j
      }
    } // i
    // Uninfect collected subfaces.
    for (i = 0; i < cavetetshlist->objects; i++) {
      parysh = (face *) fastlookup(cavetetshlist, i);
      suninfect(*parysh);
    }

    if (ivf->rejflag & 2) {
      REAL rd, cent[3];
      badface *bface;
      // Reject this point if it encroaches upon any subface.
      for (i = 0; i < cavetetshlist->objects; i++) {
        parysh = (face *) fastlookup(cavetetshlist, i);
        if (checkfac4encroach((point) parysh->sh[3], (point) parysh->sh[4], 
                              (point) parysh->sh[5], insertpt, cent, &rd)) {
          encshlist->newindex((void **) &bface);
          bface->ss = *parysh;
          bface->forg = (point) parysh->sh[3]; // Not a dad one.
          for (j = 0; j < 3; j++) bface->cent[j] = cent[j];
          bface->key = rd;
        }
      }
      if (encshlist->objects > 0) {
        insertpoint_abort(splitseg, ivf);
        ivf->iloc = (int) ENCSUBFACE;
        return 0;
      }
    }
  } // if (checksubfaceflag)

  if ((ivf->iloc == (int) OUTSIDE) && ivf->refineflag) {
    // The vertex lies outside of the domain. And it does not encroach
    //   upon any boundary segment or subface. Do not insert it.
    insertpoint_abort(splitseg, ivf);
    return 0;
  }

  if (ivf->splitbdflag) { 
    // The new point locates in surface mesh. Update the sC(p). 
    // We have already 'smarktested' the subfaces which directly intersect
    //   with p in 'caveshlist'. From them, we 'smarktest' their neighboring
    //   subfaces which are included in C(p). Do not across a segment.
    for (i = 0; i < caveshlist->objects; i++) {
      parysh = (face *) fastlookup(caveshlist, i);
      assert(smarktested(*parysh));
      checksh = *parysh;
      for (j = 0; j < 3; j++) {
        if (!isshsubseg(checksh)) {
          spivot(checksh, neighsh);
          assert(neighsh.sh != NULL);
          if (!smarktested(neighsh)) {
            stpivot(neighsh, neightet);
            if (infected(neightet)) {
              fsymself(neightet);
              if (infected(neightet)) {
                // This subface is inside C(p). 
                // Check if its diametrical circumsphere encloses 'p'.
                //   The purpose of this check is to avoid forming invalid
                //   subcavity in surface mesh.
                sign = incircle3d(sorg(neighsh), sdest(neighsh),
                                  sapex(neighsh), insertpt);
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
  } // if (ivf->splitbdflag)

  if (ivf->validflag) {
    // Validate C(p) and update it if it is not star-shaped.
    int cutcount = 0;

    if (ivf->respectbdflag) {
      // The initial cavity may include subfaces which are not on the facets
      //   being splitting. Find them and make them as boundary of C(p). 
      // Comment: We have already 'smarktested' the subfaces in sC(p). They
      //   are completely inside C(p). 
      for (i = 0; i < cavetetshlist->objects; i++) {
        parysh = (face *) fastlookup(cavetetshlist, i);
        stpivot(*parysh, neightet);
        if (infected(neightet)) {
          fsymself(neightet);
          if (infected(neightet)) {
            // Found a subface inside C(p).
            if (!smarktested(*parysh)) {
              // It is possible that this face is a boundary subface.
              // Check if it is a hull face.
              //assert(apex(neightet) != dummypoint);
              if (oppo(neightet) != dummypoint) {
                fsymself(neightet);
              }
              if (oppo(neightet) != dummypoint) {
                ori = orient3d(org(neightet), dest(neightet), apex(neightet),
                               insertpt);
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
                uninfect(neightet);
                unmarktest(neightet);
                cutcount++;
                neightet.ver = epivot[neightet.ver];
                cavebdrylist->newindex((void **) &parytet);
                *parytet = neightet;
                // Add three new faces to find new boundaries.
                for (j = 0; j < 3; j++) {
                  esym(neightet, neineitet);
                  neineitet.ver = epivot[neineitet.ver];
                  cavebdrylist->newindex((void **) &parytet);
                  *parytet = neineitet;
                  enextself(neightet);
                }
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
            neightet.ver = epivot[neightet.ver];
            cavebdrylist->newindex((void **) &parytet);
            *parytet = neightet;
            // Add three new faces to find new boundaries.
            for (j = 0; j < 3; j++) {
              esym(neightet, neineitet);
              neineitet.ver = epivot[neineitet.ver];
              cavebdrylist->newindex((void **) &parytet);
              *parytet = neineitet;
              enextself(neightet);
            }
          }
        }
      } // i
    } // if (ivf->respectbdflag)

    // Update the cavity by removing invisible faces until it is star-shaped.
    for (i = 0; i < cavebdrylist->objects; i++) {
      cavetet = (triface *) fastlookup(cavebdrylist, i);
      // 'cavetet' is an exterior tet adjacent to the cavity.
      // Check if its neighbor is inside C(p).
      fsym(*cavetet, neightet);
      if (infected(neightet)) {        
        if (apex(*cavetet) != dummypoint) {
          // It is a cavity boundary face. Check its visibility.
          if (oppo(neightet) != dummypoint) {
            ori = orient3d(org(*cavetet), dest(*cavetet), apex(*cavetet),
                           insertpt); 
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
          uninfect(neightet);
          unmarktest(neightet);
          cutcount++;
          // Add three new faces to find new boundaries.
          for (j = 0; j < 3; j++) {
            esym(neightet, neineitet);
            neineitet.ver = epivot[neineitet.ver];
            cavebdrylist->newindex((void **) &parytet);
            *parytet = neineitet;
            enextself(neightet);
          }
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
        insertpoint_abort(splitseg, ivf);
        ivf->iloc = (int) BADELEMENT;
        return 0;
      }

      if (ivf->splitbdflag) { 
        int cutshcount = 0;
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
            if ((splitsh != NULL) && (splitsh->sh != NULL)) {
              // The to-be split subface should be in sC(p).
              if (!smarktested(*splitsh)) i++;
            }
          } else if (loc == ONEDGE) {
            if ((splitseg != NULL) && (splitseg->sh != NULL)) {
              // The to-be split segment should be in sC(p).
              if (!smarktested(*splitseg)) i++;
            }
            if ((splitsh != NULL) && (splitsh->sh != NULL)) {
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
            insertpoint_abort(splitseg, ivf);
            ivf->iloc = (int) BADELEMENT;
            return 0;
          }
        } // if (cutshcount > 0)
      } // if (ivf->splitbdflag)
    } // if (cutcount > 0)

  } // if (ivf->validflag)

  if (ivf->refineflag) {
    // The new point is inserted by Delaunay refinement, i.e., it is the 
    //   circumcenter of a tetrahedron, or a subface, or a segment.
    //   Do not insert this point if the tetrahedron, or subface, or segment
    //   is not inside the final cavity.
    if (((ivf->refineflag == 1) && !infected(ivf->refinetet)) ||
        ((ivf->refineflag == 2) && !smarktested(ivf->refinesh))) {
      insertpoint_abort(splitseg, ivf);
      ivf->iloc = (int) BADELEMENT;
      return 0;
    }
  } // if (ivf->refineflag)

  if (b->plc && (loc != INSTAR)) {
    // Reject the new point if it lies too close to an existing point (b->plc),
    // or it lies inside a protecting ball of near vertex (ivf->rejflag & 4).
    // Collect the list of vertices of the initial cavity.
    if (loc == OUTSIDE) {
      pts = (point *) &(searchtet->tet[4]);
      for (i = 0; i < 3; i++) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = pts[i];
      }
    } else if (loc == INTETRAHEDRON) {
      pts = (point *) &(searchtet->tet[4]);
      for (i = 0; i < 4; i++) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = pts[i];
      } 
    } else if (loc == ONFACE) {
      pts = (point *) &(searchtet->tet[4]);
      for (i = 0; i < 3; i++) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = pts[i];
      }
      if (pts[3] != dummypoint) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = pts[3];
      }
      fsym(*searchtet, spintet);
      if (oppo(spintet) != dummypoint) {
        cavetetvertlist->newindex((void **) &parypt);
        *parypt = oppo(spintet);
      }
    } else if (loc == ONEDGE) {
      spintet = *searchtet;
      cavetetvertlist->newindex((void **) &parypt);
      *parypt = org(spintet);
      cavetetvertlist->newindex((void **) &parypt);
      *parypt = dest(spintet);
      while (1) {
        if (apex(spintet) != dummypoint) {
          cavetetvertlist->newindex((void **) &parypt);
          *parypt = apex(spintet);
        }
        fnextself(spintet);
        if (spintet.tet == searchtet->tet) break;
      }
    }

    int rejptflag = (ivf->rejflag & 4);
    REAL rd;
    pts = NULL;

    for (i = 0; i < cavetetvertlist->objects; i++) {
      parypt = (point *) fastlookup(cavetetvertlist, i);
      rd = distance(*parypt, insertpt);
      // Is the point very close to an existing point?
      if (rd < b->minedgelength) {
        pts = parypt; 
        loc = NEARVERTEX;
        break;
      }
      if (rejptflag) {
        // Is the point encroaches upon an existing point?
        if (rd < (0.5 * (*parypt)[pointmtrindex])) {
          pts = parypt;
          loc = ENCVERTEX; 
          break;
        }
      }
    }
    cavetetvertlist->restart(); // Clear the work list.

    if (pts != NULL) {
      // The point is either too close to an existing vertex (NEARVERTEX)
      //   or encroaches upon (inside the protecting ball) of that vertex.
      if (loc == NEARVERTEX) {
        if (b->nomergevertex) { // -M0/1 option.
          // In this case, we still insert this vertex. Although it is very
          //   close to an existing vertex. Give a warning, anyway.
          if (!b->quiet) {
            printf("Warning:  Two points, %d and %d, are very close.\n",
                   pointmark(insertpt), pointmark(*pts));
            printf("  Creating a very short edge (len = %g) (< %g).\n",
                   rd, b->minedgelength);
            printf("  You may try a smaller tolerance (-T) (current is %g)\n", 
                   b->epsilon);
            printf("  to avoid this warning.\n");
          }
        } else {
          insertpt[3] = rd; // Only for reporting.
          setpoint2ppt(insertpt, *pts);
          insertpoint_abort(splitseg, ivf);
          ivf->iloc = (int) loc;
          return 0;
        }
      } else { // loc == ENCVERTEX
        // The point lies inside the protection ball.
        setpoint2ppt(insertpt, *pts);
        insertpoint_abort(splitseg, ivf);
        ivf->iloc = (int) loc;
        return 0;
      }
    }
  } // if (b->plc && (loc != INSTAR))

  if (b->weighted || ivf->cdtflag || ivf->smlenflag
      ) {
    // There may be other vertices inside C(p). We need to find them.
    // Collect all vertices of C(p).
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      //assert(infected(*cavetet));
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
    // Uninfect all collected (cavity) vertices.
    for (i = 0; i < cavetetvertlist->objects; i++) {
      parypt = (point *) fastlookup(cavetetvertlist, i);
      puninfect(*parypt);
    }
    if (ivf->smlenflag) {
      REAL len;
      // Get the length of the shortest edge connecting to 'newpt'.
      parypt = (point *) fastlookup(cavetetvertlist, 0);
      ivf->smlen = distance(*parypt, insertpt);
      ivf->parentpt = *parypt;
      for (i = 1; i < cavetetvertlist->objects; i++) {
        parypt = (point *) fastlookup(cavetetvertlist, i);
        len = distance(*parypt, insertpt);
        if (len < ivf->smlen) {
          ivf->smlen = len;
          ivf->parentpt = *parypt;
        }
      } 
    }
  }


  if (ivf->cdtflag) {
    // Unmark tets.
    for (i = 0; i < caveoldtetlist->objects; i++) {
      cavetet = (triface *) fastlookup(caveoldtetlist, i);
      unmarktest(*cavetet);
    }
    for (i = 0; i < cavebdrylist->objects; i++) {
      cavetet = (triface *) fastlookup(cavebdrylist, i);
      unmarktest(*cavetet);
    }
    // Clean up arrays which are not needed.
    cavetetlist->restart();
    if (checksubsegflag) {
      cavetetseglist->restart();
    }
    if (checksubfaceflag) {
      cavetetshlist->restart();
    }
    return 1;
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
            neineitet = spintet; // An outer tet. Remember it.
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
          if (!ivf->splitbdflag) {
            checkseg = *paryseg;
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
            checksh = *parysh; // Remember this side.
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
          if (!ivf->splitbdflag) {
            checksh = *parysh;
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
  } // if (checksubfaceflag)

  // Create new tetrahedra to fill the cavity.

  for (i = 0; i < cavebdrylist->objects; i++) {
    cavetet = (triface *) fastlookup(cavebdrylist, i);
    neightet = *cavetet;
    unmarktest(neightet); // Unmark it.
    // Get the oldtet (inside the cavity).
    fsym(neightet, oldtet);
    if (apex(neightet) != dummypoint) {
      // Create a new tet in the cavity.
      maketetrahedron(&newtet);
      setorg(newtet, dest(neightet));
      setdest(newtet, org(neightet));
      setapex(newtet, apex(neightet));
      setoppo(newtet, insertpt);
    } else {
      // Create a new hull tet.
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

  // Re-use this list to save new interior cavity faces.
  cavetetlist->restart();

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
        assert(newneitet.tet[newneitet.ver & 3] == NULL);
        bond(neightet, newneitet);
        if (ivf->lawson > 1) { 
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

  if (((splitsh != NULL) && (splitsh->sh != NULL)) ||
      ((splitseg != NULL) && (splitseg->sh != NULL))) {
    // Split a subface or a segment.
    sinsertvertex(insertpt, splitsh, splitseg, ivf->sloc, ivf->sbowywat, 0);
  }

  if (checksubfaceflag) {
    if (ivf->splitbdflag) {
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
      // The Boundary recovery phase.
      // Put all new subfaces into stack for recovery.
      for (i = 0; i < caveshbdlist->objects; i++) {
        // Get an old subface at edge [a, b].
        parysh = (face *) fastlookup(caveshbdlist, i);
        spivot(*parysh, checksh); // The new subface [a, b, p].
        // Do not recover a deleted new face (degenerated).
        if (checksh.sh[3] != NULL) {
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
    if (ivf->splitbdflag) {
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
            point2tetorg(sorg(checkseg), neightet);
            finddirection(&neightet, sdest(checkseg));
            assert(dest(neightet) == sdest(checkseg));
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
      // The Boundary Recovery Phase.  
      // Queue missing segments in C(p) for recovery.
      if (splitseg != NULL) {
        // Queue two new subsegments in C(p) for recovery.
        for (i = 0; i < cavesegshlist->objects; i++) {
          paryseg = (face *) fastlookup(cavesegshlist, i);
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

  if (b->weighted
      ) {
    // Some vertices may be completed inside the cavity. They must be
    //   detected and added to recovering list.
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
        }
      }
    }
  }

  if (ivf->chkencflag & 1) {
    // Queue all segment outside C(p).
    for (i = 0; i < cavetetseglist->objects; i++) {
      paryseg = (face *) fastlookup(cavetetseglist, i);
      // Skip if it is the split segment.
      if (!sinfected(*paryseg)) {
        enqueuesubface(badsubsegs, paryseg);
      }
    }
    if (splitseg != NULL) {
      // Queue the two new subsegments inside C(p).
      for (i = 0; i < cavesegshlist->objects; i++) {
        paryseg = (face *) fastlookup(cavesegshlist, i);
        enqueuesubface(badsubsegs, paryseg);
      }
    }
  } // if (chkencflag & 1)

  if (ivf->chkencflag & 2) {
    // Queue all subfaces outside C(p).
    for (i = 0; i < cavetetshlist->objects; i++) {
      parysh = (face *) fastlookup(cavetetshlist, i);
      // Skip if it is a split subface.
      if (!sinfected(*parysh)) {
        enqueuesubface(badsubfacs, parysh);
      }
    }
    // Queue all new subfaces inside C(p).
    for (i = 0; i < caveshbdlist->objects; i++) {
      // Get an old subface at edge [a, b].
      parysh = (face *) fastlookup(caveshbdlist, i);
      spivot(*parysh, checksh); // checksh is a new subface [a, b, p].
      // Do not recover a deleted new face (degenerated).
      if (checksh.sh[3] != NULL) {
        enqueuesubface(badsubfacs, &checksh);
      }
    }
  } // if (chkencflag & 2)

  if (ivf->chkencflag & 4) {
    // Queue all new tetrahedra in C(p).
    for (i = 0; i < cavebdrylist->objects; i++) {
      cavetet = (triface *) fastlookup(cavebdrylist, i);
      enqueuetetrahedron(cavetet);
    }
  }

  // C(p) is re-meshed successfully.

  // Delete the old tets in C(p).
  for (i = 0; i < caveoldtetlist->objects; i++) {
    searchtet = (triface *) fastlookup(caveoldtetlist, i);
    if (ishulltet(*searchtet)) {
      hullsize--;
    }
    tetrahedrondealloc(searchtet->tet);
  }

  if (((splitsh != NULL) && (splitsh->sh != NULL)) ||
      ((splitseg != NULL) && (splitseg->sh != NULL))) {
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
    if ((splitseg != NULL) && (splitseg->sh != NULL)) {
      // Delete the old segment in sC(p).
      shellfacedealloc(subsegs, splitseg->sh);
    }
  }

  if (ivf->lawson) {
    for (i = 0; i < cavebdrylist->objects; i++) {
      searchtet = (triface *) fastlookup(cavebdrylist, i);
      flippush(flipstack, searchtet);
    }
    if (ivf->lawson > 1) {
      for (i = 0; i < cavetetlist->objects; i++) {
        searchtet = (triface *) fastlookup(cavetetlist, i);
        flippush(flipstack, searchtet);
      }
    }
  }


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
  
  if (b->weighted || ivf->validflag) {
    cavetetvertlist->restart();
  }
  
  if (((splitsh != NULL) && (splitsh->sh != NULL)) ||
      ((splitseg != NULL) && (splitseg->sh != NULL))) {
    caveshlist->restart();
    caveshbdlist->restart();
    cavesegshlist->restart();
  }

  return 1; // Point is inserted.
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertpoint_abort()    Abort the insertion of a new vertex.               //
//                                                                           //
// The cavity will be restored.  All working lists are cleared.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::insertpoint_abort(face *splitseg, insertvertexflags *ivf)
{
  triface *cavetet;
  face *parysh;
  int i;

  for (i = 0; i < caveoldtetlist->objects; i++) {
    cavetet = (triface *) fastlookup(caveoldtetlist, i);
    uninfect(*cavetet);
    unmarktest(*cavetet);
  }
  for (i = 0; i < cavebdrylist->objects; i++) {
    cavetet = (triface *) fastlookup(cavebdrylist, i);
    unmarktest(*cavetet); 
  }
  cavetetlist->restart();
  cavebdrylist->restart();
  caveoldtetlist->restart();
  cavetetseglist->restart();
  cavetetshlist->restart();
  if (ivf->splitbdflag) { 
    if ((splitseg != NULL) && (splitseg->sh != NULL)) {
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
}

////                                                                       ////
////                                                                       ////
//// flip_cxx /////////////////////////////////////////////////////////////////

