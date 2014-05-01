#include "../tetgen.h"
//// meshstat_cxx /////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// printfcomma()    Print a (large) number with the 'thousands separator'.   // 
//                                                                           //
// The following code was simply copied from "stackoverflow".                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::printfcomma(unsigned long n)
{
  unsigned long n2 = 0;
  int scale = 1;
  while (n >= 1000) {
    n2 = n2 + scale * (n % 1000);
    n /= 1000;
    scale *= 1000;
  }
  printf ("%ld", n);
  while (scale != 1) {
    scale /= 1000;
    n = n2 / scale;
    n2 = n2  % scale;
    printf (",%03ld", n);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkmesh()    Test the mesh for topological consistency.                 //
//                                                                           //
// If 'topoflag' is set, only check the topological connection of the mesh,  //
// i.e., do not report degenerated or inverted elements.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkmesh(int topoflag)
{
  triface tetloop, neightet, symtet;
  point pa, pb, pc, pd;
  REAL ori;
  int horrors, i;

  if (!b->quiet) {
    printf("  Checking consistency of mesh...\n");
  }

  horrors = 0;
  tetloop.ver = 0;
  // Run through the list of tetrahedra, checking each one.
  tetrahedrons->traversalinit();
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Check all four faces of the tetrahedron.
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      pa = org(tetloop);
      pb = dest(tetloop);
      pc = apex(tetloop);
      pd = oppo(tetloop);
      if (tetloop.ver == 0) {  // Only test for inversion once.
        if (!ishulltet(tetloop)) {  // Only do test if it is not a hull tet.
          if (!topoflag) {
            ori = orient3d(pa, pb, pc, pd);
            if (ori >= 0.0) {
              printf("  !! !! %s ", ori > 0.0 ? "Inverted" : "Degenerated");
              printf("  (%d, %d, %d, %d) (ori = %.17g)\n", pointmark(pa),
                     pointmark(pb), pointmark(pc), pointmark(pd), ori);
              horrors++;
            }
          }
        }
        if (infected(tetloop)) { 
          // This may be a bug. Report it.
          printf("  !! (%d, %d, %d, %d) is infected.\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          horrors++;
        }
        if (marktested(tetloop)) {
          // This may be a bug. Report it.
          printf("  !! (%d, %d, %d, %d) is marked.\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          horrors++;
        }
      }
      if (tetloop.tet[tetloop.ver] == NULL) {
        printf("  !! !! No neighbor at face (%d, %d, %d).\n", pointmark(pa),
               pointmark(pb), pointmark(pc));
        horrors++;
      } else {
        // Find the neighboring tetrahedron on this face.
        fsym(tetloop, neightet);
        // Check that the tetrahedron's neighbor knows it's a neighbor.
        fsym(neightet, symtet);
        if ((tetloop.tet != symtet.tet) || (tetloop.ver != symtet.ver)) {
          printf("  !! !! Asymmetric tetra-tetra bond:\n");
          if (tetloop.tet == symtet.tet) {
            printf("   (Right tetrahedron, wrong orientation)\n");
          }
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
                 pointmark(dest(neightet)), pointmark(apex(neightet)),
                 pointmark(oppo(neightet)));
          horrors++;
        }
        // Check if they have the same edge (the bond() operation).
        if ((org(neightet) != pb) || (dest(neightet) != pa)) {
          printf("  !! !! Wrong edge-edge bond:\n");
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
                 pointmark(dest(neightet)), pointmark(apex(neightet)),
                 pointmark(oppo(neightet)));
          horrors++;
        }
        // Check if they have the same apex.
        if (apex(neightet) != pc) {
          printf("  !! !! Wrong face-face bond:\n");
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
                 pointmark(dest(neightet)), pointmark(apex(neightet)),
                 pointmark(oppo(neightet)));
          horrors++;
        }
        // Check if they have the same opposite.
        if (oppo(neightet) == pd) {
          printf("  !! !! Two identical tetra:\n");
          printf("    First:  (%d, %d, %d, %d)\n", pointmark(pa),
                 pointmark(pb), pointmark(pc), pointmark(pd));
          printf("    Second: (%d, %d, %d, %d)\n", pointmark(org(neightet)),
                 pointmark(dest(neightet)), pointmark(apex(neightet)),
                 pointmark(oppo(neightet)));
          horrors++;
        }
      }
      if (facemarked(tetloop)) {
        // This may be a bug. Report it.
        printf("  !! tetface (%d, %d, %d) %d is marked.\n", pointmark(pa),
               pointmark(pb), pointmark(pc), pointmark(pd));
      }
    }
    // Check the six edges of this tet.
    for (i = 0; i < 6; i++) {
      tetloop.ver = edge2ver[i];
      if (edgemarked(tetloop)) {
        // This may be a bug. Report it.
        printf("  !! tetedge (%d, %d) %d, %d is marked.\n", 
               pointmark(org(tetloop)), pointmark(dest(tetloop)), 
               pointmark(apex(tetloop)), pointmark(oppo(tetloop)));
      }
    }
    tetloop.tet = alltetrahedrontraverse();
  }
  if (horrors == 0) {
    if (!b->quiet) {
      printf("  In my studied opinion, the mesh appears to be consistent.\n");
    }
  } else {
    printf("  !! !! !! !! %d %s witnessed.\n", horrors, 
           horrors > 1 ? "abnormity" : "abnormities");
  }

  return horrors;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkshells()       Test the boundary mesh for topological consistency.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkshells()
{
  triface neightet, symtet;
  face shloop, spinsh, nextsh;
  face checkseg;
  point pa, pb;
  int bakcount;
  int horrors, i;

  if (!b->quiet) {
    printf("  Checking consistency of the mesh boundary...\n");
  }
  horrors = 0;

  void **bakpathblock = subfaces->pathblock;
  void *bakpathitem = subfaces->pathitem;
  int bakpathitemsleft = subfaces->pathitemsleft;
  int bakalignbytes = subfaces->alignbytes;

  subfaces->traversalinit();
  shloop.sh = shellfacetraverse(subfaces);
  while (shloop.sh != NULL) {
    shloop.shver = 0;
    for (i = 0; i < 3; i++) {
      // Check the face ring at this edge.
      pa = sorg(shloop);
      pb = sdest(shloop);
      spinsh = shloop;
      spivot(spinsh, nextsh);
      bakcount = horrors;
      while ((nextsh.sh != NULL) && (nextsh.sh != shloop.sh)) {
        if (nextsh.sh[3] == NULL) {
          printf("  !! !! Wrong subface-subface connection (Dead subface).\n");
          printf("    First: x%lx (%d, %d, %d).\n", (uintptr_t) spinsh.sh,
                 pointmark(sorg(spinsh)), pointmark(sdest(spinsh)), 
                 pointmark(sapex(spinsh)));
          printf("    Second: x%lx (DEAD)\n", (uintptr_t) nextsh.sh);
          horrors++;
          break;
        }
        // check if they have the same edge.
        if (!(((sorg(nextsh) == pa) && (sdest(nextsh) == pb)) ||
              ((sorg(nextsh) == pb) && (sdest(nextsh) == pa)))) {
           printf("  !! !! Wrong subface-subface connection.\n");
           printf("    First: x%lx (%d, %d, %d).\n", (uintptr_t) spinsh.sh,
                  pointmark(sorg(spinsh)), pointmark(sdest(spinsh)), 
                  pointmark(sapex(spinsh)));
           printf("    Scond: x%lx (%d, %d, %d).\n", (uintptr_t) nextsh.sh,
                  pointmark(sorg(nextsh)), pointmark(sdest(nextsh)), 
                  pointmark(sapex(nextsh)));
           horrors++;
           break;
        }
        // Check they should not have the same apex.
        if (sapex(nextsh) == sapex(spinsh)) {
           printf("  !! !! Existing two duplicated subfaces.\n");
           printf("    First: x%lx (%d, %d, %d).\n", (uintptr_t) spinsh.sh,
                  pointmark(sorg(spinsh)), pointmark(sdest(spinsh)), 
                  pointmark(sapex(spinsh)));
           printf("    Scond: x%lx (%d, %d, %d).\n", (uintptr_t) nextsh.sh,
                  pointmark(sorg(nextsh)), pointmark(sdest(nextsh)), 
                  pointmark(sapex(nextsh)));
           horrors++;
           break;
        }
        spinsh = nextsh;
        spivot(spinsh, nextsh);
      }
      // Check subface-subseg bond.
      sspivot(shloop, checkseg);
      if (checkseg.sh != NULL) {
        if (checkseg.sh[3] == NULL) {
          printf("  !! !! Wrong subface-subseg connection (Dead subseg).\n");
          printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) shloop.sh,
                 pointmark(sorg(shloop)), pointmark(sdest(shloop)), 
                 pointmark(sapex(shloop)));
          printf("    Sub: x%lx (Dead)\n", (uintptr_t) checkseg.sh);
          horrors++;
        } else {
          if (!(((sorg(checkseg) == pa) && (sdest(checkseg) == pb)) ||
                ((sorg(checkseg) == pb) && (sdest(checkseg) == pa)))) {
            printf("  !! !! Wrong subface-subseg connection.\n");
            printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) shloop.sh,
                   pointmark(sorg(shloop)), pointmark(sdest(shloop)), 
                   pointmark(sapex(shloop)));
            printf("    Seg: x%lx (%d, %d).\n", (uintptr_t) checkseg.sh,
                   pointmark(sorg(checkseg)), pointmark(sdest(checkseg)));
            horrors++;
          }
        }
      }
      if (horrors > bakcount) break; // An error detected. 
      senextself(shloop);
    }
    // Check tet-subface connection.
    stpivot(shloop, neightet);
    if (neightet.tet != NULL) {
      if (neightet.tet[4] == NULL) {
        printf("  !! !! Wrong sub-to-tet connection (Dead tet)\n");
        printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) shloop.sh,
               pointmark(sorg(shloop)), pointmark(sdest(shloop)), 
               pointmark(sapex(shloop)));
        printf("    Tet: x%lx (DEAD)\n", (uintptr_t) neightet.tet);
        horrors++;
      } else {
        if (!((sorg(shloop) == org(neightet)) && 
              (sdest(shloop) == dest(neightet)))) {
          printf("  !! !! Wrong sub-to-tet connection\n");
          printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) shloop.sh,
                 pointmark(sorg(shloop)), pointmark(sdest(shloop)), 
                 pointmark(sapex(shloop)));
          printf("    Tet: x%lx (%d, %d, %d, %d).\n",
                 (uintptr_t) neightet.tet, pointmark(org(neightet)), 
                 pointmark(dest(neightet)), pointmark(apex(neightet)),
                 pointmark(oppo(neightet)));
          horrors++;
        }
        tspivot(neightet, spinsh);
        if (!((sorg(spinsh) == org(neightet)) && 
              (sdest(spinsh) == dest(neightet)))) {
          printf("  !! !! Wrong tet-sub connection.\n");
          printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) spinsh.sh,
                 pointmark(sorg(spinsh)), pointmark(sdest(spinsh)), 
                 pointmark(sapex(spinsh)));
          printf("    Tet: x%lx (%d, %d, %d, %d).\n", 
                 (uintptr_t) neightet.tet, pointmark(org(neightet)), 
                 pointmark(dest(neightet)), pointmark(apex(neightet)), 
                 pointmark(oppo(neightet)));
          horrors++;
        }
        fsym(neightet, symtet);
        tspivot(symtet, spinsh);
        if (spinsh.sh != NULL) {
          if (!((sorg(spinsh) == org(symtet)) && 
                (sdest(spinsh) == dest(symtet)))) {
            printf("  !! !! Wrong tet-sub connection.\n");
            printf("    Sub: x%lx (%d, %d, %d).\n", (uintptr_t) spinsh.sh,
                   pointmark(sorg(spinsh)), pointmark(sdest(spinsh)), 
                   pointmark(sapex(spinsh)));
            printf("    Tet: x%lx (%d, %d, %d, %d).\n", 
                   (uintptr_t) symtet.tet, pointmark(org(symtet)), 
                   pointmark(dest(symtet)), pointmark(apex(symtet)), 
                   pointmark(oppo(symtet)));
            horrors++;
          }
        } else {
          printf("  Warning: Broken tet-sub-tet connection.\n");
        }
      }
    }
    if (sinfected(shloop)) {
      // This may be a bug. report it.
      printf("  !! A infected subface: (%d, %d, %d).\n", 
             pointmark(sorg(shloop)), pointmark(sdest(shloop)), 
             pointmark(sapex(shloop)));
    }
    if (smarktested(shloop)) {
      // This may be a bug. report it.
      printf("  !! A marked subface: (%d, %d, %d).\n", pointmark(sorg(shloop)), 
             pointmark(sdest(shloop)), pointmark(sapex(shloop)));
    }
    shloop.sh = shellfacetraverse(subfaces);
  }

  if (horrors == 0) {
    if (!b->quiet) {
      printf("  Mesh boundaries connected correctly.\n");
    }
  } else {
    printf("  !! !! !! !! %d boundary connection viewed with horror.\n",
           horrors);
  }

  subfaces->pathblock = bakpathblock;
  subfaces->pathitem = bakpathitem;
  subfaces->pathitemsleft = bakpathitemsleft;
  subfaces->alignbytes = bakalignbytes;

  return horrors;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checksegments()    Check the connections between tetrahedra and segments. //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checksegments()
{
  triface tetloop, neightet, spintet;
  shellface *segs;
  face neighsh, spinsh, checksh;
  face sseg, checkseg;
  point pa, pb;
  int miscount;
  int t1ver;
  int horrors, i;


  if (!b->quiet) {
    printf("  Checking tet->seg connections...\n");
  }

  horrors = 0;
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != NULL) {
    // Loop the six edges of the tet.
    if (tetloop.tet[8] != NULL) {
      segs = (shellface *) tetloop.tet[8];
      for (i = 0; i < 6; i++) {
        sdecode(segs[i], sseg);
        if (sseg.sh != NULL) {
          // Get the edge of the tet.
          tetloop.ver = edge2ver[i];
          // Check if they are the same edge.
          pa = (point) sseg.sh[3];
          pb = (point) sseg.sh[4];
          if (!(((org(tetloop) == pa) && (dest(tetloop) == pb)) ||
                ((org(tetloop) == pb) && (dest(tetloop) == pa)))) {
            printf("  !! Wrong tet-seg connection.\n");
            printf("    Tet: x%lx (%d, %d, %d, %d) - Seg: x%lx (%d, %d).\n", 
                   (uintptr_t) tetloop.tet, pointmark(org(tetloop)),
                   pointmark(dest(tetloop)), pointmark(apex(tetloop)),
                   pointmark(oppo(tetloop)), (uintptr_t) sseg.sh,
                   pointmark(pa), pointmark(pb));
            horrors++;
          } else {
            // Loop all tets sharing at this edge.
            neightet = tetloop;
            do {
              tsspivot1(neightet, checkseg);
              if (checkseg.sh != sseg.sh) {
                printf("  !! Wrong tet->seg connection.\n");
                printf("    Tet: x%lx (%d, %d, %d, %d) - ", 
                       (uintptr_t) neightet.tet, pointmark(org(neightet)),
                       pointmark(dest(neightet)), pointmark(apex(neightet)),
                       pointmark(oppo(neightet)));
                if (checkseg.sh != NULL) {
                  printf("Seg x%lx (%d, %d).\n", (uintptr_t) checkseg.sh,
                         pointmark(sorg(checkseg)),pointmark(sdest(checkseg))); 
                } else {
                  printf("Seg: NULL.\n");
                }
                horrors++;
              }
              fnextself(neightet);
            } while (neightet.tet != tetloop.tet);
          }
          // Check the seg->tet pointer.
          sstpivot1(sseg, neightet);
          if (neightet.tet == NULL) {
            printf("  !! Wrong seg->tet connection (A NULL tet).\n");
            horrors++;
          } else {
            if (!(((org(neightet) == pa) && (dest(neightet) == pb)) ||
                ((org(neightet) == pb) && (dest(neightet) == pa)))) {
              printf("  !! Wrong seg->tet connection (Wrong edge).\n");
              printf("    Tet: x%lx (%d, %d, %d, %d) - Seg: x%lx (%d, %d).\n", 
                     (uintptr_t) neightet.tet, pointmark(org(neightet)),
                     pointmark(dest(neightet)), pointmark(apex(neightet)),
                     pointmark(oppo(neightet)), (uintptr_t) sseg.sh,
                     pointmark(pa), pointmark(pb));
              horrors++;
            }
          }
        }
      }
    }
    // Loop the six edge of this tet.
    neightet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      neightet.ver = edge2ver[i];
      if (edgemarked(neightet)) {
        // A possible bug. Report it.
        printf("  !! A marked edge: (%d, %d, %d, %d) -- x%lx %d.\n",
               pointmark(org(neightet)), pointmark(dest(neightet)),
               pointmark(apex(neightet)), pointmark(oppo(neightet)),
               (uintptr_t) neightet.tet, neightet.ver);
        // Check if all tets at the edge are marked.
        spintet = neightet;
        while (1) {
          fnextself(spintet);
          if (!edgemarked(spintet)) {
            printf("  !! !! An unmarked edge (%d, %d, %d, %d) -- x%lx %d.\n",
                   pointmark(org(spintet)), pointmark(dest(spintet)),
                   pointmark(apex(spintet)), pointmark(oppo(spintet)),
                   (uintptr_t) spintet.tet, spintet.ver);
            horrors++;
          }
          if (spintet.tet == neightet.tet) break;
        }
      }
    }
    tetloop.tet = tetrahedrontraverse();
  }

  if (!b->quiet) {
    printf("  Checking seg->tet connections...\n");
  }

  miscount = 0; // Count the number of unrecovered segments.
  subsegs->traversalinit();
  sseg.shver = 0;
  sseg.sh = shellfacetraverse(subsegs);
  while (sseg.sh != NULL) {
    pa = sorg(sseg);
    pb = sdest(sseg);
    spivot(sseg, neighsh);
    if (neighsh.sh != NULL) {      
      spinsh = neighsh;
      while (1) {
        // Check seg-subface bond.
        if (((sorg(spinsh) == pa) && (sdest(spinsh) == pb)) ||
	    ((sorg(spinsh) == pb) && (sdest(spinsh) == pa))) {
          // Keep the same rotate direction.
          //if (sorg(spinsh) != pa) {          
          //  sesymself(spinsh);
          //  printf("  !! Wrong ori at subface (%d, %d, %d) -- x%lx %d\n",
          //         pointmark(sorg(spinsh)), pointmark(sdest(spinsh)),
          //         pointmark(sapex(spinsh)), (uintptr_t) spinsh.sh,
          //         spinsh.shver);
          //  horrors++;
          //}
          stpivot(spinsh, spintet);
          if (spintet.tet != NULL) {
            // Check if all tets at this segment.
            while (1) {
              tsspivot1(spintet, checkseg);
              if (checkseg.sh == NULL) {
                printf("  !! !! No seg at tet (%d, %d, %d, %d) -- x%lx %d\n",
                       pointmark(org(spintet)), pointmark(dest(spintet)),
                       pointmark(apex(spintet)), pointmark(oppo(spintet)),
                       (uintptr_t) spintet.tet, spintet.ver);
                horrors++;
              }
              if (checkseg.sh != sseg.sh) {
                printf("  !! !! Wrong seg (%d, %d) at tet (%d, %d, %d, %d)\n",
                       pointmark(sorg(checkseg)), pointmark(sdest(checkseg)),
                       pointmark(org(spintet)), pointmark(dest(spintet)),
                       pointmark(apex(spintet)), pointmark(oppo(spintet)));
                horrors++;
              }
              fnextself(spintet);
              // Stop at the next subface.
              tspivot(spintet, checksh);
              if (checksh.sh != NULL) break;
            } // while (1)
          }
        } else { 
          printf("  !! Wrong seg-subface (%d, %d, %d) -- x%lx %d connect\n",
                 pointmark(sorg(spinsh)), pointmark(sdest(spinsh)),
                 pointmark(sapex(spinsh)), (uintptr_t) spinsh.sh,
                 spinsh.shver);
          horrors++;
          break;
        } // if pa, pb
        spivotself(spinsh);
        if (spinsh.sh == NULL) break; // A dangling segment.
        if (spinsh.sh == neighsh.sh) break;
      } // while (1)
    } // if (neighsh.sh != NULL)
    // Count the number of "un-recovered" segments.
    sstpivot1(sseg, neightet);
    if (neightet.tet == NULL) {
      miscount++;
    }
    sseg.sh = shellfacetraverse(subsegs);
  }

  if (!b->quiet) {
    printf("  Checking seg->seg connections...\n");
  }

  points->traversalinit();
  pa = pointtraverse();
  while (pa != NULL) {
    if (pointtype(pa) == FREESEGVERTEX) {
      // There should be two subsegments connected at 'pa'.
      // Get a subsegment containing 'pa'.
      sdecode(point2sh(pa), sseg);
      if ((sseg.sh == NULL) || sseg.sh[3] == NULL) {
        printf("  !! Dead point-to-seg pointer at point %d.\n",
               pointmark(pa));
        horrors++;
      } else {
        sseg.shver = 0;
        if (sorg(sseg) != pa) {
          if (sdest(sseg) != pa) {
            printf("  !! Wrong point-to-seg pointer at point %d.\n",
                   pointmark(pa));
            horrors++;
          } else {
            // Find the next subsegment at 'pa'.
            senext(sseg, checkseg);
            if ((checkseg.sh == NULL) || (checkseg.sh[3] == NULL)) {
              printf("  !! Dead seg-seg connection at point %d.\n",
                     pointmark(pa));
              horrors++;
            } else {
              spivotself(checkseg);
              checkseg.shver = 0;
              if (sorg(checkseg) != pa) {
                printf("  !! Wrong seg-seg connection at point %d.\n",
                     pointmark(pa));
                horrors++;
              }
            }
          }
        } else {
          // Find the previous subsegment at 'pa'.
          senext2(sseg, checkseg);
          if ((checkseg.sh == NULL) || (checkseg.sh[3] == NULL)) {
            printf("  !! Dead seg-seg connection at point %d.\n",
                   pointmark(pa));
            horrors++;
          } else {
            spivotself(checkseg);
            checkseg.shver = 0;
            if (sdest(checkseg) != pa) {
              printf("  !! Wrong seg-seg connection at point %d.\n",
                     pointmark(pa));
              horrors++;
            }
          }
        }
      }
    }
    pa = pointtraverse();
  }

  if (horrors == 0) {
    printf("  Segments are connected properly.\n");
  } else {
    printf("  !! !! !! !! Found %d missing connections.\n", horrors);
  }
  if (miscount > 0) {
    printf("  !! !! Found %d missing segments.\n", miscount);
  }

  return horrors;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkdelaunay()    Ensure that the mesh is (constrained) Delaunay.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkdelaunay()
{
  triface tetloop;
  triface symtet;
  face checksh;
  point pa, pb, pc, pd, pe;
  REAL sign;
  int ndcount; // Count the non-locally Delaunay faces.
  int horrors;

  if (!b->quiet) {
    printf("  Checking Delaunay property of the mesh...\n");
  }

  ndcount = 0;
  horrors = 0;
  tetloop.ver = 0;
  // Run through the list of triangles, checking each one.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Check all four faces of the tetrahedron.
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, symtet);
      // Only do test if its adjoining tet is not a hull tet or its pointer
      //   is larger (to ensure that each pair isn't tested twice).
      if (((point) symtet.tet[7] != dummypoint)&&(tetloop.tet < symtet.tet)) {
        pa = org(tetloop);
        pb = dest(tetloop);
        pc = apex(tetloop);
        pd = oppo(tetloop);
        pe = oppo(symtet);
        sign = insphere_s(pa, pb, pc, pd, pe);
        if (sign < 0.0) {
          ndcount++;
          if (checksubfaceflag) {
            tspivot(tetloop, checksh);
          }
          if (checksh.sh == NULL) {
            printf("  !! Non-locally Delaunay (%d, %d, %d) - %d, %d\n",
                   pointmark(pa), pointmark(pb), pointmark(pc), pointmark(pd),
                   pointmark(pe));
            horrors++;
          }
        }
      }
    }
    tetloop.tet = tetrahedrontraverse();
  }

  if (horrors == 0) {
    if (!b->quiet) {
      if (ndcount > 0) {
        printf("  The mesh is constrained Delaunay.\n");
      } else {
        printf("  The mesh is Delaunay.\n");
      } 
    }
  } else {
    printf("  !! !! !! !! Found %d non-Delaunay faces.\n", horrors);
  }

  return horrors;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Check if the current tetrahedralization is (constrained) regular.         //
//                                                                           //
// The parameter 'type' determines which regularity should be checked:       //
//   - 0:  check the Delaunay property.                                      //
//   - 1:  check the Delaunay property with symbolic perturbation.           //
//   - 2:  check the regular property, the weights are stored in p[3].       //
//   - 3:  check the regular property with symbolic perturbation.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkregular(int type)
{
  triface tetloop;
  triface symtet;
  face checksh;
  point p[5];
  REAL sign;
  int ndcount; // Count the non-locally Delaunay faces.
  int horrors;

  if (!b->quiet) {
    printf("  Checking %s %s property of the mesh...\n",
           (type & 2) == 0 ? "Delaunay" : "regular",
           (type & 1) == 0 ? " " : "(s)");
  }

  // Make sure orient3d(p[1], p[0], p[2], p[3]) > 0;
  //   Hence if (insphere(p[1], p[0], p[2], p[3], p[4]) > 0) means that
  //     p[4] lies inside the circumsphere of p[1], p[0], p[2], p[3].
  //   The same if orient4d(p[1], p[0], p[2], p[3], p[4]) > 0 means that
  //     p[4] lies below the oriented hyperplane passing through 
  //     p[1], p[0], p[2], p[3].

  ndcount = 0;
  horrors = 0;
  tetloop.ver = 0;
  // Run through the list of triangles, checking each one.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Check all four faces of the tetrahedron.
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, symtet);
      // Only do test if its adjoining tet is not a hull tet or its pointer
      //   is larger (to ensure that each pair isn't tested twice).
      if (((point) symtet.tet[7] != dummypoint)&&(tetloop.tet < symtet.tet)) {
        p[0] = org(tetloop);   // pa
        p[1] = dest(tetloop);  // pb
        p[2] = apex(tetloop);  // pc
        p[3] = oppo(tetloop);  // pd
        p[4] = oppo(symtet);   // pe

        if (type == 0) {
          sign = insphere(p[1], p[0], p[2], p[3], p[4]);
        } else if (type == 1) {
          sign = insphere_s(p[1], p[0], p[2], p[3], p[4]);
        } else if (type == 2) {
          sign = orient4d(p[1],    p[0],    p[2],    p[3],    p[4], 
                          p[1][3], p[0][3], p[2][3], p[3][3], p[4][3]);
        } else { // type == 3
          sign = orient4d_s(p[1],    p[0],    p[2],    p[3],    p[4], 
                            p[1][3], p[0][3], p[2][3], p[3][3], p[4][3]);
        }

        if (sign > 0.0) {
          ndcount++;
          if (checksubfaceflag) {
            tspivot(tetloop, checksh);
          }
          if (checksh.sh == NULL) {
            printf("  !! Non-locally %s (%d, %d, %d) - %d, %d\n",
                   (type & 2) == 0 ? "Delaunay" : "regular",
		   pointmark(p[0]), pointmark(p[1]), pointmark(p[2]), 
                   pointmark(p[3]), pointmark(p[4]));
            horrors++;
          }
        }
      }
    }
    tetloop.tet = tetrahedrontraverse();
  }

  if (horrors == 0) {
    if (!b->quiet) {
      if (ndcount > 0) {
        printf("  The mesh is constrained %s.\n",
               (type & 2) == 0 ? "Delaunay" : "regular");
      } else {
        printf("  The mesh is %s.\n", (type & 2) == 0 ? "Delaunay" : "regular");
      } 
    }
  } else {
    printf("  !! !! !! !! Found %d non-%s faces.\n", horrors,
           (type & 2) == 0 ? "Delaunay" : "regular");
  }

  return horrors;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// checkconforming()    Ensure that the mesh is conforming Delaunay.         //
//                                                                           //
// If 'flag' is 1, only check subsegments. If 'flag' is 2, check subfaces.   //
// If 'flag' is 3, check both subsegments and subfaces.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int tetgenmesh::checkconforming(int flag)
{
  triface searchtet, neightet, spintet;
  face shloop;
  face segloop;
  point eorg, edest, eapex, pa, pb, pc;
  REAL cent[3], radius, dist, diff, rd, len;
  bool enq;
  int encsubsegs, encsubfaces;
  int t1ver; 
  int i;

  REAL A[4][4], rhs[4], D;
  int indx[4];
  REAL elen[3];

  encsubsegs = 0;

  if (flag & 1) {
    if (!b->quiet) {
      printf("  Checking conforming property of segments...\n");
    }
    encsubsegs = 0;

    // Run through the list of subsegments, check each one.
    subsegs->traversalinit();
    segloop.sh = shellfacetraverse(subsegs);
    while (segloop.sh != (shellface *) NULL) {
      eorg = (point) segloop.sh[3];
      edest = (point) segloop.sh[4];
      radius = 0.5 * distance(eorg, edest);
      for (i = 0; i < 3; i++) cent[i] = 0.5 * (eorg[i] + edest[i]);

      enq = false;
      sstpivot1(segloop, neightet);
      if (neightet.tet != NULL) {
        spintet = neightet;
        while (1) {
          eapex= apex(spintet);
          if (eapex != dummypoint) {
            dist = distance(eapex, cent);
            diff = dist - radius;
            if (fabs(diff) / radius <= b->epsilon) diff = 0.0; // Rounding.
            if (diff < 0) {
              enq = true; break;
            }
          }
          fnextself(spintet);
          if (spintet.tet == neightet.tet) break;
        }
      }
      if (enq) {
        printf("  !! !! Non-conforming segment: (%d, %d)\n",
               pointmark(eorg), pointmark(edest));
        encsubsegs++;
      }
      segloop.sh = shellfacetraverse(subsegs);
    }

    if (encsubsegs == 0) {
      if (!b->quiet) {
        printf("  The segments are conforming Delaunay.\n");
      }
    } else {
      printf("  !! !! %d subsegments are non-conforming.\n", encsubsegs);
    }
  } // if (flag & 1)

  encsubfaces = 0;

  if (flag & 2) {
    if (!b->quiet) {
      printf("  Checking conforming property of subfaces...\n");
    }

    // Run through the list of subfaces, check each one.
    subfaces->traversalinit();
    shloop.sh = shellfacetraverse(subfaces);
    while (shloop.sh != (shellface *) NULL) {
      pa = (point) shloop.sh[3];
      pb = (point) shloop.sh[4];
      pc = (point) shloop.sh[5];

      // Compute the coefficient matrix A (3x3).
      A[0][0] = pb[0] - pa[0];
      A[0][1] = pb[1] - pa[1];
      A[0][2] = pb[2] - pa[2]; // vector V1 (pa->pb)
      A[1][0] = pc[0] - pa[0];
      A[1][1] = pc[1] - pa[1];
      A[1][2] = pc[2] - pa[2]; // vector V2 (pa->pc)
      cross(A[0], A[1], A[2]); // vector V3 (V1 X V2)

      // Compute the right hand side vector b (3x1).
      elen[0] = dot(A[0], A[0]);
      elen[1] = dot(A[1], A[1]);
      rhs[0] = 0.5 * elen[0];
      rhs[1] = 0.5 * elen[1];
      rhs[2] = 0.0;

      if (lu_decmp(A, 3, indx, &D, 0)) {
        lu_solve(A, 3, indx, rhs, 0);
        cent[0] = pa[0] + rhs[0];
        cent[1] = pa[1] + rhs[1];
        cent[2] = pa[2] + rhs[2];
        rd = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);

        // Check if this subface is encroached.
        for (i = 0; i < 2; i++) {
          stpivot(shloop, searchtet);
          if (!ishulltet(searchtet)) {
            len = distance(oppo(searchtet), cent);
            if ((fabs(len - rd) / rd) < b->epsilon) len = rd; // Rounding.
            if (len < rd) {
              printf("  !! !! Non-conforming subface: (%d, %d, %d)\n",
                     pointmark(pa), pointmark(pb), pointmark(pc));
              encsubfaces++;
              enq = true; break;
            }
          }
          sesymself(shloop);
        }
      }
      shloop.sh = shellfacetraverse(subfaces);
    }

    if (encsubfaces == 0) {
      if (!b->quiet) {
        printf("  The subfaces are conforming Delaunay.\n");
      }
    } else {
      printf("  !! !! %d subfaces are non-conforming.\n", encsubfaces);
    }
  } // if (flag & 2)

  return encsubsegs + encsubfaces;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// qualitystatistics()    Print statistics about the quality of the mesh.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::qualitystatistics()
{
  triface tetloop, neightet;
  point p[4];
  char sbuf[128];
  REAL radiusratiotable[12];
  REAL aspectratiotable[12];
  REAL A[4][4], rhs[4], D;
  REAL V[6][3], N[4][3], H[4]; // edge-vectors, face-normals, face-heights.
  REAL edgelength[6], alldihed[6], faceangle[3];
  REAL shortest, longest;
  REAL smallestvolume, biggestvolume;
  REAL smallestratio, biggestratio;
  REAL smallestdiangle, biggestdiangle;
  REAL smallestfaangle, biggestfaangle;
  REAL total_tet_vol, total_tetprism_vol;
  REAL tetvol, minaltitude;
  REAL cirradius, minheightinv; // insradius;
  REAL shortlen, longlen;
  REAL tetaspect, tetradius;
  REAL smalldiangle, bigdiangle;
  REAL smallfaangle, bigfaangle;
  unsigned long radiustable[12];
  unsigned long aspecttable[16];
  unsigned long dihedangletable[18];
  unsigned long faceangletable[18];
  int indx[4];
  int radiusindex;
  int aspectindex;
  int tendegree;
  int i, j;

  printf("Mesh quality statistics:\n\n");

  shortlen = longlen = 0.0;
  smalldiangle = bigdiangle = 0.0;
  total_tet_vol = 0.0;
  total_tetprism_vol = 0.0;

  radiusratiotable[0]  =    0.707;    radiusratiotable[1]  =     1.0;
  radiusratiotable[2]  =      1.1;    radiusratiotable[3]  =     1.2;
  radiusratiotable[4]  =      1.4;    radiusratiotable[5]  =     1.6;
  radiusratiotable[6]  =      1.8;    radiusratiotable[7]  =     2.0;
  radiusratiotable[8]  =      2.5;    radiusratiotable[9]  =     3.0;
  radiusratiotable[10] =     10.0;    radiusratiotable[11] =     0.0;

  aspectratiotable[0]  =      1.5;    aspectratiotable[1]  =     2.0;
  aspectratiotable[2]  =      2.5;    aspectratiotable[3]  =     3.0;
  aspectratiotable[4]  =      4.0;    aspectratiotable[5]  =     6.0;
  aspectratiotable[6]  =     10.0;    aspectratiotable[7]  =    15.0;
  aspectratiotable[8]  =     25.0;    aspectratiotable[9]  =    50.0;
  aspectratiotable[10] =    100.0;    aspectratiotable[11] =     0.0;
  
  for (i = 0; i < 12; i++) radiustable[i] = 0l;
  for (i = 0; i < 12; i++) aspecttable[i] = 0l;
  for (i = 0; i < 18; i++) dihedangletable[i] = 0l;
  for (i = 0; i < 18; i++) faceangletable[i] = 0l;

  minaltitude = xmax - xmin + ymax - ymin + zmax - zmin;
  minaltitude = minaltitude * minaltitude;
  shortest = minaltitude;
  longest = 0.0;
  smallestvolume = minaltitude;
  biggestvolume = 0.0;
  smallestratio = 1e+16; // minaltitude;
  biggestratio = 0.0;
  smallestdiangle = smallestfaangle = 180.0;
  biggestdiangle = biggestfaangle = 0.0;


  int attrnum = numelemattrib - 1;

  // Loop all elements, calculate quality parameters for each element.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {

    if (b->convex) {
      // Skip tets in the exterior.
      if (elemattribute(tetloop.tet, attrnum) == -1.0) {
        tetloop.tet = tetrahedrontraverse();
        continue;
      }
    }

    // Get four vertices: p0, p1, p2, p3.
    for (i = 0; i < 4; i++) p[i] = (point) tetloop.tet[4 + i];

    // Get the tet volume.
    tetvol = orient3dfast(p[1], p[0], p[2], p[3]) / 6.0;
    total_tet_vol += tetvol;
    total_tetprism_vol += tetprismvol(p[0], p[1], p[2], p[3]);

    // Calculate the largest and smallest volume.
    if (tetvol < smallestvolume) {
      smallestvolume = tetvol;
    } 
    if (tetvol > biggestvolume) {
      biggestvolume = tetvol;
    }

    // Set the edge vectors: V[0], ..., V[5]
    for (i = 0; i < 3; i++) V[0][i] = p[0][i] - p[3][i]; // V[0]: p3->p0.
    for (i = 0; i < 3; i++) V[1][i] = p[1][i] - p[3][i]; // V[1]: p3->p1.
    for (i = 0; i < 3; i++) V[2][i] = p[2][i] - p[3][i]; // V[2]: p3->p2.
    for (i = 0; i < 3; i++) V[3][i] = p[1][i] - p[0][i]; // V[3]: p0->p1.
    for (i = 0; i < 3; i++) V[4][i] = p[2][i] - p[1][i]; // V[4]: p1->p2.
    for (i = 0; i < 3; i++) V[5][i] = p[0][i] - p[2][i]; // V[5]: p2->p0.

    // Get the squares of the edge lengths.
    for (i = 0; i < 6; i++) edgelength[i] = dot(V[i], V[i]);

    // Calculate the longest and shortest edge length.
    for (i = 0; i < 6; i++) {
      if (i == 0) {
        shortlen = longlen = edgelength[i];
      } else {
        shortlen = edgelength[i] < shortlen ? edgelength[i] : shortlen;
        longlen  = edgelength[i] > longlen  ? edgelength[i] : longlen;
      }
      if (edgelength[i] > longest) {
        longest = edgelength[i];
      } 
      if (edgelength[i] < shortest) {
        shortest = edgelength[i];
      }
    }

    // Set the matrix A = [V[0], V[1], V[2]]^T.
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) A[j][i] = V[j][i];
    }

    // Decompose A just once.
    if (lu_decmp(A, 3, indx, &D, 0)) {   
      // Get the three faces normals.
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) rhs[i] = 0.0;
        rhs[j] = 1.0;  // Positive means the inside direction
        lu_solve(A, 3, indx, rhs, 0);
        for (i = 0; i < 3; i++) N[j][i] = rhs[i];
      }
      // Get the fourth face normal by summing up the first three.
      for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
      // Get the radius of the circumsphere.
      for (i = 0; i < 3; i++) rhs[i] = 0.5 * dot(V[i], V[i]);
      lu_solve(A, 3, indx, rhs, 0);
      cirradius = sqrt(dot(rhs, rhs));
      // Normalize the face normals.
      for (i = 0; i < 4; i++) {
        // H[i] is the inverse of height of its corresponding face.
        H[i] = sqrt(dot(N[i], N[i]));
        for (j = 0; j < 3; j++) N[i][j] /= H[i];
      }
      // Get the radius of the inscribed sphere.
      // insradius = 1.0 / (H[0] + H[1] + H[2] + H[3]);
      // Get the biggest H[i] (corresponding to the smallest height).
      minheightinv = H[0];
      for (i = 1; i < 3; i++) {
        if (H[i] > minheightinv) minheightinv = H[i];
      }
    } else {
      // A nearly degenerated tet.
      if (tetvol <= 0.0) {
        // assert(tetvol != 0.0);
        printf("  !! Warning:  A %s tet (%d,%d,%d,%d).\n", 
               tetvol < 0 ? "inverted" : "degenerated", pointmark(p[0]),
               pointmark(p[1]), pointmark(p[2]), pointmark(p[3]));
        // Skip it.        
        tetloop.tet = tetrahedrontraverse();
        continue;
      }
      // Calculate the four face normals.
      facenormal(p[2], p[1], p[3], N[0], 1, NULL);
      facenormal(p[0], p[2], p[3], N[1], 1, NULL);
      facenormal(p[1], p[0], p[3], N[2], 1, NULL);
      facenormal(p[0], p[1], p[2], N[3], 1, NULL);
      // Normalize the face normals.
      for (i = 0; i < 4; i++) {
        // H[i] is the twice of the area of the face.
        H[i] = sqrt(dot(N[i], N[i]));
        for (j = 0; j < 3; j++) N[i][j] /= H[i];
      }
      // Get the biggest H[i] / tetvol (corresponding to the smallest height).
      minheightinv = (H[0] / tetvol);
      for (i = 1; i < 3; i++) {
        if ((H[i] / tetvol) > minheightinv) minheightinv = (H[i] / tetvol);
      }
      // Let the circumradius to be the half of its longest edge length.
      cirradius = 0.5 * sqrt(longlen);
    }

    // Get the dihedrals (in degree) at each edges.
    j = 0;
    for (i = 1; i < 4; i++) {
      alldihed[j] = -dot(N[0], N[i]); // Edge cd, bd, bc.
      if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
      else if (alldihed[j] > 1.0) alldihed[j] = 1;
      alldihed[j] = acos(alldihed[j]) / PI * 180.0;
      j++;
    }
    for (i = 2; i < 4; i++) {
      alldihed[j] = -dot(N[1], N[i]); // Edge ad, ac.
      if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
      else if (alldihed[j] > 1.0) alldihed[j] = 1;
      alldihed[j] = acos(alldihed[j]) / PI * 180.0;
      j++;
    }
    alldihed[j] = -dot(N[2], N[3]); // Edge ab.
    if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
    else if (alldihed[j] > 1.0) alldihed[j] = 1;
    alldihed[j] = acos(alldihed[j]) / PI * 180.0;

    // Calculate the largest and smallest dihedral angles.
    for (i = 0; i < 6; i++) {
      if (i == 0) {
        smalldiangle = bigdiangle = alldihed[i];
      } else {
        smalldiangle = alldihed[i] < smalldiangle ? alldihed[i] : smalldiangle;
        bigdiangle = alldihed[i] > bigdiangle ? alldihed[i] : bigdiangle;
      }
      if (alldihed[i] < smallestdiangle) {
        smallestdiangle = alldihed[i];
      } 
      if (alldihed[i] > biggestdiangle) {
        biggestdiangle = alldihed[i];
      }
      // Accumulate the corresponding number in the dihedral angle histogram.
      if (alldihed[i] < 5.0) {
        tendegree = 0;
      } else if (alldihed[i] >= 5.0 && alldihed[i] < 10.0) {
        tendegree = 1;
      } else if (alldihed[i] >= 80.0 && alldihed[i] < 110.0) {
        tendegree = 9; // Angles between 80 to 110 degree are in one entry.
      } else if (alldihed[i] >= 170.0 && alldihed[i] < 175.0) {
        tendegree = 16;
      } else if (alldihed[i] >= 175.0) {
        tendegree = 17;
      } else {
        tendegree = (int) (alldihed[i] / 10.);
        if (alldihed[i] < 80.0) {
          tendegree++;  // In the left column.
        } else {
          tendegree--;  // In the right column.
        }
      }
      dihedangletable[tendegree]++;
    }



    // Calculate the largest and smallest face angles.
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, neightet);
      // Only do the calulation once for a face.
      if (((point) neightet.tet[7] == dummypoint) || 
          (tetloop.tet < neightet.tet)) {
        p[0] = org(tetloop);
        p[1] = dest(tetloop);
        p[2] = apex(tetloop);
        faceangle[0] = interiorangle(p[0], p[1], p[2], NULL);
        faceangle[1] = interiorangle(p[1], p[2], p[0], NULL);
        faceangle[2] = PI - (faceangle[0] + faceangle[1]);
        // Translate angles into degrees.
        for (i = 0; i < 3; i++) {
          faceangle[i] = (faceangle[i] * 180.0) / PI;
        }
        // Calculate the largest and smallest face angles.
        for (i = 0; i < 3; i++) {
          if (i == 0) {
            smallfaangle = bigfaangle = faceangle[i];
          } else {
            smallfaangle = faceangle[i] < smallfaangle ? 
              faceangle[i] : smallfaangle;
            bigfaangle = faceangle[i] > bigfaangle ? faceangle[i] : bigfaangle;
          }
          if (faceangle[i] < smallestfaangle) {
            smallestfaangle = faceangle[i];
          } 
          if (faceangle[i] > biggestfaangle) {
            biggestfaangle = faceangle[i];
          }
          tendegree = (int) (faceangle[i] / 10.);
          faceangletable[tendegree]++;
        }
      }
    }

    // Calculate aspect ratio and radius-edge ratio for this element.
    tetradius = cirradius / sqrt(shortlen);
    // tetaspect = sqrt(longlen) / (2.0 * insradius);
    tetaspect = sqrt(longlen) * minheightinv;
    // Remember the largest and smallest aspect ratio.
    if (tetaspect < smallestratio) {
      smallestratio = tetaspect;
    } 
    if (tetaspect > biggestratio) {
      biggestratio = tetaspect;
    }
    // Accumulate the corresponding number in the aspect ratio histogram.
    aspectindex = 0;
    while ((tetaspect > aspectratiotable[aspectindex]) && (aspectindex < 11)) {
      aspectindex++;
    }
    aspecttable[aspectindex]++;
    radiusindex = 0;
    while ((tetradius > radiusratiotable[radiusindex]) && (radiusindex < 11)) {
      radiusindex++;
    }
    radiustable[radiusindex]++;

    tetloop.tet = tetrahedrontraverse();
  }

  shortest = sqrt(shortest);
  longest = sqrt(longest);
  minaltitude = sqrt(minaltitude);

  printf("  Smallest volume: %16.5g   |  Largest volume: %16.5g\n",
         smallestvolume, biggestvolume);
  printf("  Shortest edge:   %16.5g   |  Longest edge:   %16.5g\n",
         shortest, longest);
  printf("  Smallest asp.ratio: %13.5g   |  Largest asp.ratio: %13.5g\n",
         smallestratio, biggestratio);
  sprintf(sbuf, "%.17g", biggestfaangle);
  if (strlen(sbuf) > 8) {
    sbuf[8] = '\0';
  }
  printf("  Smallest facangle: %14.5g   |  Largest facangle:       %s\n",
         smallestfaangle, sbuf);
  sprintf(sbuf, "%.17g", biggestdiangle);
  if (strlen(sbuf) > 8) {
    sbuf[8] = '\0';
  }
  printf("  Smallest dihedral: %14.5g   |  Largest dihedral:       %s\n\n",
         smallestdiangle, sbuf);

  printf("  Aspect ratio histogram:\n");
  printf("         < %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
         aspectratiotable[0], aspecttable[0], aspectratiotable[5],
         aspectratiotable[6], aspecttable[6]);
  for (i = 1; i < 5; i++) {
    printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
           aspectratiotable[i - 1], aspectratiotable[i], aspecttable[i],
           aspectratiotable[i + 5], aspectratiotable[i + 6],
           aspecttable[i + 6]);
  }
  printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g -            :  %8ld\n",
         aspectratiotable[4], aspectratiotable[5], aspecttable[5],
         aspectratiotable[10], aspecttable[11]);
  printf("  (A tetrahedron's aspect ratio is its longest edge length");
  printf(" divided by its\n");
  printf("    smallest side height)\n\n");

  printf("  Face angle histogram:\n");
  for (i = 0; i < 9; i++) {
    printf("    %3d - %3d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
           i * 10, i * 10 + 10, faceangletable[i],
           i * 10 + 90, i * 10 + 100, faceangletable[i + 9]);
  }
  if (minfaceang != PI) {
    printf("  Minimum input face angle is %g (degree).\n",
           minfaceang / PI * 180.0);
  }
  printf("\n");

  printf("  Dihedral angle histogram:\n");
  // Print the three two rows:
  printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
         0, 5, dihedangletable[0], 80, 110, dihedangletable[9]);
  printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
         5, 10, dihedangletable[1], 110, 120, dihedangletable[10]);
  // Print the third to seventh rows.
  for (i = 2; i < 7; i++) {
    printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
           (i - 1) * 10, (i - 1) * 10 + 10, dihedangletable[i],
           (i - 1) * 10 + 110, (i - 1) * 10 + 120, dihedangletable[i + 9]);
  }
  // Print the last two rows.
  printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
         60, 70, dihedangletable[7], 170, 175, dihedangletable[16]);
  printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
         70, 80, dihedangletable[8], 175, 180, dihedangletable[17]);
  if (minfacetdihed != PI) {
    printf("  Minimum input dihedral angle is %g (degree).\n",
           minfacetdihed / PI * 180.0);
  }
  printf("\n");

  printf("\n");
}


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// memorystatistics()    Report the memory usage.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::memorystatistics()
{
  printf("Memory usage statistics:\n\n");
 
  // Count the number of blocks of tetrahedra. 
  int tetblocks = 0;
  tetrahedrons->pathblock = tetrahedrons->firstblock;
  while (tetrahedrons->pathblock != NULL) {
    tetblocks++;
    tetrahedrons->pathblock = (void **) *(tetrahedrons->pathblock);  
  }

  // Calculate the total memory (in bytes) used by storing meshes.
  unsigned long totalmeshmemory = 0l, totalt2shmemory = 0l;
  totalmeshmemory = points->maxitems * points->itembytes +
                    tetrahedrons->maxitems * tetrahedrons->itembytes;
  if (b->plc || b->refine) {
    totalmeshmemory += (subfaces->maxitems * subfaces->itembytes +
                        subsegs->maxitems * subsegs->itembytes);
    totalt2shmemory = (tet2subpool->maxitems * tet2subpool->itembytes +
                       tet2segpool->maxitems * tet2segpool->itembytes);
  }

  unsigned long totalalgomemory = 0l;
  totalalgomemory = cavetetlist->totalmemory + cavebdrylist->totalmemory +
                    caveoldtetlist->totalmemory + 
                    flippool->maxitems * flippool->itembytes;
  if (b->plc || b->refine) {
    totalalgomemory += (subsegstack->totalmemory + subfacstack->totalmemory +
                        subvertstack->totalmemory + 
                        caveshlist->totalmemory + caveshbdlist->totalmemory +
                        cavesegshlist->totalmemory +
                        cavetetshlist->totalmemory + 
                        cavetetseglist->totalmemory +
                        caveencshlist->totalmemory +
                        caveencseglist->totalmemory +
                        cavetetvertlist->totalmemory +
                        unflipqueue->totalmemory);
  }

  printf("  Maximum number of tetrahedra:  %ld\n", tetrahedrons->maxitems);
  printf("  Maximum number of tet blocks (blocksize = %d):  %d\n",
         b->tetrahedraperblock, tetblocks);
  /*
  if (b->plc || b->refine) {
    printf("  Approximate memory for tetrahedral mesh (bytes):  %ld\n",
           totalmeshmemory);
    
    printf("  Approximate memory for extra pointers (bytes):  %ld\n",
           totalt2shmemory);
  } else {
    printf("  Approximate memory for tetrahedralization (bytes):  %ld\n",
           totalmeshmemory);
  }
  printf("  Approximate memory for algorithms (bytes):  %ld\n",
         totalalgomemory);
  printf("  Approximate memory for working arrays (bytes):  %ld\n",
         totalworkmemory);
  printf("  Approximate total used memory (bytes):  %ld\n",
         totalmeshmemory + totalt2shmemory + totalalgomemory + 
         totalworkmemory);
  */
  if (b->plc || b->refine) {
    printf("  Approximate memory for tetrahedral mesh (bytes):  ");
    printfcomma(totalmeshmemory); printf("\n");
    
    printf("  Approximate memory for extra pointers (bytes):  ");
    printfcomma(totalt2shmemory); printf("\n");
  } else {
    printf("  Approximate memory for tetrahedralization (bytes):  ");
    printfcomma(totalmeshmemory); printf("\n");
  }
  printf("  Approximate memory for algorithms (bytes):  ");
  printfcomma(totalalgomemory); printf("\n");
  printf("  Approximate memory for working arrays (bytes):  ");
  printfcomma(totalworkmemory); printf("\n");
  printf("  Approximate total used memory (bytes):  ");
  printfcomma(totalmeshmemory + totalt2shmemory + totalalgomemory + 
              totalworkmemory);
  printf("\n");

  printf("\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// statistics()    Print all sorts of cool facts.                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::statistics()
{
  long tetnumber, facenumber;

  printf("\nStatistics:\n\n");
  printf("  Input points: %d\n", in->numberofpoints);
  if (b->refine) {
    printf("  Input tetrahedra: %d\n", in->numberoftetrahedra);
  }
  if (b->plc) {
    printf("  Input facets: %d\n", in->numberoffacets);
    printf("  Input segments: %ld\n", insegments);
    printf("  Input holes: %d\n", in->numberofholes);
    printf("  Input regions: %d\n", in->numberofregions);
  }

  tetnumber = tetrahedrons->items - hullsize;
  facenumber = (tetnumber * 4l + hullsize) / 2l;

  if (b->weighted) { // -w option
    printf("\n  Mesh points: %ld\n", points->items - nonregularcount);
  } else {
    printf("\n  Mesh points: %ld\n", points->items);
  }
  printf("  Mesh tetrahedra: %ld\n", tetnumber);
  printf("  Mesh faces: %ld\n", facenumber);
  if (meshedges > 0l) {
    printf("  Mesh edges: %ld\n", meshedges);
  } else {
    if (!nonconvex) {
      long vsize = points->items - dupverts - unuverts;
      if (b->weighted) vsize -= nonregularcount;
      meshedges = vsize + facenumber - tetnumber - 1;
      printf("  Mesh edges: %ld\n", meshedges);
    }
  }

  if (b->plc || b->refine) {
    printf("  Mesh faces on facets: %ld\n", subfaces->items);
    printf("  Mesh edges on segments: %ld\n", subsegs->items);
    if (st_volref_count > 0l) {
      printf("  Steiner points inside domain: %ld\n", st_volref_count);
    }
    if (st_facref_count > 0l) {
      printf("  Steiner points on facets:  %ld\n", st_facref_count);
    }
    if (st_segref_count > 0l) {
      printf("  Steiner points on segments:  %ld\n", st_segref_count);
    }
  } else {
    printf("  Convex hull faces: %ld\n", hullsize);
    if (meshhulledges > 0l) {
      printf("  Convex hull edges: %ld\n", meshhulledges);
    }
  }
  if (b->weighted) { // -w option
    printf("  Skipped non-regular points: %ld\n", nonregularcount);
  }
  printf("\n");


  if (b->verbose > 0) {
    if (b->plc || b->refine) { // -p or -r
      if (tetrahedrons->items > 0l) {
        qualitystatistics();
      }
    }
    if (tetrahedrons->items > 0l) {
      memorystatistics();
    }
  }
}

////                                                                       ////
////                                                                       ////
//// meshstat_cxx /////////////////////////////////////////////////////////////

