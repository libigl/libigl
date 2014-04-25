#include "../tetgen.h"
//// main_cxx /////////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    The interface for users using TetGen library to       //
//                     generate tetrahedral meshes with all features.        //
//                                                                           //
// The sequence is roughly as follows.  Many of these steps can be skipped,  //
// depending on the command line switches.                                   //
//                                                                           //
// - Initialize constants and parse the command line.                        //
// - Read the vertices from a file and either                                //
//   - tetrahedralize them (no -r), or                                       //
//   - read an old mesh from files and reconstruct it (-r).                  //
// - Insert the boundary segments and facets (-p or -Y).                     //
// - Read the holes (-p), regional attributes (-pA), and regional volume     //
//   constraints (-pa).  Carve the holes and concavities, and spread the     //
//   regional attributes and volume constraints.                             //
// - Enforce the constraints on minimum quality bound (-q) and maximum       //
//   volume (-a), and a mesh size function (-m).                             //
// - Optimize the mesh wrt. specified quality measures (-O and -o).          //
// - Write the output files and print the statistics.                        //
// - Check the consistency of the mesh (-C).                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(tetgenbehavior *b, tetgenio *in, tetgenio *out,
                    tetgenio *addin, tetgenio *bgmin)
{
  tetgenmesh m;
  clock_t tv[12], ts[5]; // Timing informations (defined in time.h)
  REAL cps = (REAL) CLOCKS_PER_SEC;

  tv[0] = clock();
 
  m.b = b;
  m.in = in;
  m.addin = addin;

  if (b->metric && bgmin && (bgmin->numberofpoints > 0)) {
    m.bgm = new tetgenmesh(); // Create an empty background mesh.
    m.bgm->b = b;
    m.bgm->in = bgmin;
  }

  m.initializepools();
  m.transfernodes();

  exactinit(b->verbose, b->noexact, b->nostaticfilter,
            m.xmax - m.xmin, m.ymax - m.ymin, m.zmax - m.zmin);

  tv[1] = clock();

  if (b->refine) { // -r
    m.reconstructmesh();
  } else { // -p
    m.incrementaldelaunay(ts[0]);
  }

  tv[2] = clock();

  if (!b->quiet) {
    if (b->refine) {
      printf("Mesh reconstruction seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);
    } else {
      printf("Delaunay seconds:  %g\n", ((REAL)(tv[2]-tv[1])) / cps);
      if (b->verbose) {
        printf("  Point sorting seconds:  %g\n", ((REAL)(ts[0]-tv[1])) / cps);
      }
    }
  }

  if (b->plc && !b->refine) { // -p
    m.meshsurface();

    ts[0] = clock();

    if (!b->quiet) {
      printf("Surface mesh seconds:  %g\n", ((REAL)(ts[0]-tv[2])) / cps);
    }

    if (b->diagnose) { // -d
      m.detectinterfaces();

      ts[1] = clock();

      if (!b->quiet) {
        printf("Self-intersection seconds:  %g\n", ((REAL)(ts[1]-ts[0])) / cps);
      }

      // Only output when self-intersecting faces exist.
      if (m.subfaces->items > 0l) {
        m.outnodes(out);
        m.outsubfaces(out);
      }

      return;
    }
  }

  tv[3] = clock();

  if ((b->metric) && (m.bgm != NULL)) { // -m
    m.bgm->initializepools();
    m.bgm->transfernodes();
    m.bgm->reconstructmesh();

    ts[0] = clock();

    if (!b->quiet) {
      printf("Background mesh reconstruct seconds:  %g\n",
             ((REAL)(ts[0] - tv[3])) / cps);
    }

    if (b->metric) { // -m
      m.interpolatemeshsize();

      ts[1] = clock();

      if (!b->quiet) {
        printf("Size interpolating seconds:  %g\n",((REAL)(ts[1]-ts[0])) / cps);
      }
    }
  }

  tv[4] = clock();

  if (b->plc && !b->refine) { // -p
    if (b->nobisect) { // -Y
      m.recoverboundary(ts[0]);
    } else {
      m.constraineddelaunay(ts[0]);
    }

    ts[1] = clock();

    if (!b->quiet) {
      if (b->nobisect) {
        printf("Boundary recovery ");
      } else {
        printf("Constrained Delaunay ");
      }
      printf("seconds:  %g\n", ((REAL)(ts[1] - tv[4])) / cps);
      if (b->verbose) {
        printf("  Segment recovery seconds:  %g\n",((REAL)(ts[0]-tv[4]))/ cps);
        printf("  Facet recovery seconds:  %g\n", ((REAL)(ts[1]-ts[0])) / cps);
      }
    }

    m.carveholes();

    ts[2] = clock();

    if (!b->quiet) {
      printf("Exterior tets removal seconds:  %g\n",((REAL)(ts[2]-ts[1]))/cps);
    }

    if (b->nobisect) { // -Y
      if (m.subvertstack->objects > 0l) {
        m.suppresssteinerpoints();

        ts[3] = clock();

        if (!b->quiet) {
          printf("Steiner suppression seconds:  %g\n",
                 ((REAL)(ts[3]-ts[2]))/cps);
        }
      }
    }
  }

  tv[5] = clock();

  if (b->coarsen) { // -R
    m.meshcoarsening();
  }

  tv[6] = clock();

  if (!b->quiet) {
    if (b->coarsen) {
      printf("Mesh coarsening seconds:  %g\n", ((REAL)(tv[6] - tv[5])) / cps);
    }
  }

  if ((b->plc && b->nobisect) || b->coarsen) {
    m.recoverdelaunay();
  }

  tv[7] = clock();

  if (!b->quiet) {
    if ((b->plc && b->nobisect) || b->coarsen) {
      printf("Delaunay recovery seconds:  %g\n", ((REAL)(tv[7] - tv[6]))/cps);
    }
  }

  if ((b->plc || b->refine) && b->insertaddpoints) { // -i
    if ((addin != NULL) && (addin->numberofpoints > 0)) {
      m.insertconstrainedpoints(addin); 
    }
  }

  tv[8] = clock();

  if (!b->quiet) {
    if ((b->plc || b->refine) && b->insertaddpoints) { // -i
      if ((addin != NULL) && (addin->numberofpoints > 0)) {
        printf("Constrained points seconds:  %g\n", ((REAL)(tv[8]-tv[7]))/cps);
      }
    }
  }

  if (b->quality) {
    m.delaunayrefinement();    
  }

  tv[9] = clock();

  if (!b->quiet) {
    if (b->quality) {
      printf("Refinement seconds:  %g\n", ((REAL)(tv[9] - tv[8])) / cps);
    }
  }

  if ((b->plc || b->refine) && (b->optlevel > 0)) {
    m.optimizemesh();
  }

  tv[10] = clock();

  if (!b->quiet) {
    if ((b->plc || b->refine) && (b->optlevel > 0)) {
      printf("Optimization seconds:  %g\n", ((REAL)(tv[10] - tv[9])) / cps);
    }
  }

  if (!b->nojettison && ((m.dupverts > 0) || (m.unuverts > 0)
      || (b->refine && (in->numberofcorners == 10)))) {
    m.jettisonnodes();
  }

  if ((b->order == 2) && !b->convex) {
    m.highorder();
  }

  if (!b->quiet) {
    printf("\n");
  }

  if (out != (tetgenio *) NULL) {
    out->firstnumber = in->firstnumber;
    out->mesh_dim = in->mesh_dim;
  }

  if (b->nonodewritten || b->noiterationnum) {
    if (!b->quiet) {
      printf("NOT writing a .node file.\n");
    }
  } else {
    m.outnodes(out);
  }

  if (b->noelewritten) {
    if (!b->quiet) {
      printf("NOT writing an .ele file.\n");
    }
  } else {
    if (m.tetrahedrons->items > 0l) {
      m.outelements(out);
    }
  }

  if (b->nofacewritten) {
    if (!b->quiet) {
      printf("NOT writing an .face file.\n");
    }
  } else {
    if (b->facesout) {
      if (m.tetrahedrons->items > 0l) {
        m.outfaces(out);  // Output all faces.
      }
    } else {
      if (b->plc || b->refine) {
        if (m.subfaces->items > 0l) {
          m.outsubfaces(out); // Output boundary faces.
        }
      } else {
        if (m.tetrahedrons->items > 0l) {
          m.outhullfaces(out); // Output convex hull faces.
        }
      }
    }
  }


  if (b->nofacewritten) {
    if (!b->quiet) {
      printf("NOT writing an .edge file.\n");
    }
  } else {
    if (b->edgesout) { // -e
      m.outedges(out); // output all mesh edges. 
    } else {
      if (b->plc || b->refine) {
        m.outsubsegments(out); // output subsegments.
      }
    }
  }

  if ((b->plc || b->refine) && b->metric) { // -m
    m.outmetrics(out);
  }

  if (!out && b->plc && 
      ((b->object == tetgenbehavior::OFF) ||
       (b->object == tetgenbehavior::PLY) ||
       (b->object == tetgenbehavior::STL))) {
    m.outsmesh(b->outfilename);
  }

  if (!out && b->meditview) {
    m.outmesh2medit(b->outfilename); 
  }


  if (!out && b->vtkview) {
    m.outmesh2vtk(b->outfilename); 
  }

  if (b->neighout) {
    m.outneighbors(out);
  }

  if ((!(b->plc || b->refine)) && b->voroout) {
    m.outvoronoi(out);
  }


  tv[11] = clock();

  if (!b->quiet) {
    printf("\nOutput seconds:  %g\n", ((REAL)(tv[11] - tv[10])) / cps);
    printf("Total running seconds:  %g\n", ((REAL)(tv[11] - tv[0])) / cps);
  }

  if (b->docheck) {
    m.checkmesh(0);
    if (b->plc || b->refine) {
      m.checkshells();
      m.checksegments();
    }
    if (b->docheck > 1) {
      m.checkdelaunay();
    }
  }

  if (!b->quiet) {
    m.statistics();
  }
}

#ifndef TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// main()    The command line interface of TetGen.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])

#else // with TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    The library interface of TetGen.                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(char *switches, tetgenio *in, tetgenio *out, 
                    tetgenio *addin, tetgenio *bgmin)

#endif // not TETLIBRARY

{
  tetgenbehavior b;

#ifndef TETLIBRARY

  tetgenio in, addin, bgmin;
  
  if (!b.parse_commandline(argc, argv)) {
    terminatetetgen(NULL, 10);
  }

  // Read input files.
  if (b.refine) { // -r
    if (!in.load_tetmesh(b.infilename, (int) b.object)) {
      terminatetetgen(NULL, 10);
    }
  } else { // -p
    if (!in.load_plc(b.infilename, (int) b.object)) {
      terminatetetgen(NULL, 10);
    }
  }
  if (b.insertaddpoints) { // -i
    // Try to read a .a.node file.
    addin.load_node(b.addinfilename);
  }
  if (b.metric) { // -m
    // Try to read a background mesh in files .b.node, .b.ele.
    bgmin.load_tetmesh(b.bgmeshfilename, (int) b.object);
  }

  tetrahedralize(&b, &in, NULL, &addin, &bgmin);

  return 0;

#else // with TETLIBRARY

  if (!b.parse_commandline(switches)) {
    terminatetetgen(NULL, 10);
  }
  tetrahedralize(&b, in, out, addin, bgmin);

#endif // not TETLIBRARY
}

////                                                                       ////
////                                                                       ////
//// main_cxx /////////////////////////////////////////////////////////////////

