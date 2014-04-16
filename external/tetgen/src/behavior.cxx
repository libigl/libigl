#include "../tetgen.h"
//// behavior_cxx /////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// syntax()    Print list of command line switches.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenbehavior::syntax()
{
  printf("  tetgen [-pYrq_Aa_miO_S_T_XMwcdzfenvgkJBNEFICQVh] input_file\n");
  printf("    -p  Tetrahedralizes a piecewise linear complex (PLC).\n");
  printf("    -Y  Preserves the input surface mesh (does not modify it).\n");
  printf("    -r  Reconstructs a previously generated mesh.\n");
  printf("    -q  Refines mesh (to improve mesh quality).\n");
  printf("    -R  Mesh coarsening (to reduce the mesh elements).\n");
  printf("    -A  Assigns attributes to tetrahedra in different regions.\n");
  printf("    -a  Applies a maximum tetrahedron volume constraint.\n");
  printf("    -m  Applies a mesh sizing function.\n");
  printf("    -i  Inserts a list of additional points.\n");
  printf("    -O  Specifies the level of mesh optimization.\n");
  printf("    -S  Specifies maximum number of added points.\n");
  printf("    -T  Sets a tolerance for coplanar test (default 1e-8).\n");
  printf("    -X  Suppresses use of exact arithmetic.\n");
  printf("    -M  No merge of coplanar facets or very close vertices.\n");
  printf("    -w  Generates weighted Delaunay (regular) triangulation.\n");
  printf("    -c  Retains the convex hull of the PLC.\n");
  printf("    -d  Detects self-intersections of facets of the PLC.\n");
  printf("    -z  Numbers all output items starting from zero.\n");
  printf("    -f  Outputs all faces to .face file.\n");
  printf("    -e  Outputs all edges to .edge file.\n");
  printf("    -n  Outputs tetrahedra neighbors to .neigh file.\n");
  printf("    -v  Outputs Voronoi diagram to files.\n");
  printf("    -g  Outputs mesh to .mesh file for viewing by Medit.\n");
  printf("    -k  Outputs mesh to .vtk file for viewing by Paraview.\n");
  printf("    -J  No jettison of unused vertices from output .node file.\n");
  printf("    -B  Suppresses output of boundary information.\n");
  printf("    -N  Suppresses output of .node file.\n");
  printf("    -E  Suppresses output of .ele file.\n");
  printf("    -F  Suppresses output of .face and .edge file.\n");
  printf("    -I  Suppresses mesh iteration numbers.\n");
  printf("    -C  Checks the consistency of the final mesh.\n");
  printf("    -Q  Quiet:  No terminal output except errors.\n");
  printf("    -V  Verbose:  Detailed information, more terminal output.\n");
  printf("    -h  Help:  A brief instruction for using TetGen.\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// usage()    Print a brief instruction for using TetGen.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenbehavior::usage()
{
  printf("TetGen\n");
  printf("A Quality Tetrahedral Mesh Generator and 3D Delaunay ");
  printf("Triangulator\n");
  printf("Version 1.5\n");
  printf("November 4, 2013\n");
  printf("\n");
  printf("What Can TetGen Do?\n");
  printf("\n");
  printf("  TetGen generates Delaunay tetrahedralizations, constrained\n");
  printf("  Delaunay tetrahedralizations, and quality tetrahedral meshes.\n");
  printf("\n");
  printf("Command Line Syntax:\n");
  printf("\n");
  printf("  Below is the basic command line syntax of TetGen with a list of ");
  printf("short\n");
  printf("  descriptions. Underscores indicate that numbers may optionally\n");
  printf("  follow certain switches.  Do not leave any space between a ");
  printf("switch\n");
  printf("  and its numeric parameter.  \'input_file\' contains input data\n");
  printf("  depending on the switches you supplied which may be a ");
  printf("  piecewise\n");
  printf("  linear complex or a list of nodes.  File formats and detailed\n");
  printf("  description of command line switches are found in user's ");
  printf("manual.\n");
  printf("\n");
  syntax();
  printf("\n");
  printf("Examples of How to Use TetGen:\n");
  printf("\n");
  printf("  \'tetgen object\' reads vertices from object.node, and writes ");
  printf("their\n  Delaunay tetrahedralization to object.1.node, ");
  printf("object.1.ele\n  (tetrahedra), and object.1.face");
  printf(" (convex hull faces).\n");
  printf("\n");
  printf("  \'tetgen -p object\' reads a PLC from object.poly or object.");
  printf("smesh (and\n  possibly object.node) and writes its constrained ");
  printf("Delaunay\n  tetrahedralization to object.1.node, object.1.ele, ");
  printf("object.1.face,\n");
  printf("  (boundary faces) and object.1.edge (boundary edges).\n");
  printf("\n");
  printf("  \'tetgen -pq1.414a.1 object\' reads a PLC from object.poly or\n");
  printf("  object.smesh (and possibly object.node), generates a mesh ");
  printf("whose\n  tetrahedra have radius-edge ratio smaller than 1.414 and ");
  printf("have volume\n  of 0.1 or less, and writes the mesh to ");
  printf("object.1.node, object.1.ele,\n  object.1.face, and object.1.edge\n");
  printf("\n");
  printf("Please send bugs/comments to Hang Si <si@wias-berlin.de>\n");
  terminatetetgen(NULL, 0);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// parse_commandline()    Read the command line, identify switches, and set  //
//                        up options and file names.                         //
//                                                                           //
// 'argc' and 'argv' are the same parameters passed to the function main()   //
// of a C/C++ program. They together represent the command line user invoked //
// from an environment in which TetGen is running.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenbehavior::parse_commandline(int argc, char **argv)
{
  int startindex;
  int increment;
  int meshnumber;
  int i, j, k;
  char workstring[1024];

  // First determine the input style of the switches.
  if (argc == 0) {
    startindex = 0;                    // Switches are given without a dash.
    argc = 1;                    // For running the following for-loop once.
    commandline[0] = '\0';
  } else {
    startindex = 1;
    strcpy(commandline, argv[0]);
    strcat(commandline, " ");
  }

  for (i = startindex; i < argc; i++) {
    // Remember the command line for output.
    strcat(commandline, argv[i]);
    strcat(commandline, " ");
    if (startindex == 1) {
      // Is this string a filename?
      if (argv[i][0] != '-') {
        strncpy(infilename, argv[i], 1024 - 1);
        infilename[1024 - 1] = '\0';
        continue;                     
      }
    }
    // Parse the individual switch from the string.
    for (j = startindex; argv[i][j] != '\0'; j++) {
      if (argv[i][j] == 'p') {
        plc = 1;
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          facet_ang_tol = (REAL) strtod(workstring, (char **) NULL);
        }
      } else if (argv[i][j] == 's') {
        psc = 1;
      } else if (argv[i][j] == 'Y') {
        nobisect = 1;
        if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
          nobisect_param = (argv[i][j + 1] - '0');
          j++;
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
            addsteiner_algo = (argv[i][j + 1] - '0');
            j++;
          }
        }
      } else if (argv[i][j] == 'r') {
        refine = 1;
      } else if (argv[i][j] == 'q') {
        quality = 1;
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          minratio = (REAL) strtod(workstring, (char **) NULL);
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            mindihedral = (REAL) strtod(workstring, (char **) NULL);
          }
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            optmaxdihedral = (REAL) strtod(workstring, (char **) NULL);
          }
        }
      } else if (argv[i][j] == 'R') {
        coarsen = 1;
        if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
          coarsen_param = (argv[i][j + 1] - '0');
          j++;
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            coarsen_percent = (REAL) strtod(workstring, (char **) NULL);
          }
        }
      } else if (argv[i][j] == 'w') {
        weighted = 1;
        if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
          weighted_param = (argv[i][j + 1] - '0');
          j++;
        }
      } else if (argv[i][j] == 'b') {
        // -b(brio_threshold/brio_ratio/hilbert_limit/hilbert_order)
        brio_hilbert = 1;
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          brio_threshold = (int) strtol(workstring, (char **) &workstring, 0);
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            brio_ratio = (REAL) strtod(workstring, (char **) NULL);
          }
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.') || (argv[i][j + 1] == '-')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.') || (argv[i][j + 1] == '-')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            hilbert_limit = (int) strtol(workstring, (char **) &workstring, 0);
          }
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.') || (argv[i][j + 1] == '-')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.') || (argv[i][j + 1] == '-')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            hilbert_order = (REAL) strtod(workstring, (char **) NULL);
          }
        }
        if (brio_threshold == 0) { // -b0
          brio_hilbert = 0; // Turn off BRIO-Hilbert sorting. 
        }
        if (brio_ratio >= 1.0) { // -b/1
          no_sort = 1;
          brio_hilbert = 0; // Turn off BRIO-Hilbert sorting.
        }
      } else if (argv[i][j] == 'l') {
        incrflip = 1;
      } else if (argv[i][j] == 'L') {
        flipinsert = 1;
      } else if (argv[i][j] == 'm') {
        metric = 1;
      } else if (argv[i][j] == 'a') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          fixedvolume = 1;
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          maxvolume = (REAL) strtod(workstring, (char **) NULL);
        } else {
          varvolume = 1;
        }
      } else if (argv[i][j] == 'A') {
        regionattrib = 1;
      } else if (argv[i][j] == 'D') {
        conforming = 1;
        if ((argv[i][j + 1] >= '1') && (argv[i][j + 1] <= '3')) {
          reflevel = (argv[i][j + 1] - '1') + 1; 
          j++;
        }
      } else if (argv[i][j] == 'i') {
        insertaddpoints = 1;
      } else if (argv[i][j] == 'd') {
        diagnose = 1;
      } else if (argv[i][j] == 'c') {
        convex = 1;
      } else if (argv[i][j] == 'M') {
        nomergefacet = 1;
        nomergevertex = 1;
        if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '1')) {
          nomergefacet = (argv[i][j + 1] - '0');
          j++;
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '1')) {
            nomergevertex = (argv[i][j + 1] - '0');
            j++;
          }
        }
      } else if (argv[i][j] == 'X') {
        if (argv[i][j + 1] == '1') {
          nostaticfilter = 1;
          j++;
        } else {
          noexact = 1;
        }
      } else if (argv[i][j] == 'z') {
        zeroindex = 1;
      } else if (argv[i][j] == 'f') {
        facesout++;
      } else if (argv[i][j] == 'e') {
        edgesout++;
      } else if (argv[i][j] == 'n') {
        neighout++;
      } else if (argv[i][j] == 'v') {
        voroout = 1;
      } else if (argv[i][j] == 'g') {
        meditview = 1;
      } else if (argv[i][j] == 'k') {
        vtkview = 1;  
      } else if (argv[i][j] == 'J') {
        nojettison = 1;
      } else if (argv[i][j] == 'B') {
        nobound = 1;
      } else if (argv[i][j] == 'N') {
        nonodewritten = 1;
      } else if (argv[i][j] == 'E') {
        noelewritten = 1;
      } else if (argv[i][j] == 'F') {
        nofacewritten = 1;
      } else if (argv[i][j] == 'I') {
        noiterationnum = 1;
      } else if (argv[i][j] == 'S') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          steinerleft = (int) strtol(workstring, (char **) NULL, 0);
        }
      } else if (argv[i][j] == 'o') {
        if (argv[i][j + 1] == '2') {
          order = 2;
          j++;
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            optmaxdihedral = (REAL) strtod(workstring, (char **) NULL);
          }
        }
      } else if (argv[i][j] == 'O') {
        if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
          optlevel = (argv[i][j + 1] - '0');
          j++;
        }
        if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
          j++;
          if ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '7')) {
            optscheme = (argv[i][j + 1] - '0');
            j++;
          }
        }
      } else if (argv[i][j] == 'T') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          epsilon = (REAL) strtod(workstring, (char **) NULL);
        }
      } else if (argv[i][j] == 'R') {
        reversetetori = 1;
      } else if (argv[i][j] == 'C') {
        docheck++;
      } else if (argv[i][j] == 'Q') {
        quiet = 1;
      } else if (argv[i][j] == 'V') {
        verbose++;
      } else if (argv[i][j] == 'x') {
        if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
            (argv[i][j + 1] == '.')) {
          k = 0;
          while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                 (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                 (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
            j++;
            workstring[k] = argv[i][j];
            k++;
          }
          workstring[k] = '\0';
          tetrahedraperblock = (int) strtol(workstring, (char **) NULL, 0);
          if (tetrahedraperblock > 8188) {
            vertexperblock = tetrahedraperblock / 2;
            shellfaceperblock = vertexperblock / 2;
          } else {
            tetrahedraperblock = 8188;
          }
        }
      } else if ((argv[i][j] == 'h') || (argv[i][j] == 'H') ||
                 (argv[i][j] == '?')) {
        usage();
      } else {
        printf("Warning:  Unknown switch -%c.\n", argv[i][j]);
      }
    }
  }

  if (startindex == 0) {
    // Set a temporary filename for debugging output.
    strcpy(infilename, "tetgen-tmpfile");
  } else {
    if (infilename[0] == '\0') {
      // No input file name. Print the syntax and exit.
      syntax();
      terminatetetgen(NULL, 0);
    }
    // Recognize the object from file extension if it is available.
    if (!strcmp(&infilename[strlen(infilename) - 5], ".node")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = NODES;
    } else if (!strcmp(&infilename[strlen(infilename) - 5], ".poly")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = POLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 6], ".smesh")) {
      infilename[strlen(infilename) - 6] = '\0';
      object = POLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".off")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = OFF;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".ply")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = PLY;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".stl")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = STL;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 5], ".mesh")) {
      infilename[strlen(infilename) - 5] = '\0';
      object = MEDIT;
      if (!refine) plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".vtk")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = VTK;
      plc = 1;
    } else if (!strcmp(&infilename[strlen(infilename) - 4], ".ele")) {
      infilename[strlen(infilename) - 4] = '\0';
      object = MESH;
      refine = 1;
    }
  }

  if (nobisect && (!plc && !refine)) { // -Y
    plc = 1; // Default -p option.
  }
  if (quality && (!plc && !refine)) { // -q
    plc = 1; // Default -p option.
  }
  if (diagnose && !plc) { // -d
    plc = 1;
  }
  if (refine && !quality) { // -r only
    // Reconstruct a mesh, no mesh optimization.
    optlevel = 0;
  }
  if (insertaddpoints && (optlevel == 0)) { // with -i option
    optlevel = 2;
  }
  if (coarsen && (optlevel == 0)) { // with -R option
    optlevel = 2;
  }

  // Detect improper combinations of switches.
  if ((refine || plc) && weighted) {
    printf("Error:  Switches -w cannot use together with -p or -r.\n");
    return false;
  }

  if (convex) { // -c
    if (plc && !regionattrib) {
      // -A (region attribute) is needed for marking exterior tets (-1).
      regionattrib = 1; 
    }
  }

  // Note: -A must not used together with -r option. 
  // Be careful not to add an extra attribute to each element unless the
  //   input supports it (PLC in, but not refining a preexisting mesh).
  if (refine || !plc) {
    regionattrib = 0;
  }
  // Be careful not to allocate space for element area constraints that 
  //   will never be assigned any value (other than the default -1.0).
  if (!refine && !plc) {
    varvolume = 0;
  }
  // If '-a' or '-aa' is in use, enable '-q' option too.
  if (fixedvolume || varvolume) {
    if (quality == 0) {
      quality = 1;
      if (!plc && !refine) {
        plc = 1; // enable -p.
      }
    }
  }
  // No user-specified dihedral angle bound. Use default ones.
  if (!quality) {
    if (optmaxdihedral < 179.0) {
      if (nobisect) {  // with -Y option
        optmaxdihedral = 179.0;
      } else { // -p only
        optmaxdihedral = 179.999;
      }
    }
    if (optminsmtdihed < 179.999) {
      optminsmtdihed = 179.999;
    }
    if (optminslidihed < 179.999) {
      optminslidihed = 179.999;
    }
  }

  increment = 0;
  strcpy(workstring, infilename);
  j = 1;
  while (workstring[j] != '\0') {
    if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (workstring[j] != '\0');
  }
  if (noiterationnum) {
    strcpy(outfilename, infilename);
  } else if (increment == 0) {
    strcpy(outfilename, infilename);
    strcat(outfilename, ".1");
  } else {
    workstring[increment] = '%';
    workstring[increment + 1] = 'd';
    workstring[increment + 2] = '\0';
    sprintf(outfilename, workstring, meshnumber + 1);
  }
  // Additional input file name has the end ".a".
  strcpy(addinfilename, infilename);
  strcat(addinfilename, ".a");
  // Background filename has the form "*.b.ele", "*.b.node", ...
  strcpy(bgmeshfilename, infilename);
  strcat(bgmeshfilename, ".b");

  return true;
}

////                                                                       ////
////                                                                       ////
//// behavior_cxx /////////////////////////////////////////////////////////////

