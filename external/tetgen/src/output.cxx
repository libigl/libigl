#include "../tetgen.h"
//// output_cxx ///////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// jettisonnodes()    Jettison unused or duplicated vertices.                //
//                                                                           //
// Unused points are those input points which are outside the mesh domain or //
// have no connection (isolated) to the mesh.  Duplicated points exist for   //
// example if the input PLC is read from a .stl mesh file (marked during the //
// Delaunay tetrahedralization step. This routine remove these points from   //
// points list. All existing points are reindexed.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::jettisonnodes()
{
  point pointloop;
  bool jetflag;
  int oldidx, newidx;
  int remcount;

  if (!b->quiet) {
    printf("Jettisoning redundant points.\n");
  }

  points->traversalinit();
  pointloop = pointtraverse();
  oldidx = newidx = 0; // in->firstnumber;
  remcount = 0;
  while (pointloop != (point) NULL) {
    jetflag = (pointtype(pointloop) == DUPLICATEDVERTEX) || 
      (pointtype(pointloop) == UNUSEDVERTEX);
    if (jetflag) {
      // It is a duplicated or unused point, delete it.
      pointdealloc(pointloop);
      remcount++;
    } else {
      // Re-index it.
      setpointmark(pointloop, newidx + in->firstnumber);
      if (in->pointmarkerlist != (int *) NULL) {
        if (oldidx < in->numberofpoints) {
          // Re-index the point marker as well.
          in->pointmarkerlist[newidx] = in->pointmarkerlist[oldidx];
        }
      }
      newidx++;
    }
    oldidx++;
    pointloop = pointtraverse();
  }
  if (b->verbose) {
    printf("  %ld duplicated vertices are removed.\n", dupverts);
    printf("  %ld unused vertices are removed.\n", unuverts);
  }
  dupverts = 0l;
  unuverts = 0l;

  // The following line ensures that dead items in the pool of nodes cannot
  //   be allocated for the new created nodes. This ensures that the input
  //   nodes will occur earlier in the output files, and have lower indices.
  points->deaditemstack = (void *) NULL;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// highorder()   Create extra nodes for quadratic subparametric elements.    //
//                                                                           //
// 'highordertable' is an array (size = numberoftetrahedra * 6) for storing  //
// high-order nodes of each tetrahedron.  This routine is used only when -o2 //
// switch is used.                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::highorder()
{
  triface tetloop, worktet, spintet;
  point *extralist, *adjextralist;
  point torg, tdest, newpoint;
  int highorderindex;
  int t1ver;
  int i, j;

  if (!b->quiet) {
    printf("Adding vertices for second-order tetrahedra.\n");
  }

  // Initialize the 'highordertable'.
  highordertable = new point[tetrahedrons->items * 6];
  if (highordertable == (point *) NULL) {
    terminatetetgen(this, 1);
  }

  // This will overwrite the slot for element markers.
  highorderindex = 11;

  // The following line ensures that dead items in the pool of nodes cannot
  //   be allocated for the extra nodes associated with high order elements.
  //   This ensures that the primary nodes (at the corners of elements) will
  //   occur earlier in the output files, and have lower indices, than the
  //   extra nodes.
  points->deaditemstack = (void *) NULL;

  // Assign an entry for each tetrahedron to find its extra nodes. At the
  //   mean while, initialize all extra nodes be NULL.
  i = 0;
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    tetloop.tet[highorderindex] = (tetrahedron) &highordertable[i];
    for (j = 0; j < 6; j++) {
      highordertable[i + j] = (point) NULL;
    }
    i += 6;
    tetloop.tet = tetrahedrontraverse();
  }

  // To create a unique node on each edge. Loop over all tetrahedra, and
  //   look at the six edges of each tetrahedron.  If the extra node in
  //   the tetrahedron corresponding to this edge is NULL, create a node
  //   for this edge, at the same time, set the new node into the extra
  //   node lists of all other tetrahedra sharing this edge.  
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Get the list of extra nodes.
    extralist = (point *) tetloop.tet[highorderindex];
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      if (extralist[i] == (point) NULL) {
        // Go to the ith-edge.
        worktet.ver = edge2ver[i];
        // Create a new point in the middle of this edge.
        torg = org(worktet);
        tdest = dest(worktet);
        makepoint(&newpoint, FREEVOLVERTEX);
        for (j = 0; j < 3 + numpointattrib; j++) {
          newpoint[j] = 0.5 * (torg[j] + tdest[j]);
        }
        // Interpolate its metrics.
        for (j = 0; j < in->numberofpointmtrs; j++) {
          newpoint[pointmtrindex + j] = 
            0.5 * (torg[pointmtrindex + j] + tdest[pointmtrindex + j]);
        }
        // Set this point into all extra node lists at this edge.
        spintet = worktet;
        while (1) {
          if (!ishulltet(spintet)) {
            adjextralist = (point *) spintet.tet[highorderindex];
            adjextralist[ver2edge[spintet.ver]] = newpoint;
          }
          fnextself(spintet);
          if (spintet.tet == worktet.tet) break;
        }
      } // if (!extralist[i])
    } // i
    tetloop.tet = tetrahedrontraverse();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// numberedges()    Count the number of edges, save in "meshedges".          //
//                                                                           //
// This routine is called when '-p' or '-r', and '-E' options are used.  The //
// total number of edges depends on the genus of the input surface mesh.     //
//                                                                           //
// NOTE:  This routine must be called after outelements().  So all elements  //
// have been indexed.                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::numberedges()
{
  triface worktet, spintet;
  int ishulledge;
  int t1ver;
  int i;

  meshedges = meshhulledges = 0l;

  tetrahedrons->traversalinit();
  worktet.tet = tetrahedrontraverse();
  while (worktet.tet != NULL) {
    // Count the number of Voronoi faces. Look at the six edges of this
    //   tet. Count an edge only if this tet's index is smaller than
    //   those of other non-hull tets which share this edge.
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      ishulledge = 0;
      fnext(worktet, spintet);
      do {
        if (!ishulltet(spintet)) {          
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          ishulledge = 1;
        }
        fnextself(spintet);
      } while (spintet.tet != worktet.tet);
      // Count this edge if no adjacent tets are smaller than this tet.
      if (spintet.tet == worktet.tet) {
        meshedges++;
        if (ishulledge) meshhulledges++;
      }
    }
    worktet.tet = tetrahedrontraverse();
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outnodes()    Output the points to a .node file or a tetgenio structure.  //
//                                                                           //
// Note: each point has already been numbered on input (the first index is   //
// 'in->firstnumber').                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outnodes(tetgenio* out)
{
  FILE *outfile = NULL;
  char outnodefilename[FILENAMESIZE];
  face parentsh;
  point pointloop;
  int nextras, bmark, marker = 0, weightDT = 0; 
  int coordindex, attribindex;
  int pointnumber, firstindex;
  int index, i;

  if (out == (tetgenio *) NULL) {
    strcpy(outnodefilename, b->outfilename);
    strcat(outnodefilename, ".node");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outnodefilename);
    } else {
      printf("Writing nodes.\n");
    }
  }

  nextras = numpointattrib;
  if (b->weighted) { // -w
    if (b->weighted_param == 0) weightDT = 1; // Weighted DT.
  }

  bmark = !b->nobound && in->pointmarkerlist;

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outnodefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outnodefilename);
      terminatetetgen(this, 1);
    }
    // Number of points, number of dimensions, number of point attributes,
    //   and number of boundary markers (zero or one).
    fprintf(outfile, "%ld  %d  %d  %d\n", points->items, 3, nextras, bmark);
  } else {
    // Allocate space for 'pointlist';
    out->pointlist = new REAL[points->items * 3];
    if (out->pointlist == (REAL *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    // Allocate space for 'pointattributelist' if necessary;
    if (nextras > 0) {
      out->pointattributelist = new REAL[points->items * nextras];
      if (out->pointattributelist == (REAL *) NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    // Allocate space for 'pointmarkerlist' if necessary;
    if (bmark) {
      out->pointmarkerlist = new int[points->items];
      if (out->pointmarkerlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    if (b->psc) {
      out->pointparamlist = new tetgenio::pointparam[points->items];
      if (out->pointparamlist == NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    out->numberofpoints = points->items;
    out->numberofpointattributes = nextras;
    coordindex = 0;
    attribindex = 0;
  }
  
  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;

  points->traversalinit();
  pointloop = pointtraverse();
  pointnumber = firstindex; // in->firstnumber;
  index = 0;
  while (pointloop != (point) NULL) {
    if (bmark) {
      // Default the vertex has a zero marker.
      marker = 0;
      // Is it an input vertex?
      if (index < in->numberofpoints) {
        // Input point's marker is directly copied to output.
        marker = in->pointmarkerlist[index];       
      } else {
        if ((pointtype(pointloop) == FREESEGVERTEX) ||
            (pointtype(pointloop) == FREEFACETVERTEX)) {
          sdecode(point2sh(pointloop), parentsh);
          if (parentsh.sh != NULL) {
            marker = shellmark(parentsh);
            if (pointtype(pointloop) == FREEFACETVERTEX) {
              if (in->facetmarkerlist != NULL) {
                marker = in->facetmarkerlist[marker - 1];
              }
            }
          }
        } // if (pointtype(...))
      }
    }
    if (out == (tetgenio *) NULL) {
      // Point number, x, y and z coordinates.
      fprintf(outfile, "%4d    %.17g  %.17g  %.17g", pointnumber,
              pointloop[0], pointloop[1], pointloop[2]);
      for (i = 0; i < nextras; i++) {
        // Write an attribute.
        if ((i == 0) && weightDT) {          
          fprintf(outfile, "  %.17g", pointloop[0] * pointloop[0] +
             pointloop[1] * pointloop[1] + pointloop[2] * pointloop[2] 
             - pointloop[3 + i]);
        } else { 
          fprintf(outfile, "  %.17g", pointloop[3 + i]);
        }
      }
      if (bmark) {
        // Write the boundary marker.
        fprintf(outfile, "    %d", marker);
      }
      if (b->psc) {
        fprintf(outfile, "  %.8g  %.8g  %d", pointgeomuv(pointloop, 0),
                pointgeomuv(pointloop, 1), pointgeomtag(pointloop));
        if (pointtype(pointloop) == RIDGEVERTEX) {
          fprintf(outfile, "  0");
        } else if (pointtype(pointloop) == ACUTEVERTEX) {
          fprintf(outfile, "  0");
        } else if (pointtype(pointloop) == FREESEGVERTEX) {
          fprintf(outfile, "  1");
        } else if (pointtype(pointloop) == FREEFACETVERTEX) {
          fprintf(outfile, "  2");
        } else if (pointtype(pointloop) == FREEVOLVERTEX) {
          fprintf(outfile, "  3");
        } else {
          fprintf(outfile, "  -1"); // Unknown type.
        }
      }
      fprintf(outfile, "\n");
    } else {
      // X, y, and z coordinates.
      out->pointlist[coordindex++] = pointloop[0];
      out->pointlist[coordindex++] = pointloop[1];
      out->pointlist[coordindex++] = pointloop[2];
      // Point attributes.
      for (i = 0; i < nextras; i++) {
        // Output an attribute.
        if ((i == 0) && weightDT) {
          out->pointattributelist[attribindex++] = 
            pointloop[0] * pointloop[0] + pointloop[1] * pointloop[1] + 
            pointloop[2] * pointloop[2] - pointloop[3 + i];
        } else {
          out->pointattributelist[attribindex++] = pointloop[3 + i];
        }
      }
      if (bmark) {
        // Output the boundary marker.  
        out->pointmarkerlist[index] = marker;
      }
      if (b->psc) {
        out->pointparamlist[index].uv[0] = pointgeomuv(pointloop, 0);
        out->pointparamlist[index].uv[1] = pointgeomuv(pointloop, 1);
        out->pointparamlist[index].tag = pointgeomtag(pointloop);
        if (pointtype(pointloop) == RIDGEVERTEX) {
          out->pointparamlist[index].type = 0;
        } else if (pointtype(pointloop) == ACUTEVERTEX) {
          out->pointparamlist[index].type = 0;
        } else if (pointtype(pointloop) == FREESEGVERTEX) {
          out->pointparamlist[index].type = 1;
        } else if (pointtype(pointloop) == FREEFACETVERTEX) {
          out->pointparamlist[index].type = 2;
        } else if (pointtype(pointloop) == FREEVOLVERTEX) {
          out->pointparamlist[index].type = 3;
        } else {
          out->pointparamlist[index].type = -1; // Unknown type.
        }
      }
    }
    pointloop = pointtraverse();
    pointnumber++; 
    index++;
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outmetrics()    Output the metric to a file (*.mtr) or a tetgenio obj.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outmetrics(tetgenio* out)
{
  FILE *outfile = NULL;
  char outmtrfilename[FILENAMESIZE];
  point ptloop;
  int mtrindex;

  if (out == (tetgenio *) NULL) {
    strcpy(outmtrfilename, b->outfilename);
    strcat(outmtrfilename, ".mtr");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outmtrfilename);
    } else {
      printf("Writing metrics.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outmtrfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outmtrfilename);
      terminatetetgen(this, 3);
    }
    // Number of points, number of point metrices,
    // fprintf(outfile, "%ld  %d\n", points->items, sizeoftensor + 3);
    fprintf(outfile, "%ld  %d\n", points->items, 1);
  } else {
    // Allocate space for 'pointmtrlist' if necessary;
    // out->pointmtrlist = new REAL[points->items * (sizeoftensor + 3)];
    out->pointmtrlist = new REAL[points->items];
    if (out->pointmtrlist == (REAL *) NULL) {
      terminatetetgen(this, 1);
    }
    out->numberofpointmtrs = 1; // (sizeoftensor + 3);
    mtrindex = 0;
  }

  points->traversalinit();
  ptloop = pointtraverse();
  while (ptloop != (point) NULL) {
    if (out == (tetgenio *) NULL) {
      // for (i = 0; i < sizeoftensor; i++) {
      //   fprintf(outfile, "%-16.8e ", ptloop[pointmtrindex + i]);
      // }
      fprintf(outfile, "%-16.8e\n", ptloop[pointmtrindex]);
    } else {
      // for (i = 0; i < sizeoftensor; i++) {
      //   out->pointmtrlist[mtrindex++] = ptloop[pointmtrindex + i];
      // }
      out->pointmtrlist[mtrindex++] = ptloop[pointmtrindex];
    }
    ptloop = pointtraverse();
  }
  
  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outelements()    Output the tetrahedra to an .ele file or a tetgenio      //
//                  structure.                                               //
//                                                                           //
// This routine also indexes all tetrahedra (exclusing hull tets) (from in-> //
// firstnumber). The total number of mesh edges is counted in 'meshedges'.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outelements(tetgenio* out)
{
  FILE *outfile = NULL;
  char outelefilename[FILENAMESIZE];
  tetrahedron* tptr;
  point p1, p2, p3, p4;
  point *extralist;
  REAL *talist = NULL;
  int *tlist = NULL;
  long ntets;
  int firstindex, shift;
  int pointindex, attribindex;
  int highorderindex = 11; 
  int elementnumber;
  int eextras;
  int i;

  if (out == (tetgenio *) NULL) {
    strcpy(outelefilename, b->outfilename);
    strcat(outelefilename, ".ele");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outelefilename);
    } else {
      printf("Writing elements.\n");
    }
  }

  // The number of tets excluding hull tets.
  ntets = tetrahedrons->items - hullsize;

  eextras = numelemattrib;
  if (out == (tetgenio *) NULL) {
    outfile = fopen(outelefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outelefilename);
      terminatetetgen(this, 1);
    }
    // Number of tetras, points per tetra, attributes per tetra.
    fprintf(outfile, "%ld  %d  %d\n", ntets, b->order == 1 ? 4 : 10, eextras);
  } else {
    // Allocate memory for output tetrahedra.
    out->tetrahedronlist = new int[ntets * (b->order == 1 ? 4 : 10)];
    if (out->tetrahedronlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    // Allocate memory for output tetrahedron attributes if necessary.
    if (eextras > 0) {
      out->tetrahedronattributelist = new REAL[ntets * eextras];
      if (out->tetrahedronattributelist == (REAL *) NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    out->numberoftetrahedra = ntets;
    out->numberofcorners = b->order == 1 ? 4 : 10;
    out->numberoftetrahedronattributes = eextras;
    tlist = out->tetrahedronlist;
    talist = out->tetrahedronattributelist;
    pointindex = 0;
    attribindex = 0;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shift.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tptr = tetrahedrontraverse();
  elementnumber = firstindex; // in->firstnumber;
  while (tptr != (tetrahedron *) NULL) {
    if (!b->reversetetori) {
      p1 = (point) tptr[4];
      p2 = (point) tptr[5];
    } else {
      p1 = (point) tptr[5];
      p2 = (point) tptr[4];
    }
    p3 = (point) tptr[6];
    p4 = (point) tptr[7];
    if (out == (tetgenio *) NULL) {
      // Tetrahedron number, indices for four points.
      fprintf(outfile, "%5d   %5d %5d %5d %5d", elementnumber,
              pointmark(p1) - shift, pointmark(p2) - shift,
              pointmark(p3) - shift, pointmark(p4) - shift);
      if (b->order == 2) {
        extralist = (point *) tptr[highorderindex];
        // indices for six extra points.
        fprintf(outfile, "  %5d %5d %5d %5d %5d %5d",
          pointmark(extralist[0]) - shift, pointmark(extralist[1]) - shift,
          pointmark(extralist[2]) - shift, pointmark(extralist[3]) - shift,
          pointmark(extralist[4]) - shift, pointmark(extralist[5]) - shift);
      }
      for (i = 0; i < eextras; i++) {
        fprintf(outfile, "    %.17g", elemattribute(tptr, i));
      }
      fprintf(outfile, "\n");
    } else {
      tlist[pointindex++] = pointmark(p1) - shift;
      tlist[pointindex++] = pointmark(p2) - shift;
      tlist[pointindex++] = pointmark(p3) - shift;
      tlist[pointindex++] = pointmark(p4) - shift;
      if (b->order == 2) {
        extralist = (point *) tptr[highorderindex];
        tlist[pointindex++] = pointmark(extralist[0]) - shift;
        tlist[pointindex++] = pointmark(extralist[1]) - shift;
        tlist[pointindex++] = pointmark(extralist[2]) - shift;
        tlist[pointindex++] = pointmark(extralist[3]) - shift;
        tlist[pointindex++] = pointmark(extralist[4]) - shift;
        tlist[pointindex++] = pointmark(extralist[5]) - shift;
      }
      for (i = 0; i < eextras; i++) {
        talist[attribindex++] = elemattribute(tptr, i);
      }
    }
    // Remember the index of this element (for counting edges).
    setelemindex(tptr, elementnumber);
    tptr = tetrahedrontraverse();
    elementnumber++;
  }


  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outfaces()    Output all faces to a .face file or a tetgenio object.      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outfaces(tetgenio* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  triface tface, tsymface;
  face checkmark;
  point torg, tdest, tapex;
  long ntets, faces;
  int *elist = NULL, *emlist = NULL;
  int neigh1 = 0, neigh2 = 0;
  int faceid, marker = 0;
  int firstindex, shift;
  int facenumber;
  int index = 0;

  // For -o2 option.
  triface workface;
  point *extralist, pp[3] = {0,0,0}; 
  int highorderindex = 11; 
  int o2index = 0, i;

  if (out == (tetgenio *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  ntets = tetrahedrons->items - hullsize;
  faces = (ntets * 4l + hullsize) / 2l;

  if (out == (tetgenio *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminatetetgen(this, 1);
    }
    fprintf(outfile, "%ld  %d\n", faces, !b->nobound);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[faces * 3];
    if (out->trifacelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    if (b->order == 2) {
      out->o2facelist = new int[faces * 3];
    }
    // Allocate memory for 'trifacemarkerlist' if necessary.
    if (!b->nobound) {
      out->trifacemarkerlist = new int[faces];
      if (out->trifacemarkerlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    if (b->neighout > 1) {
      // '-nn' switch.
      out->adjtetlist = new int[faces * 2];
      if (out->adjtetlist == (int *) NULL) {
        printf("Error:  Out of memory.\n");
        terminatetetgen(this, 1);
      }
    }
    out->numberoftrifaces = faces;
    elist = out->trifacelist;
    emlist = out->trifacemarkerlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tface.tet = tetrahedrontraverse();
  facenumber = firstindex; // in->firstnumber;
  // To loop over the set of faces, loop over all tetrahedra, and look at
  //   the four faces of each one. If its adjacent tet is a hull tet,
  //   operate on the face, otherwise, operate on the face only if the
  //   current tet has a smaller index than its neighbor.
  while (tface.tet != (tetrahedron *) NULL) {
    for (tface.ver = 0; tface.ver < 4; tface.ver ++) {
      fsym(tface, tsymface);
      if (ishulltet(tsymface) || 
          (elemindex(tface.tet) < elemindex(tsymface.tet))) {
        torg = org(tface);
        tdest = dest(tface);
        tapex = apex(tface);
        if (b->order == 2) { // -o2
          // Get the three extra vertices on edges.
          extralist = (point *) (tface.tet[highorderindex]);
          // The extra vertices are on edges opposite the corners.
          enext(tface, workface);
          for (i = 0; i < 3; i++) {
            pp[i] = extralist[ver2edge[workface.ver]];
            enextself(workface);
          }
        }
        if (!b->nobound) {
          // Get the boundary marker of this face.
          if (b->plc || b->refine) { 
            // Shell face is used.
            tspivot(tface, checkmark);
            if (checkmark.sh == NULL) {
              marker = 0;  // It is an inner face. It's marker is 0.
            } else {
              if (in->facetmarkerlist) {
                // The facet marker is given, get it.
                faceid = shellmark(checkmark) - 1;
                marker = in->facetmarkerlist[faceid];
              } else {
                marker = 1; // The default marker for subface is 1.
              }
            }
          } else {
            // Shell face is not used, only distinguish outer and inner face.
            marker = (int) ishulltet(tsymface);
          }
        }
        if (b->neighout > 1) {
          // '-nn' switch. Output adjacent tets indices.
          neigh1 = elemindex(tface.tet);
          if (!ishulltet(tsymface)) {
            neigh2 = elemindex(tsymface.tet);
          } else {
            neigh2 = -1;  
          }
        }
        if (out == (tetgenio *) NULL) {
          // Face number, indices of three vertices.
          fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
                  pointmark(torg) - shift, pointmark(tdest) - shift,
                  pointmark(tapex) - shift);
          if (b->order == 2) { // -o2
            fprintf(outfile, "  %4d  %4d  %4d", pointmark(pp[0]) - shift, 
                    pointmark(pp[1]) - shift, pointmark(pp[2]) - shift);
          }
          if (!b->nobound) {
            // Output a boundary marker.
            fprintf(outfile, "  %d", marker);
          }
          if (b->neighout > 1) {
            fprintf(outfile, "    %5d  %5d", neigh1, neigh2);
          }
          fprintf(outfile, "\n");
        } else {
          // Output indices of three vertices.
          elist[index++] = pointmark(torg) - shift;
          elist[index++] = pointmark(tdest) - shift;
          elist[index++] = pointmark(tapex) - shift;
          if (b->order == 2) { // -o2
            out->o2facelist[o2index++] = pointmark(pp[0]) - shift;
            out->o2facelist[o2index++] = pointmark(pp[1]) - shift;
            out->o2facelist[o2index++] = pointmark(pp[2]) - shift;
          }
          if (!b->nobound) {
            emlist[facenumber - in->firstnumber] = marker;
          }
          if (b->neighout > 1) {
            out->adjtetlist[(facenumber - in->firstnumber) * 2]     = neigh1;
            out->adjtetlist[(facenumber - in->firstnumber) * 2 + 1] = neigh2;
          }
        }
        facenumber++;
      }
    }
    tface.tet = tetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outhullfaces()    Output hull faces to a .face file or a tetgenio object. //
//                                                                           //
// The normal of each face is pointing to the outside of the domain.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outhullfaces(tetgenio* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  triface hulltet;
  point torg, tdest, tapex;
  int *elist = NULL;
  int firstindex, shift;
  int facenumber;
  int index;

  if (out == (tetgenio *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminatetetgen(this, 1);
    }
    fprintf(outfile, "%ld  0\n", hullsize);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[hullsize * 3];
    if (out->trifacelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    out->numberoftrifaces = hullsize;
    elist = out->trifacelist;
    index = 0;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  tetrahedrons->traversalinit();
  hulltet.tet = alltetrahedrontraverse();
  facenumber = firstindex;
  while (hulltet.tet != (tetrahedron *) NULL) {
    if (ishulltet(hulltet)) {
      torg = (point) hulltet.tet[4];
      tdest = (point) hulltet.tet[5];
      tapex = (point) hulltet.tet[6];
      if (out == (tetgenio *) NULL) {
        // Face number, indices of three vertices.
        fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
                pointmark(torg) - shift, pointmark(tdest) - shift,
                pointmark(tapex) - shift);
        fprintf(outfile, "\n");
      } else {
        // Output indices of three vertices.
        elist[index++] = pointmark(torg) - shift;
        elist[index++] = pointmark(tdest) - shift;
        elist[index++] = pointmark(tapex) - shift;
      }
      facenumber++;
    }
    hulltet.tet = alltetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outsubfaces()    Output subfaces (i.e. boundary faces) to a .face file or //
//                  a tetgenio structure.                                    //
//                                                                           //
// The boundary faces are found in 'subfaces'. For listing triangle vertices //
// in the same sense for all triangles in the mesh, the direction determined //
// by right-hand rule is pointer to the inside of the volume.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outsubfaces(tetgenio* out)
{
  FILE *outfile = NULL;
  char facefilename[FILENAMESIZE];
  int *elist = NULL;
  int *emlist = NULL;
  int index = 0, index1 = 0, index2 = 0;
  triface abuttingtet;
  face faceloop;
  point torg, tdest, tapex;
  int faceid = 0, marker = 0;
  int firstindex, shift;
  int neigh1 = 0, neigh2 = 0;
  int facenumber;

  // For -o2 option.
  triface workface;
  point *extralist, pp[3] = {0,0,0}; 
  int highorderindex = 11;
  int o2index = 0, i;

  int t1ver; // used by fsymself()

  if (out == (tetgenio *) NULL) {
    strcpy(facefilename, b->outfilename);
    strcat(facefilename, ".face");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", facefilename);
    } else {
      printf("Writing faces.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(facefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", facefilename);
      terminatetetgen(this, 3);
    }
    // Number of subfaces.
    fprintf(outfile, "%ld  %d\n", subfaces->items, !b->nobound);
  } else {
    // Allocate memory for 'trifacelist'.
    out->trifacelist = new int[subfaces->items * 3];
    if (out->trifacelist == (int *) NULL) {
      terminatetetgen(this, 1);
    }
    if (b->order == 2) {
      out->o2facelist = new int[subfaces->items * 3];
    }
    if (!b->nobound) {
      // Allocate memory for 'trifacemarkerlist'.
      out->trifacemarkerlist = new int[subfaces->items];
      if (out->trifacemarkerlist == (int *) NULL) {
        terminatetetgen(this, 1);
      }
    }
    if (b->neighout > 1) {
      // '-nn' switch.
      out->adjtetlist = new int[subfaces->items * 2];
      if (out->adjtetlist == (int *) NULL) {
        terminatetetgen(this, 1);
      }
    }
    out->numberoftrifaces = subfaces->items;
    elist = out->trifacelist;
    emlist = out->trifacemarkerlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  subfaces->traversalinit();
  faceloop.sh = shellfacetraverse(subfaces);
  facenumber = firstindex; // in->firstnumber;
  while (faceloop.sh != (shellface *) NULL) {
    stpivot(faceloop, abuttingtet);
    // If there is a tetrahedron containing this subface, orient it so
    //   that the normal of this face points to inside of the volume by
    //   right-hand rule.
    if (abuttingtet.tet != NULL) {
      if (ishulltet(abuttingtet)) {
        fsymself(abuttingtet);
        assert(!ishulltet(abuttingtet));
      }
    }
    if (abuttingtet.tet != NULL) {
      torg = org(abuttingtet);
      tdest = dest(abuttingtet);
      tapex = apex(abuttingtet);
      if (b->order == 2) { // -o2
        // Get the three extra vertices on edges.
        extralist = (point *) (abuttingtet.tet[highorderindex]);
        workface = abuttingtet;
        for (i = 0; i < 3; i++) {
          pp[i] = extralist[ver2edge[workface.ver]];
          enextself(workface);
        }
      }
    } else {
      // This may happen when only a surface mesh be generated.
      torg = sorg(faceloop);
      tdest = sdest(faceloop);
      tapex = sapex(faceloop);
      if (b->order == 2) { // -o2
        // There is no extra node list available.
        pp[0] = torg;
        pp[1] = tdest;
        pp[2] = tapex;
      }
    }
    if (!b->nobound) {
      if (b->refine) { // -r option.
        if (in->trifacemarkerlist) {
          marker = shellmark(faceloop);
        } else {
          marker = 1; // Default marker for a subface is 1.
        }
      } else {
        if (in->facetmarkerlist) {
          faceid = shellmark(faceloop) - 1;
          marker = in->facetmarkerlist[faceid];
        } else {
          marker = 1; // Default marker for a subface is 1.
        }
      }
    }
    if (b->neighout > 1) {
      // '-nn' switch. Output adjacent tets indices.
      neigh1 = -1;
      neigh2 = -1;
      stpivot(faceloop, abuttingtet);
      if (abuttingtet.tet != NULL) {
        neigh1 = elemindex(abuttingtet.tet);
        fsymself(abuttingtet);
        if (!ishulltet(abuttingtet)) {
          neigh2 = elemindex(abuttingtet.tet);
        }
      }
    }
    if (out == (tetgenio *) NULL) {
      fprintf(outfile, "%5d   %4d  %4d  %4d", facenumber,
              pointmark(torg) - shift, pointmark(tdest) - shift,
              pointmark(tapex) - shift);
      if (b->order == 2) { // -o2
        fprintf(outfile, "  %4d  %4d  %4d", pointmark(pp[0]) - shift, 
                pointmark(pp[1]) - shift, pointmark(pp[2]) - shift);
      }
      if (!b->nobound) {
        fprintf(outfile, "    %d", marker);
      }
      if (b->neighout > 1) {
        fprintf(outfile, "    %5d  %5d", neigh1, neigh2);
      }
      fprintf(outfile, "\n");
    } else {
      // Output three vertices of this face;
      elist[index++] = pointmark(torg) - shift;
      elist[index++] = pointmark(tdest) - shift;
      elist[index++] = pointmark(tapex) - shift;
      if (b->order == 2) { // -o2
        out->o2facelist[o2index++] = pointmark(pp[0]) - shift;
        out->o2facelist[o2index++] = pointmark(pp[1]) - shift;
        out->o2facelist[o2index++] = pointmark(pp[2]) - shift;
      }
      if (!b->nobound) {
        emlist[index1++] = marker;
      }
      if (b->neighout > 1) {
        out->adjtetlist[index2++] = neigh1;
        out->adjtetlist[index2++] = neigh2;
      }
    }
    facenumber++;
    faceloop.sh = shellfacetraverse(subfaces);
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outedges()    Output all edges to a .edge file or a tetgenio object.      //
//                                                                           //
// Note: This routine must be called after outelements(),  so that the total //
// number of edges 'meshedges' has been counted.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outedges(tetgenio* out)
{
  FILE *outfile = NULL;
  char edgefilename[FILENAMESIZE];
  triface tetloop, worktet, spintet;
  face checkseg;
  point torg, tdest;
  int *elist = NULL, *emlist = NULL;
  int ishulledge;
  int firstindex, shift;
  int edgenumber, marker;
  int index = 0, index1 = 0, index2 = 0;
  int t1ver;
  int i;

  // For -o2 option.
  point *extralist, pp = NULL; 
  int highorderindex = 11;
  int o2index = 0;

  if (out == (tetgenio *) NULL) {
    strcpy(edgefilename, b->outfilename);
    strcat(edgefilename, ".edge");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", edgefilename);
    } else {
      printf("Writing edges.\n");
    }
  }

  if (meshedges == 0l) {
    if (nonconvex) {
      numberedges();  // Count the edges.
    } else {
      // Use Euler's characteristic to get the numbe of edges.
      // It states V - E + F - C = 1, hence E = V + F - C - 1.
      long tsize = tetrahedrons->items - hullsize;
      long fsize = (tsize * 4l + hullsize) / 2l;
      long vsize = points->items - dupverts - unuverts;
      if (b->weighted) vsize -= nonregularcount;
      meshedges = vsize + fsize - tsize - 1;
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(edgefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", edgefilename);
      terminatetetgen(this, 1);
    }
    // Write the number of edges, boundary markers (0 or 1).
    fprintf(outfile, "%ld  %d\n", meshedges, !b->nobound);
  } else {
    // Allocate memory for 'edgelist'.
    out->edgelist = new int[meshedges * 2];
    if (out->edgelist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    if (b->order == 2) { // -o2 switch
      out->o2edgelist = new int[meshedges];
    }
    if (!b->nobound) {
      out->edgemarkerlist = new int[meshedges];
    }
    if (b->neighout > 1) { // '-nn' switch.
      out->edgeadjtetlist = new int[meshedges];
    }
    out->numberofedges = meshedges;
    elist = out->edgelist;
    emlist = out->edgemarkerlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift (reduce) the output indices by 1.
  }

  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  edgenumber = firstindex; // in->firstnumber;
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi faces. 
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      ishulledge = 0;
      fnext(worktet, spintet);
      do {
        if (!ishulltet(spintet)) {
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          ishulledge = 1;
        }
        fnextself(spintet);
      } while (spintet.tet != worktet.tet);
      // Count this edge if no adjacent tets are smaller than this tet.
      if (spintet.tet == worktet.tet) {
        torg = org(worktet);
        tdest = dest(worktet);
        if (b->order == 2) { // -o2
          // Get the extra vertex on this edge.
          extralist = (point *) worktet.tet[highorderindex];
          pp = extralist[ver2edge[worktet.ver]];
        }
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, "%5d   %4d  %4d", edgenumber,
                  pointmark(torg) - shift, pointmark(tdest) - shift);
          if (b->order == 2) { // -o2
            fprintf(outfile, "  %4d", pointmark(pp) - shift);
          }
        } else {
          // Output three vertices of this face;
          elist[index++] = pointmark(torg) - shift;
          elist[index++] = pointmark(tdest) - shift;
          if (b->order == 2) { // -o2
            out->o2edgelist[o2index++] = pointmark(pp) - shift;
          }
        }
        if (!b->nobound) {
          if (b->plc || b->refine) {
            // Check if the edge is a segment.
            tsspivot1(worktet, checkseg);
            if (checkseg.sh != NULL) {
              marker = shellmark(checkseg);
              if (marker == 0) {  // Does it have no marker?
                marker = 1;  // Set the default marker for this segment.
              }
            } else {
              marker = 0;  // It's not a segment.
            }
          } else {
            // Mark it if it is a hull edge.
            marker = ishulledge ? 1 : 0;
          }
          if (out == (tetgenio *) NULL) {
            fprintf(outfile, "  %d", marker);
          } else {
            emlist[index1++] = marker;
          }
        }
        if (b->neighout > 1) { // '-nn' switch.
          if (out == (tetgenio *) NULL) {
            fprintf(outfile, "  %d", elemindex(tetloop.tet));
          } else {
            out->edgeadjtetlist[index2++] = elemindex(tetloop.tet);
          }
        }
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, "\n");
        }
        edgenumber++;
      }
    }
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outsubsegments()    Output segments to a .edge file or a structure.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outsubsegments(tetgenio* out)
{
  FILE *outfile = NULL;
  char edgefilename[FILENAMESIZE];
  int *elist = NULL;
  int index, i;
  face edgeloop;
  point torg, tdest;
  int firstindex, shift;
  int marker;
  int edgenumber;

  // For -o2 option.
  triface workface, spintet;
  point *extralist, pp = NULL; 
  int highorderindex = 11;
  int o2index = 0;

  // For -nn option.
  int neigh = -1;
  int index2 = 0;

  int t1ver; // used by fsymself()

  if (out == (tetgenio *) NULL) {
    strcpy(edgefilename, b->outfilename);
    strcat(edgefilename, ".edge");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", edgefilename);
    } else {
      printf("Writing edges.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(edgefilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", edgefilename);
      terminatetetgen(this, 3);
    }
    // Number of subsegments.
    fprintf(outfile, "%ld  1\n", subsegs->items);
  } else {
    // Allocate memory for 'edgelist'.
    out->edgelist = new int[subsegs->items * (b->order == 1 ? 2 : 3)];
    if (out->edgelist == (int *) NULL) {
      terminatetetgen(this, 1);
    }
    if (b->order == 2) {
      out->o2edgelist = new int[subsegs->items];
    }
    out->edgemarkerlist = new int[subsegs->items];
    if (out->edgemarkerlist == (int *) NULL) {
      terminatetetgen(this, 1);
    }
    if (b->neighout > 1) {
      out->edgeadjtetlist = new int[subsegs->items];
    }
    out->numberofedges = subsegs->items;
    elist = out->edgelist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }
  index = 0;
  i = 0;

  subsegs->traversalinit();
  edgeloop.sh = shellfacetraverse(subsegs);
  edgenumber = firstindex; // in->firstnumber;
  while (edgeloop.sh != (shellface *) NULL) {
    torg = sorg(edgeloop);
    tdest = sdest(edgeloop);
    if ((b->order == 2) || (b->neighout > 1)) {
      sstpivot1(edgeloop, workface);
      if (workface.tet != NULL) {
        // We must find a non-hull tet.
        if (ishulltet(workface)) {
          spintet = workface;
          while (1) {
            fnextself(spintet);
            if (!ishulltet(spintet)) break;
            if (spintet.tet == workface.tet) break;
          }
          assert(!ishulltet(spintet));
          workface = spintet;
        }
      }
    }
    if (b->order == 2) { // -o2
      // Get the extra vertex on this edge.
      if (workface.tet != NULL) {
        extralist = (point *) workface.tet[highorderindex];
        pp = extralist[ver2edge[workface.ver]];
      } else {
        pp = torg; // There is no extra node available.
      }
    }
    if (b->neighout > 1) { // -nn
      if (workface.tet != NULL) {
        neigh = elemindex(workface.tet);
      } else {
        neigh = -1;
      }
    }
    marker = shellmark(edgeloop);
    if (marker == 0) {
      marker = 1; // Default marker of a boundary edge is 1. 
    }
    if (out == (tetgenio *) NULL) {
      fprintf(outfile, "%5d   %4d  %4d", edgenumber,
              pointmark(torg) - shift, pointmark(tdest) - shift);
      if (b->order == 2) { // -o2
        fprintf(outfile, "  %4d", pointmark(pp) - shift);
      }
      fprintf(outfile, "  %d", marker);
      if (b->neighout > 1) { // -nn
        fprintf(outfile, "  %4d", neigh);
      }
      fprintf(outfile, "\n");
    } else {
      // Output three vertices of this face;
      elist[index++] = pointmark(torg) - shift;
      elist[index++] = pointmark(tdest) - shift;
      if (b->order == 2) { // -o2
        out->o2edgelist[o2index++] = pointmark(pp) - shift;
      }
      out->edgemarkerlist[i++] = marker;
      if (b->neighout > 1) { // -nn
        out->edgeadjtetlist[index2++] = neigh;
      }
    }
    edgenumber++;
    edgeloop.sh = shellfacetraverse(subsegs);
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outneighbors()    Output tet neighbors to a .neigh file or a structure.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outneighbors(tetgenio* out)
{
  FILE *outfile = NULL;
  char neighborfilename[FILENAMESIZE];
  int *nlist = NULL;
  int index = 0;
  triface tetloop, tetsym;
  int neighbori[4];
  int firstindex;
  int elementnumber;
  long ntets;

  if (out == (tetgenio *) NULL) {
    strcpy(neighborfilename, b->outfilename);
    strcat(neighborfilename, ".neigh");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", neighborfilename);
    } else {
      printf("Writing neighbors.\n");
    }
  }

  ntets = tetrahedrons->items - hullsize;

  if (out == (tetgenio *) NULL) {
    outfile = fopen(neighborfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", neighborfilename);
      terminatetetgen(this, 1);
    }
    // Number of tetrahedra, four faces per tetrahedron.
    fprintf(outfile, "%ld  %d\n", ntets, 4);
  } else {
    // Allocate memory for 'neighborlist'.
    out->neighborlist = new int[ntets * 4];
    if (out->neighborlist == (int *) NULL) {
      printf("Error:  Out of memory.\n");
      terminatetetgen(this, 1);
    }
    nlist = out->neighborlist;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;

  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  elementnumber = firstindex; // in->firstnumber;
  while (tetloop.tet != (tetrahedron *) NULL) {
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, tetsym);
      if (!ishulltet(tetsym)) {
        neighbori[tetloop.ver] = elemindex(tetsym.tet);
      } else {
        neighbori[tetloop.ver] = -1;
      }
    }
    if (out == (tetgenio *) NULL) {
      // Tetrahedra number, neighboring tetrahedron numbers.
      fprintf(outfile, "%4d    %4d  %4d  %4d  %4d\n", elementnumber,
              neighbori[0], neighbori[1], neighbori[2], neighbori[3]);
    } else {
      nlist[index++] = neighbori[0];
      nlist[index++] = neighbori[1];
      nlist[index++] = neighbori[2];
      nlist[index++] = neighbori[3];
    }
    tetloop.tet = tetrahedrontraverse();
    elementnumber++;
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outvoronoi()    Output the Voronoi diagram to .v.node, .v.edge, v.face,   //
//                 and .v.cell.                                              //
//                                                                           //
// The Voronoi diagram is the geometric dual of the Delaunay triangulation.  //
// The Voronoi vertices are the circumcenters of Delaunay tetrahedra.  Each  //
// Voronoi edge connects two Voronoi vertices at two sides of a common Dela- //
// unay face. At a face of convex hull, it becomes a ray (goto the infinity).//
// A Voronoi face is the convex hull of all Voronoi vertices around a common //
// Delaunay edge. It is a closed polygon for any internal Delaunay edge. At a//
// ridge, it is unbounded.  Each Voronoi cell is the convex hull of all Vor- //
// onoi vertices around a common Delaunay vertex. It is a polytope for any   //
// internal Delaunay vertex. It is an unbounded polyhedron for a Delaunay    //
// vertex belonging to the convex hull.                                      //
//                                                                           //
// NOTE: This routine is only used when the input is only a set of point.    //
// Comment: Special thanks to Victor Liu for finding and fixing few bugs.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outvoronoi(tetgenio* out)
{
  FILE *outfile = NULL;
  char outfilename[FILENAMESIZE];
  tetgenio::voroedge *vedge = NULL;
  tetgenio::vorofacet *vfacet = NULL;
  arraypool *tetlist, *ptlist;
  triface tetloop, worktet, spintet, firsttet;
  point pt[4], ploop, neipt;
  REAL ccent[3], infvec[3], vec1[3], vec2[3], L;
  long ntets, faces, edges;
  int *indexarray, *fidxs, *eidxs;
  int arraysize, *vertarray = NULL;
  int vpointcount, vedgecount, vfacecount, tcount;
  int ishullvert, ishullface;
  int index, shift, end1, end2;
  int i, j;

  int t1ver; // used by fsymself()

  // Output Voronoi vertices to .v.node file.
  if (out == (tetgenio *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.node");
  }

  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi vertices.\n");
    }
  }

  // Determine the first index (0 or 1).
  shift = (b->zeroindex ? 0 : in->firstnumber);

  // Each face and edge of the tetrahedral mesh will be indexed for indexing
  //   the Voronoi edges and facets. Indices of faces and edges are saved in
  //   each tetrahedron (including hull tets).

  // Allocate the total space once.
  indexarray = new int[tetrahedrons->items * 10];

  // Allocate space (10 integers) into each tetrahedron. It re-uses the slot
  //   for element markers, flags.
  i = 0;
  tetrahedrons->traversalinit();
  tetloop.tet = alltetrahedrontraverse();
  while (tetloop.tet != NULL) {
    tetloop.tet[11] = (tetrahedron) &(indexarray[i * 10]);
    i++;
    tetloop.tet = alltetrahedrontraverse();
  }

  // The number of tetrahedra (excluding hull tets) (Voronoi vertices).
  ntets = tetrahedrons->items - hullsize;
  // The number of Delaunay faces (Voronoi edges).
  faces = (4l * ntets + hullsize) / 2l;
  // The number of Delaunay edges (Voronoi faces).
  long vsize = points->items - dupverts - unuverts;
  if (b->weighted) vsize -= nonregularcount;
  edges = vsize + faces - ntets - 1;

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminatetetgen(this, 3);
    }
    // Number of voronoi points, 3 dim, no attributes, no marker.
    fprintf(outfile, "%ld  3  0  0\n", ntets);
  } else {
    // Allocate space for 'vpointlist'.
    out->numberofvpoints = (int) ntets;
    out->vpointlist = new REAL[out->numberofvpoints * 3];
    if (out->vpointlist == (REAL *) NULL) {
      terminatetetgen(this, 1);
    }
  }

  // Output Voronoi vertices (the circumcenters of tetrahedra). 
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vpointcount = 0; // The (internal) v-index always starts from 0. 
  index = 0;
  while (tetloop.tet != (tetrahedron *) NULL) {
    for (i = 0; i < 4; i++) {
      pt[i] = (point) tetloop.tet[4 + i];
      setpoint2tet(pt[i], encode(tetloop));
    }
    if (b->weighted) {
      orthosphere(pt[0], pt[1], pt[2], pt[3], pt[0][3], pt[1][3], pt[2][3], 
                  pt[3][3], ccent, NULL);
    } else {
      circumsphere(pt[0], pt[1], pt[2], pt[3], ccent, NULL);
    }
    if (out == (tetgenio *) NULL) {
      fprintf(outfile, "%4d  %16.8e %16.8e %16.8e\n", vpointcount + shift,
              ccent[0], ccent[1], ccent[2]);
    } else {
      out->vpointlist[index++] = ccent[0];
      out->vpointlist[index++] = ccent[1];
      out->vpointlist[index++] = ccent[2];
    }
    setelemindex(tetloop.tet, vpointcount);
    vpointcount++;
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi edges to .v.edge file.
  if (out == (tetgenio *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.edge");
  }
  
  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi edges.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminatetetgen(this, 3);
    }
    // Number of Voronoi edges, no marker.
    fprintf(outfile, "%ld  0\n", faces);
  } else {
    // Allocate space for 'vpointlist'.
    out->numberofvedges = (int) faces;
    out->vedgelist = new tetgenio::voroedge[out->numberofvedges];
  }

  // Output the Voronoi edges. 
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vedgecount = 0; // D-Face (V-edge) index (from zero).
  index = 0; // The Delaunay-face index.
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi edges. Look at the four faces of each
    //   tetrahedron. Count the face if the tetrahedron's index is
    //   smaller than its neighbor's or the neighbor is outside.
    end1 = elemindex(tetloop.tet);
    for (tetloop.ver = 0; tetloop.ver < 4; tetloop.ver++) {
      fsym(tetloop, worktet);
      if (ishulltet(worktet) || 
          (elemindex(tetloop.tet) < elemindex(worktet.tet))) {
        // Found a Voronoi edge. Operate on it.
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, "%4d  %4d", vedgecount + shift, end1 + shift);
        } else {
          vedge = &(out->vedgelist[index++]);
          vedge->v1 = end1 + shift;
        }
        if (!ishulltet(worktet)) {
          end2 = elemindex(worktet.tet);
        } else {
          end2 = -1;
        }
        // Note that end2 may be -1 (worktet.tet is outside).
        if (end2 == -1) {
          // Calculate the out normal of this hull face.
          pt[0] = dest(worktet);
          pt[1] = org(worktet);
          pt[2] = apex(worktet);
          for (j = 0; j < 3; j++) vec1[j] = pt[1][j] - pt[0][j];
          for (j = 0; j < 3; j++) vec2[j] = pt[2][j] - pt[0][j];
          cross(vec1, vec2, infvec);
          // Normalize it.
          L = sqrt(infvec[0] * infvec[0] + infvec[1] * infvec[1]
                   + infvec[2] * infvec[2]);
          if (L > 0) for (j = 0; j < 3; j++) infvec[j] /= L;
          if (out == (tetgenio *) NULL) {
            fprintf(outfile, " -1");
            fprintf(outfile, " %g %g %g\n", infvec[0], infvec[1], infvec[2]);
          } else {
            vedge->v2 = -1;
            vedge->vnormal[0] = infvec[0];
            vedge->vnormal[1] = infvec[1];
            vedge->vnormal[2] = infvec[2];
          }
        } else {
          if (out == (tetgenio *) NULL) {
            fprintf(outfile, " %4d\n", end2 + shift);
          } else {
            vedge->v2 = end2 + shift;
            vedge->vnormal[0] = 0.0;
            vedge->vnormal[1] = 0.0;
            vedge->vnormal[2] = 0.0;
          }
        }
        // Save the V-edge index in this tet and its neighbor.
        fidxs = (int *) (tetloop.tet[11]);
        fidxs[tetloop.ver] = vedgecount;
        fidxs = (int *) (worktet.tet[11]);
        fidxs[worktet.ver & 3] = vedgecount;
        vedgecount++;
      }
    } // tetloop.ver
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi faces to .v.face file.
  if (out == (tetgenio *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.face");
  }
  
  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi faces.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminatetetgen(this, 3);
    }
    // Number of Voronoi faces.
    fprintf(outfile, "%ld  0\n", edges);
  } else {
    out->numberofvfacets = edges;
    out->vfacetlist = new tetgenio::vorofacet[out->numberofvfacets];
    if (out->vfacetlist == (tetgenio::vorofacet *) NULL) {
      terminatetetgen(this, 1);
    }
  }

  // Output the Voronoi facets.
  tetrahedrons->traversalinit();
  tetloop.tet = tetrahedrontraverse();
  vfacecount = 0; // D-edge (V-facet) index (from zero).
  while (tetloop.tet != (tetrahedron *) NULL) {
    // Count the number of Voronoi faces. Look at the six edges of each
    //   tetrahedron. Count the edge only if the tetrahedron's index is
    //   smaller than those of all other tetrahedra that share the edge.
    worktet.tet = tetloop.tet;
    for (i = 0; i < 6; i++) {
      worktet.ver = edge2ver[i];
      // Count the number of faces at this edge. If the edge is a hull edge,
      //   the face containing dummypoint is also counted.
      //ishulledge = 0; // Is it a hull edge.
      tcount = 0;
      firsttet = worktet;
      spintet = worktet;
      while (1) {
        tcount++;
        fnextself(spintet);
        if (spintet.tet == worktet.tet) break;
        if (!ishulltet(spintet)) {
          if (elemindex(spintet.tet) < elemindex(worktet.tet)) break;
        } else {
          //ishulledge = 1;
          if (apex(spintet) == dummypoint) {
            // We make this V-edge appear in the end of the edge list.
            fnext(spintet, firsttet); 
          }
        }
      } // while (1)
      if (spintet.tet == worktet.tet) {
        // Found a Voronoi facet. Operate on it.
        pt[0] = org(worktet);
        pt[1] = dest(worktet);
        end1 = pointmark(pt[0]) - in->firstnumber; // V-cell index
        end2 = pointmark(pt[1]) - in->firstnumber;
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, "%4d  %4d %4d  %-2d ", vfacecount + shift, 
                  end1 + shift, end2 + shift, tcount);
        } else {
          vfacet = &(out->vfacetlist[vfacecount]);
          vfacet->c1 = end1 + shift;
          vfacet->c2 = end2 + shift;
          vfacet->elist = new int[tcount + 1];
          vfacet->elist[0] = tcount;
          index = 1;
        }
        // Output V-edges of this V-facet.
        spintet = firsttet; //worktet;
        while (1) {
          fidxs = (int *) (spintet.tet[11]);
          if (apex(spintet) != dummypoint) {
            vedgecount = fidxs[spintet.ver & 3];
            ishullface = 0;
          } else {
            ishullface = 1; // It's not a real face.
          }
          if (out == (tetgenio *) NULL) {
            fprintf(outfile, " %d", !ishullface ? (vedgecount + shift) : -1); 
          } else {
            vfacet->elist[index++] = !ishullface ? (vedgecount + shift) : -1;
          }
          // Save the V-facet index in this tet at this edge.
          eidxs = &(fidxs[4]);
          eidxs[ver2edge[spintet.ver]] = vfacecount;
          // Go to the next face.
          fnextself(spintet);
          if (spintet.tet == firsttet.tet) break;
        } // while (1)
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, "\n");
        }
        vfacecount++;
      } // if (spintet.tet == worktet.tet)
    } // if (i = 0; i < 6; i++)
    tetloop.tet = tetrahedrontraverse();
  }

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }

  // Output Voronoi cells to .v.cell file.
  if (out == (tetgenio *) NULL) {
    strcpy(outfilename, b->outfilename);
    strcat(outfilename, ".v.cell");
  }
  
  if (!b->quiet) {
    if (out == (tetgenio *) NULL) {
      printf("Writing %s.\n", outfilename);
    } else {
      printf("Writing Voronoi cells.\n");
    }
  }

  if (out == (tetgenio *) NULL) {
    outfile = fopen(outfilename, "w");
    if (outfile == (FILE *) NULL) {
      printf("File I/O Error:  Cannot create file %s.\n", outfilename);
      terminatetetgen(this, 3);
    }
    // Number of Voronoi cells.
    fprintf(outfile, "%ld\n", points->items - unuverts - dupverts);
  } else {
    out->numberofvcells = points->items - unuverts - dupverts;
    out->vcelllist = new int*[out->numberofvcells];
    if (out->vcelllist == (int **) NULL) {
      terminatetetgen(this, 1);
    }
  }

  // Output Voronoi cells.
  tetlist = cavetetlist;
  ptlist = cavetetvertlist;
  points->traversalinit();
  ploop = pointtraverse();
  vpointcount = 0;
  while (ploop != (point) NULL) {
    if ((pointtype(ploop) != UNUSEDVERTEX) &&
        (pointtype(ploop) != DUPLICATEDVERTEX) &&
        (pointtype(ploop) != NREGULARVERTEX)) { 
      getvertexstar(1, ploop, tetlist, ptlist, NULL);
      // Mark all vertices. Check if it is a hull vertex.
      ishullvert = 0;
      for (i = 0; i < ptlist->objects; i++) {
        neipt = * (point *) fastlookup(ptlist, i);
        if (neipt != dummypoint) {
          pinfect(neipt);
        } else {
          ishullvert = 1;
        }
      }
      tcount = (int) ptlist->objects;
      if (out == (tetgenio *) NULL) {
        fprintf(outfile, "%4d  %-2d ", vpointcount + shift, tcount);
      } else {
        arraysize = tcount;
        vertarray = new int[arraysize + 1];
        out->vcelllist[vpointcount] = vertarray;
        vertarray[0] = tcount;
        index = 1;
      }
      // List Voronoi facets bounding this cell.
      for (i = 0; i < tetlist->objects; i++) {
        worktet = * (triface *) fastlookup(tetlist, i);
        // Let 'worktet' be [a,b,c,d] where d = ploop.
        for (j = 0; j < 3; j++) {
          neipt = org(worktet); // neipt is a, or b, or c
          // Skip the dummypoint.
          if (neipt != dummypoint) {
            if (pinfected(neipt)) {
              // It's not processed yet.
              puninfect(neipt);
              // Go to the DT edge [a,d], or [b,d], or [c,d]. 
              esym(worktet, spintet);
              enextself(spintet);
              // Get the V-face dual to this edge.
              eidxs = (int *) spintet.tet[11];
              vfacecount = eidxs[4 + ver2edge[spintet.ver]];
              if (out == (tetgenio *) NULL) {
                fprintf(outfile, " %d", vfacecount + shift);
              } else {
                vertarray[index++] = vfacecount + shift;
              }
            }
          }
          enextself(worktet);
        } // j
      } // i
      if (ishullvert) {
        // Add a hull facet (-1) to the facet list.
        if (out == (tetgenio *) NULL) {
          fprintf(outfile, " -1");
        } else {
          vertarray[index++] = -1;
        }
      }
      if (out == (tetgenio *) NULL) {
        fprintf(outfile, "\n");
      }
      tetlist->restart();
      ptlist->restart();
      vpointcount++;
    }
    ploop = pointtraverse();
  }

  // Delete the space for face/edge indices.
  delete [] indexarray;

  if (out == (tetgenio *) NULL) {
    fprintf(outfile, "# Generated by %s\n", b->commandline);
    fclose(outfile);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outsmesh()    Write surface mesh to a .smesh file, which can be read and  //
//               tetrahedralized by TetGen.                                  //
//                                                                           //
// You can specify a filename (without suffix) in 'smfilename'. If you don't //
// supply a filename (let smfilename be NULL), the default name stored in    //
// 'tetgenbehavior' will be used.                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outsmesh(char* smfilename)
{
  FILE *outfile;
  char nodfilename[FILENAMESIZE];
  char smefilename[FILENAMESIZE];
  face faceloop;
  point p1, p2, p3;
  int firstindex, shift;
  int bmark;
  int faceid, marker;
  int i;

  if (smfilename != (char *) NULL && smfilename[0] != '\0') {
    strcpy(smefilename, smfilename);
  } else if (b->outfilename[0] != '\0') {
    strcpy(smefilename, b->outfilename);
  } else {
    strcpy(smefilename, "unnamed");
  }
  strcpy(nodfilename, smefilename);
  strcat(smefilename, ".smesh");
  strcat(nodfilename, ".node");

  if (!b->quiet) {
    printf("Writing %s.\n", smefilename);
  }
  outfile = fopen(smefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("File I/O Error:  Cannot create file %s.\n", smefilename);
    return;
  }

  // Determine the first index (0 or 1).
  firstindex = b->zeroindex ? 0 : in->firstnumber;
  shift = 0; // Default no shiftment.
  if ((in->firstnumber == 1) && (firstindex == 0)) {
    shift = 1; // Shift the output indices by 1.
  }

  fprintf(outfile, "# %s.  TetGen's input file.\n", smefilename);
  fprintf(outfile, "\n# part 1: node list.\n");
  fprintf(outfile, "0  3  0  0  # nodes are found in %s.\n", nodfilename);

  marker = 0; // avoid compile warning.
  bmark = !b->nobound && in->facetmarkerlist;  

  fprintf(outfile, "\n# part 2: facet list.\n");
  // Number of facets, boundary marker.
  fprintf(outfile, "%ld  %d\n", subfaces->items, bmark);
  
  subfaces->traversalinit();
  faceloop.sh = shellfacetraverse(subfaces);
  while (faceloop.sh != (shellface *) NULL) {
    p1 = sorg(faceloop);
    p2 = sdest(faceloop);
    p3 = sapex(faceloop);
    if (bmark) {
      faceid = shellmark(faceloop) - 1;
      if (faceid >= 0) { 
        marker = in->facetmarkerlist[faceid];
      } else {
        marker = 0; // This subface must be added manually later.
      }
    }
    fprintf(outfile, "3    %4d  %4d  %4d", pointmark(p1) - shift,
            pointmark(p2) - shift, pointmark(p3) - shift);
    if (bmark) {
      fprintf(outfile, "    %d", marker);
    }
    fprintf(outfile, "\n");
    faceloop.sh = shellfacetraverse(subfaces);
  }

  // Copy input holelist.
  fprintf(outfile, "\n# part 3: hole list.\n");
  fprintf(outfile, "%d\n", in->numberofholes);
  for (i = 0; i < in->numberofholes; i++) {
    fprintf(outfile, "%d  %g  %g  %g\n", i + in->firstnumber,
            in->holelist[i * 3], in->holelist[i * 3 + 1],
            in->holelist[i * 3 + 2]);
  }

  // Copy input regionlist.
  fprintf(outfile, "\n# part 4: region list.\n");
  fprintf(outfile, "%d\n", in->numberofregions);
  for (i = 0; i < in->numberofregions; i++) {
    fprintf(outfile, "%d  %g  %g  %g  %d  %g\n", i + in->firstnumber,
            in->regionlist[i * 5], in->regionlist[i * 5 + 1],
            in->regionlist[i * 5 + 2], (int) in->regionlist[i * 5 + 3],
            in->regionlist[i * 5 + 4]);
  }

  fprintf(outfile, "# Generated by %s\n", b->commandline);
  fclose(outfile);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outmesh2medit()    Write mesh to a .mesh file, which can be read and      //
//                    rendered by Medit (a free mesh viewer from INRIA).     //
//                                                                           //
// You can specify a filename (without suffix) in 'mfilename'.  If you don't //
// supply a filename (let mfilename be NULL), the default name stored in     //
// 'tetgenbehavior' will be used. The output file will have the suffix .mesh.//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outmesh2medit(char* mfilename)
{
  FILE *outfile;
  char mefilename[FILENAMESIZE];
  tetrahedron* tetptr;
  triface tface, tsymface;
  face segloop, checkmark;
  point ptloop, p1, p2, p3, p4;
  long ntets, faces;
  int pointnumber;
  int faceid, marker;
  int i;

  if (mfilename != (char *) NULL && mfilename[0] != '\0') {
    strcpy(mefilename, mfilename);
  } else if (b->outfilename[0] != '\0') {
    strcpy(mefilename, b->outfilename);
  } else {
    strcpy(mefilename, "unnamed");
  }
  strcat(mefilename, ".mesh");

  if (!b->quiet) {
    printf("Writing %s.\n", mefilename);
  }
  outfile = fopen(mefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("File I/O Error:  Cannot create file %s.\n", mefilename);
    return;
  }

  fprintf(outfile, "MeshVersionFormatted 1\n");
  fprintf(outfile, "\n");
  fprintf(outfile, "Dimension\n");
  fprintf(outfile, "3\n");
  fprintf(outfile, "\n");

  fprintf(outfile, "\n# Set of mesh vertices\n");
  fprintf(outfile, "Vertices\n");
  fprintf(outfile, "%ld\n", points->items);

  points->traversalinit();
  ptloop = pointtraverse();
  pointnumber = 1;                        // Medit need start number form 1.
  while (ptloop != (point) NULL) {
    // Point coordinates.
    fprintf(outfile, "%.17g  %.17g  %.17g", ptloop[0], ptloop[1], ptloop[2]);
    if (in->numberofpointattributes > 0) {
      // Write an attribute, ignore others if more than one.
      fprintf(outfile, "  %.17g\n", ptloop[3]);
    } else {
      fprintf(outfile, "    0\n");
    }
    setpointmark(ptloop, pointnumber);
    ptloop = pointtraverse();
    pointnumber++;
  }

  // Compute the number of faces.
  ntets = tetrahedrons->items - hullsize;
  faces = (ntets * 4l + hullsize) / 2l;

  fprintf(outfile, "\n# Set of Triangles\n");
  fprintf(outfile, "Triangles\n");
  fprintf(outfile, "%ld\n", faces);

  tetrahedrons->traversalinit();
  tface.tet = tetrahedrontraverse();
  while (tface.tet != (tetrahedron *) NULL) {
    for (tface.ver = 0; tface.ver < 4; tface.ver ++) {
      fsym(tface, tsymface);
      if (ishulltet(tsymface) || 
          (elemindex(tface.tet) < elemindex(tsymface.tet))) {
        p1 = org (tface);
        p2 = dest(tface);
        p3 = apex(tface);
        fprintf(outfile, "%5d  %5d  %5d",
                pointmark(p1), pointmark(p2), pointmark(p3));
        // Check if it is a subface.
        tspivot(tface, checkmark);
        if (checkmark.sh == NULL) {
          marker = 0;  // It is an inner face. It's marker is 0.
        } else {
          if (in->facetmarkerlist) {
            // The facet marker is given, get it.
            faceid = shellmark(checkmark) - 1;
            marker = in->facetmarkerlist[faceid];
          } else {
            marker = 1; // The default marker for subface is 1.
          }
        }
        fprintf(outfile, "    %d\n", marker);
      }
    }
    tface.tet = tetrahedrontraverse();
  }

  fprintf(outfile, "\n# Set of Tetrahedra\n");
  fprintf(outfile, "Tetrahedra\n");
  fprintf(outfile, "%ld\n", ntets);

  tetrahedrons->traversalinit();
  tetptr = tetrahedrontraverse();
  while (tetptr != (tetrahedron *) NULL) {
    if (!b->reversetetori) {
      p1 = (point) tetptr[4];
      p2 = (point) tetptr[5];
    } else {
      p1 = (point) tetptr[5];
      p2 = (point) tetptr[4];
    }
    p3 = (point) tetptr[6];
    p4 = (point) tetptr[7];
    fprintf(outfile, "%5d  %5d  %5d  %5d",
            pointmark(p1), pointmark(p2), pointmark(p3), pointmark(p4));
    if (numelemattrib > 0) {
      fprintf(outfile, "  %.17g", elemattribute(tetptr, 0));
    } else {
      fprintf(outfile, "  0");
    }
    fprintf(outfile, "\n");
    tetptr = tetrahedrontraverse();
  }

  fprintf(outfile, "\nCorners\n");
  fprintf(outfile, "%d\n", in->numberofpoints);

  for (i = 0; i < in->numberofpoints; i++) {
    fprintf(outfile, "%4d\n", i + 1);
  }

  if (b->plc || b->refine) {
    fprintf(outfile, "\nEdges\n");
    fprintf(outfile, "%ld\n", subsegs->items);

    subsegs->traversalinit();
    segloop.sh = shellfacetraverse(subsegs);
    while (segloop.sh != (shellface *) NULL) {
      p1 = sorg(segloop);
      p2 = sdest(segloop);
      fprintf(outfile, "%5d  %5d", pointmark(p1), pointmark(p2));
      marker = shellmark(segloop);
      fprintf(outfile, "    %d\n", marker);
      segloop.sh = shellfacetraverse(subsegs);
    }
  }

  fprintf(outfile, "\nEnd\n");
  fclose(outfile);
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// outmesh2vtk()    Save mesh to file in VTK Legacy format.                  //
//                                                                           //
// This function was contributed by Bryn Llyod from ETH, 2007.               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenmesh::outmesh2vtk(char* ofilename)
{
  FILE *outfile;
  char vtkfilename[FILENAMESIZE];
  point pointloop, p1, p2, p3, p4;
  tetrahedron* tptr;
  double x, y, z;
  int n1, n2, n3, n4;
  int nnodes = 4;
  int celltype = 10;

  if (b->order == 2) {
    printf("  Write VTK not implemented for order 2 elements \n");
    return;
  }

  int NEL = tetrahedrons->items - hullsize;
  int NN = points->items;

  if (ofilename != (char *) NULL && ofilename[0] != '\0') {
    strcpy(vtkfilename, ofilename);
  } else if (b->outfilename[0] != '\0') {
    strcpy(vtkfilename, b->outfilename);
  } else {
    strcpy(vtkfilename, "unnamed");
  }
  strcat(vtkfilename, ".vtk");

  if (!b->quiet) {
    printf("Writing %s.\n", vtkfilename);
  }
  outfile = fopen(vtkfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("File I/O Error:  Cannot create file %s.\n", vtkfilename);
    return;
  }

  //always write big endian
  //bool ImALittleEndian = !testIsBigEndian();

  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Unstructured Grid\n");
  fprintf(outfile, "ASCII\n"); // BINARY
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(outfile, "POINTS %d double\n", NN);

  points->traversalinit();
  pointloop = pointtraverse();
  for(int id=0; id<NN && pointloop != (point) NULL; id++){
    x = pointloop[0];
    y = pointloop[1];
    z = pointloop[2];
    fprintf(outfile, "%.17g %.17g %.17g\n", x, y, z);
    pointloop = pointtraverse();
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "CELLS %d %d\n", NEL, NEL*(4+1));
  //NEL rows, each has 1 type id + 4 node id's
 
  tetrahedrons->traversalinit();
  tptr = tetrahedrontraverse();
  //elementnumber = firstindex; // in->firstnumber;
  while (tptr != (tetrahedron *) NULL) {
    if (!b->reversetetori) {
      p1 = (point) tptr[4];
      p2 = (point) tptr[5];
    } else {
      p1 = (point) tptr[5];
      p2 = (point) tptr[4];
    }
    p3 = (point) tptr[6];
    p4 = (point) tptr[7];
    n1 = pointmark(p1) - in->firstnumber;
    n2 = pointmark(p2) - in->firstnumber;
    n3 = pointmark(p3) - in->firstnumber;
    n4 = pointmark(p4) - in->firstnumber;
    fprintf(outfile, "%d  %4d %4d %4d %4d\n", nnodes, n1, n2, n3, n4);
    tptr = tetrahedrontraverse();
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "CELL_TYPES %d\n", NEL);
  for(int tid=0; tid<NEL; tid++){
    fprintf(outfile, "%d\n", celltype);
  }
  fprintf(outfile, "\n");

  if (numelemattrib > 0) {
    // Output tetrahedra region attributes.
    fprintf(outfile, "CELL_DATA %d\n", NEL);
    fprintf(outfile, "SCALARS cell_scalars int 1\n");
    fprintf(outfile, "LOOKUP_TABLE default\n");
    tetrahedrons->traversalinit();
    tptr = tetrahedrontraverse();
    while (tptr != (tetrahedron *) NULL) {
      fprintf(outfile, "%d\n", (int) elemattribute(tptr, numelemattrib - 1));
      tptr = tetrahedrontraverse();
    }
    fprintf(outfile, "\n");
  }

  fclose(outfile);
}

////                                                                       ////
////                                                                       ////
//// output_cxx ///////////////////////////////////////////////////////////////

