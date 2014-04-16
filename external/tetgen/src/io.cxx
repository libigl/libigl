#include "../tetgen.h"
//// io_cxx ///////////////////////////////////////////////////////////////////
////                                                                       ////
////                                                                       ////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_node_call()    Read a list of points from a file.                    //
//                                                                           //
// 'infile' is the file handle contains the node list.  It may point to a    //
// .node, or .poly or .smesh file. 'markers' indicates each node contains an //
// additional marker (integer) or not. 'uvflag' indicates each node contains //
// u,v coordinates or not. It is reuqired by a PSC. 'infilename' is the name //
// of the file being read,  it is only used in error messages.               //
//                                                                           //
// The 'firstnumber' (0 or 1) is automatically determined by the number of   //
// the first index of the first point.                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_node_call(FILE* infile, int markers, int uvflag, 
                              char* infilename)
{
  char inputline[INPUTLINESIZE];
  char *stringptr;
  REAL x, y, z, attrib;
  int firstnode, currentmarker;
  int index, attribindex;
  int i, j;

  // Initialize 'pointlist', 'pointattributelist', and 'pointmarkerlist'.
  pointlist = new REAL[numberofpoints * 3];
  if (pointlist == (REAL *) NULL) {
    terminatetetgen(NULL, 1);
  }
  if (numberofpointattributes > 0) {
    pointattributelist = new REAL[numberofpoints * numberofpointattributes];
    if (pointattributelist == (REAL *) NULL) {
      terminatetetgen(NULL, 1);
    }
  }
  if (markers) {
    pointmarkerlist = new int[numberofpoints];
    if (pointmarkerlist == (int *) NULL) {
      terminatetetgen(NULL, 1);
    }
  }
  if (uvflag) {
    pointparamlist = new pointparam[numberofpoints];
    if (pointparamlist == NULL) {
      terminatetetgen(NULL, 1);
    }
  }

  // Read the point section.
  index = 0;
  attribindex = 0;
  for (i = 0; i < numberofpoints; i++) {
    stringptr = readnumberline(inputline, infile, infilename);
    if (useindex) {
      if (i == 0) {
        firstnode = (int) strtol (stringptr, &stringptr, 0);
        if ((firstnode == 0) || (firstnode == 1)) {
          firstnumber = firstnode;
        }
      }
      stringptr = findnextnumber(stringptr);
    } // if (useindex)
    if (*stringptr == '\0') {
      printf("Error:  Point %d has no x coordinate.\n", firstnumber + i);
      break;
    }
    x = (REAL) strtod(stringptr, &stringptr);
    stringptr = findnextnumber(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Point %d has no y coordinate.\n", firstnumber + i);
      break;
    }
    y = (REAL) strtod(stringptr, &stringptr);
    if (mesh_dim == 3) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d has no z coordinate.\n", firstnumber + i);
        break;
      }
      z = (REAL) strtod(stringptr, &stringptr);
    } else {
      z = 0.0; // mesh_dim == 2;
    }
    pointlist[index++] = x;
    pointlist[index++] = y;
    pointlist[index++] = z;
    // Read the point attributes.
    for (j = 0; j < numberofpointattributes; j++) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        attrib = 0.0;
      } else {
        attrib = (REAL) strtod(stringptr, &stringptr);
      }
      pointattributelist[attribindex++] = attrib;
    }
    if (markers) {
      // Read a point marker.
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        currentmarker = 0;
      } else {
        currentmarker = (int) strtol (stringptr, &stringptr, 0);
      }
      pointmarkerlist[i] = currentmarker;
    }
    if (uvflag) {
      // Read point paramteters.
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d has no uv[0].\n", firstnumber + i);
        break;
      }
      pointparamlist[i].uv[0] = (REAL) strtod(stringptr, &stringptr);
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d has no uv[1].\n", firstnumber + i);
        break;
      }
      pointparamlist[i].uv[1] = (REAL) strtod(stringptr, &stringptr);
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d has no tag.\n", firstnumber + i);
        break;
      }
      pointparamlist[i].tag = (int) strtol (stringptr, &stringptr, 0);
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d has no type.\n", firstnumber + i);
        break;
      }
      pointparamlist[i].type = (int) strtol (stringptr, &stringptr, 0);
      if ((pointparamlist[i].type < 0) || (pointparamlist[i].type > 2)) {
        printf("Error:  Point %d has an invalid type.\n", firstnumber + i);
        break;
      }
    }
  }
  if (i < numberofpoints) {
    // Failed to read points due to some error.
    delete [] pointlist;
    pointlist = (REAL *) NULL;
    if (markers) {
      delete [] pointmarkerlist;
      pointmarkerlist = (int *) NULL;
    }
    if (numberofpointattributes > 0) {
      delete [] pointattributelist;
      pointattributelist = (REAL *) NULL;
    }
    if (uvflag) {
      delete [] pointparamlist;
      pointparamlist = NULL;
    }
    numberofpoints = 0;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_node()    Load a list of points from a .node file.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_node(char* filebasename)
{
  FILE *infile;
  char innodefilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  bool okflag;
  int markers;
  int uvflag; // for psc input.

  // Assembling the actual file names we want to open.
  strcpy(innodefilename, filebasename);
  strcat(innodefilename, ".node");

  // Try to open a .node file.
  infile = fopen(innodefilename, "r");
  if (infile == (FILE *) NULL) {
    printf("  Cannot access file %s.\n", innodefilename);
    return false;
  }
  printf("Opening %s.\n", innodefilename); 

  // Set initial flags.
  mesh_dim = 3;
  numberofpointattributes = 0;  // no point attribute.
  markers = 0;  // no boundary marker.
  uvflag = 0; // no uv parameters (required by a PSC). 

  // Read the first line of the file.
  stringptr = readnumberline(inputline, infile, innodefilename);
  // Does this file contain an index column?
  stringptr = strstr(inputline, "rbox");
  if (stringptr == NULL) {
    // Read number of points, number of dimensions, number of point
    //   attributes, and number of boundary markers. 
    stringptr = inputline;
    numberofpoints = (int) strtol (stringptr, &stringptr, 0);
    stringptr = findnextnumber(stringptr);
    if (*stringptr != '\0') {
      mesh_dim = (int) strtol (stringptr, &stringptr, 0);
    }
    stringptr = findnextnumber(stringptr);
    if (*stringptr != '\0') {
      numberofpointattributes = (int) strtol (stringptr, &stringptr, 0);
    }
    stringptr = findnextnumber(stringptr);
    if (*stringptr != '\0') {
      markers = (int) strtol (stringptr, &stringptr, 0);
    }
    stringptr = findnextnumber(stringptr);
    if (*stringptr != '\0') {
      uvflag = (int) strtol (stringptr, &stringptr, 0);
    }
  } else {
    // It is a rbox (qhull) input file.
    stringptr = inputline;
    // Get the dimension.
    mesh_dim = (int) strtol (stringptr, &stringptr, 0);
    // Get the number of points.
    stringptr = readnumberline(inputline, infile, innodefilename);
    numberofpoints = (int) strtol (stringptr, &stringptr, 0);
    // There is no index column.
    useindex = 0;
  }

  // Load the list of nodes.
  okflag = load_node_call(infile, markers, uvflag, innodefilename);

  fclose(infile);
  return okflag;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_edge()    Load a list of edges from a .edge file.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_edge(char* filebasename)
{
  FILE *infile;
  char inedgefilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int markers, corner;
  int index;
  int i, j;

  strcpy(inedgefilename, filebasename);
  strcat(inedgefilename, ".edge");

  infile = fopen(inedgefilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", inedgefilename);
  } else {
    //printf("  Cannot access file %s.\n", inedgefilename);
    return false;
  }

  // Read number of boundary edges.
  stringptr = readnumberline(inputline, infile, inedgefilename);
  numberofedges = (int) strtol (stringptr, &stringptr, 0);
  if (numberofedges > 0) {
    edgelist = new int[numberofedges * 2];
    if (edgelist == (int *) NULL) {
      terminatetetgen(NULL, 1);
    }
    stringptr = findnextnumber(stringptr);
    if (*stringptr == '\0') {
      markers = 0;  // Default value.
    } else {
      markers = (int) strtol (stringptr, &stringptr, 0);
    }
    if (markers > 0) {
      edgemarkerlist = new int[numberofedges];
    }
  }

  // Read the list of edges.
  index = 0;
  for (i = 0; i < numberofedges; i++) {
    // Read edge index and the edge's two endpoints.
    stringptr = readnumberline(inputline, infile, inedgefilename);
    for (j = 0; j < 2; j++) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Edge %d is missing vertex %d in %s.\n",
               i + firstnumber, j + 1, inedgefilename);
        terminatetetgen(NULL, 1);
      }
      corner = (int) strtol(stringptr, &stringptr, 0);
      if (corner < firstnumber || corner >= numberofpoints + firstnumber) {
        printf("Error:  Edge %d has an invalid vertex index.\n",
               i + firstnumber);
        terminatetetgen(NULL, 1);
      }
      edgelist[index++] = corner;
    }
    if (numberofcorners == 10) {
      // Skip an extra vertex (generated by a previous -o2 option).
      stringptr = findnextnumber(stringptr);
    }
    // Read the edge marker if it has.
    if (markers) {
      stringptr = findnextnumber(stringptr);
      edgemarkerlist[i] = (int) strtol(stringptr, &stringptr, 0);
    }
  }

  fclose(infile);
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_face()    Load a list of faces (triangles) from a .face file.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_face(char* filebasename)
{
  FILE *infile;
  char infilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  REAL attrib;
  int markers, corner;
  int index;
  int i, j;

  strcpy(infilename, filebasename);
  strcat(infilename, ".face");

  infile = fopen(infilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", infilename);
  } else {
    return false;
  }

  // Read number of faces, boundary markers.
  stringptr = readnumberline(inputline, infile, infilename);
  numberoftrifaces = (int) strtol (stringptr, &stringptr, 0);
  stringptr = findnextnumber(stringptr);
  if (mesh_dim == 2) {
    // Skip a number.
    stringptr = findnextnumber(stringptr);
  }
  if (*stringptr == '\0') {
    markers = 0;  // Default there is no marker per face.
  } else {
    markers = (int) strtol (stringptr, &stringptr, 0);
  }
  if (numberoftrifaces > 0) {
    trifacelist = new int[numberoftrifaces * 3];
    if (trifacelist == (int *) NULL) {
      terminatetetgen(NULL, 1);
    }
    if (markers) {
      trifacemarkerlist = new int[numberoftrifaces];
      if (trifacemarkerlist == (int *) NULL) {
        terminatetetgen(NULL, 1);
      }
    }
  }

  // Read the list of faces.
  index = 0;
  for (i = 0; i < numberoftrifaces; i++) {
    // Read face index and the face's three corners.
    stringptr = readnumberline(inputline, infile, infilename);
    for (j = 0; j < 3; j++) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Face %d is missing vertex %d in %s.\n",
               i + firstnumber, j + 1, infilename);
        terminatetetgen(NULL, 1);
      }
      corner = (int) strtol(stringptr, &stringptr, 0);
      if (corner < firstnumber || corner >= numberofpoints + firstnumber) {
        printf("Error:  Face %d has an invalid vertex index.\n",
               i + firstnumber);
        terminatetetgen(NULL, 1);
      }
      trifacelist[index++] = corner;
    }
    if (numberofcorners == 10) {
      // Skip 3 extra vertices (generated by a previous -o2 option).
      for (j = 0; j < 3; j++) {
        stringptr = findnextnumber(stringptr);
      }
    }
    // Read the boundary marker if it exists.
    if (markers) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        attrib = 0.0;
      } else {
        attrib = (REAL) strtod(stringptr, &stringptr);
      }
      trifacemarkerlist[i] = (int) attrib;
    }
  }

  fclose(infile);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_tet()    Load a list of tetrahedra from a .ele file.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_tet(char* filebasename)
{
  FILE *infile;
  char infilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  REAL attrib;
  int corner;
  int index, attribindex;
  int i, j;

  strcpy(infilename, filebasename);
  strcat(infilename, ".ele");

  infile = fopen(infilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", infilename);
  } else {
    return false;
  }

  // Read number of elements, number of corners (4 or 10), number of
  //   element attributes.
  stringptr = readnumberline(inputline, infile, infilename);
  numberoftetrahedra = (int) strtol (stringptr, &stringptr, 0);
  if (numberoftetrahedra <= 0) {
    printf("Error:  Invalid number of tetrahedra.\n");
    fclose(infile);
    return false;
  }
  stringptr = findnextnumber(stringptr);
  if (*stringptr == '\0') {
    numberofcorners = 4;  // Default read 4 nodes per element.
  } else {
    numberofcorners = (int) strtol(stringptr, &stringptr, 0);
  }
  stringptr = findnextnumber(stringptr);
  if (*stringptr == '\0') {
    numberoftetrahedronattributes = 0; // Default no attribute.
  } else {
    numberoftetrahedronattributes = (int) strtol(stringptr, &stringptr, 0);
  }
  if (numberofcorners != 4 && numberofcorners != 10) {
    printf("Error:  Wrong number of corners %d (should be 4 or 10).\n", 
           numberofcorners);
    fclose(infile);
    return false;
  }

  // Allocate memory for tetrahedra.
  tetrahedronlist = new int[numberoftetrahedra * numberofcorners]; 
  if (tetrahedronlist == (int *) NULL) {
    terminatetetgen(NULL, 1);
  }
  // Allocate memory for output tetrahedron attributes if necessary.
  if (numberoftetrahedronattributes > 0) {
    tetrahedronattributelist = new REAL[numberoftetrahedra *
                                        numberoftetrahedronattributes];
    if (tetrahedronattributelist == (REAL *) NULL) {
      terminatetetgen(NULL, 1);
    }
  }

  // Read the list of tetrahedra.
  index = 0;
  attribindex = 0;
  for (i = 0; i < numberoftetrahedra; i++) {
    // Read tetrahedron index and the tetrahedron's corners.
    stringptr = readnumberline(inputline, infile, infilename);
    for (j = 0; j < numberofcorners; j++) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Tetrahedron %d is missing vertex %d in %s.\n",
               i + firstnumber, j + 1, infilename);
        terminatetetgen(NULL, 1);
      }
      corner = (int) strtol(stringptr, &stringptr, 0);
      if (corner < firstnumber || corner >= numberofpoints + firstnumber) {
        printf("Error:  Tetrahedron %d has an invalid vertex index.\n",
               i + firstnumber);
        terminatetetgen(NULL, 1);
      }
      tetrahedronlist[index++] = corner;
    }
    // Read the tetrahedron's attributes.
    for (j = 0; j < numberoftetrahedronattributes; j++) {
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        attrib = 0.0;
      } else {
        attrib = (REAL) strtod(stringptr, &stringptr);
      }
      tetrahedronattributelist[attribindex++] = attrib;
    }
  }

  fclose(infile);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_vol()    Load a list of volume constraints from a .vol file.         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_vol(char* filebasename)
{
  FILE *infile;
  char inelefilename[FILENAMESIZE];
  char infilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr; 
  REAL volume;
  int volelements;
  int i;

  strcpy(infilename, filebasename);
  strcat(infilename, ".vol");

  infile = fopen(infilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", infilename);
  } else {
    return false;
  }

  // Read number of tetrahedra.
  stringptr = readnumberline(inputline, infile, infilename);
  volelements = (int) strtol (stringptr, &stringptr, 0);
  if (volelements != numberoftetrahedra) {
    strcpy(inelefilename, filebasename);
    strcat(infilename, ".ele");
    printf("Warning:  %s and %s disagree on number of tetrahedra.\n",
           inelefilename, infilename);
    fclose(infile);
    return false;
  }

  tetrahedronvolumelist = new REAL[volelements];
  if (tetrahedronvolumelist == (REAL *) NULL) {
    terminatetetgen(NULL, 1);
  }

  // Read the list of volume constraints.
  for (i = 0; i < volelements; i++) {
    stringptr = readnumberline(inputline, infile, infilename);
    stringptr = findnextnumber(stringptr);
    if (*stringptr == '\0') {
      volume = -1.0; // No constraint on this tetrahedron.
    } else {
      volume = (REAL) strtod(stringptr, &stringptr);
    }
    tetrahedronvolumelist[i] = volume;
  }

  fclose(infile);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_var()    Load constraints applied on facets, segments, and nodes     //
//               from a .var file.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_var(char* filebasename)
{
  FILE *infile;
  char varfilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int index;
  int i;

  // Variant constraints are saved in file "filename.var".
  strcpy(varfilename, filebasename);
  strcat(varfilename, ".var");
  infile = fopen(varfilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", varfilename);
  } else {
    return false;
  }

  // Read the facet constraint section.
  stringptr = readnumberline(inputline, infile, varfilename);
  if (*stringptr != '\0') {
    numberoffacetconstraints = (int) strtol (stringptr, &stringptr, 0);
  } else {
    numberoffacetconstraints = 0;
  }
  if (numberoffacetconstraints > 0) {
    // Initialize 'facetconstraintlist'.
    facetconstraintlist = new REAL[numberoffacetconstraints * 2];
    index = 0;
    for (i = 0; i < numberoffacetconstraints; i++) {
      stringptr = readnumberline(inputline, infile, varfilename);
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  facet constraint %d has no facet marker.\n",
               firstnumber + i);
        break;
      } else {
        facetconstraintlist[index++] = (REAL) strtod(stringptr, &stringptr);
      }
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  facet constraint %d has no maximum area bound.\n",
               firstnumber + i);
        break;
      } else {
        facetconstraintlist[index++] = (REAL) strtod(stringptr, &stringptr);
      }
    }
    if (i < numberoffacetconstraints) {
      // This must be caused by an error.
      fclose(infile);
      return false;
    }
  }

  // Read the segment constraint section.
  stringptr = readnumberline(inputline, infile, varfilename);
  if (*stringptr != '\0') {
    numberofsegmentconstraints = (int) strtol (stringptr, &stringptr, 0);
  } else {
    numberofsegmentconstraints = 0;
  }
  if (numberofsegmentconstraints > 0) {
    // Initialize 'segmentconstraintlist'.
    segmentconstraintlist = new REAL[numberofsegmentconstraints * 3];
    index = 0;
    for (i = 0; i < numberofsegmentconstraints; i++) {
      stringptr = readnumberline(inputline, infile, varfilename);
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  segment constraint %d has no frist endpoint.\n",
               firstnumber + i);
        break;
      } else {
        segmentconstraintlist[index++] = (REAL) strtod(stringptr, &stringptr);
      }
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  segment constraint %d has no second endpoint.\n",
               firstnumber + i);
        break;
      } else {
        segmentconstraintlist[index++] = (REAL) strtod(stringptr, &stringptr);
      }
      stringptr = findnextnumber(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  segment constraint %d has no maximum length bound.\n",
               firstnumber + i);
        break;
      } else {
        segmentconstraintlist[index++] = (REAL) strtod(stringptr, &stringptr);
      }
    }
    if (i < numberofsegmentconstraints) {
      // This must be caused by an error.
      fclose(infile);
      return false;
    }
  }

  fclose(infile);
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_mtr()    Load a size specification map from a .mtr file.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_mtr(char* filebasename)
{
  FILE *infile;
  char mtrfilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr;
  REAL mtr;
  int ptnum;
  int mtrindex;
  int i, j;

  strcpy(mtrfilename, filebasename);
  strcat(mtrfilename, ".mtr");
  infile = fopen(mtrfilename, "r");
  if (infile != (FILE *) NULL) {
    printf("Opening %s.\n", mtrfilename);
  } else {
    return false;
  }

  // Read the number of points.
  stringptr = readnumberline(inputline, infile, mtrfilename);
  ptnum = (int) strtol (stringptr, &stringptr, 0);
  if (ptnum != numberofpoints) {
    printf("  !! Point numbers are not equal. Ignored.\n");
    fclose(infile);
    return false;
  }
  // Read the number of columns (1, 3, or 6).
  stringptr = findnextnumber(stringptr); // Skip number of points.
  if (*stringptr != '\0') {
    numberofpointmtrs = (int) strtol (stringptr, &stringptr, 0);
  }
  if (numberofpointmtrs == 0) {
    // Column number doesn't match. Set a default number (1).
    numberofpointmtrs = 1;
  }

  // Allocate space for pointmtrlist.
  pointmtrlist = new REAL[numberofpoints * numberofpointmtrs];
  if (pointmtrlist == (REAL *) NULL) {
    terminatetetgen(NULL, 1);
  }
  mtrindex = 0;
  for (i = 0; i < numberofpoints; i++) {
    // Read metrics.
    stringptr = readnumberline(inputline, infile, mtrfilename);
    for (j = 0; j < numberofpointmtrs; j++) {
      if (*stringptr == '\0') {
        printf("Error:  Metric %d is missing value #%d in %s.\n",
               i + firstnumber, j + 1, mtrfilename);
        terminatetetgen(NULL, 1);
      }
      mtr = (REAL) strtod(stringptr, &stringptr);
      pointmtrlist[mtrindex++] = mtr;
      stringptr = findnextnumber(stringptr);
    }
  }

  fclose(infile);
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_poly()    Load a PL complex from a .poly or a .smesh file.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_poly(char* filebasename)
{
  FILE *infile;
  char inpolyfilename[FILENAMESIZE];
  char insmeshfilename[FILENAMESIZE];
  char inputline[INPUTLINESIZE];
  char *stringptr, *infilename;
  int smesh, markers, uvflag, currentmarker;
  int index;
  int i, j, k;

  // Assembling the actual file names we want to open.
  strcpy(inpolyfilename, filebasename);
  strcpy(insmeshfilename, filebasename);
  strcat(inpolyfilename, ".poly");
  strcat(insmeshfilename, ".smesh");

  // First assume it is a .poly file.
  smesh = 0;
  // Try to open a .poly file.
  infile = fopen(inpolyfilename, "r");
  if (infile == (FILE *) NULL) {
    // .poly doesn't exist! Try to open a .smesh file.
    infile = fopen(insmeshfilename, "r");
    if (infile == (FILE *) NULL) {
      printf("  Cannot access file %s and %s.\n",
             inpolyfilename, insmeshfilename);
      return false;
    } else {
      printf("Opening %s.\n", insmeshfilename);
      infilename = insmeshfilename;
    }
    smesh = 1;
  } else {
    printf("Opening %s.\n", inpolyfilename);
    infilename = inpolyfilename;
  }

  // Initialize the default values.
  mesh_dim = 3;  // Three-dimensional coordinates.
  numberofpointattributes = 0;  // no point attribute.
  markers = 0;  // no boundary marker.
  uvflag = 0; // no uv parameters (required by a PSC).

  // Read number of points, number of dimensions, number of point
  //   attributes, and number of boundary markers.
  stringptr = readnumberline(inputline, infile, infilename);
  numberofpoints = (int) strtol (stringptr, &stringptr, 0);
  stringptr = findnextnumber(stringptr);
  if (*stringptr != '\0') {
    mesh_dim = (int) strtol (stringptr, &stringptr, 0);      
  }
  stringptr = findnextnumber(stringptr);
  if (*stringptr != '\0') {
    numberofpointattributes = (int) strtol (stringptr, &stringptr, 0);
  }
  stringptr = findnextnumber(stringptr);
  if (*stringptr != '\0') {
    markers = (int) strtol (stringptr, &stringptr, 0);
  }
  if (*stringptr != '\0') {
    uvflag = (int) strtol (stringptr, &stringptr, 0);
  }

  if (numberofpoints > 0) {
    // Load the list of nodes.
    if (!load_node_call(infile, markers, uvflag, infilename)) {
      fclose(infile);
      return false;
    }
  } else {
    // If the .poly or .smesh file claims there are zero points, that
    //   means the points should be read from a separate .node file.
    if (!load_node(filebasename)) {
      fclose(infile);
      return false;
    }
  }

  if ((mesh_dim != 3) && (mesh_dim != 2)) {
    printf("Input error:  TetGen only works for 2D & 3D point sets.\n");
    fclose(infile);
    return false;
  }
  if (numberofpoints < (mesh_dim + 1)) {
    printf("Input error:  TetGen needs at least %d points.\n", mesh_dim + 1);
    fclose(infile);
    return false;
  }

  facet *f;
  polygon *p;

  if (mesh_dim == 3) {

    // Read number of facets and number of boundary markers.
    stringptr = readnumberline(inputline, infile, infilename);
    if (stringptr == NULL) {
      // No facet list, return.
      fclose(infile);
      return true;
    }
    numberoffacets = (int) strtol (stringptr, &stringptr, 0);
    if (numberoffacets <= 0) {
      // No facet list, return.
      fclose(infile);
      return true;
    }
    stringptr = findnextnumber(stringptr);
    if (*stringptr == '\0') {
      markers = 0;  // no boundary marker.
    } else {
      markers = (int) strtol (stringptr, &stringptr, 0);
    }

    // Initialize the 'facetlist', 'facetmarkerlist'.
    facetlist = new facet[numberoffacets];
    if (markers == 1) {
      facetmarkerlist = new int[numberoffacets];
    }

    // Read data into 'facetlist', 'facetmarkerlist'.
    if (smesh == 0) {
      // Facets are in .poly file format.
      for (i = 1; i <= numberoffacets; i++) {
        f = &(facetlist[i - 1]);
        init(f);
        f->numberofholes = 0;
        currentmarker = 0;
        // Read number of polygons, number of holes, and a boundary marker.
        stringptr = readnumberline(inputline, infile, infilename);
        f->numberofpolygons = (int) strtol (stringptr, &stringptr, 0);
        stringptr = findnextnumber(stringptr);
        if (*stringptr != '\0') {
          f->numberofholes = (int) strtol (stringptr, &stringptr, 0);
          if (markers == 1) {
            stringptr = findnextnumber(stringptr);
            if (*stringptr != '\0') {
              currentmarker = (int) strtol(stringptr, &stringptr, 0);
            }
          }
        }
        // Initialize facetmarker if it needs.
        if (markers == 1) {
          facetmarkerlist[i - 1] = currentmarker; 
        }
        // Each facet should has at least one polygon.
        if (f->numberofpolygons <= 0) {
          printf("Error:  Wrong number of polygon in %d facet.\n", i);
          break; 
        }
        // Initialize the 'f->polygonlist'.
        f->polygonlist = new polygon[f->numberofpolygons];
        // Go through all polygons, read in their vertices.
        for (j = 1; j <= f->numberofpolygons; j++) {
          p = &(f->polygonlist[j - 1]);
          init(p);
          // Read number of vertices of this polygon.
          stringptr = readnumberline(inputline, infile, infilename);
          p->numberofvertices = (int) strtol(stringptr, &stringptr, 0);
          if (p->numberofvertices < 1) {
            printf("Error:  Wrong polygon %d in facet %d\n", j, i);
            break;
          }
          // Initialize 'p->vertexlist'.
          p->vertexlist = new int[p->numberofvertices];
          // Read all vertices of this polygon.
          for (k = 1; k <= p->numberofvertices; k++) {
            stringptr = findnextnumber(stringptr);
            if (*stringptr == '\0') {
              // Try to load another non-empty line and continue to read the
              //   rest of vertices.
              stringptr = readnumberline(inputline, infile, infilename);
              if (*stringptr == '\0') {
                printf("Error: Missing %d endpoints of polygon %d in facet %d",
                       p->numberofvertices - k, j, i);
                break;
              }
            }
            p->vertexlist[k - 1] = (int) strtol (stringptr, &stringptr, 0);
          }
        } 
        if (j <= f->numberofpolygons) {
          // This must be caused by an error. However, there're j - 1
          //   polygons have been read. Reset the 'f->numberofpolygon'.
          if (j == 1) {
            // This is the first polygon.
            delete [] f->polygonlist;
          }
          f->numberofpolygons = j - 1;
          // No hole will be read even it exists.
          f->numberofholes = 0;
          break;
        }
        // If this facet has hole pints defined, read them.
        if (f->numberofholes > 0) {
          // Initialize 'f->holelist'.
          f->holelist = new REAL[f->numberofholes * 3];
          // Read the holes' coordinates.
          index = 0;
          for (j = 1; j <= f->numberofholes; j++) {
            stringptr = readnumberline(inputline, infile, infilename);
            for (k = 1; k <= 3; k++) {
              stringptr = findnextnumber(stringptr);
              if (*stringptr == '\0') {
                printf("Error:  Hole %d in facet %d has no coordinates", j, i);
                break;
              }
              f->holelist[index++] = (REAL) strtod (stringptr, &stringptr);
            }
            if (k <= 3) {
              // This must be caused by an error.
              break;
            }
          }
          if (j <= f->numberofholes) {
            // This must be caused by an error.
            break;
          }
        }
      }
      if (i <= numberoffacets) {
        // This must be caused by an error.
        numberoffacets = i - 1;
        fclose(infile);
        return false;
      }
    } else { // poly == 0
      // Read the facets from a .smesh file.
      for (i = 1; i <= numberoffacets; i++) {
        f = &(facetlist[i - 1]);
        init(f);
        // Initialize 'f->facetlist'. In a .smesh file, each facetlist only
        //   contains exactly one polygon, no hole.
        f->numberofpolygons = 1;
        f->polygonlist = new polygon[f->numberofpolygons];
        p = &(f->polygonlist[0]);
        init(p);
        // Read number of vertices of this polygon.
        stringptr = readnumberline(inputline, infile, insmeshfilename);
        p->numberofvertices = (int) strtol (stringptr, &stringptr, 0);
        if (p->numberofvertices < 1) {
          printf("Error:  Wrong number of vertex in facet %d\n", i);
          break;
        }
        // Initialize 'p->vertexlist'.
        p->vertexlist = new int[p->numberofvertices];
        for (k = 1; k <= p->numberofvertices; k++) {
          stringptr = findnextnumber(stringptr);
          if (*stringptr == '\0') {
            // Try to load another non-empty line and continue to read the
            //   rest of vertices.
            stringptr = readnumberline(inputline, infile, infilename);
            if (*stringptr == '\0') {
              printf("Error:  Missing %d endpoints in facet %d",
                     p->numberofvertices - k, i);
              break;
            }
          }
          p->vertexlist[k - 1] = (int) strtol (stringptr, &stringptr, 0);
        }
        if (k <= p->numberofvertices) {
          // This must be caused by an error.
          break;
        }
        // Read facet's boundary marker at last.
        if (markers == 1) {
          stringptr = findnextnumber(stringptr);
          if (*stringptr == '\0') {
            currentmarker = 0;
          } else {
            currentmarker = (int) strtol(stringptr, &stringptr, 0);
          }
          facetmarkerlist[i - 1] = currentmarker;
        }
      }
      if (i <= numberoffacets) {
        // This must be caused by an error.
        numberoffacets = i - 1;
        fclose(infile);
        return false;
      }
    }

    // Read the hole section.
    stringptr = readnumberline(inputline, infile, infilename);
    if (stringptr == NULL) {
      // No hole list, return.
      fclose(infile);
      return true;
    }
    if (*stringptr != '\0') {
      numberofholes = (int) strtol (stringptr, &stringptr, 0);
    } else {
      numberofholes = 0;
    }
    if (numberofholes > 0) {
      // Initialize 'holelist'.
      holelist = new REAL[numberofholes * 3];
      for (i = 0; i < 3 * numberofholes; i += 3) {
        stringptr = readnumberline(inputline, infile, infilename);
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Hole %d has no x coord.\n", firstnumber + (i / 3));
          break;
        } else {
          holelist[i] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Hole %d has no y coord.\n", firstnumber + (i / 3));
          break;
        } else {
          holelist[i + 1] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Hole %d has no z coord.\n", firstnumber + (i / 3));
          break;
        } else {
          holelist[i + 2] = (REAL) strtod(stringptr, &stringptr);
        }
      }
      if (i < 3 * numberofholes) {
        // This must be caused by an error.
        fclose(infile);
        return false;
      }
    }

    // Read the region section.  The 'region' section is optional, if we
    //   don't reach the end-of-file, try read it in.
    stringptr = readnumberline(inputline, infile, NULL);
    if (stringptr != (char *) NULL && *stringptr != '\0') {
      numberofregions = (int) strtol (stringptr, &stringptr, 0);
    } else {
      numberofregions = 0;
    }
    if (numberofregions > 0) {
      // Initialize 'regionlist'.
      regionlist = new REAL[numberofregions * 5];
      index = 0;
      for (i = 0; i < numberofregions; i++) {
        stringptr = readnumberline(inputline, infile, infilename);
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no x coordinate.\n", firstnumber + i);
          break;
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no y coordinate.\n", firstnumber + i);
          break;
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no z coordinate.\n", firstnumber + i);
          break;
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no region attrib.\n", firstnumber + i);
          break;
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findnextnumber(stringptr);
        if (*stringptr == '\0') {
          regionlist[index] = regionlist[index - 1];
        } else {
          regionlist[index] = (REAL) strtod(stringptr, &stringptr);
        }
        index++;
      }
      if (i < numberofregions) {
        // This must be caused by an error.
        fclose(infile);
        return false;
      }
    }

  } else {

    // Read a PSLG from Triangle's poly file.
    assert(mesh_dim == 2);
    // A PSLG is a facet of a PLC.
    numberoffacets = 1;
    // Initialize the 'facetlist'.
    facetlist = new facet[numberoffacets];
    facetmarkerlist = (int *) NULL; // No facet markers.
    f = &(facetlist[0]);
    init(f);
    // Read number of segments.
    stringptr = readnumberline(inputline, infile, infilename);
    // Segments are degenerate polygons.
    f->numberofpolygons = (int) strtol (stringptr, &stringptr, 0);
    if (f->numberofpolygons > 0) {
      f->polygonlist = new polygon[f->numberofpolygons];
    }
    // Go through all segments, read in their vertices.
    for (j = 0; j < f->numberofpolygons; j++) {
      p = &(f->polygonlist[j]);
      init(p);
      // Read in a segment.
      stringptr = readnumberline(inputline, infile, infilename);
      stringptr = findnextnumber(stringptr); // Skip its index.
      p->numberofvertices = 2; // A segment always has two vertices.
      p->vertexlist = new int[p->numberofvertices];
      p->vertexlist[0] = (int) strtol (stringptr, &stringptr, 0);
      stringptr = findnextnumber(stringptr);
      p->vertexlist[1] = (int) strtol (stringptr, &stringptr, 0);
    }
    // Read number of holes.
    stringptr = readnumberline(inputline, infile, infilename);
    f->numberofholes = (int) strtol (stringptr, &stringptr, 0);
    if (f->numberofholes > 0) {
      // Initialize 'f->holelist'.
      f->holelist = new REAL[f->numberofholes * 3];
      // Read the holes' coordinates.
      for (j = 0; j < f->numberofholes; j++) {
        // Read a 2D hole point.
        stringptr = readnumberline(inputline, infile, infilename);
        stringptr = findnextnumber(stringptr); // Skip its index.
        f->holelist[j * 3] = (REAL) strtod (stringptr, &stringptr);
        stringptr = findnextnumber(stringptr);
        f->holelist[j * 3 + 1] = (REAL) strtod (stringptr, &stringptr);
        f->holelist[j * 3 + 2] = 0.0; // The z-coord.
      }
    }
    // The regions are skipped.

  }

  // End of reading poly/smesh file.
  fclose(infile);
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_off()    Load a polyhedron from a .off file.                         //
//                                                                           //
// The .off format is one of file formats of the Geomview, an interactive    //
// program for viewing and manipulating geometric objects.  More information //
// is available form: http://www.geomview.org.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_off(char* filebasename)
{
  FILE *fp;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  char infilename[FILENAMESIZE];
  char buffer[INPUTLINESIZE];
  char *bufferp;
  double *coord;
  int nverts = 0, iverts = 0;
  int nfaces = 0, ifaces = 0;
  int nedges = 0;
  int line_count = 0, i;

  // Default, the off file's index is from '0'. We check it by remembering the
  //   smallest index we found in the file. It should be either 0 or 1.
  int smallestidx = 0; 

  strncpy(infilename, filebasename, 1024 - 1);
  infilename[FILENAMESIZE - 1] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".off") != 0) {
    strcat(infilename, ".off");
  }

  if (!(fp = fopen(infilename, "r"))) {
    printf("  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    // Check section
    if (nverts == 0) {
      // Read header 
      bufferp = strstr(bufferp, "OFF");
      if (bufferp != NULL) {
        // Read mesh counts
        bufferp = findnextnumber(bufferp); // Skip field "OFF".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) 
            || (nverts == 0)) {
          printf("Syntax error reading header on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        // Allocate memory for 'tetgenio'
        if (nverts > 0) {
          numberofpoints = nverts;
          pointlist = new REAL[nverts * 3];
          smallestidx = nverts + 1; // A bigger enough number.
        }
        if (nfaces > 0) {        
          numberoffacets = nfaces;
          facetlist = new tetgenio::facet[nfaces];
        }
      }
    } else if (iverts < nverts) {
      // Read vertex coordinates
      coord = &pointlist[iverts * 3];
      for (i = 0; i < 3; i++) {
        if (*bufferp == '\0') {
          printf("Syntax error reading vertex coords on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        coord[i] = (REAL) strtod(bufferp, &bufferp);
        bufferp = findnextnumber(bufferp);
      }
      iverts++;
    } else if (ifaces < nfaces) {
      // Get next face
      f = &facetlist[ifaces];
      init(f);      
      // In .off format, each facet has one polygon, no hole.
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      p = &f->polygonlist[0];
      init(p);
      // Read the number of vertices, it should be greater than 0.
      p->numberofvertices = (int) strtol(bufferp, &bufferp, 0);
      if (p->numberofvertices == 0) {
        printf("Syntax error reading polygon on line %d in file %s\n",
               line_count, infilename);
        fclose(fp);
        return false;
      }
      // Allocate memory for face vertices
      p->vertexlist = new int[p->numberofvertices];
      for (i = 0; i < p->numberofvertices; i++) {
        bufferp = findnextnumber(bufferp);
        if (*bufferp == '\0') {
          printf("Syntax error reading polygon on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        p->vertexlist[i] = (int) strtol(bufferp, &bufferp, 0);
        // Detect the smallest index.
        if (p->vertexlist[i] < smallestidx) {
          smallestidx = p->vertexlist[i];
        }
      }
      ifaces++;
    } else {
      // Should never get here
      printf("Found extra text starting at line %d in file %s\n", line_count,
             infilename);
      break;
    }
  }

  // Close file
  fclose(fp);

  // Decide the firstnumber of the index.
  if (smallestidx == 0) {
    firstnumber = 0;  
  } else if (smallestidx == 1) {
    firstnumber = 1;
  } else {
    printf("A wrong smallest index (%d) was detected in file %s\n",
           smallestidx, infilename);
    return false;
  }

  if (iverts != nverts) {
    printf("Expected %d vertices, but read only %d vertices in file %s\n",
           nverts, iverts, infilename);
    return false;
  }
  if (ifaces != nfaces) {
    printf("Expected %d faces, but read only %d faces in file %s\n",
           nfaces, ifaces, infilename);
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_ply()    Load a polyhedron from a .ply file.                         //
//                                                                           //
// This is a simplified version of reading .ply files, which only reads the  //
// set of vertices and the set of faces. Other informations (such as color,  //
// material, texture, etc) in .ply file are ignored. Complete routines for   //
// reading and writing ,ply files are available from: http://www.cc.gatech.  //
// edu/projects/large_models/ply.html.  Except the header section, ply file  //
// format has exactly the same format for listing vertices and polygons as   //
// off file format.                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_ply(char* filebasename)
{
  FILE *fp;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  char infilename[FILENAMESIZE];
  char buffer[INPUTLINESIZE];
  char *bufferp, *str;
  double *coord;
  int endheader = 0, format = 0;
  int nverts = 0, iverts = 0;
  int nfaces = 0, ifaces = 0;
  int line_count = 0, i;

  // Default, the ply file's index is from '0'. We check it by remembering the
  //   smallest index we found in the file. It should be either 0 or 1.
  int smallestidx = 0; 

  strncpy(infilename, filebasename, FILENAMESIZE - 1);
  infilename[FILENAMESIZE - 1] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".ply") != 0) {
    strcat(infilename, ".ply");
  }

  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    if (!endheader) {
      // Find if it is the keyword "end_header".
      str = strstr(bufferp, "end_header");
      // strstr() is case sensitive.
      if (!str) str = strstr(bufferp, "End_header");
      if (!str) str = strstr(bufferp, "End_Header");
      if (str) {
        // This is the end of the header section.
        endheader = 1; 
        continue;
      }
      // Parse the number of vertices and the number of faces.
      if (nverts == 0 || nfaces == 0) {
        // Find if it si the keyword "element".
        str = strstr(bufferp, "element");
        if (!str) str = strstr(bufferp, "Element");
        if (str) {
          bufferp = findnextfield(str);
          if (*bufferp == '\0') {
            printf("Syntax error reading element type on line%d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          if (nverts == 0) {
            // Find if it is the keyword "vertex".
            str = strstr(bufferp, "vertex");
            if (!str) str = strstr(bufferp, "Vertex");
            if (str) {
              bufferp = findnextnumber(str);
              if (*bufferp == '\0') {
                printf("Syntax error reading vertex number on line");
                printf(" %d in file %s\n", line_count, infilename);
                fclose(fp);
                return false;
              }
              nverts = (int) strtol(bufferp, &bufferp, 0);
              // Allocate memory for 'tetgenio'
              if (nverts > 0) {
                numberofpoints = nverts;
                pointlist = new REAL[nverts * 3];
                smallestidx = nverts + 1; // A big enough index.
              }
            }
          }
          if (nfaces == 0) {
            // Find if it is the keyword "face".
            str = strstr(bufferp, "face");
            if (!str) str = strstr(bufferp, "Face");
            if (str) {
              bufferp = findnextnumber(str);
              if (*bufferp == '\0') {
                printf("Syntax error reading face number on line");
                printf(" %d in file %s\n", line_count, infilename);
                fclose(fp);
                return false;
              }
              nfaces = (int) strtol(bufferp, &bufferp, 0);
              // Allocate memory for 'tetgenio'
              if (nfaces > 0) {        
                numberoffacets = nfaces;
                facetlist = new tetgenio::facet[nfaces];
              }
            }
          }
        } // It is not the string "element". 
      }
      if (format == 0) {
        // Find the keyword "format".
        str = strstr(bufferp, "format");
        if (!str) str = strstr(bufferp, "Format");
        if (str) {
          format = 1;
          bufferp = findnextfield(str);
          // Find if it is the string "ascii".
          str = strstr(bufferp, "ascii");
          if (!str) str = strstr(bufferp, "ASCII");
          if (!str) {
            printf("This routine only reads ascii format of ply files.\n");
            printf("Hint: You can convert the binary to ascii format by\n");
            printf("  using the provided ply tools:\n");
            printf("  ply2ascii < %s > ascii_%s\n", infilename, infilename);
            fclose(fp);
            return false;
          }
        }
      }
    } else if (iverts < nverts) {
      // Read vertex coordinates
      coord = &pointlist[iverts * 3];
      for (i = 0; i < 3; i++) {
        if (*bufferp == '\0') {
          printf("Syntax error reading vertex coords on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        coord[i] = (REAL) strtod(bufferp, &bufferp);
        bufferp = findnextnumber(bufferp);
      }
      iverts++;
    } else if (ifaces < nfaces) {
      // Get next face
      f = &facetlist[ifaces];
      init(f);      
      // In .off format, each facet has one polygon, no hole.
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      p = &f->polygonlist[0];
      init(p);
      // Read the number of vertices, it should be greater than 0.
      p->numberofvertices = (int) strtol(bufferp, &bufferp, 0);
      if (p->numberofvertices == 0) {
        printf("Syntax error reading polygon on line %d in file %s\n",
               line_count, infilename);
        fclose(fp);
        return false;
      }
      // Allocate memory for face vertices
      p->vertexlist = new int[p->numberofvertices];
      for (i = 0; i < p->numberofvertices; i++) {
        bufferp = findnextnumber(bufferp);
        if (*bufferp == '\0') {
          printf("Syntax error reading polygon on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        p->vertexlist[i] = (int) strtol(bufferp, &bufferp, 0);
        if (p->vertexlist[i] < smallestidx) {
          smallestidx = p->vertexlist[i];
        }
      }
      ifaces++;
    } else {
      // Should never get here
      printf("Found extra text starting at line %d in file %s\n", line_count,
             infilename);
      break;
    }
  }

  // Close file
  fclose(fp);

  // Decide the firstnumber of the index.
  if (smallestidx == 0) {
    firstnumber = 0;  
  } else if (smallestidx == 1) {
    firstnumber = 1;
  } else {
    printf("A wrong smallest index (%d) was detected in file %s\n",
           smallestidx, infilename);
    return false;
  }

  if (iverts != nverts) {
    printf("Expected %d vertices, but read only %d vertices in file %s\n",
           nverts, iverts, infilename);
    return false;
  }
  if (ifaces != nfaces) {
    printf("Expected %d faces, but read only %d faces in file %s\n",
           nfaces, ifaces, infilename);
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_stl()    Load a surface mesh from a .stl file.                       //
//                                                                           //
// The .stl or stereolithography format is an ASCII or binary file used in   //
// manufacturing.  It is a list of the triangular surfaces that describe a   //
// computer generated solid model. This is the standard input for most rapid //
// prototyping machines.                                                     //
//                                                                           //
// Comment: A .stl file many contain many duplicated points.  They will be   //
// unified during the Delaunay tetrahedralization process.                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_stl(char* filebasename)
{
  FILE *fp;
  tetgenmesh::arraypool *plist;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  char infilename[FILENAMESIZE];
  char buffer[INPUTLINESIZE];
  char *bufferp, *str;
  double *coord;
  int solid = 0;
  int nverts = 0, iverts = 0;
  int nfaces = 0;
  int line_count = 0, i;

  strncpy(infilename, filebasename, FILENAMESIZE - 1);
  infilename[FILENAMESIZE - 1] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".stl") != 0) {
    strcat(infilename, ".stl");
  }

  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  // STL file has no number of points available. Use a list to read points.
  plist = new tetgenmesh::arraypool(sizeof(double) * 3, 10); 

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    // The ASCII .stl file must start with the lower case keyword solid and
    //   end with endsolid.
    if (solid == 0) {
      // Read header 
      bufferp = strstr(bufferp, "solid");
      if (bufferp != NULL) {
        solid = 1;
      }
    } else {
      // We're inside the block of the solid.
      str = bufferp;
      // Is this the end of the solid.
      bufferp = strstr(bufferp, "endsolid");
      if (bufferp != NULL) {
        solid = 0;
      } else {
        // Read the XYZ coordinates if it is a vertex.
        bufferp = str;
        bufferp = strstr(bufferp, "vertex");
        if (bufferp != NULL) {
          plist->newindex((void **) &coord);
          for (i = 0; i < 3; i++) {
            bufferp = findnextnumber(bufferp);
            if (*bufferp == '\0') {
              printf("Syntax error reading vertex coords on line %d\n",
                   line_count);
              delete plist;
              fclose(fp);
              return false;
            }
            coord[i] = (REAL) strtod(bufferp, &bufferp);
          }
        }
      }
    }
  }
  fclose(fp);

  nverts = (int) plist->objects;
  // nverts should be an integer times 3 (every 3 vertices denote a face).
  if (nverts == 0 || (nverts % 3 != 0)) {
    printf("Error:  Wrong number of vertices in file %s.\n", infilename);
    delete plist;
    return false;
  }
  numberofpoints = nverts;
  pointlist = new REAL[nverts * 3];
  for (i = 0; i < nverts; i++) {
    coord = (double *) fastlookup(plist, i);
    iverts = i * 3;
    pointlist[iverts] = (REAL) coord[0];
    pointlist[iverts + 1] = (REAL) coord[1];
    pointlist[iverts + 2] = (REAL) coord[2];
  }

  nfaces = (int) (nverts / 3);
  numberoffacets = nfaces;
  facetlist = new tetgenio::facet[nfaces];

  // Default use '1' as the array starting index.
  firstnumber = 1;
  iverts = firstnumber;
  for (i = 0; i < nfaces; i++) {
    f = &facetlist[i];
    init(f);      
    // In .stl format, each facet has one polygon, no hole.
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[1];
    p = &f->polygonlist[0];
    init(p);
    // Each polygon has three vertices.
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = iverts;
    p->vertexlist[1] = iverts + 1;
    p->vertexlist[2] = iverts + 2;
    iverts += 3;
  }

  delete plist;
  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_medit()    Load a surface mesh from a .mesh file.                    //
//                                                                           //
// The .mesh format is the file format of Medit, a user-friendly interactive //
// mesh viewer program.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_medit(char* filebasename, int istetmesh)
{
  FILE *fp;
  tetgenio::facet *tmpflist, *f;
  tetgenio::polygon *p;
  char infilename[FILENAMESIZE];
  char buffer[INPUTLINESIZE];
  char *bufferp, *str;
  double *coord;
  int *tmpfmlist;
  int dimension = 0;
  int nverts = 0;
  int nfaces = 0;
  int ntets = 0;
  int line_count = 0;
  int corners = 0; // 3 (triangle) or 4 (quad).
  int *plist;
  int i, j;

  int smallestidx = 0;

  strncpy(infilename, filebasename, FILENAMESIZE - 1);
  infilename[FILENAMESIZE - 1] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 5], ".mesh") != 0) {
    strcat(infilename, ".mesh");
  }
  
  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  while ((bufferp = readline(buffer, fp, &line_count)) != NULL) {
    if (*bufferp == '#') continue;  // A comment line is skipped.
    if (dimension == 0) {
      // Find if it is the keyword "Dimension".
      str = strstr(bufferp, "Dimension");
      if (!str) str = strstr(bufferp, "dimension");
      if (!str) str = strstr(bufferp, "DIMENSION");
      if (str) {
        // Read the dimensions
        bufferp = findnextnumber(str); // Skip field "Dimension".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        dimension = (int) strtol(bufferp, &bufferp, 0);
        if (dimension != 2 && dimension != 3) {
          printf("Unknown dimension in file on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        mesh_dim = dimension;
      }
    }
    if (nverts == 0) {
      // Find if it is the keyword "Vertices".
      str = strstr(bufferp, "Vertices");
      if (!str) str = strstr(bufferp, "vertices");
      if (!str) str = strstr(bufferp, "VERTICES");
      if (str) {
        // Read the number of vertices.
        bufferp = findnextnumber(str); // Skip field "Vertices".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        nverts = (int) strtol(bufferp, &bufferp, 0);
        // Initialize the smallest index.
        smallestidx = nverts + 1;
        // Allocate memory for 'tetgenio'
        if (nverts > 0) {
          numberofpoints = nverts;
          pointlist = new REAL[nverts * 3];
        }
        // Read the follwoing node list.
        for (i = 0; i < nverts; i++) {
          bufferp = readline(buffer, fp, &line_count);
          if (bufferp == NULL) {
            printf("Unexpected end of file on line %d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          // Read vertex coordinates
          coord = &pointlist[i * 3];
          for (j = 0; j < 3; j++) {
            if (*bufferp == '\0') {
              printf("Syntax error reading vertex coords on line");
              printf(" %d in file %s\n", line_count, infilename);
              fclose(fp);
              return false;
            }
            if ((j < 2) || (dimension == 3)) {
              coord[j] = (REAL) strtod(bufferp, &bufferp);
            } else {
              assert((j == 2) && (dimension == 2));
              coord[j] = 0.0;
            }
            bufferp = findnextnumber(bufferp);
          }
        }
        continue;
      }
    }
    if (ntets == 0) {
      // Find if it is the keyword "Tetrahedra"
      corners = 0;
      str = strstr(bufferp, "Tetrahedra");
      if (!str) str = strstr(bufferp, "tetrahedra");
      if (!str) str = strstr(bufferp, "TETRAHEDRA");
      if (str) {
        corners = 4;
      }
      if (corners == 4) {
        // Read the number of tetrahedra
        bufferp = findnextnumber(str); // Skip field "Tetrahedra".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        ntets = strtol(bufferp, &bufferp, 0);
        if (ntets > 0) {
          // It is a tetrahedral mesh.
          numberoftetrahedra = ntets;
          numberofcorners = 4;
          numberoftetrahedronattributes = 1;
          tetrahedronlist = new int[ntets * 4];
          tetrahedronattributelist = new REAL[ntets];
        }
      } // if (corners == 4)
      // Read the list of tetrahedra.
      for (i = 0; i < numberoftetrahedra; i++) {
        plist = &(tetrahedronlist[i * 4]);
        bufferp = readline(buffer, fp, &line_count);
        if (bufferp == NULL) {
          printf("Unexpected end of file on line %d in file %s\n",
                 line_count, infilename);
          fclose(fp);
          return false;
        }
        // Read the vertices of the tet.
        for (j = 0; j < corners; j++) {
          if (*bufferp == '\0') {
            printf("Syntax error reading face on line %d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          plist[j] = (int) strtol(bufferp, &bufferp, 0);
          // Remember the smallest index.
          if (plist[j] < smallestidx) smallestidx = plist[j];
          bufferp = findnextnumber(bufferp);
        }
        // Read the attribute of the tet if it exists.
        tetrahedronattributelist[i] = 0;
        if (*bufferp != '\0') {
          tetrahedronattributelist[i] = (REAL) strtol(bufferp, &bufferp, 0);
        }
      } // i
    } // Tetrahedra
    if (nfaces == 0) {
      // Find if it is the keyword "Triangles" or "Quadrilaterals".
      corners = 0;
      str = strstr(bufferp, "Triangles");
      if (!str) str = strstr(bufferp, "triangles");
      if (!str) str = strstr(bufferp, "TRIANGLES");
      if (str) {
        corners = 3;
      } else {
        str = strstr(bufferp, "Quadrilaterals");
        if (!str) str = strstr(bufferp, "quadrilaterals");
        if (!str) str = strstr(bufferp, "QUADRILATERALS");
        if (str) {
          corners = 4;
        }
      }
      if (corners == 3 || corners == 4) {
        // Read the number of triangles (or quadrilaterals).
        bufferp = findnextnumber(str); // Skip field "Triangles".
        if (*bufferp == '\0') {
          // Read a non-empty line.
          bufferp = readline(buffer, fp, &line_count);
        }
        nfaces = strtol(bufferp, &bufferp, 0);
        // Allocate memory for 'tetgenio'
        if (nfaces > 0) {
          if (!istetmesh) {
            // It is a PLC surface mesh.
            if (numberoffacets > 0) {
              // facetlist has already been allocated. Enlarge arrays.
              // This happens when the surface mesh contains mixed cells.
              tmpflist = new tetgenio::facet[numberoffacets + nfaces];
              tmpfmlist = new int[numberoffacets + nfaces];
              // Copy the data of old arrays into new arrays.
              for (i = 0; i < numberoffacets; i++) {
                f = &(tmpflist[i]);
                tetgenio::init(f);
                *f = facetlist[i];
                tmpfmlist[i] = facetmarkerlist[i];
              }
              // Release old arrays.
              delete [] facetlist;
              delete [] facetmarkerlist;
              // Remember the new arrays.
              facetlist = tmpflist;
              facetmarkerlist = tmpfmlist;
            } else {
              // This is the first time to allocate facetlist.
              facetlist = new tetgenio::facet[nfaces];
              facetmarkerlist = new int[nfaces];
            }
          } else {
            if (corners == 3) {
              // It is a surface mesh of a tetrahedral mesh.
              numberoftrifaces = nfaces;
              trifacelist = new int[nfaces * 3];
              trifacemarkerlist = new int[nfaces];
            }
          }
        } // if (nfaces > 0)
        // Read the following list of faces.
        if (!istetmesh) {
          for (i = numberoffacets; i < numberoffacets + nfaces; i++) {
            bufferp = readline(buffer, fp, &line_count);
            if (bufferp == NULL) {
              printf("Unexpected end of file on line %d in file %s\n",
                     line_count, infilename);
              fclose(fp);
              return false;
            }
            f = &facetlist[i];
            tetgenio::init(f);
            // In .mesh format, each facet has one polygon, no hole.
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[1];
            p = &f->polygonlist[0];
            tetgenio::init(p);
            p->numberofvertices = corners;
            // Allocate memory for face vertices
            p->vertexlist = new int[p->numberofvertices];
            // Read the vertices of the face.
            for (j = 0; j < corners; j++) {
              if (*bufferp == '\0') {
                printf("Syntax error reading face on line %d in file %s\n",
                       line_count, infilename);
                fclose(fp);
                return false;
              }
              p->vertexlist[j] = (int) strtol(bufferp, &bufferp, 0);
              // Remember the smallest index.
              if (p->vertexlist[j] < smallestidx) {
                smallestidx = p->vertexlist[j];
              }
              bufferp = findnextnumber(bufferp);
            }
            // Read the marker of the face if it exists.
            facetmarkerlist[i] = 0;
            if (*bufferp != '\0') {
              facetmarkerlist[i] = (int) strtol(bufferp, &bufferp, 0);
            }
          }
          // Have read in a list of triangles/quads.
          numberoffacets += nfaces;
          nfaces = 0;
        } else {
          // It is a surface mesh of a tetrahedral mesh.
          if (corners == 3) {
            for (i = 0; i < numberoftrifaces; i++) {
              plist = &(trifacelist[i * 3]);
              bufferp = readline(buffer, fp, &line_count);
              if (bufferp == NULL) {
                printf("Unexpected end of file on line %d in file %s\n",
                       line_count, infilename);
                fclose(fp);
                return false;
              }
              // Read the vertices of the face.
              for (j = 0; j < corners; j++) {
                if (*bufferp == '\0') {
                  printf("Syntax error reading face on line %d in file %s\n",
                         line_count, infilename);
                  fclose(fp);
                  return false;
                }
                plist[j] = (int) strtol(bufferp, &bufferp, 0);
                // Remember the smallest index.
                if (plist[j] < smallestidx) {
                  smallestidx = plist[j];
                }
                bufferp = findnextnumber(bufferp);
              }
              // Read the marker of the face if it exists.
              trifacemarkerlist[i] = 0;
              if (*bufferp != '\0') {
                trifacemarkerlist[i] = (int) strtol(bufferp, &bufferp, 0);
              }
            } // i
          } // if (corners == 3)
        } // if (b->refine)
      } // if (corners == 3 || corners == 4)
    }
  }

  // Close file
  fclose(fp);

  // Decide the firstnumber of the index.
  if (smallestidx == 0) {
    firstnumber = 0;  
  } else if (smallestidx == 1) {
    firstnumber = 1;
  } else {
    printf("A wrong smallest index (%d) was detected in file %s\n",
           smallestidx, infilename);
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_vtk()    Load VTK surface mesh from file (.vtk ascii or binary).     //
//                                                                           //
// This function is contributed by: Bryn Lloyd, Computer Vision Laboratory,  //
// ETH, Zuerich. May 7, 2007.                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Two inline functions used in read/write VTK files.

void swapBytes(unsigned char* var, int size)
{
  int i = 0;
  int j = size - 1;
  char c;

  while (i < j) {
    c = var[i]; var[i] = var[j]; var[j] = c;
    i++, j--;
  }
}

bool testIsBigEndian()
{
  short word = 0x4321;
  if((*(char *)& word) != 0x21)
    return true;
  else 
    return false;
}


bool tetgenio::load_vtk(char* filebasename)
{
  FILE *fp;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  char infilename[FILENAMESIZE];
  char line[INPUTLINESIZE];
  char mode[128], id[256], fmt[64];
  char *bufferp;
  double *coord;
  float _x, _y, _z;
  int nverts = 0;
  int nfaces = 0;
  int line_count = 0;
  int dummy;
  int id1, id2, id3;
  int nn = -1;
  int nn_old = -1;
  int i, j;
  bool ImALittleEndian = !testIsBigEndian();

  int smallestidx = 0;

  strncpy(infilename, filebasename, FILENAMESIZE - 1);
  infilename[FILENAMESIZE - 1] = '\0';
  if (infilename[0] == '\0') {
    printf("Error:  No filename.\n");
    return false;
  }
  if (strcmp(&infilename[strlen(infilename) - 4], ".vtk") != 0) {
    strcat(infilename, ".vtk");
  }
  if (!(fp = fopen(infilename, "r"))) {
    printf("Error:  Unable to open file %s\n", infilename);
    return false;
  }
  printf("Opening %s.\n", infilename);

  // Default uses the index starts from '0'.
  firstnumber = 0;
  strcpy(mode, "BINARY");

  while((bufferp = readline(line, fp, &line_count)) != NULL) {
    if(strlen(line) == 0) continue;
    //swallow lines beginning with a comment sign or white space
    if(line[0] == '#' || line[0]=='\n' || line[0] == 10 || line[0] == 13 || 
       line[0] == 32) continue;

    sscanf(line, "%s", id);
    if(!strcmp(id, "ASCII")) {
      strcpy(mode, "ASCII");
    }

    if(!strcmp(id, "POINTS")) {
      sscanf(line, "%s %d %s", id, &nverts, fmt);
      if (nverts > 0) {
        numberofpoints = nverts;
        pointlist = new REAL[nverts * 3];
        smallestidx = nverts + 1;
      }

      if(!strcmp(mode, "BINARY")) {
        for(i = 0; i < nverts; i++) {
          coord = &pointlist[i * 3];
          if(!strcmp(fmt, "double")) {
            fread((char*)(&(coord[0])), sizeof(double), 1, fp);
            fread((char*)(&(coord[1])), sizeof(double), 1, fp);
            fread((char*)(&(coord[2])), sizeof(double), 1, fp);
            if(ImALittleEndian){
              swapBytes((unsigned char *) &(coord[0]), sizeof(coord[0]));
              swapBytes((unsigned char *) &(coord[1]), sizeof(coord[1]));
              swapBytes((unsigned char *) &(coord[2]), sizeof(coord[2]));
            }
          } else if(!strcmp(fmt, "float")) {
            fread((char*)(&_x), sizeof(float), 1, fp);
            fread((char*)(&_y), sizeof(float), 1, fp);
            fread((char*)(&_z), sizeof(float), 1, fp);
            if(ImALittleEndian){
              swapBytes((unsigned char *) &_x, sizeof(_x));
              swapBytes((unsigned char *) &_y, sizeof(_y));
              swapBytes((unsigned char *) &_z, sizeof(_z));
            }
            coord[0] = double(_x);
            coord[1] = double(_y);
            coord[2] = double(_z);
          } else {
            printf("Error: Only float or double formats are supported!\n");
            return false;
          }
        }
      } else if(!strcmp(mode, "ASCII")) {
        for(i = 0; i < nverts; i++){
          bufferp = readline(line, fp, &line_count);
          if (bufferp == NULL) {
            printf("Unexpected end of file on line %d in file %s\n",
                   line_count, infilename);
            fclose(fp);
            return false;
          }
          // Read vertex coordinates
          coord = &pointlist[i * 3];
          for (j = 0; j < 3; j++) {
            if (*bufferp == '\0') {
              printf("Syntax error reading vertex coords on line");
              printf(" %d in file %s\n", line_count, infilename);
              fclose(fp);
              return false;
            }
            coord[j] = (REAL) strtod(bufferp, &bufferp);
            bufferp = findnextnumber(bufferp);
          }
        }
      }
      continue;
    }

    if(!strcmp(id, "POLYGONS")) {
      sscanf(line, "%s %d  %d", id, &nfaces, &dummy);
      if (nfaces > 0) {
        numberoffacets = nfaces;
        facetlist = new tetgenio::facet[nfaces];
      }

      if(!strcmp(mode, "BINARY")) {
        for(i = 0; i < nfaces; i++){
          fread((char*)(&nn), sizeof(int), 1, fp);
          if(ImALittleEndian){
            swapBytes((unsigned char *) &nn, sizeof(nn));
          }
          if (i == 0)
            nn_old = nn;
          if (nn != nn_old) {
            printf("Error:  No mixed cells are allowed.\n");
            return false;
          }

          if(nn == 3){
            fread((char*)(&id1), sizeof(int), 1, fp);
            fread((char*)(&id2), sizeof(int), 1, fp);
            fread((char*)(&id3), sizeof(int), 1, fp);
            if(ImALittleEndian){
              swapBytes((unsigned char *) &id1, sizeof(id1));
              swapBytes((unsigned char *) &id2, sizeof(id2));
              swapBytes((unsigned char *) &id3, sizeof(id3));
            }
            f = &facetlist[i];
            init(f);
            // In .off format, each facet has one polygon, no hole.
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[1];
            p = &f->polygonlist[0];
            init(p);
            // Set number of vertices
            p->numberofvertices = 3;
            // Allocate memory for face vertices
            p->vertexlist = new int[p->numberofvertices];
            p->vertexlist[0] = id1;
            p->vertexlist[1] = id2;
            p->vertexlist[2] = id3;
            // Detect the smallest index.
            for (j = 0; j < 3; j++) {
              if (p->vertexlist[j] < smallestidx) {
                smallestidx = p->vertexlist[j];
              }
            }
          } else {
            printf("Error: Only triangles are supported\n");
            return false;
          }
        }
      } else if(!strcmp(mode, "ASCII")) {
        for(i = 0; i < nfaces; i++) {
          bufferp = readline(line, fp, &line_count);
          nn = (int) strtol(bufferp, &bufferp, 0);
          if (i == 0)
            nn_old = nn;
          if (nn != nn_old) {
            printf("Error:  No mixed cells are allowed.\n");
            return false;
          }

          if (nn == 3) {
            bufferp = findnextnumber(bufferp); // Skip the first field.
            id1 = (int) strtol(bufferp, &bufferp, 0);
            bufferp = findnextnumber(bufferp);
            id2 = (int) strtol(bufferp, &bufferp, 0);
            bufferp = findnextnumber(bufferp);
            id3 = (int) strtol(bufferp, &bufferp, 0);
            f = &facetlist[i];
            init(f);
            // In .off format, each facet has one polygon, no hole.
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[1];
            p = &f->polygonlist[0];
            init(p);
            // Set number of vertices
            p->numberofvertices = 3;
            // Allocate memory for face vertices
            p->vertexlist = new int[p->numberofvertices];
            p->vertexlist[0] = id1;
            p->vertexlist[1] = id2;
            p->vertexlist[2] = id3;
            // Detect the smallest index.
            for (j = 0; j < 3; j++) {
              if (p->vertexlist[j] < smallestidx) {
                smallestidx = p->vertexlist[j];
              }
            }
          } else {
            printf("Error:  Only triangles are supported.\n");
            return false;
          }
        }
      }

      fclose(fp);

      // Decide the firstnumber of the index.
      if (smallestidx == 0) {
        firstnumber = 0;  
      } else if (smallestidx == 1) {
        firstnumber = 1;
      } else {
        printf("A wrong smallest index (%d) was detected in file %s\n",
               smallestidx, infilename);
        return false;
      }

      return true;
    }

    if(!strcmp(id,"LINES") || !strcmp(id,"CELLS")){
      printf("Warning:  load_vtk(): cannot read formats LINES, CELLS.\n");
    }
  } // while ()

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_plc()    Load a piecewise linear complex from file(s).               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_plc(char* filebasename, int object)
{
  bool success;

  if (object == (int) tetgenbehavior::NODES) {
    success = load_node(filebasename);
  } else if (object == (int) tetgenbehavior::POLY) {
    success = load_poly(filebasename);
  } else if (object == (int) tetgenbehavior::OFF) {
    success = load_off(filebasename);
  } else if (object == (int) tetgenbehavior::PLY) {
    success = load_ply(filebasename);
  } else if (object == (int) tetgenbehavior::STL) {
    success = load_stl(filebasename);
  } else if (object == (int) tetgenbehavior::MEDIT) {
    success = load_medit(filebasename, 0);
  } else if (object == (int) tetgenbehavior::VTK) {
    success = load_vtk(filebasename);
  } else {
    success = load_poly(filebasename);
  }

  if (success) {
    // Try to load the following files (.edge, .var, .mtr).
    load_edge(filebasename);
    load_var(filebasename);
    load_mtr(filebasename);
  }

  return success;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// load_mesh()    Load a tetrahedral mesh from file(s).                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool tetgenio::load_tetmesh(char* filebasename, int object)
{
  bool success;

  if (object == (int) tetgenbehavior::MEDIT) {
    success = load_medit(filebasename, 1);
  } else {
    success = load_node(filebasename);
    if (success) {
      success = load_tet(filebasename);
    }
    if (success) {
      // Try to load the following files (.face, .edge, .vol).
      load_face(filebasename);
      load_edge(filebasename);
      load_vol(filebasename);
    }
  }

  if (success) {
    // Try to load the following files (.var, .mtr).
    load_var(filebasename);
    load_mtr(filebasename);
  }

  return success;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_nodes()    Save points to a .node file.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_nodes(char* filebasename)
{
  FILE *fout;
  char outnodefilename[FILENAMESIZE];
  char outmtrfilename[FILENAMESIZE];
  int i, j;

  sprintf(outnodefilename, "%s.node", filebasename);
  printf("Saving nodes to %s\n", outnodefilename);
  fout = fopen(outnodefilename, "w");
  fprintf(fout, "%d  %d  %d  %d\n", numberofpoints, mesh_dim,
          numberofpointattributes, pointmarkerlist != NULL ? 1 : 0);
  for (i = 0; i < numberofpoints; i++) {
    if (mesh_dim == 2) {
      fprintf(fout, "%d  %.16g  %.16g", i + firstnumber, pointlist[i * 3],
              pointlist[i * 3 + 1]);
    } else {
      fprintf(fout, "%d  %.16g  %.16g  %.16g", i + firstnumber,
              pointlist[i * 3], pointlist[i * 3 + 1], pointlist[i * 3 + 2]);
    }
    for (j = 0; j < numberofpointattributes; j++) {
      fprintf(fout, "  %.16g", 
              pointattributelist[i * numberofpointattributes + j]);
    }
    if (pointmarkerlist != NULL) {
      fprintf(fout, "  %d", pointmarkerlist[i]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);

  // If the point metrics exist, output them to a .mtr file.
  if ((numberofpointmtrs > 0) && (pointmtrlist != (REAL *) NULL)) {
    sprintf(outmtrfilename, "%s.mtr", filebasename);
    printf("Saving metrics to %s\n", outmtrfilename);
    fout = fopen(outmtrfilename, "w");
    fprintf(fout, "%d  %d\n", numberofpoints, numberofpointmtrs);
    for (i = 0; i < numberofpoints; i++) {
      for (j = 0; j < numberofpointmtrs; j++) {
        fprintf(fout, "%.16g ", pointmtrlist[i * numberofpointmtrs + j]);
      }
      fprintf(fout, "\n");
    }
    fclose(fout);
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_elements()    Save elements to a .ele file.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_elements(char* filebasename)
{
  FILE *fout;
  char outelefilename[FILENAMESIZE];
  int i, j;

  sprintf(outelefilename, "%s.ele", filebasename);
  printf("Saving elements to %s\n", outelefilename);
  fout = fopen(outelefilename, "w");
  if (mesh_dim == 3) {
    fprintf(fout, "%d  %d  %d\n", numberoftetrahedra, numberofcorners,
            numberoftetrahedronattributes);
    for (i = 0; i < numberoftetrahedra; i++) {
      fprintf(fout, "%d", i + firstnumber);
      for (j = 0; j < numberofcorners; j++) {
        fprintf(fout, "  %5d", tetrahedronlist[i * numberofcorners + j]);
      }
      for (j = 0; j < numberoftetrahedronattributes; j++) {
        fprintf(fout, "  %g",
          tetrahedronattributelist[i * numberoftetrahedronattributes + j]);
      }
      fprintf(fout, "\n");
    }
  } else {
    // Save a two-dimensional mesh.
    fprintf(fout, "%d  %d  %d\n",numberoftrifaces,3,trifacemarkerlist ? 1 : 0);
    for (i = 0; i < numberoftrifaces; i++) {
      fprintf(fout, "%d", i + firstnumber);
      for (j = 0; j < 3; j++) {
        fprintf(fout, "  %5d", trifacelist[i * 3 + j]);
      }
      if (trifacemarkerlist != NULL) {
        fprintf(fout, "  %d", trifacemarkerlist[i]);
      }
      fprintf(fout, "\n");
    }
  }

  fclose(fout);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_faces()    Save faces to a .face file.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_faces(char* filebasename)
{
  FILE *fout;
  char outfacefilename[FILENAMESIZE];
  int i;

  sprintf(outfacefilename, "%s.face", filebasename);
  printf("Saving faces to %s\n", outfacefilename);
  fout = fopen(outfacefilename, "w");
  fprintf(fout, "%d  %d\n", numberoftrifaces, 
          trifacemarkerlist != NULL ? 1 : 0);
  for (i = 0; i < numberoftrifaces; i++) {
    fprintf(fout, "%d  %5d  %5d  %5d", i + firstnumber, trifacelist[i * 3],
            trifacelist[i * 3 + 1], trifacelist[i * 3 + 2]);
    if (trifacemarkerlist != NULL) {
      fprintf(fout, "  %d", trifacemarkerlist[i]);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_edges()    Save egdes to a .edge file.                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_edges(char* filebasename)
{
  FILE *fout;
  char outedgefilename[FILENAMESIZE];
  int i;

  sprintf(outedgefilename, "%s.edge", filebasename);
  printf("Saving edges to %s\n", outedgefilename);
  fout = fopen(outedgefilename, "w");
  fprintf(fout, "%d  %d\n", numberofedges, edgemarkerlist != NULL ? 1 : 0);
  for (i = 0; i < numberofedges; i++) {
    fprintf(fout, "%d  %4d  %4d", i + firstnumber, edgelist[i * 2],
            edgelist[i * 2 + 1]);
    if (edgemarkerlist != NULL) {
      fprintf(fout, "  %d", edgemarkerlist[i]);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_neighbors()    Save egdes to a .neigh file.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_neighbors(char* filebasename)
{
  FILE *fout;
  char outneighborfilename[FILENAMESIZE];
  int i;

  sprintf(outneighborfilename, "%s.neigh", filebasename);
  printf("Saving neighbors to %s\n", outneighborfilename);
  fout = fopen(outneighborfilename, "w");
  fprintf(fout, "%d  %d\n", numberoftetrahedra, mesh_dim + 1);
  for (i = 0; i < numberoftetrahedra; i++) {
    if (mesh_dim == 2) {
      fprintf(fout, "%d  %5d  %5d  %5d", i + firstnumber,  neighborlist[i * 3],
              neighborlist[i * 3 + 1], neighborlist[i * 3 + 2]);
    } else {
      fprintf(fout, "%d  %5d  %5d  %5d  %5d", i + firstnumber,
              neighborlist[i * 4], neighborlist[i * 4 + 1],
              neighborlist[i * 4 + 2], neighborlist[i * 4 + 3]);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_poly()    Save segments or facets to a .poly file.                   //
//                                                                           //
// It only save the facets, holes and regions. No .node file is saved.       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_poly(char* filebasename)
{
  FILE *fout;
  facet *f;
  polygon *p;
  char outpolyfilename[FILENAMESIZE];
  int i, j, k;

  sprintf(outpolyfilename, "%s.poly", filebasename);
  printf("Saving poly to %s\n", outpolyfilename);
  fout = fopen(outpolyfilename, "w");

  // The zero indicates that the vertices are in a separate .node file.
  //   Followed by number of dimensions, number of vertex attributes,
  //   and number of boundary markers (zero or one).
  fprintf(fout, "%d  %d  %d  %d\n", 0, mesh_dim, numberofpointattributes,
          pointmarkerlist != NULL ? 1 : 0);

  // Save segments or facets.
  if (mesh_dim == 2) {
    // Number of segments, number of boundary markers (zero or one).
    fprintf(fout, "%d  %d\n", numberofedges, edgemarkerlist != NULL ? 1 : 0);
    for (i = 0; i < numberofedges; i++) {
      fprintf(fout, "%d  %4d  %4d", i + firstnumber, edgelist[i * 2],
              edgelist[i * 2 + 1]);
      if (edgemarkerlist != NULL) {
        fprintf(fout, "  %d", edgemarkerlist[i]);
      }
      fprintf(fout, "\n");
    }
  } else {
    // Number of facets, number of boundary markers (zero or one).
    fprintf(fout, "%d  %d\n", numberoffacets, facetmarkerlist != NULL ? 1 : 0);
    for (i = 0; i < numberoffacets; i++) {
      f = &(facetlist[i]);
      fprintf(fout, "%d  %d  %d  # %d\n", f->numberofpolygons,f->numberofholes,
            facetmarkerlist != NULL ? facetmarkerlist[i] : 0, i + firstnumber);
      // Output polygons of this facet.
      for (j = 0; j < f->numberofpolygons; j++) {
        p = &(f->polygonlist[j]);
        fprintf(fout, "%d  ", p->numberofvertices);
        for (k = 0; k < p->numberofvertices; k++) {
          if (((k + 1) % 10) == 0) {
            fprintf(fout, "\n  ");
          }
          fprintf(fout, "  %d", p->vertexlist[k]);
        }
        fprintf(fout, "\n");
      }
      // Output holes of this facet.
      for (j = 0; j < f->numberofholes; j++) {
        fprintf(fout, "%d  %.12g  %.12g  %.12g\n", j + firstnumber,
           f->holelist[j * 3], f->holelist[j * 3 + 1], f->holelist[j * 3 + 2]);
      }
    }
  }

  // Save holes.
  fprintf(fout, "%d\n", numberofholes);
  for (i = 0; i < numberofholes; i++) {
    // Output x, y coordinates.
    fprintf(fout, "%d  %.12g  %.12g", i + firstnumber, holelist[i * mesh_dim],
            holelist[i * mesh_dim + 1]);
    if (mesh_dim == 3) {
      // Output z coordinate.
      fprintf(fout, "  %.12g", holelist[i * mesh_dim + 2]);
    }
    fprintf(fout, "\n");
  }

  // Save regions.
  fprintf(fout, "%d\n", numberofregions);
  for (i = 0; i < numberofregions; i++) {
    if (mesh_dim == 2) {
      // Output the index, x, y coordinates, attribute (region number)
      //   and maximum area constraint (maybe -1).
      fprintf(fout, "%d  %.12g  %.12g  %.12g  %.12g\n", i + firstnumber,
              regionlist[i * 4], regionlist[i * 4 + 1],
              regionlist[i * 4 + 2], regionlist[i * 4 + 3]);
    } else {
      // Output the index, x, y, z coordinates, attribute (region number)
      //   and maximum volume constraint (maybe -1).
      fprintf(fout, "%d  %.12g  %.12g  %.12g  %.12g  %.12g\n", i + firstnumber,
              regionlist[i * 5], regionlist[i * 5 + 1],
              regionlist[i * 5 + 2], regionlist[i * 5 + 3],
              regionlist[i * 5 + 4]);
    }
  }

  fclose(fout);  
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// save_faces2smesh()    Save triangular faces to a .smesh file.             //
//                                                                           //
// It only save the facets. No holes and regions. No .node file.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetgenio::save_faces2smesh(char* filebasename)
{
  FILE *fout;
  char outsmeshfilename[FILENAMESIZE];
  int i, j;

  sprintf(outsmeshfilename, "%s.smesh", filebasename);
  printf("Saving faces to %s\n", outsmeshfilename);
  fout = fopen(outsmeshfilename, "w");

  // The zero indicates that the vertices are in a separate .node file.
  //   Followed by number of dimensions, number of vertex attributes,
  //   and number of boundary markers (zero or one).
  fprintf(fout, "%d  %d  %d  %d\n", 0, mesh_dim, numberofpointattributes,
          pointmarkerlist != NULL ? 1 : 0);

  // Number of facets, number of boundary markers (zero or one).
  fprintf(fout, "%d  %d\n", numberoftrifaces, 
          trifacemarkerlist != NULL ? 1 : 0);

  // Output triangular facets.
  for (i = 0; i < numberoftrifaces; i++) {
    j = i * 3;
    fprintf(fout, "3  %d %d %d", trifacelist[j], trifacelist[j + 1], 
            trifacelist[j + 2]);
    if (trifacemarkerlist != NULL) {
      fprintf(fout, "  %d", trifacemarkerlist[i]);
    }
    fprintf(fout, "\n");
  }

  // No holes and regions.
  fprintf(fout, "0\n");
  fprintf(fout, "0\n");

  fclose(fout);  
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// readline()   Read a nonempty line from a file.                            //
//                                                                           //
// A line is considered "nonempty" if it contains something more than white  //
// spaces.  If a line is considered empty, it will be dropped and the next   //
// line will be read, this process ends until reaching the end-of-file or a  //
// non-empty line.  Return NULL if it is the end-of-file, otherwise, return  //
// a pointer to the first non-whitespace character of the line.              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

char* tetgenio::readline(char *string, FILE *infile, int *linenumber)
{
  char *result;

  // Search for a non-empty line.
  do {
    result = fgets(string, INPUTLINESIZE - 1, infile);
    if (linenumber) (*linenumber)++;
    if (result == (char *) NULL) {
      return (char *) NULL;
    }
    // Skip white spaces.
    while ((*result == ' ') || (*result == '\t')) result++;
    // If it's end of line, read another line and try again.
  } while ((*result == '\0') || (*result == '\r') || (*result == '\n'));
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// findnextfield()   Find the next field of a string.                        //
//                                                                           //
// Jumps past the current field by searching for whitespace or a comma, then //
// jumps past the whitespace or the comma to find the next field.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

char* tetgenio::findnextfield(char *string)
{
  char *result;

  result = string;
  // Skip the current field.  Stop upon reaching whitespace or a comma.
  while ((*result != '\0') && (*result != ' ') &&  (*result != '\t') && 
         (*result != ',') && (*result != ';')) {
    result++;
  }
  // Now skip the whitespace or the comma, stop at anything else that looks
  //   like a character, or the end of a line. 
  while ((*result == ' ') || (*result == '\t') || (*result == ',') ||
         (*result == ';')) {
    result++;
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// readnumberline()   Read a nonempty number line from a file.               //
//                                                                           //
// A line is considered "nonempty" if it contains something that looks like  //
// a number.  Comments (prefaced by `#') are ignored.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

char* tetgenio::readnumberline(char *string, FILE *infile, char *infilename)
{
  char *result;

  // Search for something that looks like a number.
  do {
    result = fgets(string, INPUTLINESIZE, infile);
    if (result == (char *) NULL) {
      return result;
    }
    // Skip anything that doesn't look like a number, a comment, 
    //   or the end of a line. 
    while ((*result != '\0') && (*result != '#')
           && (*result != '.') && (*result != '+') && (*result != '-')
           && ((*result < '0') || (*result > '9'))) {
      result++;
    }
    // If it's a comment or end of line, read another line and try again.
  } while ((*result == '#') || (*result == '\0'));
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// findnextnumber()   Find the next field of a number string.                //
//                                                                           //
// Jumps past the current field by searching for whitespace or a comma, then //
// jumps past the whitespace or the comma to find the next field that looks  //
// like a number.                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

char* tetgenio::findnextnumber(char *string)
{
  char *result;

  result = string;
  // Skip the current field.  Stop upon reaching whitespace or a comma.
  while ((*result != '\0') && (*result != '#') && (*result != ' ') && 
         (*result != '\t') && (*result != ',')) {
    result++;
  }
  // Now skip the whitespace and anything else that doesn't look like a
  //   number, a comment, or the end of a line. 
  while ((*result != '\0') && (*result != '#')
         && (*result != '.') && (*result != '+') && (*result != '-')
         && ((*result < '0') || (*result > '9'))) {
    result++;
  }
  // Check for a comment (prefixed with `#').
  if (*result == '#') {
    *result = '\0';
  }
  return result;
}

////                                                                       ////
////                                                                       ////
//// io_cxx ///////////////////////////////////////////////////////////////////

