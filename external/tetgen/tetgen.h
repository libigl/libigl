///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator         //
//                                                                           //
// Version 1.5                                                               //
// October 06, 2012                                                          //
//                                                                           //
// Copyright (C) 2002--2012                                                  //
// Hang Si                                                                   //
// Research Group: Numerical Mathematics and Scientific Computing            //
// Weierstrass Institute for Applied Analysis and Stochastics (WIAS)         //
// Mohrenstr. 39, 10117 Berlin, Germany                                      //
// Hang.Si@wias-berlin.de                                                    //
//                                                                           //
// TetGen is freely available through the website: http://www.tetgen.org.    //
//   It may be copied, modified, and redistributed for non-commercial use.   //
//   Please consult the file LICENSE for the detailed copyright notices.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef tetgenH
#define tetgenH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h> 

// To compile TetGen as a library instead of an executable program, define
//   the TETLIBRARY symbol.

// #define TETLIBRARY

// Uncomment the following line to disable assert macros. These macros were
//   inserted in the code where I hoped to catch bugs. They may slow down the
//   speed of TetGen.

// #define NDEBUG

// TetGen uses the double precision for a real number. 

#define REAL double

// The types 'intptr_t' and 'uintptr_t' are signed and unsigned integer types,
//   respectively. They are guaranteed to be the same width as a pointer.
//   They are defined in <stdint.h> by the C99 Standard. However, Microsoft 
//   Visual C++ 2003 -- 2008 (Visual C++ 7.1 - 9) doesn't ship with this header
//   file. In such case, we can define them by ourself. 
// Update (learned from Stack Overflow): Visual Studio 2010 and Visual C++ 2010
//   Express both have stdint.h

// The following piece of code was provided by Steven Johnson (MIT). Define the
//   symbol _MSC_VER if you are using Microsoft Visual C++. Moreover, define 
//   the _WIN64 symbol if you are running TetGen on Win64 systems.

#ifdef _MSC_VER // Microsoft Visual C++
#  ifdef _WIN64
     typedef __int64 intptr_t;
     typedef unsigned __int64 uintptr_t;
#  else // not _WIN64
     typedef int intptr_t;
     typedef unsigned int uintptr_t;
#  endif
#else // not Visual C++
#  include <stdint.h>
#endif

// Maximum number of characters in a file name (including the null).

#define FILENAMESIZE 1024

// Maximum number of chars in a line read from a file (including the null).

#define INPUTLINESIZE 2048

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgenio                                                                  //
//                                                                           //
// A structure for transfering data into and out of TetGen's mesh structure. //
//                                                                           //
// It holds a collection of arrays of data, i.e., points, facets, tetrahedra,//
// and so forth. It contains functions to read and write (input and output)  //
// files of TetGen as well as other supported mesh files.                    //
//                                                                           //
// Once an object of tetgenio is declared,  no array is created. One has to  //
// allocate enough memory for them. On deletion of this object, the memory   //
// occupied by these arrays needs to be freed.  The routine deinitialize()   //
// will be automatically called.  It frees the memory for an array if it is  //
// not a NULL. Note that it assumes that the memory is allocated by the C++  //
// "new" operator.  Otherwise, the user must priorily free them by theirself //
// and set the pointers to NULLs.                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenio {

public:

  // A "polygon" describes a simple polygon (no holes). It is not necessarily
  //   convex. Each polygon contains a number of corners (points) and the same
  //   number of sides (edges).  The points of the polygon must be given in
  //   either counterclockwise or clockwise order and they form a ring, so 
  //   every two consective points forms an edge of the polygon.
  typedef struct {
    int *vertexlist;
    int numberofvertices;
  } polygon;

  // A "facet" describes a polygonal region possibly with holes, edges, and 
  //   points floating in it.  Each facet consists of a list of polygons and
  //   a list of hole points (which lie strictly inside holes).
  typedef struct {
    polygon *polygonlist;
    int numberofpolygons;
    REAL *holelist;
    int numberofholes;
  } facet;

  // A "voroedge" is an edge of the Voronoi diagram. It corresponds to a
  //   Delaunay face.  Each voroedge is either a line segment connecting
  //   two Voronoi vertices or a ray starting from a Voronoi vertex to an
  //   "infinite vertex".  'v1' and 'v2' are two indices pointing to the
  //   list of Voronoi vertices. 'v1' must be non-negative, while 'v2' may
  //   be -1 if it is a ray, in this case, the unit normal of this ray is
  //   given in 'vnormal'. 
  typedef struct {
    int v1, v2;
    REAL vnormal[3];
  } voroedge;

  // A "vorofacet" is an facet of the Voronoi diagram. It corresponds to a
  //   Delaunay edge.  Each Voronoi facet is a convex polygon formed by a
  //   list of Voronoi edges, it may not be closed.  'c1' and 'c2' are two
  //   indices pointing into the list of Voronoi cells, i.e., the two cells
  //   share this facet.  'elist' is an array of indices pointing into the
  //   list of Voronoi edges, 'elist[0]' saves the number of Voronoi edges
  //   (including rays) of this facet.
  typedef struct {
    int c1, c2;
    int *elist;
  } vorofacet;

  // The periodic boundary condition group data structure.  A "pbcgroup"
  //   contains the definition of a pbc and the list of pbc point pairs.
  //   'fmark1' and 'fmark2' are the facetmarkers of the two pbc facets f1
  //   and f2, respectively. 'transmat' is the transformation matrix which
  //   maps a point in f1 into f2.  An array of pbc point pairs are saved
  //   in 'pointpairlist'. The first point pair is at indices [0] and [1],
  //   followed by remaining pairs. Two integers per pair.
  typedef struct {
    int fmark1, fmark2;
    REAL transmat[4][4];
    int numberofpointpairs;
    int *pointpairlist;
  } pbcgroup;

  // Additional parameters associated with an input (or mesh) vertex.
  //   These informations are provided by CAD libraries. 
  typedef struct {
    REAL uv[2];
    int tag;
    int type; // 0, 1, or 2.
  } pointparam;

  // A callback function for mesh refinement.
  typedef bool (* TetSizeFunc)(REAL*, REAL*, REAL*, REAL*, REAL*, REAL);

  // Callback functions for meshing PSCs.
  typedef REAL (* GetVertexParamOnEdge)(void*, int, int);
  typedef void (* GetSteinerOnEdge)(void*, int, REAL, REAL*);
  typedef void (* GetVertexParamOnFace)(void*, int, int, REAL*);
  typedef void (* GetEdgeSteinerParamOnFace)(void*, int, REAL, int, REAL*);
  typedef void (* GetSteinerOnFace)(void*, int, REAL*, REAL*);

  // Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
  int firstnumber; 

  // Dimension of the mesh (2 or 3), default is 3.
  int mesh_dim;

  // Does the lines in .node file contain index or not, default is 1.
  int useindex;

  // 'pointlist':  An array of point coordinates.  The first point's x
  //   coordinate is at index [0] and its y coordinate at index [1], its
  //   z coordinate is at index [2], followed by the coordinates of the
  //   remaining points.  Each point occupies three REALs. 
  // 'pointattributelist':  An array of point attributes.  Each point's
  //   attributes occupy 'numberofpointattributes' REALs.
  // 'pointmtrlist': An array of metric tensors at points. Each point's
  //   tensor occupies 'numberofpointmtr' REALs.
  // 'pointmarkerlist':  An array of point markers; one integer per point.
  REAL *pointlist;
  REAL *pointattributelist;
  REAL *pointmtrlist;
  int  *pointmarkerlist;
  pointparam *pointparamlist;
  int numberofpoints;
  int numberofpointattributes;
  int numberofpointmtrs;
 
  // 'elementlist':  An array of element (tetrahedron) corners.  The first 
  //   element's first corner is at index [0], followed by its other corners,
  //   followed by any other nodes if the element represents a nonlinear 
  //   element.  Each element occupies 'numberofcorners' ints.
  // 'elementattributelist':  An array of element attributes.  Each
  //   element's attributes occupy 'numberofelementattributes' REALs.
  // 'elementconstraintlist':  An array of constraints, i.e. tetrahedron's
  //   volume; one REAL per element.  Input only.
  // 'neighborlist':  An array of element neighbors; 4 ints per element. 
  //   Output only.
  int  *tetrahedronlist;
  REAL *tetrahedronattributelist;
  REAL *tetrahedronvolumelist;
  int  *neighborlist;
  int numberoftetrahedra;
  int numberofcorners;
  int numberoftetrahedronattributes;

  // 'facetlist':  An array of facets.  Each entry is a structure of facet.
  // 'facetmarkerlist':  An array of facet markers; one int per facet.
  facet *facetlist;
  int *facetmarkerlist;
  int numberoffacets;

  // 'holelist':  An array of holes (in volume).  Each hole is given by a
  //   seed (point) which lies strictly inside it. The first seed's x, y and z
  //   coordinates are at indices [0], [1] and [2], followed by the
  //   remaining seeds.  Three REALs per hole. 
  REAL *holelist;
  int numberofholes;

  // 'regionlist': An array of regions (subdomains).  Each region is given by
  //   a seed (point) which lies strictly inside it. The first seed's x, y and
  //   z coordinates are at indices [0], [1] and [2], followed by the regional
  //   attribute at index [3], followed by the maximum volume at index [4]. 
  //   Five REALs per region.
  // Note that each regional attribute is used only if you select the 'A'
  //   switch, and each volume constraint is used only if you select the
  //   'a' switch (with no number following).
  REAL *regionlist;
  int numberofregions;

  // 'facetconstraintlist':  An array of facet constraints.  Each constraint
  //   specifies a maximum area bound on the subfaces of that facet.  The
  //   first facet constraint is given by a facet marker at index [0] and its
  //   maximum area bound at index [1], followed by the remaining facet con-
  //   straints. Two REALs per facet constraint.  Note: the facet marker is
  //   actually an integer.
  REAL *facetconstraintlist;
  int numberoffacetconstraints;

  // 'segmentconstraintlist': An array of segment constraints. Each constraint 
  //   specifies a maximum length bound on the subsegments of that segment.
  //   The first constraint is given by the two endpoints of the segment at
  //   index [0] and [1], and the maximum length bound at index [2], followed
  //   by the remaining segment constraints.  Three REALs per constraint. 
  //   Note the segment endpoints are actually integers.
  REAL *segmentconstraintlist;
  int numberofsegmentconstraints;

  // 'pbcgrouplist':  An array of periodic boundary condition groups.
  pbcgroup *pbcgrouplist;
  int numberofpbcgroups;

  // 'trifacelist':  An array of face (triangle) corners.  The first face's
  //   corners are at indices [0], [1] and [2], followed by the remaining 
  //   faces.  Three ints per face.
  // 'adjtetlist':  An array of adjacent tetrahedra to the faces. The first
  //   face's two adjacent tetrahedra are at indices [0] and [1], followed by
  //   the remaining faces.  A '-1' indicates outside (no adj. tet). This list
  //   is output when '-nn' switch is used. Output only.
  // 'trifacemarkerlist':  An array of face markers; one int per face.
  int *trifacelist;
  int *adjtetlist;
  int *trifacemarkerlist;
  int numberoftrifaces;

  // `edgelist':  An array of edge endpoints.  The first edge's endpoints
  //   are at indices [0] and [1], followed by the remaining edges.  Two
  //   ints per edge.
  // `edgemarkerlist':  An array of edge markers; one int per edge.
  int *edgelist;
  int *edgemarkerlist;
  int numberofedges;

  // 'vpointlist':  An array of Voronoi vertex coordinates (like pointlist).
  // 'vedgelist':  An array of Voronoi edges.  Each entry is a 'voroedge'.
  // 'vfacetlist':  An array of Voronoi facets. Each entry is a 'vorofacet'.
  // 'vcelllist':  An array of Voronoi cells.  Each entry is an array of
  //   indices pointing into 'vfacetlist'. The 0th entry is used to store
  //   the length of this array.
  REAL *vpointlist;
  voroedge *vedgelist;
  vorofacet *vfacetlist;
  int **vcelllist;
  int numberofvpoints;
  int numberofvedges;
  int numberofvfacets;
  int numberofvcells;

  // A callback function.
  TetSizeFunc tetunsuitable;

  // Variable (and callback functions) for meshing PSCs.
  void *geomhandle;
  GetVertexParamOnEdge getvertexparamonedge;
  GetSteinerOnEdge getsteineronedge;
  GetVertexParamOnFace getvertexparamonface;
  GetEdgeSteinerParamOnFace getedgesteinerparamonface;
  GetSteinerOnFace getsteineronface;

  // Input & output routines.
  bool load_node_call(FILE* infile, int markers, int uvflag, char*);
  bool load_node(char*);
  bool load_edge(char*);
  bool load_face(char*);
  bool load_tet(char*);
  bool load_vol(char*);
  bool load_var(char*);
  bool load_mtr(char*);
  bool load_pbc(char*);
  bool load_poly(char*);
  bool load_off(char*);
  bool load_ply(char*);
  bool load_stl(char*);
  bool load_vtk(char*);
  bool load_medit(char*, int);
  bool load_plc(char*, int);
  bool load_tetmesh(char*, int);
  void save_nodes(char*);
  void save_elements(char*);
  void save_faces(char*);
  void save_edges(char*);
  void save_neighbors(char*);
  void save_poly(char*);
  void save_faces2smesh(char*);

  // Read line and parse string functions.
  char *readline(char* string, FILE* infile, int *linenumber);
  char *findnextfield(char* string);
  char *readnumberline(char* string, FILE* infile, char* infilename);
  char *findnextnumber(char* string);

  static void init(polygon* p) {
    p->vertexlist = (int *) NULL;
    p->numberofvertices = 0;
  }

  static void init(facet* f) {
    f->polygonlist = (polygon *) NULL;
    f->numberofpolygons = 0;
    f->holelist = (REAL *) NULL;
    f->numberofholes = 0;
  }

  // Initialize routine.
  void initialize()
  {
    firstnumber = 0; // Default item index is numbered from Zero.
    mesh_dim = 3; // Default mesh dimension is 3.
    useindex = 1;

    pointlist = (REAL *) NULL;
    pointattributelist = (REAL *) NULL;
    pointmtrlist = (REAL *) NULL;
    pointmarkerlist = (int *) NULL;
    pointparamlist = (pointparam *) NULL;
    numberofpoints = 0;
    numberofpointattributes = 0;
    numberofpointmtrs = 0;

    tetrahedronlist = (int *) NULL;
    tetrahedronattributelist = (REAL *) NULL;
    tetrahedronvolumelist = (REAL *) NULL;
    neighborlist = (int *) NULL;
    numberoftetrahedra = 0;
    numberofcorners = 4; // Default is 4 nodes per element.
    numberoftetrahedronattributes = 0;

    trifacelist = (int *) NULL;
    adjtetlist = (int *) NULL;
    trifacemarkerlist = (int *) NULL;
    numberoftrifaces = 0; 

    facetlist = (facet *) NULL;
    facetmarkerlist = (int *) NULL;
    numberoffacets = 0; 

    edgelist = (int *) NULL;
    edgemarkerlist = (int *) NULL;
    numberofedges = 0;

    holelist = (REAL *) NULL;
    numberofholes = 0;

    regionlist = (REAL *) NULL;
    numberofregions = 0;

    facetconstraintlist = (REAL *) NULL;
    numberoffacetconstraints = 0;
    segmentconstraintlist = (REAL *) NULL;
    numberofsegmentconstraints = 0;

    pbcgrouplist = (pbcgroup *) NULL;
    numberofpbcgroups = 0;

    vpointlist = (REAL *) NULL;
    vedgelist = (voroedge *) NULL;
    vfacetlist = (vorofacet *) NULL; 
    vcelllist = (int **) NULL; 
    numberofvpoints = 0;
    numberofvedges = 0;
    numberofvfacets = 0;
    numberofvcells = 0;

    tetunsuitable = NULL;

    geomhandle = NULL;
    getvertexparamonedge = NULL;
    getsteineronedge = NULL;
    getvertexparamonface = NULL;
    getedgesteinerparamonface = NULL;
    getsteineronface = NULL;
  }

  // Free the memory allocated in 'tetgenio'.  
  void deinitialize()
  {
    facet *f;
    polygon *p;
    pbcgroup *pg;
    int i, j;

    // Notice that this routine assumes that the memory was allocated by 
    //   C++ memory allocation operator 'new'.

    if (pointlist != (REAL *) NULL) {
      delete [] pointlist;
    }
    if (pointattributelist != (REAL *) NULL) {
      delete [] pointattributelist;
    }
    if (pointmtrlist != (REAL *) NULL) {
      delete [] pointmtrlist;
    }
    if (pointmarkerlist != (int *) NULL) {
      delete [] pointmarkerlist;
    }
    if (pointparamlist != (pointparam *) NULL) {
      delete [] pointparamlist;
    }

    if (tetrahedronlist != (int *) NULL) {
      delete [] tetrahedronlist;
    }
    if (tetrahedronattributelist != (REAL *) NULL) {
      delete [] tetrahedronattributelist;
    }
    if (tetrahedronvolumelist != (REAL *) NULL) {
      delete [] tetrahedronvolumelist;
    }
    if (neighborlist != (int *) NULL) {
      delete [] neighborlist;
    }

    if (trifacelist != (int *) NULL) {
      delete [] trifacelist;
    }
    if (adjtetlist != (int *) NULL) {
      delete [] adjtetlist;
    }
    if (trifacemarkerlist != (int *) NULL) {
      delete [] trifacemarkerlist;
    }

    if (edgelist != (int *) NULL) {
      delete [] edgelist;
    }
    if (edgemarkerlist != (int *) NULL) {
      delete [] edgemarkerlist;
    }

    if (facetlist != (facet *) NULL) {
      for (i = 0; i < numberoffacets; i++) {
        f = &facetlist[i];
        for (j = 0; j < f->numberofpolygons; j++) {
          p = &f->polygonlist[j];
          delete [] p->vertexlist;
        }
        delete [] f->polygonlist;
        if (f->holelist != (REAL *) NULL) {
          delete [] f->holelist;
        }
      }
      delete [] facetlist;
    }
    if (facetmarkerlist != (int *) NULL) {
      delete [] facetmarkerlist;
    }

    if (holelist != (REAL *) NULL) {
      delete [] holelist;
    }
    if (regionlist != (REAL *) NULL) {
      delete [] regionlist;
    }
    if (facetconstraintlist != (REAL *) NULL) {
      delete [] facetconstraintlist;
    }
    if (segmentconstraintlist != (REAL *) NULL) {
      delete [] segmentconstraintlist;
    }
    if (pbcgrouplist != (pbcgroup *) NULL) {
      for (i = 0; i < numberofpbcgroups; i++) {
        pg = &(pbcgrouplist[i]);
        if (pg->pointpairlist != (int *) NULL) {
          delete [] pg->pointpairlist;
        }
      }
      delete [] pbcgrouplist;
    }
    if (vpointlist != (REAL *) NULL) {
      delete [] vpointlist;
    }
    if (vedgelist != (voroedge *) NULL) {
      delete [] vedgelist;
    }
    if (vfacetlist != (vorofacet *) NULL) {
      for (i = 0; i < numberofvfacets; i++) {
        delete [] vfacetlist[i].elist;
      }
      delete [] vfacetlist;
    }
    if (vcelllist != (int **) NULL) {
      for (i = 0; i < numberofvcells; i++) {
        delete [] vcelllist[i];
      }
      delete [] vcelllist;
    }
  }

  // Constructor & destructor.
  tetgenio() {initialize();}
  ~tetgenio() {deinitialize();}

}; // class tetgenio

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgenbehavior                                                            //
//                                                                           //
// A structure for maintaining the switches and parameters used by TetGen's  //
// meshing algorithms.  They are specified by the command line arguments.    //
//                                                                           //
// NOTE: Some of the switches are incompatinle to each other, while some are //
// depend on others.  The routine parse_commandline() sets the switches from //
// the command line (a list of strings). Morover, it checks the consistency  //
// of the applied switches.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenbehavior {

public:

  // The list of switches of TetGen. 
  int plc;                                                         // '-p', 0.
  int psc;                                                         // '-s', 0.
  int refine;                                                      // '-r', 0.
  int quality;                                                     // '-q', 0.
  int nobisect;                                                    // '-Y', 0.
  int weighted;                                                    // '-w', 0.
  int brio_hilbert;                                                // '-b', 1.
  int incrflip;                                                    // '-l', 0.
  int flipinsert;                                                  // '-L', 0.
  int metric;                                                      // '-m', 0.
  int varvolume;                                                   // '-a', 0.
  int fixedvolume;                                                 // '-a', 0.
  int regionattrib;                                                // '-A', 0.
  int conforming;                                                  // '-D', 0.
  int insertaddpoints;                                             // '-i', 0.
  int diagnose;                                                    // '-d', 0.
  int convex;                                                      // '-c', 0.
  int nomerge;                                                     // '-M', 0.
  int noexact;                                                     // '-X', 0.
  int nostaticfilter;                                              // '-X', 0.
  int zeroindex;                                                   // '-z', 0.
  int facesout;                                                    // '-f', 0.
  int edgesout;                                                    // '-e', 0.
  int neighout;                                                    // '-n', 0.
  int voroout;                                                     // '-v', 0.
  int meditview;                                                   // '-g', 0.
  int vtkview;                                                     // '-k', 0.
  int nobound;                                                     // '-B', 0.
  int nonodewritten;                                               // '-N', 0.
  int noelewritten;                                                // '-E', 0.
  int nofacewritten;                                               // '-F', 0.
  int noiterationnum;                                              // '-I', 0.
  int nojettison;                                                  // '-J', 0.
  int reversetetori;                                               // '-R', 0.
  int docheck;                                                     // '-C', 0.
  int quiet;                                                       // '-Q', 0.
  int verbose;                                                     // '-V', 0.

  // The list of parameters of TetGen. 
  int vertexperblock;                                                 // 4092.
  int tetrahedraperblock;                                             // 8188.
  int shellfaceperblock;                                              // 4092.
  int nobisect_param;                                              // '-Y', 1.
  int weighted_param;                                              // '-w', 0.
  int hilbert_order;                                                    // -1.
  int hilbert_limit;                                                    //  8.
  int fliplinklevel;                                                    // -1.
  int flipstarsize;                                                     // -1.
  int fliplinklevelinc;                                                 //  1.
  int reflevel;                                                    // '-D', 3.
  int optlevel;                                                    // '-O', 2.
  int optscheme;                                                   // '-O', 7.
  int delmaxfliplevel;                                                   // 1.
  int order;                                                       // '-o', 1.
  int steinerleft;                                                 // '-S', 0.
  REAL facet_ang_tol;                                          // '-p', 179.9.
  REAL maxvolume;                                               // '-a', -1.0.
  REAL minratio;                                                 // '-q', 0.0.
  REAL mindihedral;                                              // '-q', 5.0.
  REAL optmaxdihedral;                                               // 165.0.
  REAL optminsmtdihed;                                               // 179.0.
  REAL optminslidihed;                                               // 179.0.  
  REAL epsilon;                                               // '-T', 1.0e-8.
  REAL minedgelength;                                                  // 0.0.

  // Strings of command line arguments and input/output file names.
  char commandline[1024];
  char infilename[1024];
  char outfilename[1024];
  char addinfilename[1024];
  char bgmeshfilename[1024];

  // The input object of TetGen. They are recognized by either the input 
  //   file extensions or by the specified options. 
  enum objecttype {NODES, POLY, OFF, PLY, STL, MEDIT, VTK, MESH} object;


  void syntax();
  void usage();

  // Command line parse routine.
  bool parse_commandline(int argc, char **argv);
  bool parse_commandline(char *switches) {
    return parse_commandline(0, &switches);
  }

  // Initialize all variables.
  tetgenbehavior()
  {
    plc = 0;
    psc = 0;
    refine = 0;
    quality = 0;
    nobisect = 0;
    metric = 0;
    weighted = 0;
    brio_hilbert = 1;
    incrflip = 0;
    flipinsert = 0;
    varvolume = 0;
    fixedvolume = 0;
    noexact = 0;
    nostaticfilter = 0;
    insertaddpoints = 0;
    regionattrib = 0;
    conforming = 0;
    diagnose = 0;
    convex = 0;
    zeroindex = 0;
    facesout = 0;
    edgesout = 0;
    neighout = 0;
    voroout = 0;
    meditview = 0;
    vtkview = 0;
    nobound = 0;
    nonodewritten = 0;
    noelewritten = 0;
    nofacewritten = 0;
    noiterationnum = 0;
    nomerge = 0;
    nojettison = 0;
    reversetetori = 0;
    docheck = 0;
    quiet = 0;
    verbose = 0;

    vertexperblock = 4092;
    tetrahedraperblock = 8188;
    shellfaceperblock = 4092;
    nobisect_param = 1;
    weighted_param = 0;
    hilbert_order = -1;
    hilbert_limit = 8;
    fliplinklevel = -1; // No limit on linklevel.
    flipstarsize = -1;  // No limit on flip star size.
    fliplinklevelinc = 1;
    reflevel = 3;
    optscheme = 7;  // 1 & 2 & 4, // min_max_dihedral.
    optlevel = 2;
    delmaxfliplevel = 1;
    order = 1;
    steinerleft = -1;
    facet_ang_tol = 179.9;
    maxvolume = -1.0;
    minratio = 2.0;
    mindihedral = 5.0; 
    optmaxdihedral = 165.00; // without -q, default is 179.0
    optminsmtdihed = 179.00; // without -q, default is 179.999
    optminslidihed = 179.00; // without -q, default is 179.999
    epsilon = 1.0e-8;
    minedgelength = 0.0;
    object = NODES;

    commandline[0] = '\0';
    infilename[0] = '\0';
    outfilename[0] = '\0';
    addinfilename[0] = '\0';
    bgmeshfilename[0] = '\0';

  }

}; // class tetgenbehavior

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Robust Geometric predicates                                               //
//                                                                           //
// Geometric predicates are simple tests of spatial relations of a set of d- //
// dimensional points, such as the orientation test and the point-in-sphere  //
// test. Each of these tests is performed by evaluating the sign of a deter- //
// minant of a matrix whose entries are the coordinates of these points.  If //
// the computation is performed by using the floating-point numbers, e.g.,   //
// the single or double numbers in C/C++, roundoff error may cause an incor- //
// rect result. This may either lead to a wrong result or eventually lead to //
// a failure of the program.                                                 //
//                                                                           //
// Various techniques are developed to avoid roundoff errors, such as exact  //
// multi-precision computations, interval arthmetics, adaptive exact arthme- //
// tics, and filtered exact arthmetics, etc. Devillers and Pion give a nice  //
// discussion and comparisons of these techniques for robustly computing the //
// Delaunay triangulations [Devillers and Pion 2002].                        //
//                                                                           //
// The following routines implemented the orientation test and the point-in- //
// sphere test use the adaptive exact floating-point arithmetics [Shewchuk   //
// 1997]. They are generously provided by Jonathan Schewchuk in the public   //
// domain, http://www.cs.cmu.edu/~quake/robust.html. The source code are in  //
// file "predicates.cxx".                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void exactinit(int, int, REAL, REAL, REAL);
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
REAL orient4d(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
              REAL ah, REAL bh, REAL ch, REAL dh, REAL eh);
void predicates_statistics(int weighted);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgenmesh                                                                //
//                                                                           //
// A structure containing the mesh data structure and the implementations of //
// tetrahedral meshing algorithms of TetGen.                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenmesh {

public:

  // Labels that signify the type of a vertex. 
  enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, RIDGEVERTEX, ACUTEVERTEX,
                 FACETVERTEX, VOLVERTEX, FREESEGVERTEX, FREEFACETVERTEX, 
                 FREEVOLVERTEX, NREGULARVERTEX, DEADVERTEX};
 
  // Labels that signify the type of a subsegment.
  enum shestype {NSHARP, SHARP, FAKESH};

  // Labels that signify the result of triangle-triangle intersection test.
  enum interresult {DISJOINT, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
                    TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE, 
                    COLLISIONFACE, ACROSSSEG, ACROSSSUB};

  // Labels that signify the result of point location.
  enum locateresult {OUTSIDE, INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX, INSTAR,
                     ENCVERTEX, ENCSEGMENT, ENCSUBFACE, NEARVERTEX,
                     NONREGULAR, BADELEMENT};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh data structure                                                       //
//                                                                           //
// A tetrahedral mesh of a 3D domain is a 3D simplicial complex T whose und- //
// erlying space is homeomorphic to the domain. T contains a 2D subcomplex S //
// which is a triangular mesh of the boundary of the domain. S contains a 1D //
// subcomplex L which is a linear mesh of the boundary of the surface. Faces //
// and edges in S and L are respectivly called subfaces and segments to dis- //
// tinguish them from others in T.                                           //
//                                                                           //
// TetGen stores the tetrahedra and vertices of T. Each tetrahedron contains //
// pointers to its vertices and adjacent tetrahedra.  Each vertex stores its //
// x-, y-, and z-coordinates. The faces and edges of T are implicitly repre- //
// sented by tetrahedra. 
//                                                                           //
// Each face of T belongs to either two tetrahedra or one tetrahedron. In    //
// the latter case, the face is an exterior boundary face of T.  TetGen adds //
// fictitious tetrahedra (one-to-one) at such faces, and connects them to an //
// "infinite vertex" (which has no geometric coordinates).  One can imagine  //
// such a vertex lies in 4D space and is visible by all exterior boundary    //
// faces.  The extended set of tetrahedra (including the infinite vertex) is //
// a tetrahedralization of a compact 3-manifold without bounday.  It has the //
// property that every face is shared by exactly two tetrahedra.             // 
//                                                                           //
// TetGen stores explicitly the subfaces and segments (which are in surface  //
// mesh S and the linear mesh L, respectively. Additional informations are   //
// stored in tetrahedra and subfaces to remember their relations.            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // The tetrahedron data structure.  It includes the following fields:
  //   - a list of four adjoining tetrahedra;
  //   - a list of four vertices;
  //   - a list of four subfaces (optional, for -p switch);
  //   - a list of  six segments (optional, for -p switch);
  //   - a list of user-defined floating-point attributes (optional);
  //   - a volume constraint (optional, for -a switch);
  //   - an integer of element marker (and flags);
  // The structure of a tetrahedron is an array of pointers.  Its actual size
  //   (the length of the array) is determined at runtime.

  typedef REAL **tetrahedron;

  // The subface data structure.  It includes the following fields:
  //   - a list of three adjoining subfaces;
  //   - a list of three vertices;
  //   - a list of three adjoining segments;
  //   - two adjoining tetrahedra;
  //   - an area constraint (optional, for -q switch);
  //   - an integer for boundary marker;
  //   - an integer for type: SHARPSEGMENT, NONSHARPSEGMENT, ...;
  //   - an integer for pbc group (optional, if in->pbcgrouplist exists);

  typedef REAL **shellface;

  // The point data structure.  It includes the following fields:
  //   - x, y and z coordinates;
  //   - a list of user-defined point attributes (optional);
  //   - u, v coordinates (optional, for -s switch);
  //   - a metric tensor (optional, for -q or -m switch);
  //   - a pointer to an adjacent tetrahedron;
  //   - a pointer to a parent (or a duplicate) point;
  //   - a pointer to an adjacent subface or segment (optional, -p switch);
  //   - a pointer to a tet in background mesh (optional, for -m switch);
  //   - an integer for boundary marker (point index);
  //   - an integer for point type (and flags).
  //   - an integer for geometry tag (optional, for -s switch).
  // The structure of a point is an array of REALs.  Its acutal size is 
  //   determined at the runtime.

  typedef REAL *point;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Ordered tetrahedra                                                        //
//                                                                           //
// The four vertices of a tetrahedron can be permuted in 24 different seque- //
// nces.  We call each sequence resulted by an even permutation an "ordered  //
// tetrahedron".  There are total 12 ordered tetrahedra.  They form a group  //
// which is isomorphic to the alternating group of 4 elements. Geometrically,//
// if we direct the three edges within a face of a tetrahedron by the count- //
// erclockwise order viewed from the opposite vertex of this face (using ei- //
// ther right-hand or left-hand rule). There are total twelve directed edges //
// in the tetrahedron. Each of them corresponds to an ordered tetrahedron.   //
//                                                                           //
// We represent an order tetrahedron by a pair (t, v), where t is a pointer  //
// to the tetrahedron and v is a four-bit integer, in the range from 0 to 11,//
// identifying the ordered version of the tetrahedron.  Assume the faces of  //
// the tetrahedron is numbered from 0 to 3, and the edges in a face is numb- //
// ered from 0 to 2.  Then the two lower bits of v encode the face number,   //
// and the two higher bits of v encode the edge number in that face.         //
//                                                                           //
// The four vertices of a tetrahedron are indexed from 0 to 3 (accodring to  //
// their storage in the data structure).  Give each face the same index as   //
// the node opposite it in the tetrahedron.  Denote the edge connecting face //
// i to face j as i/j. We number the twelve versions as follows:             //
//                                                                           //
//           |   edge 0     edge 1     edge 2                                //
//   --------|--------------------------------                               //
//    face 0 |   0 (0/1)    4 (0/3)    8 (0/2)                               //
//    face 1 |   1 (1/2)    5 (1/3)    9 (1/0)                               //
//    face 2 |   2 (2/3)    6 (2/1)   10 (2/0)                               //
//    face 3 |   3 (3/0)    7 (3/1)   11 (3/2)                               //
//                                                                           //
// Ordered triangles                                                         //
//                                                                           //
// The three vertices of a triangle can be permuted in 6 different sequences //
// which form a group isomorphic to the symmetric group of 3 elements. Each  //
// permutation of the vertices is called an ordered triangle.  The first two //
// vertices of an ordered triangle defines an directed edge. There are total //
// six directed edge in the triangle.  They can be divided into two groups,  //
// which correspond the two orientations of the triangle, respectively.      //
//                                                                           //
// We represent an ordered triangle by a pair (s, v),  where s is a pointer  //
// to the triangle and v is a three-bit integer, in the range from 0 to 5,   //
// identifying the directed edge of the triangle.  Using the first bit of v  //
// to identify the orientation, the other two bits of v identify the edge.   //
//                                                                           //
// Number the three vertices of a triangle from 0 to 2 (according to their   //
// storage in the data structure). Give each edge the same index as the node //
// opposite it in the triangle. The six versions of a triangle are:          //
//                                                                           //
//                 | edge 0   edge 1   edge 2                                //
//  ---------------|--------------------------                               //
//   ccw orieation |   0        2        4                                   //
//    cw orieation |   1        3        5                                   //
//                                                                           //
// In the following, a 'triface' is an order tetrahedron, and a 'face' is an //
// orider triangle.                                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class triface {
  public:
    tetrahedron *tet;
    int ver; // Range from 0 to 11.
    triface() : tet(0), ver(0) {}
    triface& operator=(const triface& t) {
      tet = t.tet; ver = t.ver;
      return *this;
    }
  };

  class face {
  public:
    shellface *sh;
    int shver; // Range from 0 to 5.
    face() : sh(0), shver(0) {}
    face& operator=(const face& s) {
      sh = s.sh; shver = s.shver;
      return *this;
    }
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// badface                                                                   //
//                                                                           //
// A multiple usages structure. Despite of its name, a 'badface' can be used //
// to represent the following objects:                                       //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;             //
//   - an encroached subsegment or subface;                                  //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;      //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;    //
//   - a recently flipped face (saved for undoing the flip later).           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class badface {
  public:
    triface tt; 
    face ss; 
    REAL key, cent[6];  // circumcenter or cos(dihedral angles) at 6 edges.
    point forg, fdest, fapex, foppo, noppo;
    badface *previtem, *nextitem; 
    badface() : key(0), forg(0), fdest(0), fapex(0), foppo(0), noppo(0),
      previtem(0), nextitem(0) {}
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// insertvertexflags                                                         //
//                                                                           //
// A collection of flags that pass to the routine insertvertex().            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class insertvertexflags {

  public:

    int iloc, bowywat, lawson;
    int rejflag, chkencflag;
    int sloc, sbowywat;
    int splitbdflag, validflag, respectbdflag;
    int assignmeshsize;

    // Used by Delaunay refinement.
    int refineflag; // 0, 1, 2, 3
    triface refinetet;
    face refinesh;

    insertvertexflags() {
      // All flags are initialized as 0.
      iloc = bowywat = lawson = 0;
      rejflag = chkencflag = 0;
      sloc = sbowywat = 0;
      splitbdflag = validflag = respectbdflag = 0;
      assignmeshsize = 0;

      refineflag = 0;
      refinetet.tet = NULL;
      refinesh.sh = NULL;
    }
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// flipconstraints                                                           //
//                                                                           //
// A structure of a collection of data (options and parameters) which pass   //
// to the edge flip function flipnm().                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class flipconstraints {

  public:

    point seg[2];  // A constraining edge to be recovered.
    point fac[3];  // A constraining face to be recovered.
    point remvert; // A vertex to be removed.

    // Control flags
    int unflip;  // Undo the performed flips.
    int collectnewtets; // Collect the new tets created by flips.
    int collectencsegflag;

    // Optimization flags.
    int remove_ndelaunay_edge; // Remove a non-Delaunay edge.
    REAL bak_tetprism_vol; // The value to be minimized.
    int remove_large_angle; // Remove a large dihedral angle at edge.
    REAL cosdihed_in; // The input cosine of the dihedral angle (> 0).
    REAL cosdihed_out; // The improved cosine of the dihedral angle.

    // Internal counters.
    int maxflippedlinklevelcount; // Maximal flipped link levels.
    int misfliplinklevelcount; // Number of missed flip possibilities.
    int chrismastreecount; // Number of Chrismas trees (unflippable case).
    int convexhulledgecount; // Number of convex hull edges (unflippable case).
    int encsegcount; // Number of hitted segments. 
    int rejf23count, rejf32count; // Number of rejections by checkflipeligi..

    void clearcounters() {
      maxflippedlinklevelcount = 0;
      misfliplinklevelcount = 0;
      chrismastreecount = 0;
      convexhulledgecount = 0;
      encsegcount = 0;
      rejf23count = rejf32count = 0;
    }

    flipconstraints() {
      seg[0] = NULL;
      fac[0] = NULL;
      remvert = NULL;

      unflip = 0;
      collectnewtets = 0;
      collectencsegflag = 0;

      remove_ndelaunay_edge = 0;
      bak_tetprism_vol = 0.0;
      remove_large_angle = 0;
      cosdihed_in = 0.0;
      cosdihed_out = 0.0;

      clearcounters();
    }
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// optparameters                                                             //
//                                                                           //
// Optimization options and parameters.                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class optparameters {

  public:

    // The one of goals of optimization.
    int max_min_volume;      // Maximize the minimum volume.
    int max_min_aspectratio; // Maximize the minimum aspect ratio.
    int min_max_dihedangle;  // Minimize the maxumum dihedral angle. 

    // The initial and improved value.
    REAL initval, imprval;

    int numofsearchdirs;
    REAL searchstep;
    int maxiter;  // Maximum smoothing iterations (disabled by -1).
    int smthiter; // Performed iterations.

    int expstarflag;
    int expstarcount;

    int flipflag;
    int checkencflag;

    optparameters() {
      max_min_volume = 0;
      max_min_aspectratio = 0;
      min_max_dihedangle = 0;

      initval = imprval = 0.0;

      numofsearchdirs = 10;
      searchstep = 0.01;
      maxiter = -1;   // Unlimited smoothing iterations.
      smthiter = 0;

      expstarflag = 0;
      expstarcount = 0;

      flipflag = 0;
      checkencflag = 0;
    }
  };


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Arraypool                                                                 //
//                                                                           //
// A dynamic linear array.                                                   //
// (It is from Shewchuk's Starbase.c, which is provided as part of Stellar,  //
// a program for improving tetrahedral meshes.)                              //
//                                                                           //
// Each arraypool contains an array of pointers to a number of blocks.  Each //
// block contains the same fixed number of objects.  Each index of the array //
// addesses a particular object in the pool.  The most significant bits add- //
// ress the index of the block containing the object. The less significant   //
// bits address this object within the block.                                //
//                                                                           //
// 'objectbytes' is the size of one object in blocks; 'log2objectsperblock'  //
// is the base-2 logarithm of 'objectsperblock'; 'objects' counts the number //
// of allocated objects; 'totalmemory' is the totoal memorypool in bytes.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class arraypool {

  public:

    int objectbytes;
    int objectsperblock;
    int log2objectsperblock; 
    int toparraylen;
    char **toparray;
    long objects;
    unsigned long totalmemory;

    void restart();
    void poolinit(int sizeofobject, int log2objperblk);
    char* getblock(int objectindex);
    void* lookup(int objectindex);
    int newindex(void **newptr);

    arraypool(int sizeofobject, int log2objperblk);
    ~arraypool();
  };

// fastlookup() -- A fast, unsafe operation. Return the pointer to the object
//   with a given index.  Note: The object's block must have been allocated,
//   i.e., by the function newindex().

#define fastlookup(pool, index) \
  (void *) ((pool)->toparray[(index) >> (pool)->log2objectsperblock] + \
            ((index) & ((pool)->objectsperblock - 1)) * (pool)->objectbytes)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Memorypool                                                                //
//                                                                           //
// A type used to allocate memory.                                           //
// (It is from Shewchuk's triangle.c.)                                       //
//                                                                           //
// firstblock is the first block of items. nowblock is the block from which  //
//   items are currently being allocated. nextitem points to the next slab   //
//   of free memory for an item. deaditemstack is the head of a linked list  //
//   (stack) of deallocated items that can be recycled.  unallocateditems is //
//   the number of items that remain to be allocated from nowblock.          //
//                                                                           //
// Traversal is the process of walking through the entire list of items, and //
//   is separate from allocation.  Note that a traversal will visit items on //
//   the "deaditemstack" stack as well as live items.  pathblock points to   //
//   the block currently being traversed.  pathitem points to the next item  //
//   to be traversed.  pathitemsleft is the number of items that remain to   //
//   be traversed in pathblock.                                              //
//                                                                           //
// itemwordtype is set to POINTER or FLOATINGPOINT, and is used to suggest   //
//   what sort of word the record is primarily made up of.  alignbytes       //
//   determines how new records should be aligned in memory.  itembytes and  //
//   itemwords are the length of a record in bytes (after rounding up) and   //
//   words.  itemsperblock is the number of items allocated at once in a     //
//   single block.  items is the number of currently allocated items.        //
//   maxitems is the maximum number of items that have been allocated at     //
//   once; it is the current number of items plus the number of records kept //
//   on deaditemstack.                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class memorypool {

  public:

    // Labels that signify whether a record consists primarily of pointers
    //   or of floating-point words.  Used for data alignment.
    enum wordtype {POINTER, FLOATINGPOINT};

    void **firstblock, **nowblock;
    void *nextitem;
    void *deaditemstack;
    void **pathblock;
    void *pathitem;
    wordtype itemwordtype;
    int  alignbytes;
    int  itembytes, itemwords;
    int  itemsperblock;
    long items, maxitems;
    int  unallocateditems;
    int  pathitemsleft;

    memorypool();
    memorypool(int, int, enum wordtype, int);
    ~memorypool();
    
    void poolinit(int, int, enum wordtype, int);
    void restart();
    void *alloc();
    void dealloc(void*);
    void traversalinit();
    void *traverse();
  };  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class variables                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Pointer to the input data (a set of nodes, a PLC, or a mesh).
  tetgenio *in, *addin;

  // Pointer to the switches and parameters.
  tetgenbehavior *b;

  // Pointer to a background mesh (contains size specification map).
  tetgenmesh *bgm;

  // Memorypools to store mesh elements: tetrahedra, subfaces, segments,
  //   and vertices. And memorypools for storing pointers which connect 
  //   tetrahedra and subfaces and segments.
  memorypool *tetrahedrons, *subfaces, *subsegs, *points;
  memorypool *tet2subpool, *tet2segpool;

  // Memorypools to store bad-quality (or encroached) elements.
  memorypool *badtetrahedrons, *badsubfacs, *badsubsegs;

  // A memorypool to store faces to be flipped.
  memorypool *flippool;
  // A stack of faces to be flipped.
  badface *flipstack;
  // A queue to store unflippable elements.
  arraypool *unflipqueue; 

  // Arrays used for point insertion (the Bowyer-Watson algorithm).
  arraypool *cavetetlist, *cavebdrylist, *caveoldtetlist;
  arraypool *cavetetshlist, *cavetetseglist, *cavetetvertlist;
  arraypool *caveencshlist, *caveencseglist;
  arraypool *caveshlist, *caveshbdlist, *cavesegshlist;

  // Stacks used for CDT construction and boundary recovery.
  arraypool *subsegstack, *subfacstack, *subvertstack;
  arraypool *suppsteinerptlist;

  // The infinite vertex.
  point dummypoint;

  // Two handles used for facet recovery in CDT.
  triface firsttopface, firstbotface;
  // Three points define a plane (used in formcavity()).
  point plane_pa, plane_pb, plane_pc;

  // Two arraies of encroached segments and subfaces (in mesh refinement).
  arraypool *encseglist, *encshlist;

  // Pointer to a recently visited tetrahedron, subface.
  triface recenttet;
  face recentsh;

  // PI is the ratio of a circle's circumference to its diameter.
  static REAL PI;

  // Array (size = numberoftetrahedra * 6) for storing high-order nodes of
  //   tetrahedra (only used when -o2 switch is selected).
  point *highordertable;

  // Other variables.
  REAL xmax, xmin, ymax, ymin, zmax, zmin;         // Bounding box of points.
  REAL longest;                          // The longest possible edge length.
  long hullsize;                           // Number of faces of convex hull.
  long insegments;                               // Number of input segments.
  long meshedges;                             // Number of output mesh edges.
  long meshhulledges;                           // Number of hull mesh edges.
  int steinerleft;                  // Number of Steiner points not yet used.
  int numpointattrib;                          // Number of point attributes.
  int sizeoftensor;                     // Number of REALs per metric tensor.
  int pointmtrindex;           // Index to find the metric tensor of a point.
  int pointparamindex;       // Index to find the u,v coordinates of a point.
  int point2simindex;         // Index to find a simplex adjacent to a point.
  int pointmarkindex;            // Index to find boundary marker of a point.
  int numelemattrib;                     // Number of tetrahedron attributes.
  int elemattribindex;          // Index to find attributes of a tetrahedron.
  int volumeboundindex;       // Index to find volume bound of a tetrahedron.
  int elemmarkerindex;              // Index to find marker of a tetrahedron.
  int shmarkindex;             // Index to find boundary marker of a subface.
  int areaboundindex;               // Index to find area bound of a subface.
  int checksubsegflag;   // Are there segments in the tetrahedralization yet?
  int checksubfaceflag;  // Are there subfaces in the tetrahedralization yet?
  int checkinverttetflag;       // Are there inverted (degenerated) tets yet?
  int checkconstraints;  // Are there variant (node, seg, facet) constraints?
  int nonconvex;                               // Is current mesh non-convex?
  int dupverts;                             // Are there duplicated vertices?
  int unuverts;                                 // Are there unused vertices?
  long samples;               // Number of random samples for point location.
  unsigned long randomseed;                    // Current random number seed.
  REAL cosmaxdihed, cosmindihed;    // The cosine values of max/min dihedral.
  REAL cossmtdihed;      // The cosine value of a bad dihedral tobe smoothed.
  REAL cosslidihed;      // The cosine value of the max dihedral of a sliver.
  REAL minfaceang, minfacetdihed;     // The minimum input (dihedral) angles.
  REAL sintheta_tol;                   // The tolerance for sin(small angle).
  int autofliplinklevel;    // The increasement of link levels, default is 1.
  int calc_tetprism_vol;   // Flag to calculate the tetrahedral-prism'volume.
  REAL tetprism_vol_sum;   // The total volume of tetrahedral-prisms (in 4D).

  // Algorithm statistical counters.
  int  max_hcurve_depth_count;
  long ptloc_count, ptloc_max_count;
  long insphere_sos_count, orient4d_sos_count;
  long flip14count, flip26count, flipn2ncount;
  long flip23count, flip32count, flip44count, flip22count;
  long maxbowatcavsize, totalbowatcavsize, totaldeadtets;
  long triedgcount, triedgcopcount;
  long across_face_count, across_edge_count, across_max_count;
  long fillregioncount, missingsubfacecount, crossingtetcount;
  long cavitycount, cavityexpcount, maxcavsize, maxregionsize;
  long maxcrossfacecount, maxflipsequence;
  long dbg_ignore_facecount, dbg_unflip_facecount;
  long ccent_relocate_count;
  long opt_sliver_peels;
  long r1count, r2count, r3count; 
  long maxfliplinklevel, maxflipstarsize;
  long flipstarcount, sucflipstarcount, skpflipstarcount;
  long st_segref_count, st_facref_count, st_volref_count; 
  long nonregularcount;
  long rejrefinetetcount, rejrefineshcount;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh manipulation primitives                                              //
//                                                                           //
// Mesh manipulation primitives are indeed very simple functions which take  //
// one or two handles as parameters,  perform basic operations such as "glue //
// two tetrahedra at a face",  "return the origin of a tetrahedron", "return //
// the subface adjoining at the face of a tetrahedron", and so on.           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Fast lookup tables for mesh manipulation primitives.
  static int mod12[36], mod6[18]; 
  static int orgpivot[12], destpivot[12], apexpivot[12], oppopivot[12];
  static int edgepivot[12], ver2edge[12], edge2ver[6];
  static int sorgpivot [6], sdestpivot[6], sapexpivot[6];
  static int snextpivot[6], epivot[4];

  // Primitives for tetrahedra.
  inline void decode(tetrahedron ptr, triface& t);
  inline tetrahedron encode(triface& t);
  inline tetrahedron encode2(tetrahedron* ptr, int ver);
  inline void bond(triface& t1, triface& t2);
  inline void dissolve(triface& t);
  inline void fsym(triface& t1, triface& t2);
  inline void fsymself(triface& t);
  inline void esym(triface& t1, triface& t2);
  inline void esymself(triface& t);
  inline void enext(triface& t1, triface& t2);
  inline void enextself(triface& t);
  inline void eprev(triface& t1, triface& t2);
  inline void eprevself(triface& t);
  inline void enextesym(triface& t1, triface& t2);
  inline void enextesymself(triface& t);
  inline void eprevesym(triface& t1, triface& t2);
  inline void eprevesymself(triface& t);
  inline void fnext(triface& t1, triface& t2);
  inline void fnextself(triface& t);
  inline point org (triface& t);
  inline point dest(triface& t);
  inline point apex(triface& t);
  inline point oppo(triface& t);
  inline void setorg (triface& t, point p);
  inline void setdest(triface& t, point p);
  inline void setapex(triface& t, point p);
  inline void setoppo(triface& t, point p);
  inline REAL elemattribute(tetrahedron* ptr, int attnum);
  inline void setelemattribute(tetrahedron* ptr, int attnum, REAL value);
  inline REAL volumebound(tetrahedron* ptr);
  inline void setvolumebound(tetrahedron* ptr, REAL value);
  inline int  elemindex(tetrahedron* ptr);
  inline void setelemindex(tetrahedron* ptr, int value);
  inline int  elemmarker(tetrahedron* ptr);
  inline void setelemmarker(tetrahedron* ptr, int value);
  inline void infect(triface& t);
  inline void uninfect(triface& t);
  inline bool infected(triface& t);
  inline void marktest(triface& t);
  inline void unmarktest(triface& t);
  inline bool marktested(triface& t);
  inline void markface(triface& t);
  inline void unmarkface(triface& t);
  inline bool facemarked(triface& t);
  inline void markedge(triface& t);
  inline void unmarkedge(triface& t);
  inline bool edgemarked(triface& t);
  inline void marktest2(triface& t);
  inline void unmarktest2(triface& t);
  inline bool marktest2ed(triface& t);
  inline int  elemcounter(triface& t);
  inline void setelemcounter(triface& t, int value);
  inline void increaseelemcounter(triface& t);
  inline void decreaseelemcounter(triface& t);
  inline bool ishulltet(triface& t);
  inline bool isdeadtet(triface& t);
 
  // Primitives for subfaces and subsegments.
  inline void sdecode(shellface sptr, face& s);
  inline shellface sencode(face& s);
  inline shellface sencode2(shellface *sh, int shver);
  inline void spivot(face& s1, face& s2);
  inline void spivotself(face& s);
  inline void sbond(face& s1, face& s2);
  inline void sbond1(face& s1, face& s2);
  inline void sdissolve(face& s);
  inline point sorg(face& s);
  inline point sdest(face& s);
  inline point sapex(face& s);
  inline void setsorg(face& s, point pointptr);
  inline void setsdest(face& s, point pointptr);
  inline void setsapex(face& s, point pointptr);
  inline void sesym(face& s1, face& s2);
  inline void sesymself(face& s);
  inline void senext(face& s1, face& s2);
  inline void senextself(face& s);
  inline void senext2(face& s1, face& s2);
  inline void senext2self(face& s);
  inline void sfnext(face& s1, face& s2);
  inline void sfnextself(face& s);
  inline REAL areabound(face& s);
  inline void setareabound(face& s, REAL value);
  inline int shellmark(face& s);
  inline void setshellmark(face& s, int value);
  inline enum shestype shelltype(face& s);
  inline void setshelltype(face& s, enum shestype value); 
  inline void sinfect(face& s);
  inline void suninfect(face& s);
  inline bool sinfected(face& s);
  inline void smarktest(face& s);
  inline void sunmarktest(face& s);
  inline bool smarktested(face& s);
  inline void smarktest2(face& s);
  inline void sunmarktest2(face& s);
  inline bool smarktest2ed(face& s);
  inline void smarktest3(face& s);
  inline void sunmarktest3(face& s);
  inline bool smarktest3ed(face& s);

  // Primitives for interacting tetrahedra and subfaces.
  inline void tsbond(triface& t, face& s);
  inline void tsdissolve(triface& t);
  inline void stdissolve(face& s);
  inline void tspivot(triface& t, face& s);
  inline void stpivot(face& s, triface& t);

  // Primitives for interacting tetrahedra and segments.
  inline void tssbond1(triface& t, face& seg);
  inline void sstbond1(face& s, triface& t);
  inline void tssdissolve1(triface& t);
  inline void sstdissolve1(face& s);
  inline void tsspivot1(triface& t, face& s);
  inline void sstpivot1(face& s, triface& t);

  // Primitives for interacting subfaces and segments.
  inline void ssbond(face& s, face& edge);
  inline void ssbond1(face& s, face& edge);
  inline void ssdissolve(face& s);
  inline void sspivot(face& s, face& edge);

  // Primitives for points.
  inline int  pointmark(point pt);
  inline void setpointmark(point pt, int value);
  inline enum verttype pointtype(point pt);
  inline void setpointtype(point pt, enum verttype value);
  inline int  pointgeomtag(point pt);
  inline void setpointgeomtag(point pt, int value);
  inline REAL pointgeomuv(point pt, int i);
  inline void setpointgeomuv(point pt, int i, REAL value);
  inline void pinfect(point pt);
  inline void puninfect(point pt);
  inline bool pinfected(point pt);
  inline void pmarktest(point pt);
  inline void punmarktest(point pt);
  inline bool pmarktested(point pt);
  inline void pmarktest2(point pt);
  inline void punmarktest2(point pt);
  inline bool pmarktest2ed(point pt);
  inline void pmarktest3(point pt);
  inline void punmarktest3(point pt);
  inline bool pmarktest3ed(point pt);
  inline tetrahedron point2tet(point pt);
  inline void setpoint2tet(point pt, tetrahedron value);
  inline shellface point2sh(point pt);
  inline void setpoint2sh(point pt, shellface value);
  inline point point2ppt(point pt);
  inline void setpoint2ppt(point pt, point value);
  inline tetrahedron point2bgmtet(point pt);
  inline void setpoint2bgmtet(point pt, tetrahedron value);

  // Advanced primitives.
  inline void point2tetorg(point pt, triface& t);
  inline void point2shorg(point pa, face& s);
  inline point farsorg(face& seg);
  inline point farsdest(face& seg);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Memory managment                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void tetrahedrondealloc(tetrahedron*);
  tetrahedron *tetrahedrontraverse();
  tetrahedron *alltetrahedrontraverse();
  void shellfacedealloc(memorypool*, shellface*);
  shellface *shellfacetraverse(memorypool*);
  void badfacedealloc(memorypool*, badface*);
  badface *badfacetraverse(memorypool*);
  void pointdealloc(point);
  point pointtraverse();

  void makeindex2pointmap(point*&);
  void makepoint2submap(memorypool*, int*&, face*&);
  void maketetrahedron(triface*);
  void makeshellface(memorypool*, face*);
  void makepoint(point*, enum verttype);

  void initializepools();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Geometric predicates and calculations                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Triangle-edge intersection test
  int tri_edge_2d(point, point, point, point, point, point, int, int*, int*);
  int tri_edge_tail(point, point, point, point, point, point, REAL, REAL, int, 
                    int*, int*);
  int tri_edge_test(point, point, point, point, point, point, int, int*, int*);

  // Triangle-triangle intersection test
  int tri_edge_inter_tail(point, point, point, point, point, REAL, REAL);
  int tri_tri_inter(point, point, point, point, point, point);

  // Linear algebra functions
  inline REAL dot(REAL* v1, REAL* v2);
  inline void cross(REAL* v1, REAL* v2, REAL* n);
  bool lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N);
  void lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N);

  // Geometric predicates
  REAL incircle3d(point pa, point pb, point pc, point pd);
  REAL insphere_s(REAL*, REAL*, REAL*, REAL*, REAL*);
  REAL orient4d_s(REAL*, REAL*, REAL*, REAL*, REAL*, 
                  REAL, REAL, REAL, REAL, REAL);

  // Geometric calculations
  inline REAL distance(REAL* p1, REAL* p2);
  void facenormal(point pa, point pb, point pc, REAL *n, int pivot, REAL *lav);
  REAL shortdistance(REAL* p, REAL* e1, REAL* e2);
  REAL triarea(REAL* pa, REAL* pb, REAL* pc);
  REAL interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n);
  void projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj);
  void projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj);
  REAL facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2);
  bool tetalldihedral(point, point, point, point, REAL*, REAL*, REAL*);
  void tetallnormal(point, point, point, point, REAL N[4][3], REAL* volume);
  REAL tetaspectratio(point, point, point, point);
  bool circumsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
  bool orthosphere(REAL*,REAL*,REAL*,REAL*,REAL,REAL,REAL,REAL,REAL*,REAL*);
  void planelineint(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);
  REAL tetprismvol(REAL* pa, REAL* pb, REAL* pc, REAL* pd);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Local mesh transformations                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void flippush(badface*&, triface*);

  // The elementary flips.
  void flip23(triface*, int, int, int);
  void flip32(triface*, int, int, int);
  void flip41(triface*, int, int, int);

  // A generalized edge flip.
  int flipnm(triface*, int n, int level, int, flipconstraints* fc);
  int flipnm_post(triface*, int n, int nn, int, flipconstraints* fc);

  // Incremental flips.
  long lawsonflip3d(point, int flipflag, int, int, int flipedgeflag);

  // Point insertion.
  int insertvertex(point newpt, triface *searchtet, face *splitsh, face*,
                   insertvertexflags *ivf);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Delaunay tetrahedralization                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void transfernodes();

  // Point sorting.
  int  transgc[8][3][8], tsb1mod3[8];
  void hilbert_init(int n);
  int  hilbert_split(point* vertexarray, int arraysize, int gc0, int gc1,
                     REAL, REAL, REAL, REAL, REAL, REAL);
  void hilbert_sort3(point* vertexarray, int arraysize, int e, int d,
                     REAL, REAL, REAL, REAL, REAL, REAL, int depth);

  // Point location.
  unsigned long randomnation(unsigned int choices);
  void randomsample(point searchpt, triface *searchtet);
  enum locateresult locate(point searchpt, triface *searchtet, int);

  // Incremental Delaunay construction.
  void initialdelaunay(point pa, point pb, point pc, point pd);
  void incrementaldelaunay(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Surface triangulation                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  bool calculateabovepoint(arraypool*, point*, point*, point*);
  void calculateabovepoint4(point, point, point, point);

  void flipshpush(face*);
  void flip22(face*, int, int);
  void flip31(face*, int);
  long lawsonflip();
  int sinsertvertex(point newpt, face*, face*, int iloc, int bowywat);
  int sremovevertex(point delpt, face*, face*, int lawson);

  enum locateresult slocate(point, face*, int, int, int);
  enum interresult sscoutsegment(face*, point);
  void scarveholes(int, REAL*);
  void triangulate(int, arraypool*, arraypool*, int, REAL*);

  void unifysubfaces(face*, face*);
  void unifysegments();
  void mergefacets();
  void identifypscedges(point*);
  void meshsurface();

  void interecursive(shellface** subfacearray, int arraysize, int axis,
                     REAL, REAL, REAL, REAL, REAL, REAL, int* internum);
  void detectinterfaces();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constrained Delaunay tetrahedralization                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void markacutevertices();

  void reportselfintersect(face *seg, face *shface);

  enum interresult finddirection(triface* searchtet, point endpt);
  enum interresult scoutsegment(point, point, triface*, point*, arraypool*);
  void getsteinerptonsegment(face* seg, point refpt, point steinpt);
  void delaunizesegments();

  enum interresult scoutsubface(face* searchsh, triface* searchtet);
  void formmissingregion(face*, arraypool*, arraypool*, arraypool*, arraypool*);
  int scoutcrossedge(triface& crosstet, arraypool*, arraypool*);
  bool formcavity(triface*, arraypool*, arraypool*, arraypool*, arraypool*, 
                  arraypool*, arraypool*);

  // Facet recovery by local re-tetrahedralization [Si and Gaertner'05,'11].
  void delaunizecavity(arraypool*, arraypool*, arraypool*, arraypool*, 
                       arraypool*, arraypool*);
  bool fillcavity(arraypool*, arraypool*, arraypool*, arraypool*);
  void carvecavity(arraypool*, arraypool*, arraypool*);
  void restorecavity(arraypool*, arraypool*, arraypool*);

  // Facet recovery by flips [Shewchuk'03].
  void flipcertify(triface *chkface, badface **pqueue);
  void flipinsertfacet(arraypool*, arraypool*, arraypool*, arraypool*);

  bool fillregion(arraypool* missingshs, arraypool*, arraypool* newshs);
  void refineregion();

  void constrainedfacets();  

  void constraineddelaunay(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constrained tetrahedralizations.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  int checkflipeligibility(int fliptype, point, point, point, point, point,
                           int level, int edgepivot, flipconstraints* fc);

  int removeedgebyflips(triface*, flipconstraints*);
  int removefacebyflips(triface*, flipconstraints*);

  int recoveredgebyflips(point, point, triface*, int fullsearch);
  int add_steinerpt_in_schoenhardtpoly(triface*, int, int chkencflag);
  int addsteiner4recoversegment(face*, int);
  int recoversegments(arraypool*, int fullsearch, int steinerflag);

  int recoverfacebyflips(point, point, point, face*, triface*);
  int recoversubfaces(arraypool*, int steinerflag);

  int getvertexstar(int, point searchpt, arraypool*, arraypool*, arraypool*);
  int getedge(point, point, triface*);
  int reduceedgesatvertex(point startpt, arraypool* endptlist);
  int removevertexbyflips(point steinerpt);

  int suppressssteinerpoint(point steinerpt);
  int suppresssteinerpoints();

  void recoverboundary(clock_t&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh reconstruction                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void carveholes();

  void reconstructmesh();

  int  scoutpoint(point, triface*, int randflag);
  REAL getpointmeshsize(point, triface*, int iloc);
  void interpolatemeshsize();

  void insertconstrainedpoints(tetgenio *addio);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh refinement                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void marksharpsegments();
  void decidefeaturepointsizes();

  int checkseg4encroach(point pa, point pb, point checkpt);
  int checkseg4split(face *chkseg, point&, int&);
  int splitsegment(face *splitseg, point encpt, int qflag, int chkencflag);
  void repairencsegs(int chkencflag);

  int checkfac4encroach(point, point, point, point checkpt, REAL*, REAL*);
  int checkfac4split(face *chkfac, point& encpt, int& qflag, REAL *ccent);
  int splitsubface(face *splitfac, point encpt, int qflag, REAL *ccent,
                   int chkencflag);
  void repairencfacs(int chkencflag);

  int checktet4split(triface *chktet, int& qflag, REAL *ccent);
  int splittetrahedron(triface* splittet,int qflag,REAL *ccent,int chkencflag);
  void repairbadtets(int chkencflag);
  void insertsinks();

  void delaunayrefinement();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh optimization                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void recoverdelaunay();

  int  gettetrahedron(point, point, point, point, triface *);
  long improvequalitybyflips();

  int  smoothpoint(point smtpt, arraypool*, int ccw, optparameters *opm);
  long improvequalitybysmoothing(optparameters *opm);

  int  splitsliver(triface *, REAL, int);
  long removeslivers(int);

  void optimizemesh();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh check and statistics                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Mesh validations.
  int checkmesh(int topoflag);
  int checkshells();
  int checksegments();
  int checkdelaunay();
  int checkregular(int);
  int checkconforming(int);

  //  Mesh statistics.
  void qualitystatistics();
  void statistics();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh output                                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void jettisonnodes();
  void highorder();
  void numberedges();
  void outnodes(tetgenio*);
  void outmetrics(tetgenio*);
  void outelements(tetgenio*);
  void outfaces(tetgenio*);
  void outhullfaces(tetgenio*);
  void outsubfaces(tetgenio*);
  void outedges(tetgenio*);
  void outsubsegments(tetgenio*);
  void outneighbors(tetgenio*);
  void outvoronoi(tetgenio*);
  void outsmesh(char*);
  void outmesh2medit(char*);
  void outmesh2vtk(char*);


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constructor & destructor                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  tetgenmesh()
  {
    in  = addin = NULL;
    b   = NULL;
    bgm = NULL;

    tetrahedrons = subfaces = subsegs = points = NULL;
    badtetrahedrons = badsubfacs = badsubsegs = NULL;
    tet2segpool = tet2subpool = NULL;
    flippool = NULL;

    dummypoint = NULL;
    flipstack = NULL;
    unflipqueue = NULL;

    cavetetlist = cavebdrylist = caveoldtetlist = NULL;
    cavetetshlist = cavetetseglist = cavetetvertlist = NULL;
    caveencshlist = caveencseglist = NULL;
    caveshlist = caveshbdlist = cavesegshlist = NULL;

    subsegstack = subfacstack = subvertstack = NULL;
    suppsteinerptlist = NULL;
    encseglist = encshlist = NULL;

    highordertable = NULL;

    plane_pa = plane_pb = plane_pc = (point) NULL;

    xmax = xmin = ymax = ymin = zmax = zmin = 0.0; 
    longest = 0.0;
    hullsize = 0l;
    insegments = 0l;
    meshedges = meshhulledges = 0l;
    steinerleft = -1;
    numpointattrib = 0;
    sizeoftensor = 0;
    pointmtrindex = 0;
    pointparamindex = 0;
    pointmarkindex = 0;
    point2simindex = 0;
    numelemattrib = 0;
    elemattribindex = 0;
    volumeboundindex = 0;
    shmarkindex = 0;
    areaboundindex = 0;
    checksubsegflag = 0;
    checksubfaceflag = 0;
    checkinverttetflag = 0;
    checkconstraints = 0;
    nonconvex = 0;
    dupverts = 0;
    unuverts = 0;
    samples = 0l;
    randomseed = 1l;
    minfaceang = minfacetdihed = PI;
    sintheta_tol = sin(0.001 * PI / 180.0);
    autofliplinklevel = 1;
    calc_tetprism_vol = 0;
    tetprism_vol_sum = 0.0;

    ptloc_count = ptloc_max_count = 0l;
    insphere_sos_count = orient4d_sos_count = 0l;
    flip14count = flip26count = flipn2ncount = 0l;
    flip23count = flip32count = flip44count = flip22count = 0l;
    maxbowatcavsize = totalbowatcavsize = totaldeadtets = 0l;
    triedgcount = triedgcopcount = 0l;
    across_face_count = across_edge_count = across_max_count = 0l;
    fillregioncount = missingsubfacecount = crossingtetcount = 0l;
    cavitycount = cavityexpcount = 0l;
    maxcavsize = maxregionsize = 0l;
    maxcrossfacecount = maxflipsequence = 0l;
    dbg_ignore_facecount = dbg_unflip_facecount = 0l;
    ccent_relocate_count = 0l;
    opt_sliver_peels = 0l;
    r1count = r2count = r3count = 0l;
    st_segref_count = st_facref_count = st_volref_count = 0l;
    nonregularcount = 0l;

    maxfliplinklevel = maxflipstarsize = 0l;
    flipstarcount = sucflipstarcount = skpflipstarcount = 0l;

    rejrefinetetcount = rejrefineshcount = 0l;
  } // tetgenmesh()

  ~tetgenmesh()
  {
    if (bgm != NULL) {
      delete bgm;
    }

    if (tetrahedrons != (memorypool *) NULL) {
      delete tetrahedrons;
    }
    if (subfaces != (memorypool *) NULL) {
      delete subfaces;
    }
    if (subsegs != (memorypool *) NULL) {
      delete subsegs;
    }
    if (points != (memorypool *) NULL) {
      delete points;
    }
    if (tet2segpool != NULL) {
      delete tet2segpool;
    }
    if (tet2subpool != NULL) {
      delete tet2subpool;
    }
    if (flippool != NULL) {
      delete flippool;
      delete unflipqueue;
    }
    if (dummypoint != (point) NULL) {
      delete [] dummypoint;
    }

    if (cavetetlist != NULL) {
      delete cavetetlist;
      delete cavebdrylist;
      delete caveoldtetlist;
      delete cavetetvertlist;
    }

    if (caveshlist != NULL) {
      delete caveshlist;
      delete caveshbdlist;
      delete cavesegshlist;
      delete cavetetshlist;
      delete cavetetseglist;
      delete caveencshlist;
      delete caveencseglist;
    }

    if (subsegstack != NULL) {
      delete subsegstack;
      delete subfacstack;
      delete subvertstack;
    }

    if (suppsteinerptlist != NULL) {
      delete suppsteinerptlist;
    }

    if (highordertable != NULL) {
      delete [] highordertable;
    }
  } // ~tetgenmesh()

};                                               // End of class tetgenmesh.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    Interface for using TetGen's library to generate      //
//                     Delaunay tetrahedralizations, constrained Delaunay    //
//                     tetrahedralizations, quality tetrahedral meshes.      //
//                                                                           //
// 'in' is an object of 'tetgenio' which contains a PLC you want to tetrahed-//
// ralize or a previously generated tetrahedral mesh you want to refine.  It //
// must not be a NULL. 'out' is another object of 'tetgenio' for storing the //
// generated tetrahedral mesh. It can be a NULL. If so, the output will be   //
// saved to file(s). If 'bgmin' != NULL, it contains a background mesh which //
// defines a mesh size distruction function.                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(tetgenbehavior *b, tetgenio *in, tetgenio *out, 
                    tetgenio *addin = NULL, tetgenio *bgmin = NULL);

#ifdef TETLIBRARY
void tetrahedralize(char *switches, tetgenio *in, tetgenio *out,
                    tetgenio *addin = NULL, tetgenio *bgmin = NULL);
#endif // #ifdef TETLIBRARY

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// terminatetetgen()    Terminate TetGen with a given exit code.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

inline void terminatetetgen(int x)
{
#ifdef TETLIBRARY
  throw x;
#else
  switch (x) {
  case 1: // Out of memory.
    printf("Error:  Out of memory.\n"); 
    break;
  case 2: // Encounter an internal error.
    printf("Please report this bug to Hang.Si@wias-berlin.de. Include\n");
    printf("  the message above, your input data set, and the exact\n");
    printf("  command line you used to run this program, thank you.\n");
    break;
  case 3:
    printf("A self-intersection was detected. Program stopped.\n");
    printf("Hint: use -d option to detect all self-intersections.\n"); 
    break;
  case 4:
    printf("A very small input feature was size detected. Program stopped.\n");
    printf("Hint: use -T option to set a smaller tolerance.\n");
    break;
  case 5:
    printf("Two very clsoe input facets were detected. Program stopped.\n");
    printf("Hint: use -Y option to avoid adding Steiner points in boundary.\n");
    break;
  case 10: 
    printf("An input error was detected Program stopped.\n"); 
    break;
  } // switch (x)
  exit(x);
#endif // #ifdef TETLIBRARY
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Inline functions of mesh data structures                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//
// Begin of primitives for tetrahedra
// 

// decode()  converts a pointer to an ordered tetrahedron. The version is
//   extracted from the four least significant bits of the pointer.

inline void tetgenmesh::decode(tetrahedron ptr, triface& t) {
  (t).ver = (int) ((uintptr_t) (ptr) & (uintptr_t) 15);
  (t).tet = (tetrahedron *) ((uintptr_t) (ptr) ^ (uintptr_t) (t).ver);
}

// encode()  compress an ordered tetrahedron into a single pointer.  It
//   relies on the assumption that all tetrahedra are aligned to sixteen-
//   byte boundaries, so that the last four significant bits are zero.

inline tetgenmesh::tetrahedron tetgenmesh::encode(triface& t) {
  return (tetrahedron) ((uintptr_t) (t).tet | (uintptr_t) (t).ver);
}

inline tetgenmesh::tetrahedron tetgenmesh::encode2(tetrahedron* ptr, int ver) {
  return (tetrahedron) ((uintptr_t) (ptr) | (uintptr_t) (ver));
}

// bond()  connects two adjacent tetrahedra together. t1 and t2 must refer
//   to the same face and the same edge. Note that the edge directions of
//   t1 and t2 are reversed. 
// Since an edge of t1 can be bonded to any of the three edges of t2. We 
//   choose to bond the edge of t2 which is symmetric to the 0-th edge of
//   t1, and vice versa. Now assume t1 is at i-th edge and t2 is at j-th
//   edge, where i, j in {0, 1, 2}. The edge in t2 symmetric to 0-th edge
//   of t1 is (i + j) modulo 3, and vice versa.  Since the edge number is 
//   coded in the two higher bits of the version, i.e., i, j in {0, 4, 8}.
//   Therefore the edge in t2 symmetric to 0-th edge of t1 becomes
//   (i + j) modulo 12, and vice versa.
/*
inline void tetgenmesh::bond(triface& t1, triface& t2) {
  (t1).tet[(t1).ver & 3] = encode2((t2).tet,
    ((t2).ver & 3) + mod12[((t1).ver & 12) + ((t2).ver & 12)]);
  (t2).tet[(t2).ver & 3] = encode2((t1).tet,
    ((t1).ver & 3) + mod12[((t1).ver & 12) + ((t2).ver & 12)]);
}
*/
// Comment:  The following code seems faster than the above code when
//   it is compiled with the optimization option, e.g., -O3.
inline void tetgenmesh::bond(triface& t1, triface& t2) {
  (t1).tet[(t1).ver & 3] = encode2((t2).tet,
    ((t2).ver & 3) + (((t1).ver & 12) + ((t2).ver & 12)) % 12);
  (t2).tet[(t2).ver & 3] = encode2((t1).tet,
    ((t1).ver & 3) + (((t1).ver & 12) + ((t2).ver & 12)) % 12);
}

// dissolve()  a bond (from one side).

inline void tetgenmesh::dissolve(triface& t) {
  t.tet[t.ver & 3] = NULL;
}

// fsym()  finds the adjacent tetrahedron at the same face and the same edge.

inline void tetgenmesh::fsym(triface& t1, triface& t2) {
  decode((t1).tet[(t1).ver & 3], t2);
  (t2).ver = mod12[(t2).ver + 12 - ((t1).ver & 12)];
}

inline void tetgenmesh::fsymself(triface& t) {
  int offset = 12 - ((t).ver & 12);
  decode((t).tet[(t).ver & 3], t);
  (t).ver = mod12[(t).ver + offset];
}

// enext()  finds the next edge (counterclockwise) in the same face.

inline void tetgenmesh::enext(triface& t1, triface& t2) {
  (t2).tet = (t1).tet;
  (t2).ver = mod12[(t1).ver + 4];
}

inline void tetgenmesh::enextself(triface& t) {
  (t).ver = mod12[(t).ver + 4];
}

// eprev()   finds the next edge (clockwise) in the same face.

inline void tetgenmesh::eprev(triface& t1, triface& t2) {
  (t2).tet = (t1).tet;
  (t2).ver = mod12[(t1).ver + 8];
}

inline void tetgenmesh::eprevself(triface& t) {
  (t).ver = mod12[(t).ver + 8];
}

// esym()  finds the reversed edge.  It is in the other face of the
//   same tetrahedron.

inline void tetgenmesh::esym(triface& t1, triface& t2) {
  (t2).tet = (t1).tet;
  (t2).ver = edgepivot[(t1).ver];
}

inline void tetgenmesh::esymself(triface& t) {
  (t).ver = edgepivot[(t).ver];
}

// enextesym()  finds the reversed edge of the next edge. It is in the other
//   face of the same tetrahedron. It is the combination esym() * enext(). 

inline void tetgenmesh::enextesym(triface& t1, triface& t2) {
  enext(t1, t2);
  esymself(t2);
}

inline void tetgenmesh::enextesymself(triface& t) {
  enextself(t);
  esymself(t);
}

// eprevesym()  finds the reversed edge of the previous edge.

inline void tetgenmesh::eprevesym(triface& t1, triface& t2) {
  eprev(t1, t2);
  esymself(t2);
}

inline void tetgenmesh::eprevesymself(triface& t) {
  eprevself(t);
  esymself(t);
}

// fnext()  finds the next face while rotating about an edge according to
//   a right-hand rule. The face is in the adjacent tetrahedron.  It is
//   the combination: fsym() * esym().

inline void tetgenmesh::fnext(triface& t1, triface& t2) {
  esym(t1, t2);
  fsymself(t2);
}

inline void tetgenmesh::fnextself(triface& t) {
  esymself(t);
  fsymself(t);
}


// The following primtives get or set the origin, destination, face apex,
//   or face opposite of an ordered tetrahedron.

inline tetgenmesh::point tetgenmesh::org(triface& t) {
  return (point) (t).tet[orgpivot[(t).ver]];
}

inline tetgenmesh::point tetgenmesh:: dest(triface& t) {
  return (point) (t).tet[destpivot[(t).ver]];
}

inline tetgenmesh::point tetgenmesh:: apex(triface& t) {
  return (point) (t).tet[apexpivot[(t).ver]];
}

inline tetgenmesh::point tetgenmesh:: oppo(triface& t) {
  return (point) (t).tet[oppopivot[(t).ver]];
}

inline void tetgenmesh:: setorg(triface& t, point p) {
  (t).tet[orgpivot[(t).ver]] = (tetrahedron) (p);
}

inline void tetgenmesh:: setdest(triface& t, point p) {
  (t).tet[destpivot[(t).ver]] = (tetrahedron) (p);
}

inline void tetgenmesh:: setapex(triface& t, point p) {
  (t).tet[apexpivot[(t).ver]] = (tetrahedron) (p);
}

inline void tetgenmesh:: setoppo(triface& t, point p) {
  (t).tet[oppopivot[(t).ver]] = (tetrahedron) (p);
}

#define setvertices(t, torg, tdest, tapex, toppo) \
  (t).tet[orgpivot[(t).ver]] = (tetrahedron) (torg);\
  (t).tet[destpivot[(t).ver]] = (tetrahedron) (tdest); \
  (t).tet[apexpivot[(t).ver]] = (tetrahedron) (tapex); \
  (t).tet[oppopivot[(t).ver]] = (tetrahedron) (toppo)

// Check or set a tetrahedron's attributes.

inline REAL tetgenmesh::elemattribute(tetrahedron* ptr, int attnum) {
  return ((REAL *) (ptr))[elemattribindex + attnum];
}

inline void tetgenmesh::setelemattribute(tetrahedron* ptr, int attnum, 
  REAL value) {
  ((REAL *) (ptr))[elemattribindex + attnum] = value;
}

// Check or set a tetrahedron's maximum volume bound.

inline REAL tetgenmesh::volumebound(tetrahedron* ptr) {
  return ((REAL *) (ptr))[volumeboundindex];
}

inline void tetgenmesh::setvolumebound(tetrahedron* ptr, REAL value) {
  ((REAL *) (ptr))[volumeboundindex] = value;
}

// Get or set a tetrahedron's index (only used for output).
//    These two routines use the reserved slot ptr[10].

inline int tetgenmesh::elemindex(tetrahedron* ptr) {
  int *iptr = (int *) &(ptr[10]);
  return iptr[0];
}

inline void tetgenmesh::setelemindex(tetrahedron* ptr, int value) {
  int *iptr = (int *) &(ptr[10]);
  iptr[0] = value;
}

// Get or set a tetrahedron's marker. 
//   Set 'value = 0' cleans all the face/edge flags.

inline int tetgenmesh::elemmarker(tetrahedron* ptr) {
  return ((int *) (ptr))[elemmarkerindex];
}

inline void tetgenmesh::setelemmarker(tetrahedron* ptr, int value) {
  ((int *) (ptr))[elemmarkerindex] = value;
}

// infect(), infected(), uninfect() -- primitives to flag or unflag a
//   tetrahedron. The last bit of the element marker is flagged (1)
//   or unflagged (0).

inline void tetgenmesh::infect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= 1;
}

inline void tetgenmesh::uninfect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~1;
}

// Test a tetrahedron for viral infection.

inline bool tetgenmesh::infected(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 1) != 0;
}

// marktest(), marktested(), unmarktest() -- primitives to flag or unflag a
//   tetrahedron.  Use the second lowerest bit of the element marker.

inline void tetgenmesh::marktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= 2;
}

inline void tetgenmesh::unmarktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~2;
}
    
inline bool tetgenmesh::marktested(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 2) != 0;
}

// markface(), unmarkface(), facemarked() -- primitives to flag or unflag a
//   face of a tetrahedron.  From the last 3rd to 6th bits are used for
//   face markers, e.g., the last third bit corresponds to loc = 0. 

inline void tetgenmesh::markface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (4 << (t.ver & 3));
}

inline void tetgenmesh::unmarkface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(4 << (t.ver & 3));
}

inline bool tetgenmesh::facemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (4 << (t.ver & 3))) != 0;
}

// markedge(), unmarkedge(), edgemarked() -- primitives to flag or unflag an
//   edge of a tetrahedron.  From the last 7th to 12th bits are used for
//   edge markers, e.g., the last 7th bit corresponds to the 0th edge, etc. 
//   Remark: The last 7th bit is marked by 2^6 = 64.

inline void tetgenmesh::markedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) (64 << ver2edge[(t).ver]);
}

inline void tetgenmesh::unmarkedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) (64 << ver2edge[(t).ver]);
}

inline bool tetgenmesh::edgemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 
           (int) (64 << ver2edge[(t).ver])) != 0;
}

// marktest2(), unmarktest2(), marktest2ed() -- primitives to flag and unflag
//   a tetrahedron. The 13th bit (2^12 = 4096) is used for this flag.

inline void tetgenmesh::marktest2(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) (4096);
}

inline void tetgenmesh::unmarktest2(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) (4096);
}

inline bool tetgenmesh::marktest2ed(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (int) (4096)) != 0;
}

// elemcounter(), setelemcounter() -- primitives to read or ser a (samll)
//   integer counter in this tet. It is saved from the 16th bit. On 32 bit
//   system, the range of the counter is [0, 2^15 = 32768]. 

inline int tetgenmesh::elemcounter(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex]) >> 16;
}

inline void tetgenmesh::setelemcounter(triface& t, int value) {
  int c = ((int *) (t.tet))[elemmarkerindex];
  // Clear the old counter while keep the other flags.
  c &= 65535; // sum_{i=0^15} 2^i
  c |= (value << 16);
  ((int *) (t.tet))[elemmarkerindex] = c;
}

inline void tetgenmesh::increaseelemcounter(triface& t) {
  int c = elemcounter(t);
  setelemcounter(t, c + 1);
}

inline void tetgenmesh::decreaseelemcounter(triface& t) {
  int c = elemcounter(t);
  assert(c > 0); // Never get a negative counter.
  setelemcounter(t, c - 1);
}

// ishulltet()  tests if t is a hull tetrahedron.

inline bool tetgenmesh::ishulltet(triface& t) {
  return (point) (t).tet[7] == dummypoint;
}

// isdeadtet()  tests if t is a tetrahedron is dead.

inline bool tetgenmesh::isdeadtet(triface& t) {
  return ((t.tet == NULL) || (t.tet[4] == NULL));
}

//
// End of primitives for tetrahedra
//

//
// Begin of primitives for subfaces/subsegments
//

// Each subface contains three pointers to its neighboring subfaces, with
//   edge versions.  To save memory, both information are kept in a single
//   pointer. To make this possible, all subfaces are aligned to eight-byte
//   boundaries, so that the last three bits of each pointer are zeros. An
//   edge version (in the range 0 to 5) is compressed into the last three
//   bits of each pointer by 'sencode()'.  'sdecode()' decodes a pointer,
//   extracting an edge version and a pointer to the beginning of a subface.

inline void tetgenmesh::sdecode(shellface sptr, face& s) {
  s.shver = (int) ((uintptr_t) (sptr) & (uintptr_t) 7);
  s.sh = (shellface *) ((uintptr_t) (sptr) ^ (uintptr_t) (s.shver));
}

inline tetgenmesh::shellface tetgenmesh::sencode(face& s) {
  return (shellface) ((uintptr_t) s.sh | (uintptr_t) s.shver);
}

inline tetgenmesh::shellface tetgenmesh::sencode2(shellface *sh, int shver) {
  return (shellface) ((uintptr_t) sh | (uintptr_t) shver);
}

// sbond() bonds two subfaces (s1) and (s2) together. s1 and s2 must refer
//   to the same edge. No requirement is needed on their orientations.

inline void tetgenmesh::sbond(face& s1, face& s2) 
{
  s1.sh[s1.shver >> 1] = sencode(s2);
  s2.sh[s2.shver >> 1] = sencode(s1);
}

// sbond1() bonds s1 <== s2, i.e., after bonding, s1 is pointing to s2,
//   but s2 is not pointing to s1.  s1 and s2 must refer to the same edge.
//   No requirement is needed on their orientations.

inline void tetgenmesh::sbond1(face& s1, face& s2) 
{
  s1.sh[s1.shver >> 1] = sencode(s2);
}

// Dissolve a subface bond (from one side).  Note that the other subface
//   will still think it's connected to this subface.

inline void tetgenmesh::sdissolve(face& s)
{
  s.sh[s.shver >> 1] = NULL;
}

// spivot() finds the adjacent subface (s2) for a given subface (s1).
//   s1 and s2 share at the same edge.

inline void tetgenmesh::spivot(face& s1, face& s2) 
{
  shellface sptr = s1.sh[s1.shver >> 1];
  sdecode(sptr, s2);
}

inline void tetgenmesh::spivotself(face& s) 
{
  shellface sptr = s.sh[s.shver >> 1];
  sdecode(sptr, s);
}

// These primitives determine or set the origin, destination, or apex
//   of a subface with respect to the edge version.

inline tetgenmesh::point tetgenmesh::sorg(face& s) 
{
  return (point) s.sh[sorgpivot[s.shver]];
}

inline tetgenmesh::point tetgenmesh::sdest(face& s) 
{
  return (point) s.sh[sdestpivot[s.shver]];
}

inline tetgenmesh::point tetgenmesh::sapex(face& s) 
{
  return (point) s.sh[sapexpivot[s.shver]];
}

inline void tetgenmesh::setsorg(face& s, point pointptr) 
{
  s.sh[sorgpivot[s.shver]] = (shellface) pointptr;
}

inline void tetgenmesh::setsdest(face& s, point pointptr) 
{
  s.sh[sdestpivot[s.shver]] = (shellface) pointptr;
}

inline void tetgenmesh::setsapex(face& s, point pointptr) 
{
  s.sh[sapexpivot[s.shver]] = (shellface) pointptr;
}

#define setshvertices(s, pa, pb, pc)\
  setsorg(s, pa);\
  setsdest(s, pb);\
  setsapex(s, pc)

// sesym()  reserves the direction of the lead edge.

inline void tetgenmesh::sesym(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = (s1.shver ^ 1);  // Inverse the last bit.
}

inline void tetgenmesh::sesymself(face& s) 
{
  s.shver ^= 1;
}

// senext()  finds the next edge (counterclockwise) in the same orientaion
//   of this face.

inline void tetgenmesh::senext(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = snextpivot[s1.shver];
}

inline void tetgenmesh::senextself(face& s) 
{
  s.shver = snextpivot[s.shver];
}

inline void tetgenmesh::senext2(face& s1, face& s2) 
{
  s2.sh = s1.sh;
  s2.shver = snextpivot[snextpivot[s1.shver]];
}

inline void tetgenmesh::senext2self(face& s) 
{
  s.shver = snextpivot[snextpivot[s.shver]];
}

// sfnext()  finds the next face (s2) in the same face ring of s1.
//           s2 and s1 have the same edge orientation.
// s2 is found through the following determinations.
//   If the edge of s1 is not a segment, then s2 = spivot(s1).
//   Otherwise, suppose the segment's 0th version is [a,b]. 
//   To find the next face in the face ring we have two cases:
//   (1) s1 is edge [a,b], then s2 = spivot(s1).
//   (2) s1 is edge [b,a], then s1 = spivot(s2).
//   In the case (2), we need to travese in the face ring of [a,b] to
//   get s2.
// Comment: The correctness of this function is guaranteed by the
//   surface mesh data structure, i.e., all subfaces at the face ring
//   of [a,b] have the same edge orientation as [a,b].

inline void tetgenmesh::sfnext(face& s1, face& s2)
{
  face seg, s3;

  spivot(s1, s2);

  if (s2.sh != NULL) {
    sspivot(s1, seg);
    if (seg.sh != NULL) {
      seg.shver = 0;
      if (sorg(s1) != sorg(seg)) {      
        while (1) {
          spivot(s2, s3);
          if (s3.sh == s1.sh) break;
          s2 = s3;
        }
        sesymself(s2);
      }
    } else {
      if (sorg(s2) != sorg(s1)) {
        sesymself(s2);
      }
    }
  }
}

inline void tetgenmesh::sfnextself(face& s)
{
  face seg, s2, s3;

  spivot(s, s2);

  if (s2.sh != NULL) {
    sspivot(s, seg);
    if (seg.sh != NULL) {
      seg.shver = 0;
      if (sorg(s) != sorg(seg)) {      
        while (1) {
          spivot(s2, s3);
          if (s3.sh == s.sh) break;
          s2 = s3;
        }
        sesymself(s2);
      }
    } else {
      if (sorg(s2) != sorg(s)) {
        sesymself(s2);
      }
    }
  }

  s = s2;
}

// Check or set a subface's maximum area bound.

inline REAL tetgenmesh::areabound(face& s) 
{
  return ((REAL *) (s.sh))[areaboundindex];
}

inline void tetgenmesh::setareabound(face& s, REAL value) 
{
  ((REAL *) (s.sh))[areaboundindex] = value;
}

// These two primitives read or set a shell marker.  Shell markers are used
//   to hold user boundary information.

inline int tetgenmesh::shellmark(face& s) 
{
  return ((int *) (s.sh))[shmarkindex];
}

inline void tetgenmesh::setshellmark(face& s, int value) 
{
  ((int *) (s.sh))[shmarkindex] = value;
}


// These two primitives set or read the type of the subface or subsegment.

inline enum tetgenmesh::shestype tetgenmesh::shelltype(face& s) 
{
  return (enum shestype) ((((int *) (s.sh))[shmarkindex + 1]) >> 8);
}

inline void tetgenmesh::setshelltype(face& s, enum shestype value) 
{
  ((int *) (s.sh))[shmarkindex + 1] = ((int) value << 8) +
    ((((int *) ((s).sh))[shmarkindex + 1]) & 255);
}


// sinfect(), sinfected(), suninfect() -- primitives to flag or unflag a
//   subface. The last bit of ((int *) ((s).sh))[shmarkindex+1] is flaged.

inline void tetgenmesh::sinfect(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *) ((s).sh))[shmarkindex+1] | (int) 1);
}

inline void tetgenmesh::suninfect(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *) ((s).sh))[shmarkindex+1] & ~(int) 1);
}

// Test a subface for viral infection.

inline bool tetgenmesh::sinfected(face& s) 
{
  return (((int *) ((s).sh))[shmarkindex+1] & (int) 1) != 0;
}

// smarktest(), smarktested(), sunmarktest() -- primitives to flag or unflag
//   a subface.The last 2nd bit of the integer is flaged.

inline void tetgenmesh::smarktest(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 2);
}

inline void tetgenmesh::sunmarktest(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)2);
}

inline bool tetgenmesh::smarktested(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 2) != 0);
}

// smarktest2(), smarktest2ed(), sunmarktest2() -- primitives to flag or 
//   unflag a subface. The last 3rd bit of the integer is flaged.

inline void tetgenmesh::smarktest2(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 4);
}

inline void tetgenmesh::sunmarktest2(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)4);
}

inline bool tetgenmesh::smarktest2ed(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 4) != 0);
}

// The last 4th bit of ((int *) ((s).sh))[shmarkindex+1] is flaged.

inline void tetgenmesh::smarktest3(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] | (int) 8);
}

inline void tetgenmesh::sunmarktest3(face& s) 
{
  ((int *) ((s).sh))[shmarkindex+1] = 
    (((int *)((s).sh))[shmarkindex+1] & ~(int)8);
}

inline bool tetgenmesh::smarktest3ed(face& s) 
{
  return ((((int *) ((s).sh))[shmarkindex+1] & (int) 8) != 0);
}

//
// End of primitives for subfaces/subsegments
//

//
// Begin of primitives for interacting between tetrahedra and subfaces
//
// tsbond() bond a tetrahedron (t) and a subface (s) together.
// Note that t and s must be the same face and the same edge. Moreover,
//   t and s have the same orientation. 
// Since the edge number in t and in s can be any number in {0,1,2}. We bond
//   the edge in s which corresponds to t's 0th edge, and vice versa.

inline void tetgenmesh::tsbond(triface& t, face& s) 
{
  int soffset, toffset, ver;

  if ((t).tet[9] == NULL) {
    // Allocate space for this tet.
    (t).tet[9] = (tetrahedron) tet2subpool->alloc();
    // NULL all fields in this space.
    for (int i = 0; i < 4; i++) {
      ((shellface *) (t).tet[9])[i] = NULL;
    }
  }

  assert(org(t) == sorg(s)); // FOR DEBUG

  if (((s).shver & 1) == 0) {
    // t and s have the same orientation.
    soffset = mod6[6 - (((t).ver & 12) >> 1)]; // {0,2,4}
    toffset = mod12[12 - (((s).shver & 6) << 1)]; // {0,4,8}
  } else {
    // t and s have revsered orientations.
    soffset = (((t).ver & 12) >> 1); // {0,2,4}
    toffset = (((s).shver & 6) << 1); // {0,4,8}
  }

  // Bond t <== s.
  ver = ((s).shver & 1) + mod6[((s).shver & 6) + soffset];
  ((shellface *) (t).tet[9])[(t).ver & 3] = sencode2((s).sh, ver);
  // Bond s <== t.
  ver = ((t).ver & 3) + mod12[((t).ver & 12) + toffset];
  s.sh[9 + ((s).shver & 1)] = (shellface) encode2((t).tet, ver);
}

// tspivot() finds a subface (s) abutting on the given tetrahdera (t).
//   Return s.sh = NULL if there is no subface at t. Otherwise, return
//   the subface s, and s and t must be at the same edge wth the same
//   orientation.

inline void tetgenmesh::tspivot(triface& t, face& s) 
{
  int soffset;

  if ((t).tet[9] == NULL) {
    (s).sh = NULL;
    return;
  }

  // Get the attached subface s.
  sdecode(((shellface *) (t).tet[9])[(t).ver & 3], (s));

  // Set the right edge in s.
  if (((s).shver & 1) == 0) {
    soffset = (((t).ver & 12) >> 1); // {0,2,4}
  } else {
    soffset = mod6[6 - (((t).ver & 12) >> 1)]; // {0,2,4}
  }
  (s).shver = ((s).shver & 1) + mod6[((s).shver & 6) + soffset];
}

// stpivot() finds a tetrahedron (t) abutting a given subface (s).
//   Return the t (if it exists) with the same edge and the same
//   orientation of s.

inline void tetgenmesh::stpivot(face& s, triface& t) 
{
  int toffset;

  decode((tetrahedron) s.sh[9 + (s.shver & 1)], t);

  if ((t).tet == NULL) {
    return;
  }

  if (((s).shver & 1) == 0) {
    toffset = (((s).shver & 6) << 1); // {0,4,8}
  } else {
    toffset = mod12[12 - (((s).shver & 6) << 1)]; // {0,4,8}
  }
  (t).ver = ((t).ver & 3) + mod12[((t).ver & 12) + toffset];
}

// tsdissolve() dissolve a bond (from the tetrahedron side).

inline void tetgenmesh::tsdissolve(triface& t) 
{
  if ((t).tet[9] != NULL) {
    ((shellface *) (t).tet[9])[(t).ver & 3] = NULL;
  }
}

// stdissolve() dissolve a bond (from the subface side).

inline void tetgenmesh::stdissolve(face& s) 
{
  (s).sh[9] = NULL;
  (s).sh[10] = NULL;
}

//
// End of primitives for interacting between tetrahedra and subfaces
//

//
// Begin of primitives for interacting between subfaces and subsegs
//

// ssbond() bond a subface to a subsegment.

inline void tetgenmesh::ssbond(face& s, face& edge) 
{
  s.sh[6 + (s.shver >> 1)] = sencode(edge);
  edge.sh[0] = sencode(s);
}

inline void tetgenmesh::ssbond1(face& s, face& edge) 
{
  s.sh[6 + (s.shver >> 1)] = sencode(edge);
  //edge.sh[0] = sencode(s);
}

// ssdisolve() dissolve a bond (from the subface side)

inline void tetgenmesh::ssdissolve(face& s) 
{
  s.sh[6 + (s.shver >> 1)] = NULL;
}

// sspivot() finds a subsegment abutting a subface.

inline void tetgenmesh::sspivot(face& s, face& edge) 
{
  shellface sptr = (shellface) s.sh[6 + (s.shver >> 1)];
  sdecode(sptr, edge);
}

//
// End of primitives for interacting between subfaces and subsegs
//

//
// Begin of primitives for interacting between tet and subsegs.
//

inline void tetgenmesh::tssbond1(triface& t, face& s)
{
  if ((t).tet[8] == NULL) {
    // Allocate space for this tet.
    (t).tet[8] = (tetrahedron) tet2segpool->alloc();
    // NULL all fields in this space.
    for (int i = 0; i < 6; i++) {
      ((shellface *) (t).tet[8])[i] = NULL;
    }
  }
  ((shellface *) (t).tet[8])[ver2edge[(t).ver]] = sencode((s)); 
}

inline void tetgenmesh::sstbond1(face& s, triface& t) 
{
  ((tetrahedron *) (s).sh)[9] = encode(t);
}

inline void tetgenmesh::tssdissolve1(triface& t)
{
  if ((t).tet[8] != NULL) {
    ((shellface *) (t).tet[8])[ver2edge[(t).ver]] = NULL;
  }
}

inline void tetgenmesh::sstdissolve1(face& s) 
{
  ((tetrahedron *) (s).sh)[9] = NULL;
}

inline void tetgenmesh::tsspivot1(triface& t, face& s)
{
  if ((t).tet[8] != NULL) {
    sdecode(((shellface *) (t).tet[8])[ver2edge[(t).ver]], s);
  } else {
    (s).sh = NULL;
  }
}

inline void tetgenmesh::sstpivot1(face& s, triface& t) 
{
  decode((tetrahedron) s.sh[9], t);
}

//
// End of primitives for interacting between tet and subsegs.
//

//
// Begin of primitives for points
//

inline int tetgenmesh::pointmark(point pt) { 
  return ((int *) (pt))[pointmarkindex]; 
}

inline void tetgenmesh::setpointmark(point pt, int value) {
  ((int *) (pt))[pointmarkindex] = value;
}


// These two primitives set and read the type of the point.

inline enum tetgenmesh::verttype tetgenmesh::pointtype(point pt) {
  return (enum verttype) (((int *) (pt))[pointmarkindex + 1] >> (int) 8);
}

inline void tetgenmesh::setpointtype(point pt, enum verttype value) {
  ((int *) (pt))[pointmarkindex + 1] = 
    ((int) value << 8) + (((int *) (pt))[pointmarkindex + 1] & (int) 255);
}

// Read and set the geometry tag of the point (used by -s option).

inline int tetgenmesh::pointgeomtag(point pt) { 
  return ((int *) (pt))[pointmarkindex + 2]; 
}

inline void tetgenmesh::setpointgeomtag(point pt, int value) {
  ((int *) (pt))[pointmarkindex + 2] = value;
}

// Read and set the u,v coordinates of the point (used by -s option).

inline REAL tetgenmesh::pointgeomuv(point pt, int i) {
  return pt[pointparamindex + i];
}

inline void tetgenmesh::setpointgeomuv(point pt, int i, REAL value) {
  pt[pointparamindex + i] = value;
}

// pinfect(), puninfect(), pinfected() -- primitives to flag or unflag
//   a point. The last bit of the integer '[pointindex+1]' is flaged.

inline void tetgenmesh::pinfect(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 1;
}

inline void tetgenmesh::puninfect(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 1;
}

inline bool tetgenmesh::pinfected(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 1) != 0;
}

// pmarktest(), punmarktest(), pmarktested() -- primitives to mark or unmark
//   a point. 

inline void tetgenmesh::pmarktest(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 2;
}

inline void tetgenmesh::punmarktest(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 2;
}

inline bool tetgenmesh::pmarktested(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 2) != 0;
}

// pmarktest2(), ...

inline void tetgenmesh::pmarktest2(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 4;
}

inline void tetgenmesh::punmarktest2(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 4;
}

inline bool tetgenmesh::pmarktest2ed(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 4) != 0;
}

// pmarktest3(), ...

inline void tetgenmesh::pmarktest3(point pt) {
  ((int *) (pt))[pointmarkindex + 1] |= (int) 8;
}

inline void tetgenmesh::punmarktest3(point pt) {
  ((int *) (pt))[pointmarkindex + 1] &= ~(int) 8;
}

inline bool tetgenmesh::pmarktest3ed(point pt) {
  return (((int *) (pt))[pointmarkindex + 1] & (int) 8) != 0;
}


// These following primitives set and read a pointer to a tetrahedron
//   a subface/subsegment, a point, or a tet of background mesh.

inline tetgenmesh::tetrahedron tetgenmesh::point2tet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex];
}

inline void tetgenmesh::setpoint2tet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex] = value;
}

inline tetgenmesh::point tetgenmesh::point2ppt(point pt) {
  return (point) ((tetrahedron *) (pt))[point2simindex + 1];
}

inline void tetgenmesh::setpoint2ppt(point pt, point value) {
  ((tetrahedron *) (pt))[point2simindex + 1] = (tetrahedron) value;
}

inline tetgenmesh::shellface tetgenmesh::point2sh(point pt) {
  return (shellface) ((tetrahedron *) (pt))[point2simindex + 2];
}

inline void tetgenmesh::setpoint2sh(point pt, shellface value) {
  ((tetrahedron *) (pt))[point2simindex + 2] = (tetrahedron) value;
}


inline tetgenmesh::tetrahedron tetgenmesh::point2bgmtet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex + 3];
}

inline void tetgenmesh::setpoint2bgmtet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex + 3] = value;
}


// point2tetorg()    Get the tetrahedron whose origin is the point.

inline void tetgenmesh::point2tetorg(point pa, triface& searchtet)
{
  decode(point2tet(pa), searchtet);
  if ((point) searchtet.tet[4] == pa) {
    searchtet.ver = 11;
  } else if ((point) searchtet.tet[5] == pa) {
    searchtet.ver = 3;
  } else if ((point) searchtet.tet[6] == pa) {
    searchtet.ver = 7;
  } else {
    assert((point) searchtet.tet[7] == pa); // SELF_CHECK
    searchtet.ver = 0;
  }
}

// point2shorg()    Get the subface/segment whose origin is the point.

inline void tetgenmesh::point2shorg(point pa, face& searchsh)
{
  sdecode(point2sh(pa), searchsh);
  if ((point) searchsh.sh[3] == pa) {
    searchsh.shver = 0;
  } else if ((point) searchsh.sh[4] == pa) {
    searchsh.shver = (searchsh.sh[5] != NULL ? 2 : 1); 
  } else {
    assert((point) searchsh.sh[5] == pa); // SELF_CHECK
    searchsh.shver = 4;
  }
}

// farsorg()    Return the origin of the subsegment.
// farsdest()   Return the destination of the subsegment.

inline tetgenmesh::point tetgenmesh::farsorg(face& s)
{
  face travesh, neighsh;

  travesh = s;
  while (1) {
    senext2(travesh, neighsh);
    spivotself(neighsh); 
    if (neighsh.sh == NULL) break;
    if (sorg(neighsh) != sorg(travesh)) sesymself(neighsh);
    assert(sorg(neighsh) == sorg(travesh)); // SELF_CHECK
    senext2(neighsh, travesh); 
  }
  return sorg(travesh);
}

inline tetgenmesh::point tetgenmesh::farsdest(face& s) 
{
  face travesh, neighsh;

  travesh = s;
  while (1) {
    senext(travesh, neighsh);
    spivotself(neighsh); 
    if (neighsh.sh == NULL) break;
    if (sdest(neighsh) != sdest(travesh)) sesymself(neighsh);
    assert(sdest(neighsh) == sdest(travesh)); // SELF_CHECK
    senext(neighsh, travesh); 
  }
  return sdest(travesh);
}

//
// End of primitives for points
//

// dot() returns the dot product: v1 dot v2.

inline REAL tetgenmesh::dot(REAL* v1, REAL* v2) 
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.

inline void tetgenmesh::cross(REAL* v1, REAL* v2, REAL* n) 
{
  n[0] =   v1[1] * v2[2] - v2[1] * v1[2];
  n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
  n[2] =   v1[0] * v2[1] - v2[0] * v1[1];
}

// distance() computs the Euclidean distance between two points.

inline REAL tetgenmesh::distance(REAL* p1, REAL* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}

// Linear algebra operators.

#define NORM2(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))

#define DIST(p1, p2) \
  sqrt(NORM2((p2)[0] - (p1)[0], (p2)[1] - (p1)[1], (p2)[2] - (p1)[2]))

#define DOT(v1, v2) \
  ((v1)[0] * (v2)[0] + (v1)[1] * (v2)[1] + (v1)[2] * (v2)[2])

#define CROSS(v1, v2, n) \
  (n)[0] =   (v1)[1] * (v2)[2] - (v2)[1] * (v1)[2];\
  (n)[1] = -((v1)[0] * (v2)[2] - (v2)[0] * (v1)[2]);\
  (n)[2] =   (v1)[0] * (v2)[1] - (v2)[0] * (v1)[1]

#define SETVECTOR3(V, a0, a1, a2) (V)[0] = (a0); (V)[1] = (a1); (V)[2] = (a2)

#define SWAP2(a0, a1, tmp) (tmp) = (a0); (a0) = (a1); (a1) = (tmp)

#endif // #ifndef tetgenH

