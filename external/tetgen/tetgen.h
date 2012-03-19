///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator         //
//                                                                           //
// Version 1.4                                                               //
// September 6, December 13, 2010                                            //
// January 19, 2011                                                          //
//                                                                           //
// Copyright (C) 2002--2011                                                  //
// Hang Si                                                                   //
// Research Group: Numerical Mathematics and Scientific Computing            //
// Weierstrass Institute for Applied Analysis and Stochastics (WIAS)         //
// Mohrenstr. 39, 10117 Berlin, Germany                                      //
// si@wias-berlin.de                                                         //
//                                                                           //
// TetGen is freely available through the website: http://tetgen.berlios.de. //
//   It may be copied, modified, and redistributed for non-commercial use.   //
//   Please consult the file LICENSE for the detailed copyright notices.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen is a library to generate tetrahedral meshes for 3D domains.  It's  //
//   main purpose is to generate suitable tetrahedral meshes for numerical   //
//   simulations using finite element and finite volume methods.             //
//                                                                           //
// TetGen incorporates a suit of geometrical and mesh generation algorithms. //
//   A brief description of algorithms used in TetGen is found in the first  //
//   section of the user's manual.  References are given for users who are   //
//   interesting in these approaches. The main references are given below:   //
//                                                                           //
//   The efficient Delaunay tetrahedralization algorithm is: H. Edelsbrunner //
//   and N. R. Shah, "Incremental Topological Flipping Works for Regular     //
//   Triangulations". Algorithmica 15: 223--241, 1996.                       //
//                                                                           //
//   The constrained Delaunay tetrahedralization algorithm is described in:  //
//   H. Si and K. Gaertner,  "Meshing Piecewise Linear Complexes by Constr-  //
//   ained Delaunay Tetrahedralizations".  In Proceeding of the 14th Inter-  //
//   national Meshing Roundtable. September 2005.                            //
//                                                                           //
//   The mesh refinement algorithm is from:  Hang Si, "Adaptive Tetrahedral  //
//   Mesh Generation by Constrained Delaunay Refinement".  International     //
//   Journal for Numerical Methods in Engineering, 75(7): 856--880, 2008.    //
//                                                                           //
// The mesh data structure of TetGen is a combination of two types of mesh   //
//   data structures.  The tetrahedron-based mesh data structure introduced  //
//   by Shewchuk is eligible for tetrahedralization algorithms. The triangle //
//   -edge data structure developed by Muecke is adopted for representing    //
//   boundary elements: subfaces and subsegments.                            //
//                                                                           //
//   J. R. Shewchuk, "Delaunay Refinement Mesh Generation". PhD thesis,      //
//   Carnegie Mellon University, Pittsburgh, PA, 1997.                       //
//                                                                           //
//   E. P. Muecke, "Shapes and Implementations in Three-Dimensional          //
//   Geometry". PhD thesis, Univ. of Illinois, Urbana, Illinois, 1993.       //
//                                                                           //
// The research of mesh generation is definitly on the move. Many State-of-  //
//   the-art algorithms need implementing and evaluating. I heartily welcome //
//   any new algorithm especially for generating quality conforming Delaunay //
//   meshes and anisotropic conforming Delaunay meshes.                      //
//                                                                           //
// TetGen is supported by the "pdelib" project of Weierstrass Institute for  //
//   Applied Analysis and Stochastics (WIAS) in Berlin.  It is a collection  //
//   of software components for solving non-linear partial differential      //
//   equations including 2D and 3D mesh generators, sparse matrix solvers,   //
//   and scientific visualization tools, etc.  For more information please   //
//   visit: http://www.wias-berlin.de/software/pdelib.                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgen.h                                                                  //
//                                                                           //
// Header file of the TetGen library. Also is the user-level header file.    //
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

// The types 'intptr_t' and 'uintptr_t' are signed and unsigned integer types,
//   respectively. They are guaranteed to be the same width as a pointer.
//   They are defined in <stdint.h> by the C99 Standard.
//   However, Microsoft Visual C++ doesn't ship with this header file yet. We
//   need to define them. (Thanks to Steven G. Johnson from MIT for the 
//   following piece of code.) 

// Define the _MSC_VER symbol if you are using Microsoft Visual C++.

// #define _MSC_VER

// Define the _WIN64 symbol if you are running TetGen on Win64.

// #define _WIN64

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

// To compile TetGen as a library instead of an executable program, define
//   the TETLIBRARY symbol.

// #define TETLIBRARY

// Uncomment the following line to disable assert macros. These macros are
//   inserted in places where I hope to catch bugs.

// #define NDEBUG

// To insert lots of self-checks for internal errors, define the SELF_CHECK
//   symbol.  This will slow down the program a bit. 

// #define SELF_CHECK

// For single precision ( which will save some memory and reduce paging ),
//   define the symbol SINGLE by using the -DSINGLE compiler switch or by
//   writing "#define SINGLE" below.
//
// For double precision ( which will allow you to refine meshes to a smaller
//   edge length), leave SINGLE undefined.

// #define SINGLE

#ifdef SINGLE
  #define REAL float
#else
  #define REAL double
#endif 	// not defined SINGLE

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen Library Overview                                                   //
//                                                                           //
// TetGen library is comprised by several data types and global functions.   //
//                                                                           //
// There are three main data types: tetgenio, tetgenbehavior, and tetgenmesh.//
// Tetgenio is used to pass data into and out of TetGen library; tetgenbeha- //
// vior keeps the runtime options and thus controls the behaviors of TetGen; //
// tetgenmesh, the biggest data type I've ever defined, contains mesh data   //
// structures and mesh traversing and transformation operators.  The meshing //
// algorithms are implemented on top of it.  These data types are defined as //
// C++ classes.                                                              //
//                                                                           //
// There are few global functions. tetrahedralize() is provided for calling  //
// TetGen from another program. Two functions: orient3d() and insphere() are //
// incorporated from a public C code provided by Shewchuk.  They performing  //
// exact geometrical tests.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class tetgenio                                                            //
//                                                                           //
// The interface for passing data into and out of the library of TetGen.     //
//                                                                           //
// The tetgenio data structure is actually a collection of arrays of points, //
// facets, tetrahedra, and so forth.  The library will read and write these  //
// arrays according to the options specified in tetgenbehavior structure.    //
//                                                                           //
// If you want to program with the library of TetGen, it's necessary for you //
// to understand this data type,while the other two structures can be hidden //
// through calling the global function "tetrahedralize()". Each array corre- //
// sponds to a list of data in the file formats of TetGen.  It is necessary  //
// to understand TetGen's input/output file formats (see user's manual).     //
//                                                                           //
// Once an object of tetgenio is declared,  no array is created. One has to  //
// allocate enough memory for them, e.g., use the "new" operator in C++. On  //
// deletion of the object, the memory occupied by these arrays needs to be   //
// freed.  Routine deinitialize() will be automatically called. It will de-  //
// allocate the memory for an array if it is not a NULL. However, it assumes //
// that the memory is allocated by the C++ "new" operator. If you use malloc //
// (), you should free() them and set the pointers to NULLs before reaching  //
// deinitialize().                                                           //
//                                                                           //
// tetgenio ontains routines for reading and writing TetGen's files, i.e.,   //
// .node, .poly, .smesh, .ele, .face, and .edge files.  Both the library of  //
// TetGen and TetView use these routines to process input files.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenio {

  public:

  // Maximum number of characters in a file name (including the null).
  enum {FILENAMESIZE = 1024};

  // Maxi. numbers of chars in a line read from a file (incl. the null).
  enum {INPUTLINESIZE = 1024};

  // The polygon data structure.  A "polygon" describes a simple polygon
  //   (no holes). It is not necessarily convex.  Each polygon contains a
  //   number of corners (points) and the same number of sides (edges).
  // Note that the points of the polygon must be given in either counter-
  //   clockwise or clockwise order and they form a ring, so every two
  //   consective points forms an edge of the polygon.
  typedef struct {
    int *vertexlist;
    int numberofvertices;
  } polygon;

  static void init(polygon* p) {
    p->vertexlist = (int *) NULL;
    p->numberofvertices = 0;
  }

  // The facet data structure.  A "facet" describes a facet. Each facet is 
  //   a polygonal region possibly with holes, edges, and points in it.
  typedef struct {
    polygon *polygonlist;
    int numberofpolygons;
    REAL *holelist;
    int numberofholes;
  } facet;

  static void init(facet* f) {
    f->polygonlist = (polygon *) NULL;
    f->numberofpolygons = 0;
    f->holelist = (REAL *) NULL;
    f->numberofholes = 0;
  }

  // A 'voroedge' is an edge of the Voronoi diagram. It corresponds to a
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

  // A 'vorofacet' is an facet of the Voronoi diagram. It corresponds to a
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

  // A callback function for mesh refinement.
  typedef bool (* TetSizeFunc)(REAL*, REAL*, REAL*, REAL*, REAL*, REAL);

  // Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
  int firstnumber; 

  // Dimension of the mesh (2 or 3), default is 3.
  int mesh_dim;

  // Does the lines in .node file contain index or not, default is TRUE.
  bool useindex;

  // 'pointlist':  An array of point coordinates.  The first point's x
  //   coordinate is at index [0] and its y coordinate at index [1], its
  //   z coordinate is at index [2], followed by the coordinates of the
  //   remaining points.  Each point occupies three REALs. 
  // 'pointattributelist':  An array of point attributes.  Each point's
  //   attributes occupy 'numberofpointattributes' REALs.
  // 'pointmtrlist': An array of metric tensors at points. Each point's
  //   tensor occupies 'numberofpointmtr' REALs.
  // `pointmarkerlist':  An array of point markers; one int per point.
  REAL *pointlist;
  REAL *pointattributelist;
  REAL *pointmtrlist;
  int *pointmarkerlist;
  int numberofpoints;
  int numberofpointattributes;
  int numberofpointmtrs;
 
  // `elementlist':  An array of element (triangle or tetrahedron) corners.
  //   The first element's first corner is at index [0], followed by its
  //   other corners in counterclockwise order, followed by any other
  //   nodes if the element represents a nonlinear element.  Each element
  //   occupies `numberofcorners' ints.
  // `elementattributelist':  An array of element attributes.  Each
  //   element's attributes occupy `numberofelementattributes' REALs.
  // `elementconstraintlist':  An array of constraints, i.e. triangle's
  //   area or tetrahedron's volume; one REAL per element.  Input only.
  // `neighborlist':  An array of element neighbors; 3 or 4 ints per
  //   element.  Output only.
  int *tetrahedronlist;
  REAL *tetrahedronattributelist;
  REAL *tetrahedronvolumelist;
  int *neighborlist;
  int numberoftetrahedra;
  int numberofcorners;
  int numberoftetrahedronattributes;

  // `facetlist':  An array of facets.  Each entry is a structure of facet.
  // `facetmarkerlist':  An array of facet markers; one int per facet.
  facet *facetlist;
  int *facetmarkerlist;
  int numberoffacets;

  // `holelist':  An array of holes.  The first hole's x, y and z
  //   coordinates  are at indices [0], [1] and [2], followed by the
  //   remaining holes. Three REALs per hole. 
  REAL *holelist;
  int numberofholes;

  // `regionlist': An array of regional attributes and volume constraints.
  //   The first constraint's x, y and z coordinates are at indices [0],
  //   [1] and [2], followed by the regional attribute at index [3], foll-
  //   owed by the maximum volume at index [4]. Five REALs per constraint.
  // Note that each regional attribute is used only if you select the `A'
  //   switch, and each volume constraint is used only if you select the
  //   `a' switch (with no number following).
  REAL *regionlist;
  int numberofregions;

  // `facetconstraintlist': An array of facet maximal area constraints.
  //   Two REALs per constraint. The first (at index [0]) is the facet
  //   marker (cast it to int), the second (at index [1]) is its maximum
  //   area bound.
  REAL *facetconstraintlist;
  int numberoffacetconstraints;

  // `segmentconstraintlist': An array of segment max. length constraints.
  //   Three REALs per constraint. The first two (at indcies [0] and [1]) 
  //   are the indices of the endpoints of the segment, the third (at index
  //   [2]) is its maximum length bound.
  REAL *segmentconstraintlist;
  int numberofsegmentconstraints;

  // 'pbcgrouplist':  An array of periodic boundary condition groups.
  pbcgroup *pbcgrouplist;
  int numberofpbcgroups;

  // `trifacelist':  An array of triangular face endpoints.  The first
  //   face's endpoints are at indices [0], [1] and [2], followed by the
  //   remaining faces.  Three ints per face.
  // `adjtetlist':  An array of adjacent tetrahedra to the faces of
  //   trifacelist. Each face has at most two adjacent tets, the first
  //   face's adjacent tets are at [0], [1]. Two ints per face. A '-1'
  //   indicates outside (no adj. tet). This list is output when '-nn'
  //   switch is used.
  // `trifacemarkerlist':  An array of face markers; one int per face.
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

  // Input & output routines.
  bool load_node_call(FILE* infile, int markers, char* nodefilename);
  bool load_node(char* filebasename);
  bool load_var(char*);
  bool load_mtr(char*);
  bool load_poly(char*);
  bool load_pbc(char*);
  bool load_off(char*);
  bool load_ply(char*);
  bool load_stl(char*);
  bool load_medit(char*);
  bool load_vtk(char*);
  bool load_plc(char*, int);
  bool load_tetmesh(char*);
  void save_nodes(char*);
  void save_elements(char*);
  void save_faces(char*);
  void save_edges(char*);
  void save_neighbors(char*);
  void save_poly(char*);

  // Read line and parse string functions.
  char *readline(char* string, FILE* infile, int *linenumber);
  char *findnextfield(char* string);
  char *readnumberline(char* string, FILE* infile, char* infilename);
  char *findnextnumber(char* string);

  // Initialize routine.
  void initialize()
  {
    firstnumber = 0; // Default item index is numbered from Zero.
    mesh_dim = 3; // Default mesh dimension is 3.
    useindex = true;

    pointlist = (REAL *) NULL;
    pointattributelist = (REAL *) NULL;
    pointmtrlist = (REAL *) NULL;
    pointmarkerlist = (int *) NULL;
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
  }

  // Free the memory allocated in 'tetgenio'.  
  void deinitialize()
  {
    facet *f;
    polygon *p;
    pbcgroup *pg;
    int i, j;

    // This routine assumes that the memory was allocated by 
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
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class tetgenbehavior                                                      //
//                                                                           //
// The object holding a collection of options controlling TetGen's behavior. //
// See "command line switches" in User's manual.                             //
//                                                                           //
// parse_commandline() provides an simple interface to set the vaules of the //
// variables.  It accepts the standard parameters (e.g., 'argc' and 'argv')  //
// that pass to C/C++ main() function. Alternatively a string which contains //
// the command line options can be used as its parameter.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenbehavior {

  public:

  // Labels define the objects which are acceptable by TetGen. They are 
  //   recognized by the file extensions.
  //   - NODES, a list of nodes (.node); 
  //   - POLY, a piecewise linear complex (.poly or .smesh); 
  //   - OFF, a polyhedron (.off, Geomview's file format); 
  //   - PLY, a polyhedron (.ply, file format from gatech);
  //   - STL, a surface mesh (.stl, stereolithography format);
  //   - MEDIT, a surface mesh (.mesh, Medit's file format); 
  //   - MESH, a tetrahedral mesh (.ele).
  //   If no extension is available, the imposed commandline switch
  //   (-p or -r) implies the object. 

  enum objecttype {NONE, NODES, POLY, OFF, PLY, STL, MEDIT, VTK, MESH};

  // Variables of command line switches. Each variable corresponds to a
  //   switch and will be initialized.

  int plc;                                                 // '-p' switch, 0.
  int quality;                                             // '-q' switch, 0.
  int refine;                                              // '-r' switch, 0.
  int coarse;                                              // '-R' switch, 0.
  int metric;                                              // '-m' switch, 0.
  int varvolume;                            // '-a' switch without number, 0.
  int fixedvolume;                             // '-a' switch with number, 0.
  int insertaddpoints;                                     // '-i' switch, 0.
  int regionattrib;                                        // '-A' switch, 0.
  int conformdel;                                          // '-D' switch, 0.
  int diagnose;                                            // '-d' switch, 0.
  int zeroindex;                                           // '-z' switch, 0.
  int btree;                                                        // -u, 1.
  int max_btreenode_size;                            // number after -u, 100.
  int optlevel;                     // number specified after '-s' switch, 3.
  int optpasses;                   // number specified after '-ss' switch, 3.
  int order;                // element order, specified after '-o' switch, 1.
  int facesout;                                            // '-f' switch, 0.
  int edgesout;                                            // '-e' switch, 0.
  int neighout;                                            // '-n' switch, 0.
  int voroout;                                             // '-v',switch, 0.
  int meditview;                                           // '-g' switch, 0.
  int gidview;                                             // '-G' switch, 0.
  int geomview;                                            // '-O' switch, 0.
  int vtkview;                                             // '-K' switch, 0.
  int nobound;                                             // '-B' switch, 0.
  int nonodewritten;                                       // '-N' switch, 0.
  int noelewritten;                                        // '-E' switch, 0.
  int nofacewritten;                                       // '-F' switch, 0.
  int noiterationnum;                                      // '-I' switch, 0.
  int nomerge;                                             // '-M',switch, 0.
  int nobisect;             // count of how often '-Y' switch is selected, 0.
  int noflip;                        // do not perform flips. '-X' switch. 0.
  int nojettison;        // do not jettison redundants nodes. '-J' switch. 0.
  int steiner;                                // number after '-S' switch. 0.
  int fliprepair;                                          // '-X' switch, 1.
  int offcenter;                                           // '-R' switch, 0.
  int docheck;                                             // '-C' switch, 0.
  int quiet;                                               // '-Q' switch, 0.
  int verbose;              // count of how often '-V' switch is selected, 0.
  int useshelles;               // '-p', '-r', '-q', '-d', or '-R' switch, 0.
  int maxflipedgelinksize;        // The maximum flippable edge link size 10.
  REAL minratio;                            // number after '-q' switch, 2.0.
  REAL goodratio;                  // number calculated from 'minratio', 0.0.
  REAL minangle;                                // minimum angle bound, 20.0.
  REAL goodangle;                         // cosine squared of minangle, 0.0.
  REAL maxvolume;                          // number after '-a' switch, -1.0.
  REAL mindihedral;                        // number after '-qq' switch, 5.0.
  REAL maxdihedral;                     // number after '-qqq' switch, 165.0.
  REAL alpha1;                          // number after '-m' switch, sqrt(2).
  REAL alpha2;                             // number after '-mm' switch, 1.0.
  REAL alpha3;                            // number after '-mmm' switch, 0.6.
  REAL epsilon;                          // number after '-T' switch, 1.0e-8.
  REAL epsilon2;                        // number after '-TT' switch, 1.0e-5.
  enum objecttype object;            // determined by -p, or -r switch. NONE.

  // Variables used to save command line switches and in/out file names.
  char commandline[1024];
  char infilename[1024];
  char outfilename[1024];
  char addinfilename[1024];
  char bgmeshfilename[1024];

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
    quality = 0;
    refine = 0;
    coarse = 0;
    metric = 0;
    minratio = 2.0;
    goodratio = 0.0;
    minangle = 20.0;
    goodangle = 0.0;
    maxdihedral = 165.0;
    mindihedral = 5.0;
    varvolume = 0;
    fixedvolume = 0;
    maxvolume = -1.0;
    regionattrib = 0;
    insertaddpoints = 0;
    diagnose = 0;
    offcenter = 0;
    conformdel = 0;
    alpha1 = sqrt(2.0);
    alpha2 = 1.0;
    alpha3 = 0.6;
    zeroindex = 0;
    btree = 1;
    max_btreenode_size = 100;
    facesout = 0;
    edgesout = 0;
    neighout = 0;
    voroout = 0;
    meditview = 0;
    gidview = 0;
    geomview = 0;
    vtkview = 0;
    optlevel = 3;
    optpasses = 3;
    order = 1;
    nojettison = 0;
    nobound = 0;
    nonodewritten = 0;
    noelewritten = 0;
    nofacewritten = 0;
    noiterationnum = 0;
    nobisect = 0;
    noflip = 0;
    steiner = -1;
    fliprepair = 1;
    nomerge = 0;
    docheck = 0;
    quiet = 0;
    verbose = 0;
    useshelles = 0;
    maxflipedgelinksize = 10;
    epsilon = 1.0e-8;
    epsilon2 = 1.0e-5;
    object = NONE;

    commandline[0] = '\0';
    infilename[0] = '\0';
    outfilename[0] = '\0';
    addinfilename[0] = '\0';
    bgmeshfilename[0] = '\0';
  }
  
  ~tetgenbehavior() 
  {
  }
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class tetgenmesh                                                          //
//                                                                           //
// The object to store, generate, and refine a tetrahedral mesh.             //
//                                                                           //
// It implements the mesh data structures and functions to create and update //
// a tetrahedral mesh according to the specified options.                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenmesh {

  public:

  // Maximum number of characters in a file name (including the null).
  enum {FILENAMESIZE = 1024};

  // For efficiency, a variety of data structures are allocated in bulk.
  //   The following constants determine how many of each structure is
  //   allocated at once.
  enum {VERPERBLOCK = 4092, SUBPERBLOCK = 4092, ELEPERBLOCK = 8188};

  // Used for the point location scheme of Mucke, Saias, and Zhu, to
  //   decide how large a random sample of tetrahedra to inspect.
  enum {SAMPLEFACTOR = 11};

  // Labels that signify two edge rings of a triangle (see Muecke's thesis).
  enum {CCW = 0, CW = 1};

  // Labels that signify whether a record consists primarily of pointers
  //   or of floating-point words.  Used for data alignment.
  enum wordtype {POINTER, FLOATINGPOINT};

  // Labels that signify the type of a vertex. 
  enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, NACUTEVERTEX, ACUTEVERTEX,
    FREESEGVERTEX, FREESUBVERTEX, FREEVOLVERTEX, DEADVERTEX = -32768};
 
  // Labels that signify the type of a subface/subsegment.
  enum shestype {NSHARP, SHARP};

  // Labels that signify the type of flips can be applied on a face.
  enum fliptype {T23, T32, T22, T44, N32, N40, FORBIDDENFACE, FORBIDDENEDGE};

  // Labels that signify the result of triangle-triangle intersection test.
  enum interresult {DISJOINT, INTERSECT, SHAREVERTEX, SHAREEDGE, SHAREFACE,
    TOUCHEDGE, TOUCHFACE, INTERVERT, INTEREDGE, INTERFACE, INTERTET,
    TRIEDGEINT, EDGETRIINT, COLLISIONFACE, INTERSUBSEG, INTERSUBFACE,
    BELOWHULL2};

  // Labels that signify the result of point location.
  enum locateresult {INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX, OUTSIDE,
    ENCSEGMENT};

  // Labels that signify the result of vertex insertion. 
  enum insertsiteresult {SUCCESSINTET, SUCCESSONFACE, SUCCESSONEDGE,
    DUPLICATEPOINT, OUTSIDEPOINT};

  // Labels that signify the result of direction finding. 
  enum finddirectionresult {ACROSSEDGE, ACROSSFACE, LEFTCOLLINEAR,
    RIGHTCOLLINEAR, TOPCOLLINEAR, BELOWHULL};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh elements                                                             //
//                                                                           //
// There are four types of mesh elements: tetrahedra, subfaces, subsegments, //
// and points,  where subfaces and subsegments are triangles and edges which //
// appear on boundaries.  A tetrahedralization of a 3D point set comprises   //
// tetrahedra and points;  a surface mesh of a 3D domain comprises subfaces  //
// subsegments and points.  The elements of all the four types consist of a  //
// tetrahedral mesh of a 3D domain.  However, TetGen uses three data types:  //
// 'tetrahedron', 'shellface', and 'point'. A 'tetrahedron' is a tetrahedron;//
// while a 'shellface' can be either a subface or a subsegment; and a 'point'//
// is a point.  These three data types, linked by pointers comprise a mesh.  //
//                                                                           //
// A tetrahedron primarily consists of a list of 4 pointers to its corners,  //
// a list of 4 pointers to its adjoining tetrahedra, a list of 4 pointers to //
// its adjoining subfaces (when subfaces are needed). Optinoally, (depending //
// on the selected switches), it may contain an arbitrary number of user-    //
// defined floating-point attributes,  an optional maximum volume constraint //
// (for -a switch), and a pointer to a list of high-order nodes (-o2 switch).//
// Since the size of a tetrahedron is not determined until running time.     //
//                                                                           //
// The data structure of tetrahedron also stores the geometrical information.//
// Let t be a tetrahedron, v0, v1, v2, and v3 be the 4 nodes corresponding   //
// to the order of their storage in t.  v3 always has a negative orientation //
// with respect to v0, v1, v2 (ie,, v3 lies above the oriented plane passes  //
// through v0, v1, v2). Let the 4 faces of t be f0, f1, f2, and f3. Vertices //
// of each face are stipulated as follows: f0 (v0, v1, v2), f1 (v0, v3, v1), //
// f2 (v1, v3, v2), f3 (v2, v3, v0).                                         //
//                                                                           //
// A subface has 3 pointers to vertices, 3 pointers to adjoining subfaces, 3 //
// pointers to adjoining subsegments, 2 pointers to adjoining tetrahedra, a  //
// boundary marker(an integer). Like a tetrahedron, the pointers to vertices,//
// subfaces, and subsegments are ordered in a way that indicates their geom- //
// etric relation.  Let s be a subface, v0, v1 and v2 be the 3 nodes corres- //
// ponding to the order of their storage in s,  e0, e1 and e2 be the 3 edges,//
// then we have: e0 (v0, v1), e1 (v1, v2), e2 (v2, v0).                      //
//                                                                           //
// A subsegment has exactly the same data fields as a subface has, but only  //
// uses some of them. It has 2 pointers to its endpoints, 2 pointers to its  //
// adjoining (and collinear) subsegments, a pointer to a subface containing  //
// it (there may exist any number of subfaces having it, choose one of them  //
// arbitrarily). The geometric relation between its endpoints and adjoining  //
// subsegments is kept with respect to the storing order of its endpoints.   //
//                                                                           //
// The data structure of point is relatively simple.  A point is a list of   //
// floating-point numbers, starting with the x, y, and z coords, followed by //
// an arbitrary number of optional user-defined floating-point attributes,   //
// an integer boundary marker, an integer for the point type, and a pointer  //
// to a tetrahedron (used for speeding up point location).                   //
//                                                                           //
// For a tetrahedron on a boundary (or a hull) of the mesh, some or all of   //
// the adjoining tetrahedra may not be present. For an interior tetrahedron, //
// often no neighboring subfaces are present,  Such absent tetrahedra and    //
// subfaces are never represented by the NULL pointers; they are represented //
// by two special records: `dummytet', the tetrahedron fills "outer space",  //
// and `dummysh',  the vacuous subfaces which are omnipresent.               //
//                                                                           //
// Tetrahedra and adjoining subfaces are glued together through the pointers //
// saved in each data fields of them. Subfaces and adjoining subsegments are //
// connected in the same fashion.  However, there are no pointers directly   //
// gluing tetrahedra and adjoining subsegments.  For the purpose of saving   //
// space, the connections between tetrahedra and subsegments are entirely    //
// mediated through subfaces.  The following part explains how subfaces are  //
// connected in TetGen.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Subface-subface and subface-subsegment connections                        //
//                                                                           //
// Adjoining subfaces sharing a common edge are connected in such a way that //
// they form a face ring around the edge. It is indeed a single linked list  //
// which is cyclic, e.g., one can start from any subface in it and traverse  //
// back. When the edge is not a subsegment, the ring only has two coplanar   //
// subfaces which are pointing to each other. Otherwise, the face ring may   //
// have any number of subfaces (and are not all coplanar).                   //
//                                                                           //
// How is the face ring formed?  Let s be a subsegment, f is one of subfaces //
// containing s as an edge.  The direction of s is stipulated from its first //
// endpoint to its second (according to their storage in s). Once the dir of //
// s is determined, the other two edges of f are oriented to follow this dir.//
// The "directional normal" N_f is a vector formed from any point in f and a //
// points orthogonally above f.                                              //
//                                                                           //
// The face ring of s is a cyclic ordered set of subfaces containing s, i.e.,//
// F(s) = {f1, f2, ..., fn}, n >= 1.  Where the order is defined as follows: //
// let fi, fj be two faces in F(s), the "normal-angle", NAngle(i,j) (range   //
// from 0 to 360 degree) is the angle between the N_fi and N_fj;  then fi is //
// in front of fj (or symbolically, fi < fj) if there exists another fk in   //
// F(s), and NAangle(k, i) < NAngle(k, j).  The face ring of s is: f1 < f2 < //
// ... < fn < f1.                                                            //
//                                                                           //
// The easiest way to imagine how a face ring is formed is to use the right- //
// hand rule.  Make a fist using your right hand with the thumb pointing to  //
// the direction of the subsegment. The face ring is connected following the //
// direction of your fingers.                                                //
//                                                                           //
// The subface and subsegment are also connected through pointers stored in  //
// their own data fields.  Every subface has a pointer to its adjoining sub- //
// segment. However, a subsegment only has one pointer to a subface which is //
// containing it. Such subface can be chosen arbitrarily, other subfaces are //
// found through the face ring.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // The tetrahedron data structure.  Fields of a tetrahedron contains:
  //   - a list of four adjoining tetrahedra;
  //   - a list of four vertices;
  //   - a list of four subfaces (optional, used for -p switch);
  //   - a list of user-defined floating-point attributes (optional);
  //   - a volume constraint (optional, used for -a switch);
  //   - an integer of element marker (optional, used for -n switch);
  //   - a pointer to a list of high-ordered nodes (optional, -o2 switch);

  typedef REAL **tetrahedron;

  // The shellface data structure.  Fields of a shellface contains:
  //   - a list of three adjoining subfaces;
  //   - a list of three vertices;
  //   - a list of two adjoining tetrahedra;
  //   - a list of three adjoining subsegments;
  //   - a pointer to a badface containing it (used for -q);
  //   - an area constraint (optional, used for -q);
  //   - an integer for boundary marker;
  //   - an integer for type: SHARPSEGMENT, NONSHARPSEGMENT, ...;
  //   - an integer for pbc group (optional, if in->pbcgrouplist exists);

  typedef REAL **shellface;

  // The point data structure.  It is actually an array of REALs:
  //   - x, y and z coordinates;
  //   - a list of user-defined point attributes (optional);
  //   - a list of REALs of a user-defined metric tensor (optional);
  //   - a pointer to a simplex (tet, tri, edge, or vertex);
  //   - a pointer to a parent (or duplicate) point;
  //   - a pointer to a tet in background mesh (optional);
  //   - a pointer to another pbc point (optional);
  //   - an integer for boundary marker;
  //   - an integer for verttype: INPUTVERTEX, FREEVERTEX, ...;

  typedef REAL *point;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh handles                                                              //
//                                                                           //
// Two special data types, 'triface' and 'face' are defined for maintaining  //
// and updating meshes. They are like pointers (or handles), which allow you //
// to hold one particular part of the mesh, i.e., a tetrahedron, a triangle, //
// an edge and a vertex.  However, these data types do not themselves store  //
// any part of the mesh. The mesh is made of the data types defined above.   //
//                                                                           //
// Muecke's "triangle-edge" data structure is the prototype for these data   //
// types.  It allows a universal representation for every tetrahedron,       //
// triangle, edge and vertex.  For understanding the following descriptions  //
// of these handle data structures,  readers are required to read both the   //
// introduction and implementation detail of "triangle-edge" data structure  //
// in Muecke's thesis.                                                       //
//                                                                           //
// A 'triface' represents a face of a tetrahedron and an oriented edge of    //
// the face simultaneously.  It has a pointer 'tet' to a tetrahedron, an     //
// integer 'loc' (range from 0 to 3) as the face index, and an integer 'ver' //
// (range from 0 to 5) as the edge version. A face of the tetrahedron can be //
// uniquly determined by the pair (tet, loc), and an oriented edge of this   //
// face can be uniquly determined by the triple (tet, loc, ver).  Therefore, //
// different usages of one triface are possible.  If we only use the pair    //
// (tet, loc), it refers to a face, and if we add the 'ver' additionally to  //
// the pair, it is an oriented edge of this face.                            //
//                                                                           //
// A 'face' represents a subface and an oriented edge of it simultaneously.  //
// It has a pointer 'sh' to a subface, an integer 'shver'(range from 0 to 5) //
// as the edge version.  The pair (sh, shver) determines a unique oriented   //
// edge of this subface.  A 'face' is also used to represent a subsegment,   //
// in this case, 'sh' points to the subsegment, and 'shver' indicates the    //
// one of two orientations of this subsegment, hence, it only can be 0 or 1. //
//                                                                           //
// Mesh navigation and updating are accomplished through a set of mesh       //
// manipulation primitives which operate on trifaces and faces.  They are    //
// introduced below.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class triface {

    public:

    tetrahedron* tet;
    int loc, ver;

    // Constructors;
    triface() : tet(0), loc(0), ver(0) {}
    // Operators;
    triface& operator=(const triface& t) {
      tet = t.tet; loc = t.loc; ver = t.ver;
      return *this;
    }
    bool operator==(triface& t) {
      return tet == t.tet && loc == t.loc && ver == t.ver;
    }
    bool operator!=(triface& t) {
      return tet != t.tet || loc != t.loc || ver != t.ver;
    }
  };

  class face {

    public:

    shellface *sh;
    int shver;

    // Constructors;
    face() : sh(0), shver(0) {}
    // Operators;
    face& operator=(const face& s) {
      sh = s.sh; shver = s.shver;
      return *this;
    }
    bool operator==(face& s) {return (sh == s.sh) && (shver == s.shver);}
    bool operator!=(face& s) {return (sh != s.sh) || (shver != s.shver);}
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The badface structure                                                     //
//                                                                           //
// A multiple usages structure. Despite of its name, a 'badface' can be used //
// to represent the following objects:                                       //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;             //
//   - an encroached subsegment or subface;                                  //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;      //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;    //
//   - a degenerate tetrahedron (see routine checkdegetet()).                //
//   - a recently flipped face (saved for undoing the flip later).           //
//                                                                           //
// It has the following fields:  'tt' holds a tetrahedron; 'ss' holds a sub- //
// segment or subface; 'cent' is the circumcent of 'tt' or 'ss', 'key' is a  //
// special value depending on the use, it can be either the square of the    //
// radius-edge ratio of 'tt' or the flipped type of 'tt';  'forg', 'fdest',  //
// 'fapex', and 'foppo' are vertices saved for checking the object in 'tt'   //
// or 'ss' is still the same when it was stored; 'noppo' is the fifth vertex //
// of a degenerate point set.  'previtem' and 'nextitem' implement a double  //
// link for managing many basfaces.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  struct badface {
    triface tt; 
    face ss; 
    REAL key;
    REAL cent[3];
    point forg, fdest, fapex, foppo;
    point noppo;
    struct badface *previtem, *nextitem; 
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Elementary flip data structure                                            //
//                                                                           //
// A data structure to record three types of elementary flips, which are     //
// 2-to-3, 3-to-2, and 2-to-2 flips.                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class elemflip {

    public:

    enum fliptype ft; // ft \in {T23, T32, T22}.
    point pset1[3];
    point pset2[3];

    elemflip() {
      ft = T23; // Default.
      pset1[0] = pset1[1] = pset1[2] = (point) NULL;
      pset2[0] = pset2[1] = pset2[2] = (point) NULL;
    }

  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The pbcdata structure                                                     //
//                                                                           //
// A pbcdata stores data of a periodic boundary condition defined on a pair  //
// of facets or segments. Let f1 and f2 define a pbcgroup. 'fmark' saves the //
// facet markers of f1 and f2;  'ss' contains two subfaces belong to f1 and  //
// f2, respectively.  Let s1 and s2 define a segment pbcgroup. 'segid' are   //
// the segment ids of s1 and s2; 'ss' contains two segments belong to s1 and //
// s2, respectively. 'transmat' are two transformation matrices. transmat[0] //
// transforms a point of f1 (or s1) into a point of f2 (or s2),  transmat[1] //
// does the inverse.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  struct pbcdata {
    int fmark[2];
    int segid[2];
    face ss[2];
    REAL transmat[2][4][4];
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Fast lookup tables for mesh manipulation primitives.                      //
//                                                                           //
// Mesh manipulation primitives (given below) are basic operations on mesh   //
// data structures. They answer basic queries on mesh handles, such as "what //
// is the origin (or destination, or apex) of the face?", "what is the next  //
// (or previous) edge in the edge ring?", and "what is the next face in the  //
// face ring?", and so on.                                                   //
//                                                                           //
// The implementation of teste basic queries can take advangtage of the fact //
// that the mesh data structures additionally store geometric informations.  //
// For example, we have ordered the 4 vertices (from 0 to 3) and the 4 faces //
// (from 0 to 3) of a tetrahedron,  and for each face of the tetrahedron, a  //
// sequence of vertices has stipulated,  therefore the origin of any face of //
// the tetrahedron can be quickly determined by a table 'locver2org', which  //
// takes the index of the face and the edge version as inputs.  A list of    //
// fast lookup tables are defined below. They're just like global variables. //
// These tables are initialized at the runtime.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // For enext() primitive, uses 'ver' as the index. 
  static int ve[6];

  // For org(), dest() and apex() primitives, uses 'ver' as the index.
  static int vo[6], vd[6], va[6];

  // For org(), dest() and apex() primitives, uses 'loc' as the first
  //   index and 'ver' as the second index.
  static int locver2org[4][6];
  static int locver2dest[4][6];
  static int locver2apex[4][6];

  // For oppo() primitives, uses 'loc' as the index.
  static int loc2oppo[4];

  // For fnext() primitives, uses 'loc' as the first index and 'ver' as
  //   the second index,  returns an array containing a new 'loc' and a
  //   new 'ver'. Note: Only valid for 'ver' equals one of {0, 2, 4}.
  static int locver2nextf[4][6][2];

  // The edge number (from 0 to 5) of a tet is defined as follows:
  static int locver2edge[4][6];
  static int edge2locver[6][2];

  // The map from a given face ('loc') to the other three faces in the tet.
  //   and the map from a given face's edge ('loc', 'ver') to other two
  //   faces in the tet opposite to this edge. (used in speeding the Bowyer-
  //   Watson cavity construction).
  static int locpivot[4][3];
  static int locverpivot[4][6][2];

  // For enumerating three edges of a triangle.
  static int plus1mod3[3];
  static int minus1mod3[3];

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh manipulation primitives                                              //
//                                                                           //
// A serial of mesh operations such as topological maintenance,  navigation, //
// local modification, etc.,  is accomplished through a set of mesh manipul- //
// ation primitives. These primitives are indeed very simple functions which //
// take one or two handles ('triface's and 'face's) as parameters,  perform  //
// basic operations such as "glue two tetrahedra at a face",  "return the    //
// origin of a tetrahedron", "return the subface adjoining at the face of a  //
// tetrahedron", and so on.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Primitives for tetrahedra.
  inline void decode(tetrahedron ptr, triface& t);
  inline tetrahedron encode(triface& t);
  inline void sym(triface& t1, triface& t2);
  inline void symself(triface& t);
  inline void bond(triface& t1, triface& t2);
  inline void dissolve(triface& t);
  inline point org(triface& t);
  inline point dest(triface& t);
  inline point apex(triface& t);
  inline point oppo(triface& t);
  inline void setorg(triface& t, point pointptr);
  inline void setdest(triface& t, point pointptr);
  inline void setapex(triface& t, point pointptr);
  inline void setoppo(triface& t, point pointptr);
  inline void esym(triface& t1, triface& t2);
  inline void esymself(triface& t);
  inline void enext(triface& t1, triface& t2);
  inline void enextself(triface& t);
  inline void enext2(triface& t1, triface& t2);
  inline void enext2self(triface& t);
  inline bool fnext(triface& t1, triface& t2);
  inline bool fnextself(triface& t);
  inline void symedge(triface& t1, triface& t2);
  inline void symedgeself(triface& t);
  inline void tfnext(triface& t1, triface& t2);
  inline void tfnextself(triface& t);
  inline void enextfnext(triface& t1, triface& t2);
  inline void enextfnextself(triface& t);
  inline void enext2fnext(triface& t1, triface& t2);
  inline void enext2fnextself(triface& t);
  inline REAL elemattribute(tetrahedron* ptr, int attnum);
  inline void setelemattribute(tetrahedron* ptr, int attnum, REAL value);
  inline REAL volumebound(tetrahedron* ptr);
  inline void setvolumebound(tetrahedron* ptr, REAL value);
  inline int getelemmarker(tetrahedron* ptr);
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
 
  // Primitives for subfaces and subsegments.
  inline void sdecode(shellface sptr, face& s);
  inline shellface sencode(face& s);
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
  inline void sfnext(face&, face&);
  inline void sfnextself(face&);
  inline badface* shell2badface(face& s);
  inline void setshell2badface(face& s, badface* value);
  inline REAL areabound(face& s);
  inline void setareabound(face& s, REAL value);
  inline int shellmark(face& s);
  inline void setshellmark(face& s, int value);
  inline enum shestype shelltype(face& s);
  inline void setshelltype(face& s, enum shestype value); 
  inline int shellpbcgroup(face& s);
  inline void setshellpbcgroup(face& s, int value);
  inline void sinfect(face& s);
  inline void suninfect(face& s);
  inline bool sinfected(face& s);

  // Primitives for interacting tetrahedra and subfaces.
  inline void tspivot(triface& t, face& s);
  inline void stpivot(face& s, triface& t);
  inline void tsbond(triface& t, face& s);
  inline void tsdissolve(triface& t);
  inline void stdissolve(face& s);

  // Primitives for interacting subfaces and subsegs.
  inline void sspivot(face& s, face& edge);
  inline void ssbond(face& s, face& edge);
  inline void ssdissolve(face& s);

  inline void tsspivot1(triface& t, face& seg);
  inline void tssbond1(triface& t, face& seg);
  inline void tssdissolve1(triface& t);

  // Primitives for points.
  inline int  pointmark(point pt);
  inline void setpointmark(point pt, int value);
  inline enum verttype pointtype(point pt);
  inline void setpointtype(point pt, enum verttype value);
  inline void pinfect(point pt);
  inline void puninfect(point pt);
  inline bool pinfected(point pt);
  inline tetrahedron point2tet(point pt);
  inline void setpoint2tet(point pt, tetrahedron value);
  inline shellface point2sh(point pt);
  inline void setpoint2sh(point pt, shellface value);
  inline shellface point2seg(point pt);
  inline void setpoint2seg(point pt, shellface value);
  inline point point2ppt(point pt);
  inline void setpoint2ppt(point pt, point value);
  inline tetrahedron point2bgmtet(point pt);
  inline void setpoint2bgmtet(point pt, tetrahedron value);
  inline point point2pbcpt(point pt);
  inline void setpoint2pbcpt(point pt, point value);

  // Advanced primitives.
  inline void adjustedgering(triface& t, int direction);
  inline void adjustedgering(face& s, int direction);
  inline bool isdead(triface* t);
  inline bool isdead(face* s);
  inline bool isfacehaspoint(triface* t, point testpoint);
  inline bool isfacehaspoint(face* t, point testpoint);
  inline bool isfacehasedge(face* s, point tend1, point tend2);
  inline bool issymexist(triface* t);
  void getnextsface(face*, face*);
  void tsspivot(triface*, face*);
  void sstpivot(face*, triface*);
  void point2tetorg(point, triface&);
  void point2shorg(point, face&);
  void point2segorg(point, face&);
  bool findorg(triface* t, point dorg);
  bool findorg(face* s, point dorg);
  void findedge(triface* t, point eorg, point edest);
  void findedge(face* s, point eorg, point edest);
  void getonextseg(face* s, face* lseg);
  void getseghasorg(face* sseg, point dorg);
  point getsubsegfarorg(face* sseg);
  point getsubsegfardest(face* sseg);
  void printtet(triface*);
  void printsh(face*);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Arraypool                                                                 //
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


// A function: int cmp(const T &, const T &),  is said to realize a
//   linear order on the type T if there is a linear order <= on T such
//   that for all x and y in T satisfy the following relation:
//                 -1  if x < y.
//   comp(x, y) =   0  if x is equivalent to y.
//                 +1  if x > y.
// A 'compfunc' is a pointer to a linear-order function. 

  typedef int (*compfunc) (const void *, const void *);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// List                                                                      //
//                                                                           //
// An array of items with automatically reallocation of memory.              //
//                                                                           //
// 'base' is the starting address of the array.  'itembytes' is the size of  //
//   each item in byte.                                                      //
//                                                                           //
// 'items' is the number of items stored in list.  'maxitems' indicates how  //
//   many items can be stored in this list. 'expandsize' is the increasing   //
//   size (items) when the list is full.                                     //
//                                                                           //
// The index of list always starts from zero, i.e., for a list L contains    //
//   n elements, the first element is L[0], and the last element is L[n-1].  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class list {

    public:

    char *base;
    int  itembytes;
    int  items, maxitems, expandsize;
    compfunc comp;

    list(int itbytes, compfunc pcomp, int mitems = 256, int exsize = 128) {
      listinit(itbytes, pcomp, mitems, exsize);
    }
    ~list() { free(base); }

    void *operator[](int i) { return (void *) (base + i * itembytes); }

    void listinit(int itbytes, compfunc pcomp, int mitems, int exsize);
    void setcomp(compfunc compf) { comp = compf; }    
    void clear() { items = 0; }
    int  len() { return items; }
    void *append(void* appitem);
    void *insert(int pos, void* insitem);
    void del(int pos, int order);
    int  hasitem(void* checkitem);
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Memorypool                                                                //
//                                                                           //
// A type used to allocate memory.                                           //
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
// Queue                                                                     //
//                                                                           //
// A 'queue' is a FIFO data structure.                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  class queue : public memorypool {

    public:

    void **head, **tail;
    int  linkitembytes;
    int  linkitems; // Not count 'head' and 'tail'.

     queue(int bytecount, int itemcount = 256) {
       linkitembytes = bytecount;
       poolinit(bytecount + sizeof(void *), itemcount, POINTER, 0);
       head = (void **) alloc();
       tail = (void **) alloc();
       *head = (void *) tail;
       *tail = NULL;
       linkitems = 0;
     }

     void clear() {
       // Reset the pool.
       restart();
       // Initialize all variables.
       head = (void **) alloc();
       tail = (void **) alloc();
       *head = (void *) tail;
       *tail = NULL;
       linkitems = 0;
     }

     long len() { return linkitems; }
     bool empty() { return linkitems == 0; }

     void *push(void* newitem) {
       void **newnode = tail;
       if (newitem != (void *) NULL) {
         memcpy((void *)(newnode + 1), newitem, linkitembytes);
       }
       tail = (void **) alloc();
       *tail = NULL;
       *newnode = (void *) tail;
       linkitems++;
       return (void *)(newnode + 1);
     }

     void *pop() {
       if (linkitems > 0) {
         void **deadnode = (void **) *head;
         *head = *deadnode;
         dealloc((void *) deadnode);
         linkitems--;
         return (void *)(deadnode + 1);
       } else {
         return NULL;
       }
     }
  };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Memory managment routines                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void dummyinit(int, int);
  void initializepools();
  void tetrahedrondealloc(tetrahedron*);
  tetrahedron *tetrahedrontraverse();
  void shellfacedealloc(memorypool*, shellface*);
  shellface *shellfacetraverse(memorypool*);
  void badfacedealloc(memorypool*, badface*);
  badface *badfacetraverse(memorypool*);
  void pointdealloc(point);
  point pointtraverse();
  void maketetrahedron(triface*);
  void makeshellface(memorypool*, face*);
  void makepoint(point*);

  void makepoint2tetmap();
  void makepoint2segmap();
  void makeindex2pointmap(point*&);
  void makesegmentmap(int*&, shellface**&);
  void makesubfacemap(int*&, shellface**&);
  void maketetrahedronmap(int*&, tetrahedron**&);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Geometric functions                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // PI is the ratio of a circle's circumference to its diameter.
  static REAL PI;

  // Triangle-triangle intersection test
  enum interresult edge_vert_col_inter(REAL*, REAL*, REAL*);
  enum interresult edge_edge_cop_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
  enum interresult tri_vert_cop_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
  enum interresult tri_edge_cop_inter(REAL*, REAL*, REAL*,REAL*,REAL*,REAL*);
  enum interresult tri_edge_inter_tail(REAL*, REAL*, REAL*, REAL*, REAL*,
                                        REAL, REAL);
  enum interresult tri_edge_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
  enum interresult tri_tri_inter(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);
  int tri_edge_2d(point, point, point, point, point, point, int, int*, int*);
  int tri_edge_test(point, point, point, point, point, point, int, int*, int*);

  // Geometric tests
  REAL incircle3d(point pa, point pb, point pc, point pd);
  REAL insphere_s(REAL*, REAL*, REAL*, REAL*, REAL*);
  bool iscollinear(REAL*, REAL*, REAL*, REAL eps);
  bool iscoplanar(REAL*, REAL*, REAL*, REAL*, REAL vol6, REAL eps);
  bool iscospheric(REAL*, REAL*, REAL*, REAL*, REAL*, REAL vol24, REAL eps);

  // Linear algebra functions
  inline REAL dot(REAL* v1, REAL* v2);
  inline void cross(REAL* v1, REAL* v2, REAL* n);
  bool lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N);
  void lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N);

  // Geometric calculations
  inline REAL distance(REAL* p1, REAL* p2);
  REAL shortdistance(REAL* p, REAL* e1, REAL* e2);
  REAL shortdistance(REAL* p, REAL* e1, REAL* e2, REAL* e3);
  REAL interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n);
  void projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj);
  void projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj);
  void facenormal(REAL* pa, REAL* pb, REAL* pc, REAL* n, REAL* nlen);
  void facenormal2(point pa, point pb, point pc, REAL *n, int pivot);
  void edgeorthonormal(REAL* e1, REAL* e2, REAL* op, REAL* n);
  REAL facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2);
  void tetalldihedral(point, point, point, point, REAL*, REAL*, REAL*);
  void tetallnormal(point, point, point, point, REAL N[4][3], REAL* volume);
  REAL tetaspectratio(point, point, point, point);
  bool circumsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
  void inscribedsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
  void rotatepoint(REAL* p, REAL rotangle, REAL* p1, REAL* p2);
  void planelineint(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);

  // Point location routines.
  unsigned long randomnation(unsigned int choices);
  REAL distance2(tetrahedron* tetptr, point p);
  void randomsample(point searchpt, triface *searchtet);
  enum locateresult locate(point searchpt, triface* searchtet);
  enum locateresult locate2(point searchpt, triface* searchtet, arraypool*);
  enum locateresult preciselocate(point searchpt, triface* searchtet, long);
  enum locateresult adjustlocate(point, triface*, enum locateresult, REAL);
  enum locateresult hullwalk(point searchpt, triface* hulltet);
  enum locateresult locatesub(point searchpt, face* searchsh, int, REAL);
  enum locateresult adjustlocatesub(point, face*, enum locateresult, REAL);
  enum locateresult locateseg(point searchpt, face* searchseg);
  enum locateresult adjustlocateseg(point, face*, enum locateresult, REAL);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh update functions                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void enqueueflipface(triface&, queue*);
  void enqueueflipedge(face&, queue*);
  void flip23(triface*, queue*);
  void flip32(triface*, queue*);
  void flip22(triface*, queue*);
  void flip22sub(face*, queue*);
  long lawson3d(queue* flipqueue);
  long lawson(queue* flipqueue);

  bool removetetbypeeloff(triface *striptet, triface*);
  bool removefacebyflip23(REAL *key, triface*, triface*, queue*);
  bool removeedgebyflip22(REAL *key, int, triface*, queue*);
  bool removeedgebyflip32(REAL *key, triface*, triface*, queue*);
  bool removeedgebytranNM(REAL*,int,triface*,triface*,point,point,queue*);
  bool removeedgebycombNM(REAL*,int,triface*,int*,triface*,triface*,queue*);

  void splittetrahedron(point, triface*, queue*);
  void splittetface(point, triface*, queue*);
  void splitsubface(point, face*, queue*);
  bool splittetedge(point, triface*, queue*);
  void splitsubedge(point, face*, queue*);

  void formstarpolyhedron(point pt, list* tetlist, list* verlist, bool);
  void formbowatcavitysub(point, face*, list*, list*);
  void formbowatcavityquad(point, list*, list*);
  void formbowatcavitysegquad(point, list*, list*);
  void formbowatcavity(point bp, face* bpseg, face* bpsh, int* n, int* nmax,
                       list** sublists, list** subceillists, list** tetlists,
                       list** ceillists);
  void releasebowatcavity(face*, int, list**, list**, list**, list**);
  bool validatebowatcavityquad(point bp, list* ceillist, REAL maxcosd);
  void updatebowatcavityquad(list* tetlist, list* ceillist);
  void updatebowatcavitysub(list* sublist, list* subceillist, int* cutcount);
  bool trimbowatcavity(point bp, face* bpseg, int n, list** sublists,
                       list** subceillists, list** tetlists,list** ceillists,
                       REAL maxcosd);
  void bowatinsertsite(point bp, face* splitseg, int n, list** sublists,
                       list** subceillists, list** tetlists, list** ceillists,
                       list* verlist, queue* flipque, bool chkencseg, 
                       bool chkencsub, bool chkbadtet);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Delaunay tetrahedralization functions                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Point sorting routines.
  void btree_sort(point*, int, int, REAL, REAL, REAL, REAL, REAL, REAL, int);
  void btree_insert(point insertpt);
  void btree_search(point searchpt, triface* searchtet);
  void ordervertices(point* vertexarray, int arraysize);

  enum locateresult insertvertexbw(point insertpt, triface *searchtet, 
                                   bool bwflag, bool visflag, 
                                   bool noencsegflag, bool noencsubflag);
  bool unifypoint(point testpt, triface*, enum locateresult, REAL);
  bool incrflipdelaunay(triface*, point*, long, bool, bool, REAL, queue*);
  long delaunizevertices();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Surface triangulation functions                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  enum locateresult sinsertvertex(point insertpt, face *splitsh,face *splitseg,
                                  bool bwflag, bool cflag);
  void formstarpolygon(point pt, list* trilist, list* verlist);
  void getfacetabovepoint(face* facetsh);
  bool incrflipdelaunaysub(int shmark, REAL eps, list*, int, REAL*, queue*);
  enum finddirectionresult finddirectionsub(face* searchsh, point tend);
  void insertsubseg(face* tri);
  bool scoutsegmentsub(face* searchsh, point tend);
  void flipedgerecursive(face* flipedge, queue* flipqueue);
  void constrainededge(face* startsh, point tend, queue* flipqueue);
  void recoversegment(point tstart, point tend, queue* flipqueue);
  void infecthullsub(memorypool* viri);
  void plaguesub(memorypool* viri);
  void carveholessub(int holes, REAL* holelist, memorypool* viri);
  void triangulate(int shmark, REAL eps, list* ptlist, list* conlist,int holes,
                   REAL* holelist, memorypool* viri, queue*);
  void retrievenewsubs(list* newshlist, bool removeseg);
  void unifysegments();
  void assignsegmentmarkers();
  void mergefacets(queue* flipqueue);
  long meshsurface();

  // Detect intersecting facets of PLC.
  void interecursive(shellface** subfacearray, int arraysize, int axis,
                     REAL bxmin, REAL bxmax, REAL bymin, REAL bymax,
                     REAL bzmin, REAL bzmax, int* internum);
  void detectinterfaces(); 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Constrained Delaunay tetrahedralization functions                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Segment recovery routines.
  void markacutevertices(REAL acuteangle);
  enum finddirectionresult finddirection(triface* searchtet, point, long);
  enum interresult finddirection2(triface* searchtet, point);
  enum interresult finddirection3(triface* searchtet, point);
  enum interresult scoutsegment2(face*, triface*, point*);
  void getsegmentsplitpoint2(face* sseg, point refpt, REAL* vt);
  void getsegmentsplitpoint3(face* sseg, point refpt, REAL* vt);
  void delaunizesegments2();

  // Facets recovery routines.
  enum interresult scoutsubface(face* ssub, triface* searchtet, int);
  enum interresult scoutcrosstet(face* ssub, triface* searchtet, arraypool*);
  void recoversubfacebyflips(face* pssub, triface* crossface, arraypool*);
  void formcavity(face*, arraypool*, arraypool*, arraypool*, arraypool*, 
                  arraypool*, arraypool*, arraypool*);
  bool delaunizecavity(arraypool*, arraypool*, arraypool*, arraypool*,
                       arraypool*, arraypool*);
  bool fillcavity(arraypool*, arraypool*, arraypool*, arraypool*);
  void carvecavity(arraypool*, arraypool*, arraypool*);
  void restorecavity(arraypool*, arraypool*, arraypool*);
  void splitsubedge(point, face*, arraypool*, arraypool*);
  void constrainedfacets2();

  void formskeleton(clock_t&);

  // Carving out holes and concavities routines.
  void infecthull(memorypool *viri);
  void plague(memorypool *viri);
  void regionplague(memorypool *viri, REAL attribute, REAL volume);
  void removeholetets(memorypool *viri);
  void assignregionattribs();
  void carveholes();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Steiner points removal functions                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void initializecavity(list* floorlist, list* ceillist, list* frontlist,
                        list* ptlist, list* gluelist);
  bool delaunizecavvertices(triface*, list*, list*, list*, queue*);
  void retrievenewtets(list* newtetlist);
  void insertauxsubface(triface* front, triface* idfront);
  bool scoutfront(triface* front, triface* idfront);
  void gluefronts(triface* front, triface* front1, list* gluetetlist,
                  list* glueshlist);
  bool identifyfronts(list* frontlist,list* misfrontlist,list* gluetetlist,
                      list* glueshlist);
  void detachauxsubfaces(list* newtetlist);
  bool carvecavity(list* newtetlist, list* outtetlist, list* gluetetlist,
                   queue* flipque);

  void replacepolygonsubs(list* oldshlist, list* newshlist);
  void orientnewsubs(list* newshlist, face* orientsh, REAL* norm);
  bool registerelemflip(enum fliptype ft, point pa1, point pb1, point pc1,
                        point pa2, point pb2, point pc2);
  bool check4fixededge(point pa, point pb);
  bool removeedgebyflips(triface* remedge, int*);
  bool removefacebyflips(triface* remface, int*);
  bool recoveredgebyflips(triface* searchtet, point pb, int*);
  bool recoverfacebyflips(triface* front, int*);
  bool constrainedcavity(triface* oldtet, list* floorlist, list* ceillist,
                         list* ptlist, list* frontlist, list* misfrontlist,
                         list* newtetlist, list* gluetetlist, list* glueshlist,
                         queue* flipque);
  bool findrelocatepoint2(point sp, point np, REAL* n, list*, list*);
  bool relocatepoint(point steinpt, triface* oldtet, list*, list*, queue*);
  bool findcollapseedge(point suppt, point* conpt, list* oldtetlist, list*);
  void collapseedge(point suppt, point conpt, list* oldtetlist, list*);
  void deallocfaketets(list* frontlist);
  void restorepolyhedron(list* oldtetlist);
  bool suppressfacetpoint(face* supsh, list* frontlist, list* misfrontlist,
                          list* ptlist, list* conlist, memorypool* viri,
                          queue* flipque, bool noreloc, bool optflag);
  bool suppresssegpoint(face* supseg, list* spinshlist, list* newsegshlist,
                        list* frontlist, list* misfrontlist, list* ptlist,
                        list* conlist, memorypool* viri, queue* flipque,
                        bool noreloc, bool optflag);
  bool suppressvolpoint(triface* suptet, list* frontlist, list* misfrontlist,
                        list* ptlist, queue* flipque, bool optflag);
  void removesteiners2();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh rebuild functions                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void transfernodes();
  long reconstructmesh();
  void insertconstrainedpoints(tetgenio *addio);
  bool p1interpolatebgm(point pt, triface* bgmtet, long *scount);
  void interpolatesizemap();
  void duplicatebgmesh();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh refinement functions                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  void marksharpsegments(REAL sharpangle);
  void decidefeaturepointsizes();
  void enqueueencsub(face* ss, point encpt, int quenumber, REAL* cent);
  badface* dequeueencsub(int* quenumber);
  void enqueuebadtet(triface* tt, REAL key, REAL* cent);
  badface* topbadtetra();
  void dequeuebadtet();
  bool checkseg4encroach(face* testseg, point testpt, point*, bool enqflag);
  bool checksub4encroach(face* testsub, point testpt, bool enqflag);
  bool checktet4badqual(triface* testtet, bool enqflag);
  bool acceptsegpt(point segpt, point refpt, face* splitseg);
  bool acceptfacpt(point facpt, list* subceillist, list* verlist);
  bool acceptvolpt(point volpt, list* ceillist, list* verlist);
  void getsplitpoint(point e1, point e2, point refpt, point newpt);
  void setnewpointsize(point newpt, point e1, point e2);
  bool splitencseg(point, face*, list*, list*, list*,queue*,bool,bool,bool);
  bool tallencsegs(point testpt, int n, list** ceillists);
  bool tallencsubs(point testpt, int n, list** ceillists);
  void tallbadtetrahedrons();
  void repairencsegs(bool chkencsub, bool chkbadtet);
  void repairencsubs(bool chkbadtet);
  void repairbadtets();
  void enforcequality();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh optimization routines                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  bool checktet4ill(triface* testtet, bool enqflag);
  bool checktet4opt(triface* testtet, bool enqflag);
  bool removeedge(badface* remedge, bool optflag);
  bool smoothpoint(point smthpt, point, point, list*, bool, REAL*);
  bool smoothsliver(badface* remedge, list *starlist);
  bool splitsliver(badface* remedge, list *tetlist, list *ceillist);
  void tallslivers(bool optflag);
  void optimizemesh2(bool optflag);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh output functions                                                     //
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
  void outmesh2gid(char*);
  void outmesh2off(char*);
  void outmesh2vtk(char*);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh check functions                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  int checkmesh();
  int checkshells();
  int checksegments();
  int checkdelaunay(REAL, queue*);
  void checkconforming();
  void algorithmicstatistics();
  void qualitystatistics();
  void statistics();

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Debug functions                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
  /*
  void ptet(triface* t);
  void psh(face* s);
  int pteti(int i, int j, int k, int l);
  void pface(int i, int j, int k);
  bool pedge(int i, int j);
  int psubface(int i, int j, int k);
  void psubseg(int i, int j);
  int pmark(point p);
  void pvert(point p);
  int pverti(int i);
  REAL test_orient3d(int i, int j, int k, int l);
  REAL test_insphere(int i, int j, int k, int l, int m);
  REAL test_insphere_s(int i, int j, int k, int l, int m);
  void print_tetarray(arraypool* tetarray);
  void print_tetlist(list* tetlist);
  void print_facearray(arraypool* facearray);
  void print_facelist(list* facelist);
  void print_subfacearray(arraypool* subfacearray);
  void print_subfacelist(list* subfacelist);
  void dump_facetof(face* pssub);
  void print_fliptetlist(triface *fliptet);
  void print_deaditemstack(void* deaditemstack);
  int check_deaditemstack(void* deaditemstack, uintptr_t addr);
  void print_abtetlist(triface *abtetlist, int len);
  int checkpoint2tetmap();
  int checkpoint2submap();
  int checkpoint2segmap();
  */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class variables                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  // Pointer to the input data (a set of nodes, a PLC, or a mesh).
  tetgenio *in;

  // Pointer to the options (and filenames).
  tetgenbehavior *b;

  // Pointer to a background mesh (contains size specification map).
  tetgenmesh *bgm;

  // Variables used to allocate and access memory for tetrahedra, subfaces
  //   subsegments, points, encroached subfaces, encroached subsegments,
  //   bad-quality tetrahedra, and so on.
  memorypool *tetrahedrons;
  memorypool *subfaces;
  memorypool *subsegs;
  memorypool *points;
  memorypool *badsubsegs;
  memorypool *badsubfaces;
  memorypool *badtetrahedrons;
  memorypool *tet2segpool, *tet2subpool;

  // Pointer to the 'tetrahedron' that occupies all of "outer space".
  tetrahedron *dummytet;
  tetrahedron *dummytetbase; // Keep base address so we can free it later.

  // Pointer to the omnipresent subface.  Referenced by any tetrahedron,
  //   or subface that isn't connected to a subface at that location.
  shellface *dummysh;
  shellface *dummyshbase;    // Keep base address so we can free it later.

  // Entry to find the binary tree nodes (-u option).
  arraypool *btreenode_list;
  // The maximum size of a btree node (number after -u option) is
  int max_btreenode_size; // <= b->max_btreenode_size.
  // The maximum btree depth (for bookkeeping).
  int max_btree_depth; 

  // Arrays used by Bowyer-Watson algorithm.
  arraypool *cavetetlist, *cavebdrylist, *caveoldtetlist;
  arraypool *caveshlist, *caveshbdlist;
  // Stacks used by the boundary recovery algorithm.
  arraypool *subsegstack, *subfacstack;

  // Two handles used in constrained facet recovery.
  triface firsttopface, firstbotface;

  // An array for registering elementary flips.
  arraypool *elemfliplist;

  // An array of fixed edges for facet recovering by flips.
  arraypool *fixededgelist;

  // A point above the plane in which the facet currently being used lies.
  //   It is used as a reference point for orient3d().
  point *facetabovepointarray, abovepoint, dummypoint;

  // Array (size = numberoftetrahedra * 6) for storing high-order nodes of
  //   tetrahedra (only used when -o2 switch is selected).
  point *highordertable;

  // Arrays for storing and searching pbc data. 'subpbcgrouptable', (size
  //   is numberofpbcgroups) for pbcgroup of subfaces. 'segpbcgrouptable',
  //   a list for pbcgroup of segments. Because a segment can have several
  //   pbcgroup incident on it, its size is unknown on input, it will be
  //   found in 'createsegpbcgrouptable()'.
  pbcdata *subpbcgrouptable;
  list *segpbcgrouptable;
  // A map for searching the pbcgroups of a given segment. 'idx2segpglist'
  //   (size = number of input segments + 1), and 'segpglist'.  
  int *idx2segpglist, *segpglist;

  // Queues that maintain the bad (badly-shaped or too large) tetrahedra.
  //   The tails are pointers to the pointers that have to be filled in to
  //   enqueue an item.  The queues are ordered from 63 (highest priority)
  //   to 0 (lowest priority).
  badface *subquefront[3], **subquetail[3];
  badface *tetquefront[64], *tetquetail[64];
  int nextnonemptyq[64];
  int firstnonemptyq, recentq;

  // Pointer to a recently visited tetrahedron. Improves point location
  //   if proximate points are inserted sequentially.
  triface recenttet;

  REAL xmax, xmin, ymax, ymin, zmax, zmin;         // Bounding box of points.
  REAL longest;                          // The longest possible edge length.
  REAL lengthlimit;                     // The limiting length of a new edge.
  long hullsize;                           // Number of faces of convex hull.
  long insegments;                               // Number of input segments.
  long meshedges;                             // Number of output mesh edges.
  int steinerleft;                  // Number of Steiner points not yet used.
  int sizeoftensor;                     // Number of REALs per metric tensor.
  int pointmtrindex;           // Index to find the metric tensor of a point.
  int point2simindex;         // Index to find a simplex adjacent to a point.
  int pointmarkindex;            // Index to find boundary marker of a point.
  int point2pbcptindex;              // Index to find a pbc point to a point.
  int highorderindex;    // Index to find extra nodes for highorder elements.
  int elemattribindex;          // Index to find attributes of a tetrahedron.
  int volumeboundindex;       // Index to find volume bound of a tetrahedron.
  int elemmarkerindex;              // Index to find marker of a tetrahedron.
  int shmarkindex;             // Index to find boundary marker of a subface.
  int areaboundindex;               // Index to find area bound of a subface.
  int checksubfaces;                   // Are there subfaces in the mesh yet?
  int checksubsegs;                     // Are there subsegs in the mesh yet?
  int checkpbcs;                   // Are there periodic boundary conditions?
  int varconstraint;     // Are there variant (node, seg, facet) constraints?
  int nonconvex;                               // Is current mesh non-convex?
  int dupverts;                             // Are there duplicated vertices?
  int unuverts;                                 // Are there unused vertices?
  int relverts;                          // The number of relocated vertices.
  int suprelverts;            // The number of suppressed relocated vertices.
  int collapverts;             // The number of collapsed relocated vertices.
  int unsupverts;                     // The number of unsuppressed vertices.
  int smoothsegverts;                     // The number of smoothed vertices.
  int jettisoninverts;            // The number of jettisoned input vertices.
  long samples;               // Number of random samples for point location.
  unsigned long randomseed;                    // Current random number seed.
  REAL macheps;                                       // The machine epsilon.
  REAL cosmaxdihed, cosmindihed;    // The cosine values of max/min dihedral.
  REAL minfaceang, minfacetdihed;     // The minimum input (dihedral) angles.
  int maxcavfaces, maxcavverts;            // The size of the largest cavity.
  bool b_steinerflag;

  // Algorithm statistical counters.
  long ptloc_count, ptloc_max_count;
  long orient3dcount;
  long inspherecount, insphere_sos_count;
  long flip14count, flip26count, flipn2ncount;
  long flip22count;
  long inserthullcount;
  long maxbowatcavsize, totalbowatcavsize, totaldeadtets;
  long across_face_count, across_edge_count, across_max_count;
  long maxcavsize, maxregionsize;
  long ndelaunayedgecount, cavityexpcount;
  long opt_tet_peels, opt_face_flips, opt_edge_flips;

  long abovecount;                     // Number of abovepoints calculation.
  long bowatvolcount, bowatsubcount, bowatsegcount;       // Bowyer-Watsons.
  long updvolcount, updsubcount, updsegcount;   // Bow-Wat cavities updates.
  long failvolcount, failsubcount, failsegcount;           // Bow-Wat fails.
  long outbowatcircumcount;    // Number of circumcenters outside Bowat-cav.
  long r1count, r2count, r3count;        // Numbers of edge splitting rules.
  long cdtenforcesegpts;                // Number of CDT enforcement points.
  long rejsegpts, rejsubpts, rejtetpts;        // Number of rejected points.
  long optcount[10];            // Numbers of various optimizing operations.
  long flip23s, flip32s, flip22s, flip44s;     // Number of flips performed.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class constructor & destructor                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

  tetgenmesh()
  {
    bgm = (tetgenmesh *) NULL;
    in = (tetgenio *) NULL;
    b = (tetgenbehavior *) NULL;

    tetrahedrons = (memorypool *) NULL;
    subfaces = (memorypool *) NULL;
    subsegs = (memorypool *) NULL;
    points = (memorypool *) NULL;
    badsubsegs = (memorypool *) NULL;
    badsubfaces = (memorypool *) NULL;
    badtetrahedrons = (memorypool *) NULL;
    tet2segpool = NULL;
    tet2subpool = NULL;

    dummytet = (tetrahedron *) NULL;
    dummytetbase = (tetrahedron *) NULL;
    dummysh = (shellface *) NULL;
    dummyshbase = (shellface *) NULL;

    facetabovepointarray = (point *) NULL;
    abovepoint = (point) NULL;
    dummypoint = NULL;
    btreenode_list = (arraypool *) NULL;
    highordertable = (point *) NULL;
    subpbcgrouptable = (pbcdata *) NULL;
    segpbcgrouptable = (list *) NULL;
    idx2segpglist = (int *) NULL;
    segpglist = (int *) NULL;

    cavetetlist = NULL;
    cavebdrylist = NULL;
    caveoldtetlist = NULL;
    caveshlist = caveshbdlist = NULL;
    subsegstack = subfacstack = NULL;

    elemfliplist = (arraypool *) NULL;
    fixededgelist = (arraypool *) NULL;

    xmax = xmin = ymax = ymin = zmax = zmin = 0.0; 
    longest = 0.0;
    hullsize = 0l;
    insegments = 0l;
    meshedges = 0l;
    pointmtrindex = 0;
    pointmarkindex = 0;
    point2simindex = 0;
    point2pbcptindex = 0;
    highorderindex = 0;
    elemattribindex = 0;
    volumeboundindex = 0;
    shmarkindex = 0;
    areaboundindex = 0;
    checksubfaces = 0;
    checksubsegs = 0;
    checkpbcs = 0;
    varconstraint = 0;
    nonconvex = 0;
    dupverts = 0;
    unuverts = 0;
    relverts = 0;
    suprelverts = 0;
    collapverts = 0;
    unsupverts = 0;
    jettisoninverts = 0;
    samples = 0l;
    randomseed = 1l;
    macheps = 0.0;
    minfaceang = minfacetdihed = PI;
    b_steinerflag = false;

    ptloc_count = ptloc_max_count = 0l;
    orient3dcount = 0l;
    inspherecount = insphere_sos_count = 0l;
    flip14count = flip26count = flipn2ncount = 0l;
    flip22count = 0l;
    inserthullcount = 0l;
    maxbowatcavsize = totalbowatcavsize = totaldeadtets = 0l;
    across_face_count = across_edge_count = across_max_count = 0l;
    maxcavsize = maxregionsize = 0l;
    ndelaunayedgecount = cavityexpcount = 0l;
    opt_tet_peels = opt_face_flips = opt_edge_flips = 0l;

    maxcavfaces = maxcavverts = 0;
    abovecount = 0l;
    bowatvolcount = bowatsubcount = bowatsegcount = 0l;
    updvolcount = updsubcount = updsegcount = 0l;
    outbowatcircumcount = 0l;
    failvolcount = failsubcount = failsegcount = 0l;
    r1count = r2count = r3count = 0l;
    cdtenforcesegpts = 0l;
    rejsegpts = rejsubpts = rejtetpts = 0l;
    flip23s = flip32s = flip22s = flip44s = 0l;
  } // tetgenmesh()
    
  ~tetgenmesh()
  {
    bgm = (tetgenmesh *) NULL;
    in = (tetgenio *) NULL;
    b = (tetgenbehavior *) NULL;

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
    if (dummytetbase != (tetrahedron *) NULL) {
      delete [] dummytetbase;
    }
    if (dummyshbase != (shellface *) NULL) {
      delete [] dummyshbase;
    }
    if (facetabovepointarray != (point *) NULL) {
      delete [] facetabovepointarray;
    }
    if (dummypoint != NULL) {
      delete [] dummypoint;
    }
    if (highordertable != (point *) NULL) {
      delete [] highordertable;
    }
    if (subpbcgrouptable != (pbcdata *) NULL) {
      delete [] subpbcgrouptable;
    }
    if (segpbcgrouptable != (list *) NULL) {
      delete segpbcgrouptable;
      delete [] idx2segpglist;
      delete [] segpglist;
    }

    if (cavetetlist != NULL) {
      delete cavetetlist;
      delete cavebdrylist;
      delete caveoldtetlist;
    }
    if (subsegstack != NULL) {
      delete subsegstack;
    }
    if (subfacstack != NULL) {
      delete subfacstack;
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
    printf("  Please report this bug to sihang@mail.berlios.de. Include\n");
    printf("    the message above, your input data set, and the exact\n");
    printf("     command line you used to run this program, thank you.\n");
    break;
  default:
    printf("Program stopped.\n"); 
  } // switch (x)
  exit(x);
#endif // #ifdef TETLIBRARY
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Geometric predicates                                                      //
//                                                                           //
// Return one of the values +1, 0, and -1 on basic geometric questions such  //
// as the orientation of point sets, in-circle, and in-sphere tests.  They   //
// are basic units for implmenting geometric algorithms.  TetGen uses two 3D //
// geometric predicates: the orientation and in-sphere tests.                //
//                                                                           //
// Orientation test:  let a, b, c be a sequence of 3 non-collinear points in //
// R^3.  They defines a unique hypeplane H.  Let H+ and H- be the two spaces //
// separated by H, which are defined as follows (using the left-hand rule):  //
// make a fist using your left hand in such a way that your fingers follow   //
// the order of a, b and c, then your thumb is pointing to H+.  Given any    //
// point d in R^3, the orientation test returns +1 if d lies in H+, -1 if d  //
// lies in H-, or 0 if d lies on H.                                          //
//                                                                           //
// In-sphere test:  let a, b, c, d be 4 non-coplanar points in R^3.  They    //
// defines a unique circumsphere S.  Given any point e in R^3, the in-sphere //
// test returns +1 if e lies inside S, or -1 if e lies outside S, or 0 if e  //
// lies on S.                                                                //
//                                                                           //
// The following routines use arbitrary precision floating-point arithmetic. //
// They are provided by J. R. Schewchuk in public domain (http://www.cs.cmu. //
// edu/~quake/robust.html). The source code are in "predicates.cxx".         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL exactinit();
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Inline functions of mesh data structures                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Some macros for convenience

#define Div2  >> 1
#define Mod2  & 01

// NOTE: These bit operators should only be used in macros below.

// Get orient(Range from 0 to 2) from face version(Range from 0 to 5).

#define Orient(V)   ((V) Div2)

// Determine edge ring(0 or 1) from face version(Range from 0 to 5).

#define EdgeRing(V) ((V) Mod2)

//
// Begin of primitives for tetrahedra
// 

// Each tetrahedron contains four pointers to its neighboring tetrahedra,
//   with face indices.  To save memory, both information are kept in a
//   single pointer. To make this possible, all tetrahedra are aligned to
//   eight-byte boundaries, so that the last three bits of each pointer are
//   zeros. A face index (in the range 0 to 3) is compressed into the last
//   two bits of each pointer by the function 'encode()'.  The function
//   'decode()' decodes a pointer, extracting a face index and a pointer to
//   the beginning of a tetrahedron.

inline void tetgenmesh::decode(tetrahedron ptr, triface& t) {
  t.loc = (int) ((uintptr_t) (ptr) & (uintptr_t) 3);
  t.tet = (tetrahedron *) ((uintptr_t) (ptr) & ~(uintptr_t) 7);
}

inline tetgenmesh::tetrahedron tetgenmesh::encode(triface& t) {
  return (tetrahedron) ((uintptr_t) t.tet | (uintptr_t) t.loc);
}

// sym() finds the abutting tetrahedron on the same face.

inline void tetgenmesh::sym(triface& t1, triface& t2) {
  tetrahedron ptr = t1.tet[t1.loc];
  decode(ptr, t2);
}

inline void tetgenmesh::symself(triface& t) {
  tetrahedron ptr = t.tet[t.loc];
  decode(ptr, t);
}

// Bond two tetrahedra together at their faces.

inline void tetgenmesh::bond(triface& t1, triface& t2) {
  t1.tet[t1.loc] = encode(t2);
  t2.tet[t2.loc] = encode(t1);
}

// Dissolve a bond (from one side).  Note that the other tetrahedron will
//   still think it is connected to this tetrahedron.  Usually, however,
//   the other tetrahedron is being deleted entirely, or bonded to another
//   tetrahedron, so it doesn't matter.

inline void tetgenmesh::dissolve(triface& t) {
  t.tet[t.loc] = (tetrahedron) dummytet;
}

// These primitives determine or set the origin, destination, apex or
//   opposition of a tetrahedron with respect to 'loc' and 'ver'.

inline tetgenmesh::point tetgenmesh::org(triface& t) {
  return (point) t.tet[locver2org[t.loc][t.ver] + 4];
}

inline tetgenmesh::point tetgenmesh::dest(triface& t) {
  return (point) t.tet[locver2dest[t.loc][t.ver] + 4];
}

inline tetgenmesh::point tetgenmesh::apex(triface& t) {
  return (point) t.tet[locver2apex[t.loc][t.ver] + 4];
}

inline tetgenmesh::point tetgenmesh::oppo(triface& t) {
  return (point) t.tet[loc2oppo[t.loc] + 4];
}

inline void tetgenmesh::setorg(triface& t, point pointptr) {
  t.tet[locver2org[t.loc][t.ver] + 4] = (tetrahedron) pointptr;
}

inline void tetgenmesh::setdest(triface& t, point pointptr) {
  t.tet[locver2dest[t.loc][t.ver] + 4] = (tetrahedron) pointptr;
}

inline void tetgenmesh::setapex(triface& t, point pointptr) {
  t.tet[locver2apex[t.loc][t.ver] + 4] = (tetrahedron) pointptr;
}

inline void tetgenmesh::setoppo(triface& t, point pointptr) {
  t.tet[loc2oppo[t.loc] + 4] = (tetrahedron) pointptr;
}

// These primitives were drived from Mucke's triangle-edge data structure
//   to change face-edge relation in a tetrahedron (esym, enext and enext2)
//   or between two tetrahedra (fnext).

// If e0 = e(i, j), e1 = e(j, i), that is e0 and e1 are the two directions
//   of the same undirected edge of a face. e0.sym() = e1 and vice versa.

inline void tetgenmesh::esym(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.loc = t1.loc;
  t2.ver = t1.ver + (EdgeRing(t1.ver) ? -1 : 1);
}

inline void tetgenmesh::esymself(triface& t) {
  t.ver += (EdgeRing(t.ver) ? -1 : 1);
}

// If e0 and e1 are both in the same edge ring of a face, e1 = e0.enext().

inline void tetgenmesh::enext(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.loc = t1.loc;
  t2.ver = ve[t1.ver];
}

inline void tetgenmesh::enextself(triface& t) {
  t.ver = ve[t.ver];
}

// enext2() is equal to e2 = e0.enext().enext()

inline void tetgenmesh::enext2(triface& t1, triface& t2) {
  t2.tet = t1.tet;
  t2.loc = t1.loc;
  t2.ver = ve[ve[t1.ver]];
}

inline void tetgenmesh::enext2self(triface& t) {
  t.ver = ve[ve[t.ver]];
}

// If f0 and f1 are both in the same face ring of a face, f1 = f0.fnext().
//   If f1 exists, return true. Otherwise, return false, i.e., f0 is a
//   boundary or hull face.

inline bool tetgenmesh::fnext(triface& t1, triface& t2) 
{
  // Get the next face.
  t2.loc = locver2nextf[t1.loc][t1.ver][0];
  // Is the next face in the same tet?
  if (t2.loc != -1) {
    // It's in the same tet. Get the edge version.
    t2.ver = locver2nextf[t1.loc][t1.ver][1];
    t2.tet = t1.tet;
  } else {
    // The next face is in the neigbhour of 't1'.
    sym(t1, t2);
    if (t2.tet != dummytet) {
      // Find the corresponding edge in t2.
      point torg;
      int tloc, tver, i;
      t2.ver = 0;
      torg = org(t1);
      for (i = 0; (i < 3) && (org(t2) != torg); i++) {
        enextself(t2);
      }
      // Go to the next face in t2.
      tloc = t2.loc;
      tver = t2.ver;
      t2.loc = locver2nextf[tloc][tver][0];
      t2.ver = locver2nextf[tloc][tver][1];
    }
  }
  return t2.tet != dummytet;
}

inline bool tetgenmesh::fnextself(triface& t1) 
{
  triface t2;

  // Get the next face.
  t2.loc = locver2nextf[t1.loc][t1.ver][0];
  // Is the next face in the same tet?
  if (t2.loc != -1) {
    // It's in the same tet. Get the edge version.
    t2.ver = locver2nextf[t1.loc][t1.ver][1];
    t1.loc = t2.loc;
    t1.ver = t2.ver;
  } else {
    // The next face is in the neigbhour of 't1'.
    sym(t1, t2);
    if (t2.tet != dummytet) {
      // Find the corresponding edge in t2.
      point torg;
      int i;
      t2.ver = 0;
      torg = org(t1);
      for (i = 0; (i < 3) && (org(t2) != torg); i++) {
        enextself(t2);
      }
      t1.loc = locver2nextf[t2.loc][t2.ver][0];
      t1.ver = locver2nextf[t2.loc][t2.ver][1];
      t1.tet = t2.tet;
    }
  }
  return t2.tet != dummytet;
}

// Given a face t1, find the face f2 in the adjacent tet. If t2 is not
//   a dummytet, then t1 and t2 refer to the same edge. Moreover, t2's
//   edge must be in 0th edge ring, e.g., t2.ver is one of {0, 2, 4}.
//   No matter what edge version t1 is.

inline void tetgenmesh::symedge(triface& t1, triface& t2)
{
  decode(t1.tet[t1.loc], t2);
  if (t2.tet != dummytet) {
    // Search the edge of t1 in t2.
    point tapex = apex(t1);
    if ((point) (t2.tet[locver2apex[t2.loc][0] + 4]) == tapex) {
      t2.ver = 0;
    } else if ((point) (t2.tet[locver2apex[t2.loc][2] + 4]) == tapex) {
      t2.ver = 2;
    } else {
      assert((point) (t2.tet[locver2apex[t2.loc][4] + 4]) == tapex);
      t2.ver = 4;
    }
  }
}

inline void tetgenmesh::symedgeself(triface& t)
{
  tetrahedron ptr;
  point tapex;

  ptr = t.tet[t.loc];
  tapex = apex(t);

  decode(ptr, t);
  if (t.tet != dummytet) {
    // Search the edge of t1 in t2.
    if ((point) (t.tet[locver2apex[t.loc][0] + 4]) == tapex) {
      t.ver = 0;
    } else if ((point) (t.tet[locver2apex[t.loc][2] + 4]) == tapex) {
      t.ver = 2;
    } else {
      assert((point) (t.tet[locver2apex[t.loc][4] + 4]) == tapex);
      t.ver = 4;
    }
  }
}

// Given a face t1, find the next face t2 in the face ring, t1 and t2
//   are in two different tetrahedra. If the next face is a hull face,
//   t2 is dummytet.

inline void tetgenmesh::tfnext(triface& t1, triface& t2)
{
  int *iptr;

  if ((t1.ver & 1) == 0) {
    t2.tet = t1.tet;
    iptr = locver2nextf[t1.loc][t1.ver];
    t2.loc = iptr[0];
    t2.ver = iptr[1];
    symedgeself(t2);  // t2.tet may be dummytet.
  } else {
    symedge(t1, t2);
    if (t2.tet != dummytet) {
      iptr = locver2nextf[t2.loc][t2.ver];
      t2.loc = iptr[0];
      t2.ver = iptr[1];
    }
  }
}

inline void tetgenmesh::tfnextself(triface& t)
{
  int *iptr;

  if ((t.ver & 1) == 0) {
    iptr = locver2nextf[t.loc][t.ver];
    t.loc = iptr[0];
    t.ver = iptr[1];
    symedgeself(t); // t.tet may be dummytet.
  } else {
    symedgeself(t);
    if (t.tet != dummytet) {
      iptr = locver2nextf[t.loc][t.ver];
      t.loc = iptr[0];
      t.ver = iptr[1];
    }
  }
}

// enextfnext() and enext2fnext() are combination primitives of enext(),
//   enext2() and fnext().

inline void tetgenmesh::enextfnext(triface& t1, triface& t2) {
  enext(t1, t2);
  fnextself(t2);
}

inline void tetgenmesh::enextfnextself(triface& t) {
  enextself(t);
  fnextself(t);
}

inline void tetgenmesh::enext2fnext(triface& t1, triface& t2) {
  enext2(t1, t2);
  fnextself(t2);
}

inline void tetgenmesh::enext2fnextself(triface& t) {
  enext2self(t);
  fnextself(t);
}

// Check or set a tetrahedron's attributes.

inline REAL tetgenmesh::elemattribute(tetrahedron* ptr, int attnum) {
  return ((REAL *) (ptr))[elemattribindex + attnum];
}

inline void tetgenmesh::
setelemattribute(tetrahedron* ptr, int attnum, REAL value){
  ((REAL *) (ptr))[elemattribindex + attnum] = value;
}

// Check or set a tetrahedron's maximum volume bound.

inline REAL tetgenmesh::volumebound(tetrahedron* ptr) {
  return ((REAL *) (ptr))[volumeboundindex];
}

inline void tetgenmesh::setvolumebound(tetrahedron* ptr, REAL value) {
  ((REAL *) (ptr))[volumeboundindex] = value;
}

// Check or set a tetrahedron's marker.

inline int tetgenmesh::getelemmarker(tetrahedron* ptr) {
  return ((int *) (ptr))[elemmarkerindex];
}

inline void tetgenmesh::setelemmarker(tetrahedron* ptr, int value) {
  ((int *) (ptr))[elemmarkerindex] = value;
}

// infect(), infected(), uninfect() -- primitives to flag or unflag a
//   tetrahedron. The last bit of the element marker is flagged (1)
//   or unflagged (0).

inline void tetgenmesh::infect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) 1;
}

inline void tetgenmesh::uninfect(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) 1;
}

// Test a tetrahedron for viral infection.

inline bool tetgenmesh::infected(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (int) 1) != 0;
}

// marktest(), marktested(), unmarktest() -- primitives to flag or unflag a
//   tetrahedron.  The last second bit of the element marker is marked (1)
//   or unmarked (0).
// One needs them in forming Bowyer-Watson cavity, to mark a tetrahedron if
//   it has been checked (for Delaunay case) so later check can be avoided.

inline void tetgenmesh::marktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) 2;
}

inline void tetgenmesh::unmarktest(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) 2;
}
    
inline bool tetgenmesh::marktested(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (int) 2) != 0;
}

// markface(), unmarkface(), facemarked() -- primitives to flag or unflag a
//   face of a tetrahedron.  From the last 3rd to 6th bits are used for
//   face markers, e.g., the last third bit corresponds to loc = 0. 
// One use of the face marker is in flip algorithm. Each queued face (check
//   for locally Delaunay) is marked.

inline void tetgenmesh::markface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= (int) (4<<(t).loc);
}

inline void tetgenmesh::unmarkface(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= ~(int) (4<<(t).loc);
}

inline bool tetgenmesh::facemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & (int) (4<<(t).loc)) != 0;
}

// markedge(), unmarkedge(), edgemarked() -- primitives to flag or unflag an
//   edge of a tetrahedron.  From the last 7th to 12th bits are used for
//   edge markers, e.g., the last 7th bit corresponds to the 0th edge, etc. 
// Remark: The last 7th bit is marked by 2^6 = 64.

inline void tetgenmesh::markedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] |= 
    (int) (64<<locver2edge[(t).loc][(t).ver]);
}

inline void tetgenmesh::unmarkedge(triface& t) {
  ((int *) (t.tet))[elemmarkerindex] &= 
    ~(int) (64<<locver2edge[(t).loc][(t).ver]);
}

inline bool tetgenmesh::edgemarked(triface& t) {
  return (((int *) (t.tet))[elemmarkerindex] & 
            (int) (64<<locver2edge[(t).loc][(t).ver])) != 0;
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
  s.sh = (shellface *) ((uintptr_t) (sptr) & ~ (uintptr_t) 7);
}

inline tetgenmesh::shellface tetgenmesh::sencode(face& s) {
  return (shellface) ((uintptr_t) s.sh | (uintptr_t) s.shver);
}

// spivot() finds the other subface (from this subface) that shares the
//   same edge.

inline void tetgenmesh::spivot(face& s1, face& s2) {
  shellface sptr = s1.sh[Orient(s1.shver)];
  sdecode(sptr, s2);
}

inline void tetgenmesh::spivotself(face& s) {
  shellface sptr = s.sh[Orient(s.shver)];
  sdecode(sptr, s);
}

// sbond() bonds two subfaces together, i.e., after bonding, both faces
//   are pointing to each other.

inline void tetgenmesh::sbond(face& s1, face& s2) {
  s1.sh[Orient(s1.shver)] = sencode(s2);
  s2.sh[Orient(s2.shver)] = sencode(s1);
}

// sbond1() only bonds s2 to s1, i.e., after bonding, s1 is pointing to s2,
//   but s2 is not pointing to s1.

inline void tetgenmesh::sbond1(face& s1, face& s2) {
  s1.sh[Orient(s1.shver)] = sencode(s2);
}

// Dissolve a subface bond (from one side).  Note that the other subface
//   will still think it's connected to this subface.

inline void tetgenmesh::sdissolve(face& s) {
  s.sh[Orient(s.shver)] = (shellface) dummysh;
}

// These primitives determine or set the origin, destination, or apex
//   of a subface with respect to the edge version.

inline tetgenmesh::point tetgenmesh::sorg(face& s) {
  return (point) s.sh[3 + vo[s.shver]];
}

inline tetgenmesh::point tetgenmesh::sdest(face& s) {
  return (point) s.sh[3 + vd[s.shver]];
}

inline tetgenmesh::point tetgenmesh::sapex(face& s) {
  return (point) s.sh[3 + va[s.shver]];
}

inline void tetgenmesh::setsorg(face& s, point pointptr) {
  s.sh[3 + vo[s.shver]] = (shellface) pointptr;
}

inline void tetgenmesh::setsdest(face& s, point pointptr) {
  s.sh[3 + vd[s.shver]] = (shellface) pointptr;
}

inline void tetgenmesh::setsapex(face& s, point pointptr) {
  s.sh[3 + va[s.shver]] = (shellface) pointptr;
}

// These primitives were drived from Mucke[2]'s triangle-edge data structure
//   to change face-edge relation in a subface (sesym, senext and senext2).

inline void tetgenmesh::sesym(face& s1, face& s2) {
  s2.sh = s1.sh;
  s2.shver = s1.shver + (EdgeRing(s1.shver) ? -1 : 1);
}

inline void tetgenmesh::sesymself(face& s) {
  s.shver += (EdgeRing(s.shver) ? -1 : 1);
}

inline void tetgenmesh::senext(face& s1, face& s2) {
  s2.sh = s1.sh;
  s2.shver = ve[s1.shver];
}

inline void tetgenmesh::senextself(face& s) { 
  s.shver = ve[s.shver]; 
}

inline void tetgenmesh::senext2(face& s1, face& s2) {
  s2.sh = s1.sh;
  s2.shver = ve[ve[s1.shver]];
}

inline void tetgenmesh::senext2self(face& s) {
  s.shver = ve[ve[s.shver]];
}

// If f0 and f1 are both in the same face ring, then f1 = f0.fnext(),

inline void tetgenmesh::sfnext(face& s1, face& s2) {
  getnextsface(&s1, &s2);
}

inline void tetgenmesh::sfnextself(face& s) {
  getnextsface(&s, NULL);
}

// These primitives read or set a pointer of the badface structure.  The
//   pointer is stored sh[11].

inline tetgenmesh::badface* tetgenmesh::shell2badface(face& s) {
  return (badface*) s.sh[11];
}

inline void tetgenmesh::setshell2badface(face& s, badface* value) {
  s.sh[11] = (shellface) value;
}

// Check or set a subface's maximum area bound.

inline REAL tetgenmesh::areabound(face& s) {
  return ((REAL *) (s.sh))[areaboundindex];
}

inline void tetgenmesh::setareabound(face& s, REAL value) {
  ((REAL *) (s.sh))[areaboundindex] = value;
}

// These two primitives read or set a shell marker.  Shell markers are used
//   to hold user boundary information.
// The last two bits of the int ((int *) ((s).sh))[shmarkindex] are used
//   by sinfect() and smarktest().

inline int tetgenmesh::shellmark(face& s) {
  return (((int *) ((s).sh))[shmarkindex]) >> (int) 2; 
  // return ((int *) (s.sh))[shmarkindex];
}

inline void tetgenmesh::setshellmark(face& s, int value) {
  ((int *) ((s).sh))[shmarkindex] = (value << (int) 2) + 
    ((((int *) ((s).sh))[shmarkindex]) & (int) 3);
  // ((int *) (s.sh))[shmarkindex] = value;
}

// These two primitives set or read the type of the subface or subsegment.

inline enum tetgenmesh::shestype tetgenmesh::shelltype(face& s) {
  return (enum shestype) ((int *) (s.sh))[shmarkindex + 1];
}

inline void tetgenmesh::setshelltype(face& s, enum shestype value) {
  ((int *) (s.sh))[shmarkindex + 1] = (int) value;
}

// These two primitives set or read the pbc group of the subface.

inline int tetgenmesh::shellpbcgroup(face& s) {
  return ((int *) (s.sh))[shmarkindex + 2];
}

inline void tetgenmesh::setshellpbcgroup(face& s, int value) {
  ((int *) (s.sh))[shmarkindex + 2] = value;
}

// sinfect(), sinfected(), suninfect() -- primitives to flag or unflag a
//   subface. The last bit of ((int *) ((s).sh))[shmarkindex] is flaged.

inline void tetgenmesh::sinfect(face& s) {
  ((int *) ((s).sh))[shmarkindex] = 
    (((int *) ((s).sh))[shmarkindex] | (int) 1);
  // s.sh[6] = (shellface) ((unsigned long) s.sh[6] | (unsigned long) 4l);
}

inline void tetgenmesh::suninfect(face& s) {
  ((int *) ((s).sh))[shmarkindex] = 
    (((int *) ((s).sh))[shmarkindex] & ~(int) 1);
  // s.sh[6] = (shellface)((unsigned long) s.sh[6] & ~(unsigned long) 4l);
}

// Test a subface for viral infection.

inline bool tetgenmesh::sinfected(face& s) {
  return (((int *) ((s).sh))[shmarkindex] & (int) 1) != 0;
}

// smarktest(), smarktested(), sunmarktest() -- primitives to flag or unflag
//   a subface. The last 2nd bit of ((int *) ((s).sh))[shmarkindex] is flaged.

#define smarktest(s) \
  ((int *) ((s).sh))[shmarkindex] = (((int *)((s).sh))[shmarkindex] | (int) 2)

#define sunmarktest(s) \
  ((int *) ((s).sh))[shmarkindex] = (((int *)((s).sh))[shmarkindex] & ~(int) 2)

#define smarktested(s) ((((int *) ((s).sh))[shmarkindex] & (int) 2) != 0)

//
// End of primitives for subfaces/subsegments
//

//
// Begin of primitives for interacting between tetrahedra and subfaces
//

// tspivot() finds a subface abutting on this tetrahdera.

inline void tetgenmesh::tspivot(triface& t, face& s) {
  if ((t).tet[9] != NULL) {
    sdecode(((shellface *) (t).tet[9])[(t).loc], s);
  } else {
    (s).sh = dummysh;
  }
  //shellface sptr = (shellface) t.tet[8 + t.loc];
  //sdecode(sptr, s);
}

// stpivot() finds a tetrahedron abutting a subface.

inline void tetgenmesh::stpivot(face& s, triface& t) {
  tetrahedron ptr = (tetrahedron) s.sh[6 + EdgeRing(s.shver)];
  decode(ptr, t);
}

// tsbond() bond a tetrahedron to a subface.

inline void tetgenmesh::tsbond(triface& t, face& s) {
  if ((t).tet[9] == NULL) {
    // Allocate space for this tet.
    (t).tet[9] = (tetrahedron) tet2subpool->alloc();
    // NULL all fields in this space.
    for (int i = 0; i < 4; i++) {
      ((shellface *) (t).tet[9])[i] = (shellface) dummysh;
    }
  }
  // Bond t <==> s.
  ((shellface *) (t).tet[9])[(t).loc] = sencode(s);
  //t.tet[8 + t.loc] = (tetrahedron) sencode(s);
  s.sh[6 + EdgeRing(s.shver)] = (shellface) encode(t);
}

// tsdissolve() dissolve a bond (from the tetrahedron side).

inline void tetgenmesh::tsdissolve(triface& t) {
  if ((t).tet[9] != NULL) {
    ((shellface *) (t).tet[9])[(t).loc] = (shellface) dummysh;
  }
  // t.tet[8 + t.loc] = (tetrahedron) dummysh;
}

// stdissolve() dissolve a bond (from the subface side).

inline void tetgenmesh::stdissolve(face& s) {
  s.sh[6 + EdgeRing(s.shver)] = (shellface) dummytet;
}

//
// End of primitives for interacting between tetrahedra and subfaces
//

//
// Begin of primitives for interacting between subfaces and subsegs
//

// sspivot() finds a subsegment abutting a subface.

inline void tetgenmesh::sspivot(face& s, face& edge) {
  shellface sptr = (shellface) s.sh[8 + Orient(s.shver)];
  sdecode(sptr, edge);
}

// ssbond() bond a subface to a subsegment.

inline void tetgenmesh::ssbond(face& s, face& edge) {
  s.sh[8 + Orient(s.shver)] = sencode(edge);
  edge.sh[0] = sencode(s);
}

// ssdisolve() dissolve a bond (from the subface side)

inline void tetgenmesh::ssdissolve(face& s) {
  s.sh[8 + Orient(s.shver)] = (shellface) dummysh;
}

//
// End of primitives for interacting between subfaces and subsegs
//

//
// Begin of primitives for interacting between tet and subsegs.
//

inline void tetgenmesh::tsspivot1(triface& t, face& s)
{
  if ((t).tet[8] != NULL) {
    sdecode(((shellface *) (t).tet[8])[locver2edge[(t).loc][(t).ver]], s);
  } else {
    (s).sh = dummysh;
  }
  // shellface sptr = (shellface) t.tet[8 + locver2edge[t.loc][t.ver]];
  // sdecode(sptr, seg);
}

// Only bond/dissolve at tet's side, but not vice versa.

inline void tetgenmesh::tssbond1(triface& t, face& s)
{
  if ((t).tet[8] == NULL) {
    // Allocate space for this tet.
    (t).tet[8] = (tetrahedron) tet2segpool->alloc();
    // NULL all fields in this space.
    for (int i = 0; i < 6; i++) {
      ((shellface *) (t).tet[8])[i] = (shellface) dummysh;
    }
  }
  // Bond the segment.
  ((shellface *) (t).tet[8])[locver2edge[(t).loc][(t).ver]] = sencode((s));
  // t.tet[8 + locver2edge[t.loc][t.ver]] = (tetrahedron) sencode(seg);
}

inline void tetgenmesh::tssdissolve1(triface& t)
{
  if ((t).tet[8] != NULL) {
    ((shellface *) (t).tet[8])[locver2edge[(t).loc][(t).ver]] 
      = (shellface) dummysh;
  }
  // t.tet[8 + locver2edge[t.loc][t.ver]] = (tetrahedron) dummysh;
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
// The last significant bit of this integer is used by pinfect/puninfect.

inline enum tetgenmesh::verttype tetgenmesh::pointtype(point pt) {
  return (enum verttype) (((int *) (pt))[pointmarkindex + 1] >> (int) 1);
}

inline void tetgenmesh::setpointtype(point pt, enum verttype value) {
  ((int *) (pt))[pointmarkindex + 1] = 
    ((int) value << 1) + (((int *) (pt))[pointmarkindex + 1] & (int) 1);
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

// These following primitives set and read a pointer to a tetrahedron
//   a subface/subsegment, a point, or a tet of background mesh.

inline tetgenmesh::tetrahedron tetgenmesh::point2tet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex];
}

inline void tetgenmesh::setpoint2tet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex] = value;
}

inline tetgenmesh::shellface tetgenmesh::point2sh(point pt) {
  return (shellface) ((tetrahedron *) (pt))[point2simindex + 1];
}

inline void tetgenmesh::setpoint2sh(point pt, shellface value) {
  ((tetrahedron *) (pt))[point2simindex + 1] = (tetrahedron) value;
}

inline tetgenmesh::shellface tetgenmesh::point2seg(point pt) {
  return (shellface) ((tetrahedron *) (pt))[point2simindex + 2];
}

inline void tetgenmesh::setpoint2seg(point pt, shellface value) {
  ((tetrahedron *) (pt))[point2simindex + 2] = (tetrahedron) value;
}

inline tetgenmesh::point tetgenmesh::point2ppt(point pt) {
  return (point) ((tetrahedron *) (pt))[point2simindex + 3];
}

inline void tetgenmesh::setpoint2ppt(point pt, point value) {
  ((tetrahedron *) (pt))[point2simindex + 3] = (tetrahedron) value;
}

inline tetgenmesh::tetrahedron tetgenmesh::point2bgmtet(point pt) {
  return ((tetrahedron *) (pt))[point2simindex + 4];
}

inline void tetgenmesh::setpoint2bgmtet(point pt, tetrahedron value) {
  ((tetrahedron *) (pt))[point2simindex + 4] = value;
}

// These primitives set and read a pointer to its pbc point.

inline tetgenmesh::point tetgenmesh::point2pbcpt(point pt) {
  return (point) ((tetrahedron *) (pt))[point2pbcptindex];
}

inline void tetgenmesh::setpoint2pbcpt(point pt, point value) {
  ((tetrahedron *) (pt))[point2pbcptindex] = (tetrahedron) value;
}

//
// End of primitives for points
//

//
// Begin of advanced primitives
//

// adjustedgering() adjusts the edge version so that it belongs to the
//   indicated edge ring.  The 'direction' only can be 0(CCW) or 1(CW).
//   If the edge is not in the wanted edge ring, reverse it.

inline void tetgenmesh::adjustedgering(triface& t, int direction) {
  if (EdgeRing(t.ver) != direction) {
    esymself(t);
  }
}

inline void tetgenmesh::adjustedgering(face& s, int direction) {
  if (EdgeRing(s.shver) != direction) {
    sesymself(s);
  }
}

// isdead() returns TRUE if the tetrahedron or subface has been dealloced.

inline bool tetgenmesh::isdead(triface* t) {
  if (t->tet == (tetrahedron *) NULL) return true;
  else return t->tet[4] == (tetrahedron) NULL;
}

inline bool tetgenmesh::isdead(face* s) {
  if (s->sh == (shellface *) NULL) return true;
  else return s->sh[3] == (shellface) NULL;
}

// isfacehaspoint() returns TRUE if the 'testpoint' is one of the vertices
//   of the tetface 't' subface 's'.

inline bool tetgenmesh::isfacehaspoint(triface* t, point testpoint) {
  return ((org(*t) == testpoint) || (dest(*t) == testpoint) ||
          (apex(*t) == testpoint));
}

inline bool tetgenmesh::isfacehaspoint(face* s, point testpoint) {
  return (s->sh[3] == (shellface) testpoint) || 
         (s->sh[4] == (shellface) testpoint) ||
         (s->sh[5] == (shellface) testpoint);
}

// isfacehasedge() returns TRUE if the edge (given by its two endpoints) is
//   one of the three edges of the subface 's'.

inline bool tetgenmesh::isfacehasedge(face* s, point tend1, point tend2) {
  return (isfacehaspoint(s, tend1) && isfacehaspoint(s, tend2));
}

// issymexist() returns TRUE if the adjoining tetrahedron is not 'duumytet'.

inline bool tetgenmesh::issymexist(triface* t) {
  tetrahedron *ptr = (tetrahedron *) 
    ((unsigned long)(t->tet[t->loc]) & ~(unsigned long)7l);
  return ptr != dummytet;
}

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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Two inline functions used in read/write VTK files.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

inline void swapBytes(unsigned char* var, int size)
{
  int i = 0;
  int j = size - 1;
  char c;

  while (i < j) {
    c = var[i]; var[i] = var[j]; var[j] = c;
    i++, j--;
  }
}

inline bool testIsBigEndian()
{
  short word = 0x4321;
  if((*(char *)& word) != 0x21)
    return true;
  else 
    return false;
}

#endif // #ifndef tetgenH
