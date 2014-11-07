/*****************************************************************************/
/*                                                                           */
/*      888888888        ,o,                          / 888                  */
/*         888    88o88o  "    o8888o  88o8888o o88888o 888  o88888o         */
/*         888    888    888       88b 888  888 888 888 888 d888  88b        */
/*         888    888    888  o88^o888 888  888 "88888" 888 8888oo888        */
/*         888    888    888 C888  888 888  888  /      888 q888             */
/*         888    888    888  "88o^888 888  888 Cb      888  "88oooo"        */
/*                                              "8oo8D                       */
/*                                                                           */
/*  A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.      */
/*  (triangle.c)                                                             */
/*                                                                           */
/*  Version 1.6                                                              */
/*  July 28, 2005                                                            */
/*                                                                           */
/*  Copyright 1993, 1995, 1997, 1998, 2002, 2005                             */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*  This program may be freely redistributed under the condition that the    */
/*    copyright notices (including this entire header and the copyright      */
/*    notice printed when the `-h' switch is selected) are not removed, and  */
/*    no compensation is received.  Private, research, and institutional     */
/*    use is free.  You may distribute modified versions of this code UNDER  */
/*    THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE   */
/*    SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE   */
/*    AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR    */
/*    NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution of this code as    */
/*    part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT  */
/*    WITH THE AUTHOR.  (If you are not directly supplying this code to a    */
/*    customer, and you are instead telling them how they can obtain it for  */
/*    free, then you are not required to make any arrangement with me.)      */
/*                                                                           */
/*  Hypertext instructions for Triangle are available on the Web at          */
/*                                                                           */
/*      http://www.cs.cmu.edu/~quake/triangle.html                           */
/*                                                                           */
/*  Disclaimer:  Neither I nor Carnegie Mellon warrant this code in any way  */
/*    whatsoever.  This code is provided "as-is".  Use at your own risk.     */
/*                                                                           */
/*  Some of the references listed below are marked with an asterisk.  [*]    */
/*    These references are available for downloading from the Web page       */
/*                                                                           */
/*      http://www.cs.cmu.edu/~quake/triangle.research.html                  */
/*                                                                           */
/*  Three papers discussing aspects of Triangle are available.  A short      */
/*    overview appears in "Triangle:  Engineering a 2D Quality Mesh          */
/*    Generator and Delaunay Triangulator," in Applied Computational         */
/*    Geometry:  Towards Geometric Engineering, Ming C. Lin and Dinesh       */
/*    Manocha, editors, Lecture Notes in Computer Science volume 1148,       */
/*    pages 203-222, Springer-Verlag, Berlin, May 1996 (from the First ACM   */
/*    Workshop on Applied Computational Geometry).  [*]                      */
/*                                                                           */
/*    The algorithms are discussed in the greatest detail in "Delaunay       */
/*    Refinement Algorithms for Triangular Mesh Generation," Computational   */
/*    Geometry:  Theory and Applications 22(1-3):21-74, May 2002.  [*]       */
/*                                                                           */
/*    More detail about the data structures may be found in my dissertation: */
/*    "Delaunay Refinement Mesh Generation," Ph.D. thesis, Technical Report  */
/*    CMU-CS-97-137, School of Computer Science, Carnegie Mellon University, */
/*    Pittsburgh, Pennsylvania, 18 May 1997.  [*]                            */
/*                                                                           */
/*  Triangle was created as part of the Quake Project in the School of       */
/*    Computer Science at Carnegie Mellon University.  For further           */
/*    information, see Hesheng Bao, Jacobo Bielak, Omar Ghattas, Loukas F.   */
/*    Kallivokas, David R. O'Hallaron, Jonathan R. Shewchuk, and Jifeng Xu,  */
/*    "Large-scale Simulation of Elastic Wave Propagation in Heterogeneous   */
/*    Media on Parallel Computers," Computer Methods in Applied Mechanics    */
/*    and Engineering 152(1-2):85-102, 22 January 1998.                      */
/*                                                                           */
/*  Triangle's Delaunay refinement algorithm for quality mesh generation is  */
/*    a hybrid of one due to Jim Ruppert, "A Delaunay Refinement Algorithm   */
/*    for Quality 2-Dimensional Mesh Generation," Journal of Algorithms      */
/*    18(3):548-585, May 1995 [*], and one due to L. Paul Chew, "Guaranteed- */
/*    Quality Mesh Generation for Curved Surfaces," Proceedings of the Ninth */
/*    Annual Symposium on Computational Geometry (San Diego, California),    */
/*    pages 274-280, Association for Computing Machinery, May 1993,          */
/*    http://portal.acm.org/citation.cfm?id=161150 .                         */
/*                                                                           */
/*  The Delaunay refinement algorithm has been modified so that it meshes    */
/*    domains with small input angles well, as described in Gary L. Miller,  */
/*    Steven E. Pav, and Noel J. Walkington, "When and Why Ruppert's         */
/*    Algorithm Works," Twelfth International Meshing Roundtable, pages      */
/*    91-102, Sandia National Laboratories, September 2003.  [*]             */
/*                                                                           */
/*  My implementation of the divide-and-conquer and incremental Delaunay     */
/*    triangulation algorithms follows closely the presentation of Guibas    */
/*    and Stolfi, even though I use a triangle-based data structure instead  */
/*    of their quad-edge data structure.  (In fact, I originally implemented */
/*    Triangle using the quad-edge data structure, but the switch to a       */
/*    triangle-based data structure sped Triangle by a factor of two.)  The  */
/*    mesh manipulation primitives and the two aforementioned Delaunay       */
/*    triangulation algorithms are described by Leonidas J. Guibas and Jorge */
/*    Stolfi, "Primitives for the Manipulation of General Subdivisions and   */
/*    the Computation of Voronoi Diagrams," ACM Transactions on Graphics     */
/*    4(2):74-123, April 1985, http://portal.acm.org/citation.cfm?id=282923 .*/
/*                                                                           */
/*  Their O(n log n) divide-and-conquer algorithm is adapted from Der-Tsai   */
/*    Lee and Bruce J. Schachter, "Two Algorithms for Constructing the       */
/*    Delaunay Triangulation," International Journal of Computer and         */
/*    Information Science 9(3):219-242, 1980.  Triangle's improvement of the */
/*    divide-and-conquer algorithm by alternating between vertical and       */
/*    horizontal cuts was introduced by Rex A. Dwyer, "A Faster Divide-and-  */
/*    Conquer Algorithm for Constructing Delaunay Triangulations,"           */
/*    Algorithmica 2(2):137-151, 1987.                                       */
/*                                                                           */
/*  The incremental insertion algorithm was first proposed by C. L. Lawson,  */
/*    "Software for C1 Surface Interpolation," in Mathematical Software III, */
/*    John R. Rice, editor, Academic Press, New York, pp. 161-194, 1977.     */
/*    For point location, I use the algorithm of Ernst P. Mucke, Isaac       */
/*    Saias, and Binhai Zhu, "Fast Randomized Point Location Without         */
/*    Preprocessing in Two- and Three-Dimensional Delaunay Triangulations,"  */
/*    Proceedings of the Twelfth Annual Symposium on Computational Geometry, */
/*    ACM, May 1996.  [*]  If I were to randomize the order of vertex        */
/*    insertion (I currently don't bother), their result combined with the   */
/*    result of Kenneth L. Clarkson and Peter W. Shor, "Applications of      */
/*    Random Sampling in Computational Geometry II," Discrete &              */
/*    Computational Geometry 4(1):387-421, 1989, would yield an expected     */
/*    O(n^{4/3}) bound on running time.                                      */
/*                                                                           */
/*  The O(n log n) sweepline Delaunay triangulation algorithm is taken from  */
/*    Steven Fortune, "A Sweepline Algorithm for Voronoi Diagrams",          */
/*    Algorithmica 2(2):153-174, 1987.  A random sample of edges on the      */
/*    boundary of the triangulation are maintained in a splay tree for the   */
/*    purpose of point location.  Splay trees are described by Daniel        */
/*    Dominic Sleator and Robert Endre Tarjan, "Self-Adjusting Binary Search */
/*    Trees," Journal of the ACM 32(3):652-686, July 1985,                   */
/*    http://portal.acm.org/citation.cfm?id=3835 .                           */
/*                                                                           */
/*  The algorithms for exact computation of the signs of determinants are    */
/*    described in Jonathan Richard Shewchuk, "Adaptive Precision Floating-  */
/*    Point Arithmetic and Fast Robust Geometric Predicates," Discrete &     */
/*    Computational Geometry 18(3):305-363, October 1997.  (Also available   */
/*    as Technical Report CMU-CS-96-140, School of Computer Science,         */
/*    Carnegie Mellon University, Pittsburgh, Pennsylvania, May 1996.)  [*]  */
/*    An abbreviated version appears as Jonathan Richard Shewchuk, "Robust   */
/*    Adaptive Floating-Point Geometric Predicates," Proceedings of the      */
/*    Twelfth Annual Symposium on Computational Geometry, ACM, May 1996. [*] */
/*    Many of the ideas for my exact arithmetic routines originate with      */
/*    Douglas M. Priest, "Algorithms for Arbitrary Precision Floating Point  */
/*    Arithmetic," Tenth Symposium on Computer Arithmetic, pp. 132-143, IEEE */
/*    Computer Society Press, 1991.  [*]  Many of the ideas for the correct  */
/*    evaluation of the signs of determinants are taken from Steven Fortune  */
/*    and Christopher J. Van Wyk, "Efficient Exact Arithmetic for Computa-   */
/*    tional Geometry," Proceedings of the Ninth Annual Symposium on         */
/*    Computational Geometry, ACM, pp. 163-172, May 1993, and from Steven    */
/*    Fortune, "Numerical Stability of Algorithms for 2D Delaunay Triangu-   */
/*    lations," International Journal of Computational Geometry & Applica-   */
/*    tions 5(1-2):193-213, March-June 1995.                                 */
/*                                                                           */
/*  The method of inserting new vertices off-center (not precisely at the    */
/*    circumcenter of every poor-quality triangle) is from Alper Ungor,      */
/*    "Off-centers:  A New Type of Steiner Points for Computing Size-Optimal */
/*    Quality-Guaranteed Delaunay Triangulations," Proceedings of LATIN      */
/*    2004 (Buenos Aires, Argentina), April 2004.                            */
/*                                                                           */
/*  For definitions of and results involving Delaunay triangulations,        */
/*    constrained and conforming versions thereof, and other aspects of      */
/*    triangular mesh generation, see the excellent survey by Marshall Bern  */
/*    and David Eppstein, "Mesh Generation and Optimal Triangulation," in    */
/*    Computing and Euclidean Geometry, Ding-Zhu Du and Frank Hwang,         */
/*    editors, World Scientific, Singapore, pp. 23-90, 1992.  [*]            */
/*                                                                           */
/*  The time for incrementally adding PSLG (planar straight line graph)      */
/*    segments to create a constrained Delaunay triangulation is probably    */
/*    O(t^2) per segment in the worst case and O(t) per segment in the       */
/*    common case, where t is the number of triangles that intersect the     */
/*    segment before it is inserted.  This doesn't count point location,     */
/*    which can be much more expensive.  I could improve this to O(d log d)  */
/*    time, but d is usually quite small, so it's not worth the bother.      */
/*    (This note does not apply when the -s switch is used, invoking a       */
/*    different method is used to insert segments.)                          */
/*                                                                           */
/*  The time for deleting a vertex from a Delaunay triangulation is O(d^2)   */
/*    in the worst case and O(d) in the common case, where d is the degree   */
/*    of the vertex being deleted.  I could improve this to O(d log d) time, */
/*    but d is usually quite small, so it's not worth the bother.            */
/*                                                                           */
/*  Ruppert's Delaunay refinement algorithm typically generates triangles    */
/*    at a linear rate (constant time per triangle) after the initial        */
/*    triangulation is formed.  There may be pathological cases where        */
/*    quadratic time is required, but these never arise in practice.         */
/*                                                                           */
/*  The geometric predicates (circumcenter calculations, segment             */
/*    intersection formulae, etc.) appear in my "Lecture Notes on Geometric  */
/*    Robustness" at http://www.cs.berkeley.edu/~jrs/mesh .                  */
/*                                                                           */
/*  If you make any improvements to this code, please please please let me   */
/*    know, so that I may obtain the improvements.  Even if you don't change */
/*    the code, I'd still love to hear what it's being used for.             */
/*                                                                           */
/*****************************************************************************/

/* For single precision (which will save some memory and reduce paging),     */
/*   define the symbol SINGLE by using the -DSINGLE compiler switch or by    */
/*   writing "#define SINGLE" below.                                         */
/*                                                                           */
/* For double precision (which will allow you to refine meshes to a smaller  */
/*   edge length), leave SINGLE undefined.                                   */
/*                                                                           */
/* Double precision uses more memory, but improves the resolution of the     */
/*   meshes you can generate with Triangle.  It also reduces the likelihood  */
/*   of a floating exception due to overflow.  Finally, it is much faster    */
/*   than single precision on 64-bit architectures like the DEC Alpha.  I    */
/*   recommend double precision unless you want to generate a mesh for which */
/*   you do not have enough memory.                                          */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

/* If yours is not a Unix system, define the NO_TIMER compiler switch to     */
/*   remove the Unix-specific timing code.                                   */

/* #define NO_TIMER */

/* To insert lots of self-checks for internal errors, define the SELF_CHECK  */
/*   symbol.  This will slow down the program significantly.  It is best to  */
/*   define the symbol using the -DSELF_CHECK compiler switch, but you could */
/*   write "#define SELF_CHECK" below.  If you are modifying this code, I    */
/*   recommend you turn self-checks on until your work is debugged.          */

/* #define SELF_CHECK */

/* To compile Triangle as a callable object library (triangle.o), define the */
/*   TRILIBRARY symbol.  Read the file triangle.h for details on how to call */
/*   the procedure triangulate() that results.                               */

/* #define TRILIBRARY */

/* It is possible to generate a smaller version of Triangle using one or     */
/*   both of the following symbols.  Define the REDUCED symbol to eliminate  */
/*   all features that are primarily of research interest; specifically, the */
/*   -i, -F, -s, and -C switches.  Define the CDT_ONLY symbol to eliminate   */
/*   all meshing algorithms above and beyond constrained Delaunay            */
/*   triangulation; specifically, the -r, -q, -a, -u, -D, -S, and -s         */
/*   switches.  These reductions are most likely to be useful when           */
/*   generating an object library (triangle.o) by defining the TRILIBRARY    */
/*   symbol.                                                                 */

/* #define REDUCED */
/* #define CDT_ONLY */

/* On some machines, my exact arithmetic routines might be defeated by the   */
/*   use of internal extended precision floating-point registers.  The best  */
/*   way to solve this problem is to set the floating-point registers to use */
/*   single or double precision internally.  On 80x86 processors, this may   */
/*   be accomplished by setting the CPU86 symbol for the Microsoft C         */
/*   compiler, or the LINUX symbol for the gcc compiler running on Linux.    */
/*                                                                           */
/* An inferior solution is to declare certain values as `volatile', thus     */
/*   forcing them to be stored to memory and rounded off.  Unfortunately,    */
/*   this solution might slow Triangle down quite a bit.  To use volatile    */
/*   values, write "#define INEXACT volatile" below.  Normally, however,     */
/*   INEXACT should be defined to be nothing.  ("#define INEXACT".)          */
/*                                                                           */
/* For more discussion, see http://www.cs.cmu.edu/~quake/robust.pc.html .    */
/*   For yet more discussion, see Section 5 of my paper, "Adaptive Precision */
/*   Floating-Point Arithmetic and Fast Robust Geometric Predicates" (also   */
/*   available as Section 6.6 of my dissertation).                           */

/* #define CPU86 */
/* #define LINUX */

#define INEXACT /* Nothing */
/* #define INEXACT volatile */

/* Maximum number of characters in a file name (including the null).         */

#define FILENAMESIZE 2048

/* Maximum number of characters in a line read from a file (including the    */
/*   null).                                                                  */

#define INPUTLINESIZE 1024

/* For efficiency, a variety of data structures are allocated in bulk.  The  */
/*   following constants determine how many of each structure is allocated   */
/*   at once.                                                                */

#define TRIPERBLOCK 4092           /* Number of triangles allocated at once. */
#define SUBSEGPERBLOCK 508       /* Number of subsegments allocated at once. */
#define VERTEXPERBLOCK 4092         /* Number of vertices allocated at once. */
#define VIRUSPERBLOCK 1020   /* Number of virus triangles allocated at once. */
/* Number of encroached subsegments allocated at once. */
#define BADSUBSEGPERBLOCK 252
/* Number of skinny triangles allocated at once. */
#define BADTRIPERBLOCK 4092
/* Number of flipped triangles allocated at once. */
#define FLIPSTACKERPERBLOCK 252
/* Number of splay tree nodes allocated at once. */
#define SPLAYNODEPERBLOCK 508

/* The vertex types.   A DEADVERTEX has been deleted entirely.  An           */
/*   UNDEADVERTEX is not part of the mesh, but is written to the output      */
/*   .node file and affects the node indexing in the other output files.     */

#define INPUTVERTEX 0
#define SEGMENTVERTEX 1
#define FREEVERTEX 2
#define DEADVERTEX -32768
#define UNDEADVERTEX -32767

/* The next line is used to outsmart some very stupid compilers.  If your    */
/*   compiler is smarter, feel free to replace the "int" with "void".        */
/*   Not that it matters.                                                    */

#define VOID int

/* Two constants for algorithms based on random sampling.  Both constants    */
/*   have been chosen empirically to optimize their respective algorithms.   */

/* Used for the point location scheme of Mucke, Saias, and Zhu, to decide    */
/*   how large a random sample of triangles to inspect.                      */

#define SAMPLEFACTOR 11

/* Used in Fortune's sweepline Delaunay algorithm to determine what fraction */
/*   of boundary edges should be maintained in the splay tree for point      */
/*   location on the front.                                                  */

#define SAMPLERATE 10

/* A number that speaks for itself, every kissable digit.                    */

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

/* Another fave.                                                             */

#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732

/* And here's one for those of you who are intimidated by math.              */

#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NO_TIMER
#include <sys/time.h>
#endif /* not NO_TIMER */
#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX */
#ifdef TRILIBRARY
#include "triangle.h"
#endif /* TRILIBRARY */

/* A few forward declarations.                                               */

#ifndef TRILIBRARY
char *readline();
char *findfield();
#endif /* not TRILIBRARY */

/* Labels that signify the result of point location.  The result of a        */
/*   search indicates that the point falls in the interior of a triangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */

enum locateresult {INTRIANGLE, ONEDGE, ONVERTEX, OUTSIDE};

/* Labels that signify the result of vertex insertion.  The result indicates */
/*   that the vertex was inserted with complete success, was inserted but    */
/*   encroaches upon a subsegment, was not inserted because it lies on a     */
/*   segment, or was not inserted because another vertex occupies the same   */
/*   location.                                                               */

enum insertvertexresult {SUCCESSFULVERTEX, ENCROACHINGVERTEX, VIOLATINGVERTEX,
                         DUPLICATEVERTEX};

/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction triangle, along the left edge of the direction triangle,  */
/*   or along the right edge of the direction triangle.                      */

enum finddirectionresult {WITHIN, LEFTCOLLINEAR, RIGHTCOLLINEAR};

/*****************************************************************************/
/*                                                                           */
/*  The basic mesh data structures                                           */
/*                                                                           */
/*  There are three:  vertices, triangles, and subsegments (abbreviated      */
/*  `subseg').  These three data structures, linked by pointers, comprise    */
/*  the mesh.  A vertex simply represents a mesh vertex and its properties.  */
/*  A triangle is a triangle.  A subsegment is a special data structure used */
/*  to represent an impenetrable edge of the mesh (perhaps on the outer      */
/*  boundary, on the boundary of a hole, or part of an internal boundary     */
/*  separating two triangulated regions).  Subsegments represent boundaries, */
/*  defined by the user, that triangles may not lie across.                  */
/*                                                                           */
/*  A triangle consists of a list of three vertices, a list of three         */
/*  adjoining triangles, a list of three adjoining subsegments (when         */
/*  segments exist), an arbitrary number of optional user-defined            */
/*  floating-point attributes, and an optional area constraint.  The latter  */
/*  is an upper bound on the permissible area of each triangle in a region,  */
/*  used for mesh refinement.                                                */
/*                                                                           */
/*  For a triangle on a boundary of the mesh, some or all of the neighboring */
/*  triangles may not be present.  For a triangle in the interior of the     */
/*  mesh, often no neighboring subsegments are present.  Such absent         */
/*  triangles and subsegments are never represented by NULL pointers; they   */
/*  are represented by two special records:  `dummytri', the triangle that   */
/*  fills "outer space", and `dummysub', the omnipresent subsegment.         */
/*  `dummytri' and `dummysub' are used for several reasons; for instance,    */
/*  they can be dereferenced and their contents examined without violating   */
/*  protected memory.                                                        */
/*                                                                           */
/*  However, it is important to understand that a triangle includes other    */
/*  information as well.  The pointers to adjoining vertices, triangles, and */
/*  subsegments are ordered in a way that indicates their geometric relation */
/*  to each other.  Furthermore, each of these pointers contains orientation */
/*  information.  Each pointer to an adjoining triangle indicates which face */
/*  of that triangle is contacted.  Similarly, each pointer to an adjoining  */
/*  subsegment indicates which side of that subsegment is contacted, and how */
/*  the subsegment is oriented relative to the triangle.                     */
/*                                                                           */
/*  The data structure representing a subsegment may be thought to be        */
/*  abutting the edge of one or two triangle data structures:  either        */
/*  sandwiched between two triangles, or resting against one triangle on an  */
/*  exterior boundary or hole boundary.                                      */
/*                                                                           */
/*  A subsegment consists of a list of four vertices--the vertices of the    */
/*  subsegment, and the vertices of the segment it is a part of--a list of   */
/*  two adjoining subsegments, and a list of two adjoining triangles.  One   */
/*  of the two adjoining triangles may not be present (though there should   */
/*  always be one), and neighboring subsegments might not be present.        */
/*  Subsegments also store a user-defined integer "boundary marker".         */
/*  Typically, this integer is used to indicate what boundary conditions are */
/*  to be applied at that location in a finite element simulation.           */
/*                                                                           */
/*  Like triangles, subsegments maintain information about the relative      */
/*  orientation of neighboring objects.                                      */
/*                                                                           */
/*  Vertices are relatively simple.  A vertex is a list of floating-point    */
/*  numbers, starting with the x, and y coordinates, followed by an          */
/*  arbitrary number of optional user-defined floating-point attributes,     */
/*  followed by an integer boundary marker.  During the segment insertion    */
/*  phase, there is also a pointer from each vertex to a triangle that may   */
/*  contain it.  Each pointer is not always correct, but when one is, it     */
/*  speeds up segment insertion.  These pointers are assigned values once    */
/*  at the beginning of the segment insertion phase, and are not used or     */
/*  updated except during this phase.  Edge flipping during segment          */
/*  insertion will render some of them incorrect.  Hence, don't rely upon    */
/*  them for anything.                                                       */
/*                                                                           */
/*  Other than the exception mentioned above, vertices have no information   */
/*  about what triangles, subfacets, or subsegments they are linked to.      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Handles                                                                  */
/*                                                                           */
/*  The oriented triangle (`otri') and oriented subsegment (`osub') data     */
/*  structures defined below do not themselves store any part of the mesh.   */
/*  The mesh itself is made of `triangle's, `subseg's, and `vertex's.        */
/*                                                                           */
/*  Oriented triangles and oriented subsegments will usually be referred to  */
/*  as "handles."  A handle is essentially a pointer into the mesh; it       */
/*  allows you to "hold" one particular part of the mesh.  Handles are used  */
/*  to specify the regions in which one is traversing and modifying the mesh.*/
/*  A single `triangle' may be held by many handles, or none at all.  (The   */
/*  latter case is not a memory leak, because the triangle is still          */
/*  connected to other triangles in the mesh.)                               */
/*                                                                           */
/*  An `otri' is a handle that holds a triangle.  It holds a specific edge   */
/*  of the triangle.  An `osub' is a handle that holds a subsegment.  It     */
/*  holds either the left or right side of the subsegment.                   */
/*                                                                           */
/*  Navigation about the mesh is accomplished through a set of mesh          */
/*  manipulation primitives, further below.  Many of these primitives take   */
/*  a handle and produce a new handle that holds the mesh near the first     */
/*  handle.  Other primitives take two handles and glue the corresponding    */
/*  parts of the mesh together.  The orientation of the handles is           */
/*  important.  For instance, when two triangles are glued together by the   */
/*  bond() primitive, they are glued at the edges on which the handles lie.  */
/*                                                                           */
/*  Because vertices have no information about which triangles they are      */
/*  attached to, I commonly represent a vertex by use of a handle whose      */
/*  origin is the vertex.  A single handle can simultaneously represent a    */
/*  triangle, an edge, and a vertex.                                         */
/*                                                                           */
/*****************************************************************************/

/* The triangle data structure.  Each triangle contains three pointers to    */
/*   adjoining triangles, plus three pointers to vertices, plus three        */
/*   pointers to subsegments (declared below; these pointers are usually     */
/*   `dummysub').  It may or may not also contain user-defined attributes    */
/*   and/or a floating-point "area constraint."  It may also contain extra   */
/*   pointers for nodes, when the user asks for high-order elements.         */
/*   Because the size and structure of a `triangle' is not decided until     */
/*   runtime, I haven't simply declared the type `triangle' as a struct.     */

typedef REAL **triangle;            /* Really:  typedef triangle *triangle   */

/* An oriented triangle:  includes a pointer to a triangle and orientation.  */
/*   The orientation denotes an edge of the triangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge always points    */
/*   counterclockwise about the corresponding triangle.                      */

struct otri {
  triangle *tri;
  int orient;                                         /* Ranges from 0 to 2. */
};

/* The subsegment data structure.  Each subsegment contains two pointers to  */
/*   adjoining subsegments, plus four pointers to vertices, plus two         */
/*   pointers to adjoining triangles, plus one boundary marker, plus one     */
/*   segment number.                                                         */

typedef REAL **subseg;                  /* Really:  typedef subseg *subseg   */

/* An oriented subsegment:  includes a pointer to a subsegment and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */

struct osub {
  subseg *ss;
  int ssorient;                                       /* Ranges from 0 to 1. */
};

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */

typedef REAL *vertex;

/* A queue used to store encroached subsegments.  Each subsegment's vertices */
/*   are stored so that we can check whether a subsegment is still the same. */

struct badsubseg {
  subseg encsubseg;                             /* An encroached subsegment. */
  vertex subsegorg, subsegdest;                         /* Its two vertices. */
};

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the triangle.  Each triangle's vertices are    */
/*   stored so that one can check whether a triangle is still the same.      */

struct badtriang {
  triangle poortri;                       /* A skinny or too-large triangle. */
  REAL key;                             /* cos^2 of smallest (apical) angle. */
  vertex triangorg, triangdest, triangapex;           /* Its three vertices. */
  struct badtriang *nexttriang;             /* Pointer to next bad triangle. */
};

/* A stack of triangles flipped during the most recent vertex insertion.     */
/*   The stack is used to undo the vertex insertion if the vertex encroaches */
/*   upon a subsegment.                                                      */

struct flipstacker {
  triangle flippedtri;                       /* A recently flipped triangle. */
  struct flipstacker *prevflip;               /* Previous flip in the stack. */
};

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not point directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `vertex' or to a triangle (in encoded format, so     */
/*   that an orientation is included).  In the latter case, the origin of    */
/*   the oriented triangle is the apex of a "circle event" of the sweepline  */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */

struct event {
  REAL xkey, ykey;                              /* Coordinates of the event. */
  VOID *eventptr;      /* Can be a vertex or the location of a circle event. */
  int heapposition;              /* Marks this event's position in the heap. */
};

/* A node in the splay tree.  Each node holds an oriented ghost triangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a triangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a triangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */

struct splaynode {
  struct otri keyedge;                     /* Lprev of an edge on the front. */
  vertex keydest;           /* Used to verify that splay node is still live. */
  struct splaynode *lchild, *rchild;              /* Children in splay tree. */
};

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* alignbytes determines how new records should be aligned in memory.        */
/*   itembytes is the length of a record in bytes (after rounding up).       */
/*   itemsperblock is the number of items allocated at once in a single      */
/*   block.  itemsfirstblock is the number of items in the first block,      */
/*   which can vary from the others.  items is the number of currently       */
/*   allocated items.  maxitems is the maximum number of items that have     */
/*   been allocated at once; it is the current number of items plus the      */
/*   number of records kept on deaditemstack.                                */

struct memorypool {
  VOID **firstblock, **nowblock;
  VOID *nextitem;
  VOID *deaditemstack;
  VOID **pathblock;
  VOID *pathitem;
  int alignbytes;
  int itembytes;
  int itemsperblock;
  int itemsfirstblock;
  long items, maxitems;
  int unallocateditems;
  int pathitemsleft;
};


/* Global constants.                                                         */

REAL splitter;       /* Used to split REAL factors for exact multiplication. */
REAL epsilon;                             /* Floating-point machine epsilon. */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;

/* Random number seed is not constant, but I've made it global anyway.       */

unsigned long randomseed;                     /* Current random number seed. */


/* Mesh data structure.  Triangle operates on only one mesh, but the mesh    */
/*   structure is used (instead of global variables) to allow reentrancy.    */

struct mesh {

/* Variables used to allocate memory for triangles, subsegments, vertices,   */
/*   viri (triangles being eaten), encroached segments, bad (skinny or too   */
/*   large) triangles, and splay tree nodes.                                 */

  struct memorypool triangles;
  struct memorypool subsegs;
  struct memorypool vertices;
  struct memorypool viri;
  struct memorypool badsubsegs;
  struct memorypool badtriangles;
  struct memorypool flipstackers;
  struct memorypool splaynodes;

/* Variables that maintain the bad triangle queues.  The queues are          */
/*   ordered from 4095 (highest priority) to 0 (lowest priority).            */

  struct badtriang *queuefront[4096];
  struct badtriang *queuetail[4096];
  int nextnonemptyq[4096];
  int firstnonemptyq;

/* Variable that maintains the stack of recently flipped triangles.          */

  struct flipstacker *lastflip;

/* Other variables. */

  REAL xmin, xmax, ymin, ymax;                            /* x and y bounds. */
  REAL xminextreme;      /* Nonexistent x value used as a flag in sweepline. */
  int invertices;                               /* Number of input vertices. */
  int inelements;                              /* Number of input triangles. */
  int insegments;                               /* Number of input segments. */
  int holes;                                       /* Number of input holes. */
  int regions;                                   /* Number of input regions. */
  int undeads;    /* Number of input vertices that don't appear in the mesh. */
  long edges;                                     /* Number of output edges. */
  int mesh_dim;                                /* Dimension (ought to be 2). */
  int nextras;                           /* Number of attributes per vertex. */
  int eextras;                         /* Number of attributes per triangle. */
  long hullsize;                          /* Number of edges in convex hull. */
  int steinerleft;                 /* Number of Steiner points not yet used. */
  int vertexmarkindex;         /* Index to find boundary marker of a vertex. */
  int vertex2triindex;     /* Index to find a triangle adjacent to a vertex. */
  int highorderindex;  /* Index to find extra nodes for high-order elements. */
  int elemattribindex;            /* Index to find attributes of a triangle. */
  int areaboundindex;             /* Index to find area bound of a triangle. */
  int checksegments;         /* Are there segments in the triangulation yet? */
  int checkquality;                  /* Has quality triangulation begun yet? */
  int readnodefile;                           /* Has a .node file been read? */
  long samples;              /* Number of random samples for point location. */

  long incirclecount;                 /* Number of incircle tests performed. */
  long counterclockcount;     /* Number of counterclockwise tests performed. */
  long orient3dcount;           /* Number of 3D orientation tests performed. */
  long hyperbolacount;      /* Number of right-of-hyperbola tests performed. */
  long circumcentercount;  /* Number of circumcenter calculations performed. */
  long circletopcount;       /* Number of circle top calculations performed. */

/* Triangular bounding box vertices.                                         */

  vertex infvertex1, infvertex2, infvertex3;

/* Pointer to the `triangle' that occupies all of "outer space."             */

  triangle *dummytri;
  triangle *dummytribase;    /* Keep base address so we can free() it later. */

/* Pointer to the omnipresent subsegment.  Referenced by any triangle or     */
/*   subsegment that isn't really connected to a subsegment at that          */
/*   location.                                                               */

  subseg *dummysub;
  subseg *dummysubbase;      /* Keep base address so we can free() it later. */

/* Pointer to a recently visited triangle.  Improves point location if       */
/*   proximate vertices are inserted sequentially.                           */

  struct otri recenttri;

};                                                  /* End of `struct mesh'. */


/* Data structure for command line switches and file names.  This structure  */
/*   is used (instead of global variables) to allow reentrancy.              */

struct behavior {

/* Switches for the triangulator.                                            */
/*   poly: -p switch.  refine: -r switch.                                    */
/*   quality: -q switch.                                                     */
/*     minangle: minimum angle bound, specified after -q switch.             */
/*     goodangle: cosine squared of minangle.                                */
/*     offconstant: constant used to place off-center Steiner points.        */
/*   vararea: -a switch without number.                                      */
/*   fixedarea: -a switch with number.                                       */
/*     maxarea: maximum area bound, specified after -a switch.               */
/*   usertest: -u switch.                                                    */
/*   regionattrib: -A switch.  convex: -c switch.                            */
/*   weighted: 1 for -w switch, 2 for -W switch.  jettison: -j switch        */
/*   firstnumber: inverse of -z switch.  All items are numbered starting     */
/*     from `firstnumber'.                                                   */
/*   edgesout: -e switch.  voronoi: -v switch.                               */
/*   neighbors: -n switch.  geomview: -g switch.                             */
/*   nobound: -B switch.  nopolywritten: -P switch.                          */
/*   nonodewritten: -N switch.  noelewritten: -E switch.                     */
/*   noiterationnum: -I switch.  noholes: -O switch.                         */
/*   noexact: -X switch.                                                     */
/*   order: element order, specified after -o switch.                        */
/*   nobisect: count of how often -Y switch is selected.                     */
/*   steiner: maximum number of Steiner points, specified after -S switch.   */
/*   incremental: -i switch.  sweepline: -F switch.                          */
/*   dwyer: inverse of -l switch.                                            */
/*   splitseg: -s switch.                                                    */
/*   conformdel: -D switch.  docheck: -C switch.                             */
/*   quiet: -Q switch.  verbose: count of how often -V switch is selected.   */
/*   usesegments: -p, -r, -q, or -c switch; determines whether segments are  */
/*     used at all.                                                          */
/*                                                                           */
/* Read the instructions to find out the meaning of these switches.          */

  int poly, refine, quality, vararea, fixedarea, usertest;
  int regionattrib, convex, weighted, jettison;
  int firstnumber;
  int edgesout, voronoi, neighbors, geomview;
  int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
  int noholes, noexact, conformdel;
  int incremental, sweepline, dwyer;
  int splitseg;
  int docheck;
  int quiet, verbose;
  int usesegments;
  int order;
  int nobisect;
  int steiner;
  REAL minangle, goodangle, offconstant;
  REAL maxarea;

/* Variables for file names.                                                 */

#ifndef TRILIBRARY
  char innodefilename[FILENAMESIZE];
  char inelefilename[FILENAMESIZE];
  char inpolyfilename[FILENAMESIZE];
  char areafilename[FILENAMESIZE];
  char outnodefilename[FILENAMESIZE];
  char outelefilename[FILENAMESIZE];
  char outpolyfilename[FILENAMESIZE];
  char edgefilename[FILENAMESIZE];
  char vnodefilename[FILENAMESIZE];
  char vedgefilename[FILENAMESIZE];
  char neighborfilename[FILENAMESIZE];
  char offfilename[FILENAMESIZE];
#endif /* not TRILIBRARY */

};                                              /* End of `struct behavior'. */


/*****************************************************************************/
/*                                                                           */
/*  Mesh manipulation primitives.  Each triangle contains three pointers to  */
/*  other triangles, with orientations.  Each pointer points not to the      */
/*  first byte of a triangle, but to one of the first three bytes of a       */
/*  triangle.  It is necessary to extract both the triangle itself and the   */
/*  orientation.  To save memory, I keep both pieces of information in one   */
/*  pointer.  To make this possible, I assume that all triangles are aligned */
/*  to four-byte boundaries.  The decode() routine below decodes a pointer,  */
/*  extracting an orientation (in the range 0 to 2) and a pointer to the     */
/*  beginning of a triangle.  The encode() routine compresses a pointer to a */
/*  triangle and an orientation into a single pointer.  My assumptions that  */
/*  triangles are four-byte-aligned and that the `unsigned long' type is     */
/*  long enough to hold a pointer are two of the few kludges in this program.*/
/*                                                                           */
/*  Subsegments are manipulated similarly.  A pointer to a subsegment        */
/*  carries both an address and an orientation in the range 0 to 1.          */
/*                                                                           */
/*  The other primitives take an oriented triangle or oriented subsegment,   */
/*  and return an oriented triangle or oriented subsegment or vertex; or     */
/*  they change the connections in the data structure.                       */
/*                                                                           */
/*  Below, triangles and subsegments are denoted by their vertices.  The     */
/*  triangle abc has origin (org) a, destination (dest) b, and apex (apex)   */
/*  c.  These vertices occur in counterclockwise order about the triangle.   */
/*  The handle abc may simultaneously denote vertex a, edge ab, and triangle */
/*  abc.                                                                     */
/*                                                                           */
/*  Similarly, the subsegment ab has origin (sorg) a and destination (sdest) */
/*  b.  If ab is thought to be directed upward (with b directly above a),    */
/*  then the handle ab is thought to grasp the right side of ab, and may     */
/*  simultaneously denote vertex a and edge ab.                              */
/*                                                                           */
/*  An asterisk (*) denotes a vertex whose identity is unknown.              */
/*                                                                           */
/*  Given this notation, a partial list of mesh manipulation primitives      */
/*  follows.                                                                 */
/*                                                                           */
/*                                                                           */
/*  For triangles:                                                           */
/*                                                                           */
/*  sym:  Find the abutting triangle; same edge.                             */
/*  sym(abc) -> ba*                                                          */
/*                                                                           */
/*  lnext:  Find the next edge (counterclockwise) of a triangle.             */
/*  lnext(abc) -> bca                                                        */
/*                                                                           */
/*  lprev:  Find the previous edge (clockwise) of a triangle.                */
/*  lprev(abc) -> cab                                                        */
/*                                                                           */
/*  onext:  Find the next edge counterclockwise with the same origin.        */
/*  onext(abc) -> ac*                                                        */
/*                                                                           */
/*  oprev:  Find the next edge clockwise with the same origin.               */
/*  oprev(abc) -> a*b                                                        */
/*                                                                           */
/*  dnext:  Find the next edge counterclockwise with the same destination.   */
/*  dnext(abc) -> *ba                                                        */
/*                                                                           */
/*  dprev:  Find the next edge clockwise with the same destination.          */
/*  dprev(abc) -> cb*                                                        */
/*                                                                           */
/*  rnext:  Find the next edge (counterclockwise) of the adjacent triangle.  */
/*  rnext(abc) -> *a*                                                        */
/*                                                                           */
/*  rprev:  Find the previous edge (clockwise) of the adjacent triangle.     */
/*  rprev(abc) -> b**                                                        */
/*                                                                           */
/*  org:  Origin          dest:  Destination          apex:  Apex            */
/*  org(abc) -> a         dest(abc) -> b              apex(abc) -> c         */
/*                                                                           */
/*  bond:  Bond two triangles together at the resepective handles.           */
/*  bond(abc, bad)                                                           */
/*                                                                           */
/*                                                                           */
/*  For subsegments:                                                         */
/*                                                                           */
/*  ssym:  Reverse the orientation of a subsegment.                          */
/*  ssym(ab) -> ba                                                           */
/*                                                                           */
/*  spivot:  Find adjoining subsegment with the same origin.                 */
/*  spivot(ab) -> a*                                                         */
/*                                                                           */
/*  snext:  Find next subsegment in sequence.                                */
/*  snext(ab) -> b*                                                          */
/*                                                                           */
/*  sorg:  Origin                      sdest:  Destination                   */
/*  sorg(ab) -> a                      sdest(ab) -> b                        */
/*                                                                           */
/*  sbond:  Bond two subsegments together at the respective origins.         */
/*  sbond(ab, ac)                                                            */
/*                                                                           */
/*                                                                           */
/*  For interacting tetrahedra and subfacets:                                */
/*                                                                           */
/*  tspivot:  Find a subsegment abutting a triangle.                         */
/*  tspivot(abc) -> ba                                                       */
/*                                                                           */
/*  stpivot:  Find a triangle abutting a subsegment.                         */
/*  stpivot(ab) -> ba*                                                       */
/*                                                                           */
/*  tsbond:  Bond a triangle to a subsegment.                                */
/*  tsbond(abc, ba)                                                          */
/*                                                                           */
/*****************************************************************************/

/********* Mesh manipulation primitives begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/* Fast lookup arrays to speed some of the mesh manipulation primitives.     */

int plus1mod3[3] = {1, 2, 0};
int minus1mod3[3] = {2, 0, 1};

/********* Primitives for triangles                                  *********/
/*                                                                           */
/*                                                                           */

/* decode() converts a pointer to an oriented triangle.  The orientation is  */
/*   extracted from the two least significant bits of the pointer.           */

#define decode(ptr, otri)                                                     \
  (otri).orient = (int) ((unsigned long) (ptr) & (unsigned long) 3l);         \
  (otri).tri = (triangle *)                                                   \
                  ((unsigned long) (ptr) ^ (unsigned long) (otri).orient)

/* encode() compresses an oriented triangle into a single pointer.  It       */
/*   relies on the assumption that all triangles are aligned to four-byte    */
/*   boundaries, so the two least significant bits of (otri).tri are zero.   */

#define encode(otri)                                                          \
  (triangle) ((unsigned long) (otri).tri | (unsigned long) (otri).orient)

/* The following handle manipulation primitives are all described by Guibas  */
/*   and Stolfi.  However, Guibas and Stolfi use an edge-based data          */
/*   structure, whereas I use a triangle-based data structure.               */

/* sym() finds the abutting triangle, on the same edge.  Note that the edge  */
/*   direction is necessarily reversed, because the handle specified by an   */
/*   oriented triangle is directed counterclockwise around the triangle.     */

#define sym(otri1, otri2)                                                     \
  ptr = (otri1).tri[(otri1).orient];                                          \
  decode(ptr, otri2);

#define symself(otri)                                                         \
  ptr = (otri).tri[(otri).orient];                                            \
  decode(ptr, otri);

/* lnext() finds the next edge (counterclockwise) of a triangle.             */

#define lnext(otri1, otri2)                                                   \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = plus1mod3[(otri1).orient]

#define lnextself(otri)                                                       \
  (otri).orient = plus1mod3[(otri).orient]

/* lprev() finds the previous edge (clockwise) of a triangle.                */

#define lprev(otri1, otri2)                                                   \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = minus1mod3[(otri1).orient]

#define lprevself(otri)                                                       \
  (otri).orient = minus1mod3[(otri).orient]

/* onext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same origin in the counterclockwise direction.  This */
/*   edge is part of a different triangle.                                   */

#define onext(otri1, otri2)                                                   \
  lprev(otri1, otri2);                                                        \
  symself(otri2);

#define onextself(otri)                                                       \
  lprevself(otri);                                                            \
  symself(otri);

/* oprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same origin in the clockwise direction.  This edge is part of  */
/*   a different triangle.                                                   */

#define oprev(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lnextself(otri2);

#define oprevself(otri)                                                       \
  symself(otri);                                                              \
  lnextself(otri);

/* dnext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same destination in the counterclockwise direction.  */
/*   This edge is part of a different triangle.                              */

#define dnext(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lprevself(otri2);

#define dnextself(otri)                                                       \
  symself(otri);                                                              \
  lprevself(otri);

/* dprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same destination in the clockwise direction.  This edge is     */
/*   part of a different triangle.                                           */

#define dprev(otri1, otri2)                                                   \
  lnext(otri1, otri2);                                                        \
  symself(otri2);

#define dprevself(otri)                                                       \
  lnextself(otri);                                                            \
  symself(otri);

/* rnext() moves one edge counterclockwise about the adjacent triangle.      */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rnext(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lnextself(otri2);                                                           \
  symself(otri2);

#define rnextself(otri)                                                       \
  symself(otri);                                                              \
  lnextself(otri);                                                            \
  symself(otri);

/* rprev() moves one edge clockwise about the adjacent triangle.             */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rprev(otri1, otri2)                                                   \
  sym(otri1, otri2);                                                          \
  lprevself(otri2);                                                           \
  symself(otri2);

#define rprevself(otri)                                                       \
  symself(otri);                                                              \
  lprevself(otri);                                                            \
  symself(otri);

/* These primitives determine or set the origin, destination, or apex of a   */
/* triangle.                                                                 */

#define org(otri, vertexptr)                                                  \
  vertexptr = (vertex) (otri).tri[plus1mod3[(otri).orient] + 3]

#define dest(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[minus1mod3[(otri).orient] + 3]

#define apex(otri, vertexptr)                                                 \
  vertexptr = (vertex) (otri).tri[(otri).orient + 3]

#define setorg(otri, vertexptr)                                               \
  (otri).tri[plus1mod3[(otri).orient] + 3] = (triangle) vertexptr

#define setdest(otri, vertexptr)                                              \
  (otri).tri[minus1mod3[(otri).orient] + 3] = (triangle) vertexptr

#define setapex(otri, vertexptr)                                              \
  (otri).tri[(otri).orient + 3] = (triangle) vertexptr

/* Bond two triangles together.                                              */

#define bond(otri1, otri2)                                                    \
  (otri1).tri[(otri1).orient] = encode(otri2);                                \
  (otri2).tri[(otri2).orient] = encode(otri1)

/* Dissolve a bond (from one side).  Note that the other triangle will still */
/*   think it's connected to this triangle.  Usually, however, the other     */
/*   triangle is being deleted entirely, or bonded to another triangle, so   */
/*   it doesn't matter.                                                      */

#define dissolve(otri)                                                        \
  (otri).tri[(otri).orient] = (triangle) m->dummytri

/* Copy an oriented triangle.                                                */

#define otricopy(otri1, otri2)                                                \
  (otri2).tri = (otri1).tri;                                                  \
  (otri2).orient = (otri1).orient

/* Test for equality of oriented triangles.                                  */

#define otriequal(otri1, otri2)                                               \
  (((otri1).tri == (otri2).tri) &&                                            \
   ((otri1).orient == (otri2).orient))

/* Primitives to infect or cure a triangle with the virus.  These rely on    */
/*   the assumption that all subsegments are aligned to four-byte boundaries.*/

#define infect(otri)                                                          \
  (otri).tri[6] = (triangle)                                                  \
                    ((unsigned long) (otri).tri[6] | (unsigned long) 2l)

#define uninfect(otri)                                                        \
  (otri).tri[6] = (triangle)                                                  \
                    ((unsigned long) (otri).tri[6] & ~ (unsigned long) 2l)

/* Test a triangle for viral infection.                                      */

#define infected(otri)                                                        \
  (((unsigned long) (otri).tri[6] & (unsigned long) 2l) != 0l)

/* Check or set a triangle's attributes.                                     */

#define elemattribute(otri, attnum)                                           \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)]

#define setelemattribute(otri, attnum, value)                                 \
  ((REAL *) (otri).tri)[m->elemattribindex + (attnum)] = value

/* Check or set a triangle's maximum area bound.                             */

#define areabound(otri)  ((REAL *) (otri).tri)[m->areaboundindex]

#define setareabound(otri, value)                                             \
  ((REAL *) (otri).tri)[m->areaboundindex] = value

/* Check or set a triangle's deallocation.  Its second pointer is set to     */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its fourth pointer (its first vertex)    */
/*   is set to NULL in case a `badtriang' structure points to it.            */

#define deadtri(tria)  ((tria)[1] == (triangle) NULL)

#define killtri(tria)                                                         \
  (tria)[1] = (triangle) NULL;                                                \
  (tria)[3] = (triangle) NULL

/********* Primitives for subsegments                                *********/
/*                                                                           */
/*                                                                           */

/* sdecode() converts a pointer to an oriented subsegment.  The orientation  */
/*   is extracted from the least significant bit of the pointer.  The two    */
/*   least significant bits (one for orientation, one for viral infection)   */
/*   are masked out to produce the real pointer.                             */

#define sdecode(sptr, osub)                                                   \
  (osub).ssorient = (int) ((unsigned long) (sptr) & (unsigned long) 1l);      \
  (osub).ss = (subseg *)                                                      \
              ((unsigned long) (sptr) & ~ (unsigned long) 3l)

/* sencode() compresses an oriented subsegment into a single pointer.  It    */
/*   relies on the assumption that all subsegments are aligned to two-byte   */
/*   boundaries, so the least significant bit of (osub).ss is zero.          */

#define sencode(osub)                                                         \
  (subseg) ((unsigned long) (osub).ss | (unsigned long) (osub).ssorient)

/* ssym() toggles the orientation of a subsegment.                           */

#define ssym(osub1, osub2)                                                    \
  (osub2).ss = (osub1).ss;                                                    \
  (osub2).ssorient = 1 - (osub1).ssorient

#define ssymself(osub)                                                        \
  (osub).ssorient = 1 - (osub).ssorient

/* spivot() finds the other subsegment (from the same segment) that shares   */
/*   the same origin.                                                        */

#define spivot(osub1, osub2)                                                  \
  sptr = (osub1).ss[(osub1).ssorient];                                        \
  sdecode(sptr, osub2)

#define spivotself(osub)                                                      \
  sptr = (osub).ss[(osub).ssorient];                                          \
  sdecode(sptr, osub)

/* snext() finds the next subsegment (from the same segment) in sequence;    */
/*   one whose origin is the input subsegment's destination.                 */

#define snext(osub1, osub2)                                                   \
  sptr = (osub1).ss[1 - (osub1).ssorient];                                    \
  sdecode(sptr, osub2)

#define snextself(osub)                                                       \
  sptr = (osub).ss[1 - (osub).ssorient];                                      \
  sdecode(sptr, osub)

/* These primitives determine or set the origin or destination of a          */
/*   subsegment or the segment that includes it.                             */

#define sorg(osub, vertexptr)                                                 \
  vertexptr = (vertex) (osub).ss[2 + (osub).ssorient]

#define sdest(osub, vertexptr)                                                \
  vertexptr = (vertex) (osub).ss[3 - (osub).ssorient]

#define setsorg(osub, vertexptr)                                              \
  (osub).ss[2 + (osub).ssorient] = (subseg) vertexptr

#define setsdest(osub, vertexptr)                                             \
  (osub).ss[3 - (osub).ssorient] = (subseg) vertexptr

#define segorg(osub, vertexptr)                                               \
  vertexptr = (vertex) (osub).ss[4 + (osub).ssorient]

#define segdest(osub, vertexptr)                                              \
  vertexptr = (vertex) (osub).ss[5 - (osub).ssorient]

#define setsegorg(osub, vertexptr)                                            \
  (osub).ss[4 + (osub).ssorient] = (subseg) vertexptr

#define setsegdest(osub, vertexptr)                                           \
  (osub).ss[5 - (osub).ssorient] = (subseg) vertexptr

/* These primitives read or set a boundary marker.  Boundary markers are     */
/*   used to hold user-defined tags for setting boundary conditions in       */
/*   finite element solvers.                                                 */

#define mark(osub)  (* (int *) ((osub).ss + 8))

#define setmark(osub, value)                                                  \
  * (int *) ((osub).ss + 8) = value

/* Bond two subsegments together.                                            */

#define sbond(osub1, osub2)                                                   \
  (osub1).ss[(osub1).ssorient] = sencode(osub2);                              \
  (osub2).ss[(osub2).ssorient] = sencode(osub1)

/* Dissolve a subsegment bond (from one side).  Note that the other          */
/*   subsegment will still think it's connected to this subsegment.          */

#define sdissolve(osub)                                                       \
  (osub).ss[(osub).ssorient] = (subseg) m->dummysub

/* Copy a subsegment.                                                        */

#define subsegcopy(osub1, osub2)                                              \
  (osub2).ss = (osub1).ss;                                                    \
  (osub2).ssorient = (osub1).ssorient

/* Test for equality of subsegments.                                         */

#define subsegequal(osub1, osub2)                                             \
  (((osub1).ss == (osub2).ss) &&                                              \
   ((osub1).ssorient == (osub2).ssorient))

/* Check or set a subsegment's deallocation.  Its second pointer is set to   */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its third pointer (its first vertex)     */
/*   is set to NULL in case a `badsubseg' structure points to it.            */

#define deadsubseg(sub)  ((sub)[1] == (subseg) NULL)

#define killsubseg(sub)                                                       \
  (sub)[1] = (subseg) NULL;                                                   \
  (sub)[2] = (subseg) NULL

/********* Primitives for interacting triangles and subsegments      *********/
/*                                                                           */
/*                                                                           */

/* tspivot() finds a subsegment abutting a triangle.                         */

#define tspivot(otri, osub)                                                   \
  sptr = (subseg) (otri).tri[6 + (otri).orient];                              \
  sdecode(sptr, osub)

/* stpivot() finds a triangle abutting a subsegment.  It requires that the   */
/*   variable `ptr' of type `triangle' be defined.                           */

#define stpivot(osub, otri)                                                   \
  ptr = (triangle) (osub).ss[6 + (osub).ssorient];                            \
  decode(ptr, otri)

/* Bond a triangle to a subsegment.                                          */

#define tsbond(otri, osub)                                                    \
  (otri).tri[6 + (otri).orient] = (triangle) sencode(osub);                   \
  (osub).ss[6 + (osub).ssorient] = (subseg) encode(otri)

/* Dissolve a bond (from the triangle side).                                 */

#define tsdissolve(otri)                                                      \
  (otri).tri[6 + (otri).orient] = (triangle) m->dummysub

/* Dissolve a bond (from the subsegment side).                               */

#define stdissolve(osub)                                                      \
  (osub).ss[6 + (osub).ssorient] = (subseg) m->dummytri

/********* Primitives for vertices                                   *********/
/*                                                                           */
/*                                                                           */

#define vertexmark(vx)  ((int *) (vx))[m->vertexmarkindex]

#define setvertexmark(vx, value)                                              \
  ((int *) (vx))[m->vertexmarkindex] = value

#define vertextype(vx)  ((int *) (vx))[m->vertexmarkindex + 1]

#define setvertextype(vx, value)                                              \
  ((int *) (vx))[m->vertexmarkindex + 1] = value

#define vertex2tri(vx)  ((triangle *) (vx))[m->vertex2triindex]

#define setvertex2tri(vx, value)                                              \
  ((triangle *) (vx))[m->vertex2triindex] = value

/**                                                                         **/
/**                                                                         **/
/********* Mesh manipulation primitives end here                     *********/

/********* User-defined triangle evaluation routine begins here      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  triunsuitable()   Determine if a triangle is unsuitable, and thus must   */
/*                    be further refined.                                    */
/*                                                                           */
/*  You may write your own procedure that decides whether or not a selected  */
/*  triangle is too big (and needs to be refined).  There are two ways to do */
/*  this.                                                                    */
/*                                                                           */
/*  (1)  Modify the procedure `triunsuitable' below, then recompile          */
/*  Triangle.                                                                */
/*                                                                           */
/*  (2)  Define the symbol EXTERNAL_TEST (either by adding the definition    */
/*  to this file, or by using the appropriate compiler switch).  This way,   */
/*  you can compile triangle.c separately from your test.  Write your own    */
/*  `triunsuitable' procedure in a separate C file (using the same prototype */
/*  as below).  Compile it and link the object code with triangle.o.         */
/*                                                                           */
/*  This procedure returns 1 if the triangle is too large and should be      */
/*  refined; 0 otherwise.                                                    */
/*                                                                           */
/*****************************************************************************/

#ifdef EXTERNAL_TEST

int triunsuitable();

#else /* not EXTERNAL_TEST */

#ifdef ANSI_DECLARATORS
int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area)
#else /* not ANSI_DECLARATORS */
int triunsuitable(triorg, tridest, triapex, area)
vertex triorg;                              /* The triangle's origin vertex. */
vertex tridest;                        /* The triangle's destination vertex. */
vertex triapex;                               /* The triangle's apex vertex. */
REAL area;                                      /* The area of the triangle. */
#endif /* not ANSI_DECLARATORS */

{
  REAL dxoa, dxda, dxod;
  REAL dyoa, dyda, dyod;
  REAL oalen, dalen, odlen;
  REAL maxlen;

  dxoa = triorg[0] - triapex[0];
  dyoa = triorg[1] - triapex[1];
  dxda = tridest[0] - triapex[0];
  dyda = tridest[1] - triapex[1];
  dxod = triorg[0] - tridest[0];
  dyod = triorg[1] - tridest[1];
  /* Find the squares of the lengths of the triangle's three edges. */
  oalen = dxoa * dxoa + dyoa * dyoa;
  dalen = dxda * dxda + dyda * dyda;
  odlen = dxod * dxod + dyod * dyod;
  /* Find the square of the length of the longest edge. */
  maxlen = (dalen > oalen) ? dalen : oalen;
  maxlen = (odlen > maxlen) ? odlen : maxlen;

  if (maxlen > 0.05 * (triorg[0] * triorg[0] + triorg[1] * triorg[1]) + 0.02) {
    return 1;
  } else {
    return 0;
  }
}

#endif /* not EXTERNAL_TEST */

/**                                                                         **/
/**                                                                         **/
/********* User-defined triangle evaluation routine ends here        *********/

/********* Memory allocation and program exit wrappers begin here    *********/
/**                                                                         **/
/**                                                                         **/

#ifdef ANSI_DECLARATORS
void triexit(int status)
#else /* not ANSI_DECLARATORS */
void triexit(status)
int status;
#endif /* not ANSI_DECLARATORS */

{
  exit(status);
}

#ifdef ANSI_DECLARATORS
VOID *trimalloc(int size)
#else /* not ANSI_DECLARATORS */
VOID *trimalloc(size)
int size;
#endif /* not ANSI_DECLARATORS */

{
  VOID *memptr;

  memptr = (VOID *) malloc((unsigned int) size);
  if (memptr == (VOID *) NULL) {
    printf("Error:  Out of memory.\n");
    triexit(1);
  }
  return(memptr);
}

#ifdef ANSI_DECLARATORS
void trifree(void *memptr)
#else /* not ANSI_DECLARATORS */
void trifree(memptr)
VOID *memptr;
#endif /* not ANSI_DECLARATORS */
{
  free(memptr);
}

/**                                                                         **/
/**                                                                         **/
/********* Memory allocation and program exit wrappers end here      *********/

/********* User interaction routines begin here                      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  syntax()   Print list of command line switches.                          */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

void syntax()
{
#ifdef CDT_ONLY
#ifdef REDUCED
  printf("triangle [-pAcjevngBPNEIOXzo_lQVh] input_file\n");
#else /* not REDUCED */
  printf("triangle [-pAcjevngBPNEIOXzo_iFlCQVh] input_file\n");
#endif /* not REDUCED */
#else /* not CDT_ONLY */
#ifdef REDUCED
  printf("triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__lQVh] input_file\n");
#else /* not REDUCED */
  printf("triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file\n");
#endif /* not REDUCED */
#endif /* not CDT_ONLY */

  printf("    -p  Triangulates a Planar Straight Line Graph (.poly file).\n");
#ifndef CDT_ONLY
  printf("    -r  Refines a previously generated mesh.\n");
  printf(
    "    -q  Quality mesh generation.  A minimum angle may be specified.\n");
  printf("    -a  Applies a maximum triangle area constraint.\n");
  printf("    -u  Applies a user-defined triangle constraint.\n");
#endif /* not CDT_ONLY */
  printf(
    "    -A  Applies attributes to identify triangles in certain regions.\n");
  printf("    -c  Encloses the convex hull with segments.\n");
#ifndef CDT_ONLY
  printf("    -D  Conforming Delaunay:  all triangles are truly Delaunay.\n");
#endif /* not CDT_ONLY */
/*
  printf("    -w  Weighted Delaunay triangulation.\n");
  printf("    -W  Regular triangulation (lower hull of a height field).\n");
*/
  printf("    -j  Jettison unused vertices from output .node file.\n");
  printf("    -e  Generates an edge list.\n");
  printf("    -v  Generates a Voronoi diagram.\n");
  printf("    -n  Generates a list of triangle neighbors.\n");
  printf("    -g  Generates an .off file for Geomview.\n");
  printf("    -B  Suppresses output of boundary information.\n");
  printf("    -P  Suppresses output of .poly file.\n");
  printf("    -N  Suppresses output of .node file.\n");
  printf("    -E  Suppresses output of .ele file.\n");
  printf("    -I  Suppresses mesh iteration numbers.\n");
  printf("    -O  Ignores holes in .poly file.\n");
  printf("    -X  Suppresses use of exact arithmetic.\n");
  printf("    -z  Numbers all items starting from zero (rather than one).\n");
  printf("    -o2 Generates second-order subparametric elements.\n");
#ifndef CDT_ONLY
  printf("    -Y  Suppresses boundary segment splitting.\n");
  printf("    -S  Specifies maximum number of added Steiner points.\n");
#endif /* not CDT_ONLY */
#ifndef REDUCED
  printf("    -i  Uses incremental method, rather than divide-and-conquer.\n");
  printf("    -F  Uses Fortune's sweepline algorithm, rather than d-and-c.\n");
#endif /* not REDUCED */
  printf("    -l  Uses vertical cuts only, rather than alternating cuts.\n");
#ifndef REDUCED
#ifndef CDT_ONLY
  printf(
    "    -s  Force segments into mesh by splitting (instead of using CDT).\n");
#endif /* not CDT_ONLY */
  printf("    -C  Check consistency of final mesh.\n");
#endif /* not REDUCED */
  printf("    -Q  Quiet:  No terminal output except errors.\n");
  printf("    -V  Verbose:  Detailed information on what I'm doing.\n");
  printf("    -h  Help:  Detailed instructions for Triangle.\n");
  triexit(0);
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  info()   Print out complete instructions.                                */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

void info()
{
  printf("Triangle\n");
  printf(
"A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.\n");
  printf("Version 1.6\n\n");
  printf(
"Copyright 1993, 1995, 1997, 1998, 2002, 2005 Jonathan Richard Shewchuk\n");
  printf("2360 Woolsey #H / Berkeley, California 94705-1927\n");
  printf("Bugs/comments to jrs@cs.berkeley.edu\n");
  printf(
"Created as part of the Quake project (tools for earthquake simulation).\n");
  printf(
"Supported in part by NSF Grant CMS-9318163 and an NSERC 1967 Scholarship.\n");
  printf("There is no warranty whatsoever.  Use at your own risk.\n");
#ifdef SINGLE
  printf("This executable is compiled for single precision arithmetic.\n\n\n");
#else /* not SINGLE */
  printf("This executable is compiled for double precision arithmetic.\n\n\n");
#endif /* not SINGLE */
  printf(
"Triangle generates exact Delaunay triangulations, constrained Delaunay\n");
  printf(
"triangulations, conforming Delaunay triangulations, Voronoi diagrams, and\n");
  printf(
"high-quality triangular meshes.  The latter can be generated with no small\n"
);
  printf(
"or large angles, and are thus suitable for finite element analysis.  If no\n"
);
  printf(
"command line switch is specified, your .node input file is read, and the\n");
  printf(
"Delaunay triangulation is returned in .node and .ele output files.  The\n");
  printf("command syntax is:\n\n");
  printf("triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file\n\n");
  printf(
"Underscores indicate that numbers may optionally follow certain switches.\n");
  printf(
"Do not leave any space between a switch and its numeric parameter.\n");
  printf(
"input_file must be a file with extension .node, or extension .poly if the\n");
  printf(
"-p switch is used.  If -r is used, you must supply .node and .ele files,\n");
  printf(
"and possibly a .poly file and an .area file as well.  The formats of these\n"
);
  printf("files are described below.\n\n");
  printf("Command Line Switches:\n\n");
  printf(
"    -p  Reads a Planar Straight Line Graph (.poly file), which can specify\n"
);
  printf(
"        vertices, segments, holes, regional attributes, and regional area\n");
  printf(
"        constraints.  Generates a constrained Delaunay triangulation (CDT)\n"
);
  printf(
"        fitting the input; or, if -s, -q, -a, or -u is used, a conforming\n");
  printf(
"        constrained Delaunay triangulation (CCDT).  If you want a truly\n");
  printf(
"        Delaunay (not just constrained Delaunay) triangulation, use -D as\n");
  printf(
"        well.  When -p is not used, Triangle reads a .node file by default.\n"
);
  printf(
"    -r  Refines a previously generated mesh.  The mesh is read from a .node\n"
);
  printf(
"        file and an .ele file.  If -p is also used, a .poly file is read\n");
  printf(
"        and used to constrain segments in the mesh.  If -a is also used\n");
  printf(
"        (with no number following), an .area file is read and used to\n");
  printf(
"        impose area constraints on the mesh.  Further details on refinement\n"
);
  printf("        appear below.\n");
  printf(
"    -q  Quality mesh generation by Delaunay refinement (a hybrid of Paul\n");
  printf(
"        Chew's and Jim Ruppert's algorithms).  Adds vertices to the mesh to\n"
);
  printf(
"        ensure that all angles are between 20 and 140 degrees.  An\n");
  printf(
"        alternative bound on the minimum angle, replacing 20 degrees, may\n");
  printf(
"        be specified after the `q'.  The specified angle may include a\n");
  printf(
"        decimal point, but not exponential notation.  Note that a bound of\n"
);
  printf(
"        theta degrees on the smallest angle also implies a bound of\n");
  printf(
"        (180 - 2 theta) on the largest angle.  If the minimum angle is 28.6\n"
);
  printf(
"        degrees or smaller, Triangle is mathematically guaranteed to\n");
  printf(
"        terminate (assuming infinite precision arithmetic--Triangle may\n");
  printf(
"        fail to terminate if you run out of precision).  In practice,\n");
  printf(
"        Triangle often succeeds for minimum angles up to 34 degrees.  For\n");
  printf(
"        some meshes, however, you might need to reduce the minimum angle to\n"
);
  printf(
"        avoid problems associated with insufficient floating-point\n");
  printf("        precision.\n");
  printf(
"    -a  Imposes a maximum triangle area.  If a number follows the `a', no\n");
  printf(
"        triangle is generated whose area is larger than that number.  If no\n"
);
  printf(
"        number is specified, an .area file (if -r is used) or .poly file\n");
  printf(
"        (if -r is not used) specifies a set of maximum area constraints.\n");
  printf(
"        An .area file contains a separate area constraint for each\n");
  printf(
"        triangle, and is useful for refining a finite element mesh based on\n"
);
  printf(
"        a posteriori error estimates.  A .poly file can optionally contain\n"
);
  printf(
"        an area constraint for each segment-bounded region, thereby\n");
  printf(
"        controlling triangle densities in a first triangulation of a PSLG.\n"
);
  printf(
"        You can impose both a fixed area constraint and a varying area\n");
  printf(
"        constraint by invoking the -a switch twice, once with and once\n");
  printf(
"        without a number following.  Each area specified may include a\n");
  printf("        decimal point.\n");
  printf(
"    -u  Imposes a user-defined constraint on triangle size.  There are two\n"
);
  printf(
"        ways to use this feature.  One is to edit the triunsuitable()\n");
  printf(
"        procedure in triangle.c to encode any constraint you like, then\n");
  printf(
"        recompile Triangle.  The other is to compile triangle.c with the\n");
  printf(
"        EXTERNAL_TEST symbol set (compiler switch -DEXTERNAL_TEST), then\n");
  printf(
"        link Triangle with a separate object file that implements\n");
  printf(
"        triunsuitable().  In either case, the -u switch causes the user-\n");
  printf("        defined test to be applied to every triangle.\n");
  printf(
"    -A  Assigns an additional floating-point attribute to each triangle\n");
  printf(
"        that identifies what segment-bounded region each triangle belongs\n");
  printf(
"        to.  Attributes are assigned to regions by the .poly file.  If a\n");
  printf(
"        region is not explicitly marked by the .poly file, triangles in\n");
  printf(
"        that region are assigned an attribute of zero.  The -A switch has\n");
  printf(
"        an effect only when the -p switch is used and the -r switch is not.\n"
);
  printf(
"    -c  Creates segments on the convex hull of the triangulation.  If you\n");
  printf(
"        are triangulating a vertex set, this switch causes a .poly file to\n"
);
  printf(
"        be written, containing all edges of the convex hull.  If you are\n");
  printf(
"        triangulating a PSLG, this switch specifies that the whole convex\n");
  printf(
"        hull of the PSLG should be triangulated, regardless of what\n");
  printf(
"        segments the PSLG has.  If you do not use this switch when\n");
  printf(
"        triangulating a PSLG, Triangle assumes that you have identified the\n"
);
  printf(
"        region to be triangulated by surrounding it with segments of the\n");
  printf(
"        input PSLG.  Beware:  if you are not careful, this switch can cause\n"
);
  printf(
"        the introduction of an extremely thin angle between a PSLG segment\n"
);
  printf(
"        and a convex hull segment, which can cause overrefinement (and\n");
  printf(
"        possibly failure if Triangle runs out of precision).  If you are\n");
  printf(
"        refining a mesh, the -c switch works differently:  it causes a\n");
  printf(
"        .poly file to be written containing the boundary edges of the mesh\n"
);
  printf("        (useful if no .poly file was read).\n");
  printf(
"    -D  Conforming Delaunay triangulation:  use this switch if you want to\n"
);
  printf(
"        ensure that all the triangles in the mesh are Delaunay, and not\n");
  printf(
"        merely constrained Delaunay; or if you want to ensure that all the\n"
);
  printf(
"        Voronoi vertices lie within the triangulation.  (Some finite volume\n"
);
  printf(
"        methods have this requirement.)  This switch invokes Ruppert's\n");
  printf(
"        original algorithm, which splits every subsegment whose diametral\n");
  printf(
"        circle is encroached.  It usually increases the number of vertices\n"
);
  printf("        and triangles.\n");
  printf(
"    -j  Jettisons vertices that are not part of the final triangulation\n");
  printf(
"        from the output .node file.  By default, Triangle copies all\n");
  printf(
"        vertices in the input .node file to the output .node file, in the\n");
  printf(
"        same order, so their indices do not change.  The -j switch prevents\n"
);
  printf(
"        duplicated input vertices, or vertices `eaten' by holes, from\n");
  printf(
"        appearing in the output .node file.  Thus, if two input vertices\n");
  printf(
"        have exactly the same coordinates, only the first appears in the\n");
  printf(
"        output.  If any vertices are jettisoned, the vertex numbering in\n");
  printf(
"        the output .node file differs from that of the input .node file.\n");
  printf(
"    -e  Outputs (to an .edge file) a list of edges of the triangulation.\n");
  printf(
"    -v  Outputs the Voronoi diagram associated with the triangulation.\n");
  printf(
"        Does not attempt to detect degeneracies, so some Voronoi vertices\n");
  printf(
"        may be duplicated.  See the discussion of Voronoi diagrams below.\n");
  printf(
"    -n  Outputs (to a .neigh file) a list of triangles neighboring each\n");
  printf("        triangle.\n");
  printf(
"    -g  Outputs the mesh to an Object File Format (.off) file, suitable for\n"
);
  printf("        viewing with the Geometry Center's Geomview package.\n");
  printf(
"    -B  No boundary markers in the output .node, .poly, and .edge output\n");
  printf(
"        files.  See the detailed discussion of boundary markers below.\n");
  printf(
"    -P  No output .poly file.  Saves disk space, but you lose the ability\n");
  printf(
"        to maintain constraining segments on later refinements of the mesh.\n"
);
  printf("    -N  No output .node file.\n");
  printf("    -E  No output .ele file.\n");
  printf(
"    -I  No iteration numbers.  Suppresses the output of .node and .poly\n");
  printf(
"        files, so your input files won't be overwritten.  (If your input is\n"
);
  printf(
"        a .poly file only, a .node file is written.)  Cannot be used with\n");
  printf(
"        the -r switch, because that would overwrite your input .ele file.\n");
  printf(
"        Shouldn't be used with the -q, -a, -u, or -s switch if you are\n");
  printf(
"        using a .node file for input, because no .node file is written, so\n"
);
  printf("        there is no record of any added Steiner points.\n");
  printf("    -O  No holes.  Ignores the holes in the .poly file.\n");
  printf(
"    -X  No exact arithmetic.  Normally, Triangle uses exact floating-point\n"
);
  printf(
"        arithmetic for certain tests if it thinks the inexact tests are not\n"
);
  printf(
"        accurate enough.  Exact arithmetic ensures the robustness of the\n");
  printf(
"        triangulation algorithms, despite floating-point roundoff error.\n");
  printf(
"        Disabling exact arithmetic with the -X switch causes a small\n");
  printf(
"        improvement in speed and creates the possibility that Triangle will\n"
);
  printf("        fail to produce a valid mesh.  Not recommended.\n");
  printf(
"    -z  Numbers all items starting from zero (rather than one).  Note that\n"
);
  printf(
"        this switch is normally overridden by the value used to number the\n"
);
  printf(
"        first vertex of the input .node or .poly file.  However, this\n");
  printf(
"        switch is useful when calling Triangle from another program.\n");
  printf(
"    -o2 Generates second-order subparametric elements with six nodes each.\n"
);
  printf(
"    -Y  No new vertices on the boundary.  This switch is useful when the\n");
  printf(
"        mesh boundary must be preserved so that it conforms to some\n");
  printf(
"        adjacent mesh.  Be forewarned that you will probably sacrifice much\n"
);
  printf(
"        of the quality of the mesh; Triangle will try, but the resulting\n");
  printf(
"        mesh may contain poorly shaped triangles.  Works well if all the\n");
  printf(
"        boundary vertices are closely spaced.  Specify this switch twice\n");
  printf(
"        (`-YY') to prevent all segment splitting, including internal\n");
  printf("        boundaries.\n");
  printf(
"    -S  Specifies the maximum number of Steiner points (vertices that are\n");
  printf(
"        not in the input, but are added to meet the constraints on minimum\n"
);
  printf(
"        angle and maximum area).  The default is to allow an unlimited\n");
  printf(
"        number.  If you specify this switch with no number after it,\n");
  printf(
"        the limit is set to zero.  Triangle always adds vertices at segment\n"
);
  printf(
"        intersections, even if it needs to use more vertices than the limit\n"
);
  printf(
"        you set.  When Triangle inserts segments by splitting (-s), it\n");
  printf(
"        always adds enough vertices to ensure that all the segments of the\n"
);
  printf("        PLSG are recovered, ignoring the limit if necessary.\n");
  printf(
"    -i  Uses an incremental rather than a divide-and-conquer algorithm to\n");
  printf(
"        construct a Delaunay triangulation.  Try it if the divide-and-\n");
  printf("        conquer algorithm fails.\n");
  printf(
"    -F  Uses Steven Fortune's sweepline algorithm to construct a Delaunay\n");
  printf(
"        triangulation.  Warning:  does not use exact arithmetic for all\n");
  printf("        calculations.  An exact result is not guaranteed.\n");
  printf(
"    -l  Uses only vertical cuts in the divide-and-conquer algorithm.  By\n");
  printf(
"        default, Triangle alternates between vertical and horizontal cuts,\n"
);
  printf(
"        which usually improve the speed except with vertex sets that are\n");
  printf(
"        small or short and wide.  This switch is primarily of theoretical\n");
  printf("        interest.\n");
  printf(
"    -s  Specifies that segments should be forced into the triangulation by\n"
);
  printf(
"        recursively splitting them at their midpoints, rather than by\n");
  printf(
"        generating a constrained Delaunay triangulation.  Segment splitting\n"
);
  printf(
"        is true to Ruppert's original algorithm, but can create needlessly\n"
);
  printf(
"        small triangles.  This switch is primarily of theoretical interest.\n"
);
  printf(
"    -C  Check the consistency of the final mesh.  Uses exact arithmetic for\n"
);
  printf(
"        checking, even if the -X switch is used.  Useful if you suspect\n");
  printf("        Triangle is buggy.\n");
  printf(
"    -Q  Quiet:  Suppresses all explanation of what Triangle is doing,\n");
  printf("        unless an error occurs.\n");
  printf(
"    -V  Verbose:  Gives detailed information about what Triangle is doing.\n"
);
  printf(
"        Add more `V's for increasing amount of detail.  `-V' is most\n");
  printf(
"        useful; itgives information on algorithmic progress and much more\n");
  printf(
"        detailed statistics.  `-VV' gives vertex-by-vertex details, and\n");
  printf(
"        prints so much that Triangle runs much more slowly.  `-VVVV' gives\n"
);
  printf("        information only a debugger could love.\n");
  printf("    -h  Help:  Displays these instructions.\n");
  printf("\n");
  printf("Definitions:\n");
  printf("\n");
  printf(
"  A Delaunay triangulation of a vertex set is a triangulation whose\n");
  printf(
"  vertices are the vertex set, that covers the convex hull of the vertex\n");
  printf(
"  set.  A Delaunay triangulation has the property that no vertex lies\n");
  printf(
"  inside the circumscribing circle (circle that passes through all three\n");
  printf("  vertices) of any triangle in the triangulation.\n\n");
  printf(
"  A Voronoi diagram of a vertex set is a subdivision of the plane into\n");
  printf(
"  polygonal cells (some of which may be unbounded, meaning infinitely\n");
  printf(
"  large), where each cell is the set of points in the plane that are closer\n"
);
  printf(
"  to some input vertex than to any other input vertex.  The Voronoi diagram\n"
);
  printf("  is a geometric dual of the Delaunay triangulation.\n\n");
  printf(
"  A Planar Straight Line Graph (PSLG) is a set of vertices and segments.\n");
  printf(
"  Segments are simply edges, whose endpoints are all vertices in the PSLG.\n"
);
  printf(
"  Segments may intersect each other only at their endpoints.  The file\n");
  printf("  format for PSLGs (.poly files) is described below.\n\n");
  printf(
"  A constrained Delaunay triangulation (CDT) of a PSLG is similar to a\n");
  printf(
"  Delaunay triangulation, but each PSLG segment is present as a single edge\n"
);
  printf(
"  of the CDT.  (A constrained Delaunay triangulation is not truly a\n");
  printf(
"  Delaunay triangulation, because some of its triangles might not be\n");
  printf(
"  Delaunay.)  By definition, a CDT does not have any vertices other than\n");
  printf(
"  those specified in the input PSLG.  Depending on context, a CDT might\n");
  printf(
"  cover the convex hull of the PSLG, or it might cover only a segment-\n");
  printf("  bounded region (e.g. a polygon).\n\n");
  printf(
"  A conforming Delaunay triangulation of a PSLG is a triangulation in which\n"
);
  printf(
"  each triangle is truly Delaunay, and each PSLG segment is represented by\n"
);
  printf(
"  a linear contiguous sequence of edges of the triangulation.  New vertices\n"
);
  printf(
"  (not part of the PSLG) may appear, and each input segment may have been\n");
  printf(
"  subdivided into shorter edges (subsegments) by these additional vertices.\n"
);
  printf(
"  The new vertices are frequently necessary to maintain the Delaunay\n");
  printf("  property while ensuring that every segment is represented.\n\n");
  printf(
"  A conforming constrained Delaunay triangulation (CCDT) of a PSLG is a\n");
  printf(
"  triangulation of a PSLG whose triangles are constrained Delaunay.  New\n");
  printf("  vertices may appear, and input segments may be subdivided into\n");
  printf(
"  subsegments, but not to guarantee that segments are respected; rather, to\n"
);
  printf(
"  improve the quality of the triangles.  The high-quality meshes produced\n");
  printf(
"  by the -q switch are usually CCDTs, but can be made conforming Delaunay\n");
  printf("  with the -D switch.\n\n");
  printf("File Formats:\n\n");
  printf(
"  All files may contain comments prefixed by the character '#'.  Vertices,\n"
);
  printf(
"  triangles, edges, holes, and maximum area constraints must be numbered\n");
  printf(
"  consecutively, starting from either 1 or 0.  Whichever you choose, all\n");
  printf(
"  input files must be consistent; if the vertices are numbered from 1, so\n");
  printf(
"  must be all other objects.  Triangle automatically detects your choice\n");
  printf(
"  while reading the .node (or .poly) file.  (When calling Triangle from\n");
  printf(
"  another program, use the -z switch if you wish to number objects from\n");
  printf("  zero.)  Examples of these file formats are given below.\n\n");
  printf("  .node files:\n");
  printf(
"    First line:  <# of vertices> <dimension (must be 2)> <# of attributes>\n"
);
  printf(
"                                           <# of boundary markers (0 or 1)>\n"
);
  printf(
"    Remaining lines:  <vertex #> <x> <y> [attributes] [boundary marker]\n");
  printf("\n");
  printf(
"    The attributes, which are typically floating-point values of physical\n");
  printf(
"    quantities (such as mass or conductivity) associated with the nodes of\n"
);
  printf(
"    a finite element mesh, are copied unchanged to the output mesh.  If -q,\n"
);
  printf(
"    -a, -u, -D, or -s is selected, each new Steiner point added to the mesh\n"
);
  printf("    has attributes assigned to it by linear interpolation.\n\n");
  printf(
"    If the fourth entry of the first line is `1', the last column of the\n");
  printf(
"    remainder of the file is assumed to contain boundary markers.  Boundary\n"
);
  printf(
"    markers are used to identify boundary vertices and vertices resting on\n"
);
  printf(
"    PSLG segments; a complete description appears in a section below.  The\n"
);
  printf(
"    .node file produced by Triangle contains boundary markers in the last\n");
  printf("    column unless they are suppressed by the -B switch.\n\n");
  printf("  .ele files:\n");
  printf(
"    First line:  <# of triangles> <nodes per triangle> <# of attributes>\n");
  printf(
"    Remaining lines:  <triangle #> <node> <node> <node> ... [attributes]\n");
  printf("\n");
  printf(
"    Nodes are indices into the corresponding .node file.  The first three\n");
  printf(
"    nodes are the corner vertices, and are listed in counterclockwise order\n"
);
  printf(
"    around each triangle.  (The remaining nodes, if any, depend on the type\n"
);
  printf("    of finite element used.)\n\n");
  printf(
"    The attributes are just like those of .node files.  Because there is no\n"
);
  printf(
"    simple mapping from input to output triangles, Triangle attempts to\n");
  printf(
"    interpolate attributes, and may cause a lot of diffusion of attributes\n"
);
  printf(
"    among nearby triangles as the triangulation is refined.  Attributes do\n"
);
  printf("    not diffuse across segments, so attributes used to identify\n");
  printf("    segment-bounded regions remain intact.\n\n");
  printf(
"    In .ele files produced by Triangle, each triangular element has three\n");
  printf(
"    nodes (vertices) unless the -o2 switch is used, in which case\n");
  printf(
"    subparametric quadratic elements with six nodes each are generated.\n");
  printf(
"    The first three nodes are the corners in counterclockwise order, and\n");
  printf(
"    the fourth, fifth, and sixth nodes lie on the midpoints of the edges\n");
  printf(
"    opposite the first, second, and third vertices, respectively.\n");
  printf("\n");
  printf("  .poly files:\n");
  printf(
"    First line:  <# of vertices> <dimension (must be 2)> <# of attributes>\n"
);
  printf(
"                                           <# of boundary markers (0 or 1)>\n"
);
  printf(
"    Following lines:  <vertex #> <x> <y> [attributes] [boundary marker]\n");
  printf("    One line:  <# of segments> <# of boundary markers (0 or 1)>\n");
  printf(
"    Following lines:  <segment #> <endpoint> <endpoint> [boundary marker]\n");
  printf("    One line:  <# of holes>\n");
  printf("    Following lines:  <hole #> <x> <y>\n");
  printf(
"    Optional line:  <# of regional attributes and/or area constraints>\n");
  printf(
"    Optional following lines:  <region #> <x> <y> <attribute> <max area>\n");
  printf("\n");
  printf(
"    A .poly file represents a PSLG, as well as some additional information.\n"
);
  printf(
"    The first section lists all the vertices, and is identical to the\n");
  printf(
"    format of .node files.  <# of vertices> may be set to zero to indicate\n"
);
  printf(
"    that the vertices are listed in a separate .node file; .poly files\n");
  printf(
"    produced by Triangle always have this format.  A vertex set represented\n"
);
  printf(
"    this way has the advantage that it may easily be triangulated with or\n");
  printf(
"    without segments (depending on whether the -p switch is invoked).\n");
  printf("\n");
  printf(
"    The second section lists the segments.  Segments are edges whose\n");
  printf(
"    presence in the triangulation is enforced.  (Depending on the choice of\n"
);
  printf(
"    switches, segment might be subdivided into smaller edges).  Each\n");
  printf(
"    segment is specified by listing the indices of its two endpoints.  This\n"
);
  printf(
"    means that you must include its endpoints in the vertex list.  Each\n");
  printf("    segment, like each point, may have a boundary marker.\n\n");
  printf(
"    If -q, -a, -u, and -s are not selected, Triangle produces a constrained\n"
);
  printf(
"    Delaunay triangulation (CDT), in which each segment appears as a single\n"
);
  printf(
"    edge in the triangulation.  If -q, -a, -u, or -s is selected, Triangle\n"
);
  printf(
"    produces a conforming constrained Delaunay triangulation (CCDT), in\n");
  printf(
"    which segments may be subdivided into smaller edges.  If -D is\n");
  printf(
"    selected, Triangle produces a conforming Delaunay triangulation, so\n");
  printf(
"    that every triangle is Delaunay, and not just constrained Delaunay.\n");
  printf("\n");
  printf(
"    The third section lists holes (and concavities, if -c is selected) in\n");
  printf(
"    the triangulation.  Holes are specified by identifying a point inside\n");
  printf(
"    each hole.  After the triangulation is formed, Triangle creates holes\n");
  printf(
"    by eating triangles, spreading out from each hole point until its\n");
  printf(
"    progress is blocked by segments in the PSLG.  You must be careful to\n");
  printf(
"    enclose each hole in segments, or your whole triangulation might be\n");
  printf(
"    eaten away.  If the two triangles abutting a segment are eaten, the\n");
  printf(
"    segment itself is also eaten.  Do not place a hole directly on a\n");
  printf("    segment; if you do, Triangle chooses one side of the segment\n");
  printf("    arbitrarily.\n\n");
  printf(
"    The optional fourth section lists regional attributes (to be assigned\n");
  printf(
"    to all triangles in a region) and regional constraints on the maximum\n");
  printf(
"    triangle area.  Triangle reads this section only if the -A switch is\n");
  printf(
"    used or the -a switch is used without a number following it, and the -r\n"
);
  printf(
"    switch is not used.  Regional attributes and area constraints are\n");
  printf(
"    propagated in the same manner as holes:  you specify a point for each\n");
  printf(
"    attribute and/or constraint, and the attribute and/or constraint\n");
  printf(
"    affects the whole region (bounded by segments) containing the point.\n");
  printf(
"    If two values are written on a line after the x and y coordinate, the\n");
  printf(
"    first such value is assumed to be a regional attribute (but is only\n");
  printf(
"    applied if the -A switch is selected), and the second value is assumed\n"
);
  printf(
"    to be a regional area constraint (but is only applied if the -a switch\n"
);
  printf(
"    is selected).  You may specify just one value after the coordinates,\n");
  printf(
"    which can serve as both an attribute and an area constraint, depending\n"
);
  printf(
"    on the choice of switches.  If you are using the -A and -a switches\n");
  printf(
"    simultaneously and wish to assign an attribute to some region without\n");
  printf("    imposing an area constraint, use a negative maximum area.\n\n");
  printf(
"    When a triangulation is created from a .poly file, you must either\n");
  printf(
"    enclose the entire region to be triangulated in PSLG segments, or\n");
  printf(
"    use the -c switch, which automatically creates extra segments that\n");
  printf(
"    enclose the convex hull of the PSLG.  If you do not use the -c switch,\n"
);
  printf(
"    Triangle eats all triangles that are not enclosed by segments; if you\n");
  printf(
"    are not careful, your whole triangulation may be eaten away.  If you do\n"
);
  printf(
"    use the -c switch, you can still produce concavities by the appropriate\n"
);
  printf(
"    placement of holes just inside the boundary of the convex hull.\n");
  printf("\n");
  printf(
"    An ideal PSLG has no intersecting segments, nor any vertices that lie\n");
  printf(
"    upon segments (except, of course, the endpoints of each segment).  You\n"
);
  printf(
"    aren't required to make your .poly files ideal, but you should be aware\n"
);
  printf(
"    of what can go wrong.  Segment intersections are relatively safe--\n");
  printf(
"    Triangle calculates the intersection points for you and adds them to\n");
  printf(
"    the triangulation--as long as your machine's floating-point precision\n");
  printf(
"    doesn't become a problem.  You are tempting the fates if you have three\n"
);
  printf(
"    segments that cross at the same location, and expect Triangle to figure\n"
);
  printf(
"    out where the intersection point is.  Thanks to floating-point roundoff\n"
);
  printf(
"    error, Triangle will probably decide that the three segments intersect\n"
);
  printf(
"    at three different points, and you will find a minuscule triangle in\n");
  printf(
"    your output--unless Triangle tries to refine the tiny triangle, uses\n");
  printf(
"    up the last bit of machine precision, and fails to terminate at all.\n");
  printf(
"    You're better off putting the intersection point in the input files,\n");
  printf(
"    and manually breaking up each segment into two.  Similarly, if you\n");
  printf(
"    place a vertex at the middle of a segment, and hope that Triangle will\n"
);
  printf(
"    break up the segment at that vertex, you might get lucky.  On the other\n"
);
  printf(
"    hand, Triangle might decide that the vertex doesn't lie precisely on\n");
  printf(
"    the segment, and you'll have a needle-sharp triangle in your output--or\n"
);
  printf("    a lot of tiny triangles if you're generating a quality mesh.\n");
  printf("\n");
  printf(
"    When Triangle reads a .poly file, it also writes a .poly file, which\n");
  printf(
"    includes all the subsegments--the edges that are parts of input\n");
  printf(
"    segments.  If the -c switch is used, the output .poly file also\n");
  printf(
"    includes all of the edges on the convex hull.  Hence, the output .poly\n"
);
  printf(
"    file is useful for finding edges associated with input segments and for\n"
);
  printf(
"    setting boundary conditions in finite element simulations.  Moreover,\n");
  printf(
"    you will need the output .poly file if you plan to refine the output\n");
  printf(
"    mesh, and don't want segments to be missing in later triangulations.\n");
  printf("\n");
  printf("  .area files:\n");
  printf("    First line:  <# of triangles>\n");
  printf("    Following lines:  <triangle #> <maximum area>\n");
  printf("\n");
  printf(
"    An .area file associates with each triangle a maximum area that is used\n"
);
  printf(
"    for mesh refinement.  As with other file formats, every triangle must\n");
  printf(
"    be represented, and the triangles must be numbered consecutively.  A\n");
  printf(
"    triangle may be left unconstrained by assigning it a negative maximum\n");
  printf("    area.\n\n");
  printf("  .edge files:\n");
  printf("    First line:  <# of edges> <# of boundary markers (0 or 1)>\n");
  printf(
"    Following lines:  <edge #> <endpoint> <endpoint> [boundary marker]\n");
  printf("\n");
  printf(
"    Endpoints are indices into the corresponding .node file.  Triangle can\n"
);
  printf(
"    produce .edge files (use the -e switch), but cannot read them.  The\n");
  printf(
"    optional column of boundary markers is suppressed by the -B switch.\n");
  printf("\n");
  printf(
"    In Voronoi diagrams, one also finds a special kind of edge that is an\n");
  printf(
"    infinite ray with only one endpoint.  For these edges, a different\n");
  printf("    format is used:\n\n");
  printf("        <edge #> <endpoint> -1 <direction x> <direction y>\n\n");
  printf(
"    The `direction' is a floating-point vector that indicates the direction\n"
);
  printf("    of the infinite ray.\n\n");
  printf("  .neigh files:\n");
  printf(
"    First line:  <# of triangles> <# of neighbors per triangle (always 3)>\n"
);
  printf(
"    Following lines:  <triangle #> <neighbor> <neighbor> <neighbor>\n");
  printf("\n");
  printf(
"    Neighbors are indices into the corresponding .ele file.  An index of -1\n"
);
  printf(
"    indicates no neighbor (because the triangle is on an exterior\n");
  printf(
"    boundary).  The first neighbor of triangle i is opposite the first\n");
  printf("    corner of triangle i, and so on.\n\n");
  printf(
"    Triangle can produce .neigh files (use the -n switch), but cannot read\n"
);
  printf("    them.\n\n");
  printf("Boundary Markers:\n\n");
  printf(
"  Boundary markers are tags used mainly to identify which output vertices\n");
  printf(
"  and edges are associated with which PSLG segment, and to identify which\n");
  printf(
"  vertices and edges occur on a boundary of the triangulation.  A common\n");
  printf(
"  use is to determine where boundary conditions should be applied to a\n");
  printf(
"  finite element mesh.  You can prevent boundary markers from being written\n"
);
  printf("  into files produced by Triangle by using the -B switch.\n\n");
  printf(
"  The boundary marker associated with each segment in an output .poly file\n"
);
  printf("  and each edge in an output .edge file is chosen as follows:\n");
  printf(
"    - If an output edge is part or all of a PSLG segment with a nonzero\n");
  printf(
"      boundary marker, then the edge is assigned the same marker.\n");
  printf(
"    - Otherwise, if the edge lies on a boundary of the triangulation\n");
  printf(
"      (even the boundary of a hole), then the edge is assigned the marker\n");
  printf("      one (1).\n");
  printf("    - Otherwise, the edge is assigned the marker zero (0).\n");
  printf(
"  The boundary marker associated with each vertex in an output .node file\n");
  printf("  is chosen as follows:\n");
  printf(
"    - If a vertex is assigned a nonzero boundary marker in the input file,\n"
);
  printf(
"      then it is assigned the same marker in the output .node file.\n");
  printf(
"    - Otherwise, if the vertex lies on a PSLG segment (even if it is an\n");
  printf(
"      endpoint of the segment) with a nonzero boundary marker, then the\n");
  printf(
"      vertex is assigned the same marker.  If the vertex lies on several\n");
  printf("      such segments, one of the markers is chosen arbitrarily.\n");
  printf(
"    - Otherwise, if the vertex occurs on a boundary of the triangulation,\n");
  printf("      then the vertex is assigned the marker one (1).\n");
  printf("    - Otherwise, the vertex is assigned the marker zero (0).\n");
  printf("\n");
  printf(
"  If you want Triangle to determine for you which vertices and edges are on\n"
);
  printf(
"  the boundary, assign them the boundary marker zero (or use no markers at\n"
);
  printf(
"  all) in your input files.  In the output files, all boundary vertices,\n");
  printf("  edges, and segments will be assigned the value one.\n\n");
  printf("Triangulation Iteration Numbers:\n\n");
  printf(
"  Because Triangle can read and refine its own triangulations, input\n");
  printf(
"  and output files have iteration numbers.  For instance, Triangle might\n");
  printf(
"  read the files mesh.3.node, mesh.3.ele, and mesh.3.poly, refine the\n");
  printf(
"  triangulation, and output the files mesh.4.node, mesh.4.ele, and\n");
  printf("  mesh.4.poly.  Files with no iteration number are treated as if\n");
  printf(
"  their iteration number is zero; hence, Triangle might read the file\n");
  printf(
"  points.node, triangulate it, and produce the files points.1.node and\n");
  printf("  points.1.ele.\n\n");
  printf(
"  Iteration numbers allow you to create a sequence of successively finer\n");
  printf(
"  meshes suitable for multigrid methods.  They also allow you to produce a\n"
);
  printf(
"  sequence of meshes using error estimate-driven mesh refinement.\n");
  printf("\n");
  printf(
"  If you're not using refinement or quality meshing, and you don't like\n");
  printf(
"  iteration numbers, use the -I switch to disable them.  This switch also\n");
  printf(
"  disables output of .node and .poly files to prevent your input files from\n"
);
  printf(
"  being overwritten.  (If the input is a .poly file that contains its own\n");
  printf(
"  points, a .node file is written.  This can be quite convenient for\n");
  printf("  computing CDTs or quality meshes.)\n\n");
  printf("Examples of How to Use Triangle:\n\n");
  printf(
"  `triangle dots' reads vertices from dots.node, and writes their Delaunay\n"
);
  printf(
"  triangulation to dots.1.node and dots.1.ele.  (dots.1.node is identical\n");
  printf(
"  to dots.node.)  `triangle -I dots' writes the triangulation to dots.ele\n");
  printf(
"  instead.  (No additional .node file is needed, so none is written.)\n");
  printf("\n");
  printf(
"  `triangle -pe object.1' reads a PSLG from object.1.poly (and possibly\n");
  printf(
"  object.1.node, if the vertices are omitted from object.1.poly) and writes\n"
);
  printf(
"  its constrained Delaunay triangulation to object.2.node and object.2.ele.\n"
);
  printf(
"  The segments are copied to object.2.poly, and all edges are written to\n");
  printf("  object.2.edge.\n\n");
  printf(
"  `triangle -pq31.5a.1 object' reads a PSLG from object.poly (and possibly\n"
);
  printf(
"  object.node), generates a mesh whose angles are all between 31.5 and 117\n"
);
  printf(
"  degrees and whose triangles all have areas of 0.1 or less, and writes the\n"
);
  printf(
"  mesh to object.1.node and object.1.ele.  Each segment may be broken up\n");
  printf("  into multiple subsegments; these are written to object.1.poly.\n");
  printf("\n");
  printf(
"  Here is a sample file `box.poly' describing a square with a square hole:\n"
);
  printf("\n");
  printf(
"    # A box with eight vertices in 2D, no attributes, one boundary marker.\n"
);
  printf("    8 2 0 1\n");
  printf("     # Outer box has these vertices:\n");
  printf("     1   0 0   0\n");
  printf("     2   0 3   0\n");
  printf("     3   3 0   0\n");
  printf("     4   3 3   33     # A special marker for this vertex.\n");
  printf("     # Inner square has these vertices:\n");
  printf("     5   1 1   0\n");
  printf("     6   1 2   0\n");
  printf("     7   2 1   0\n");
  printf("     8   2 2   0\n");
  printf("    # Five segments with boundary markers.\n");
  printf("    5 1\n");
  printf("     1   1 2   5      # Left side of outer box.\n");
  printf("     # Square hole has these segments:\n");
  printf("     2   5 7   0\n");
  printf("     3   7 8   0\n");
  printf("     4   8 6   10\n");
  printf("     5   6 5   0\n");
  printf("    # One hole in the middle of the inner square.\n");
  printf("    1\n");
  printf("     1   1.5 1.5\n");
  printf("\n");
  printf(
"  Note that some segments are missing from the outer square, so you must\n");
  printf(
"  use the `-c' switch.  After `triangle -pqc box.poly', here is the output\n"
);
  printf(
"  file `box.1.node', with twelve vertices.  The last four vertices were\n");
  printf(
"  added to meet the angle constraint.  Vertices 1, 2, and 9 have markers\n");
  printf(
"  from segment 1.  Vertices 6 and 8 have markers from segment 4.  All the\n");
  printf(
"  other vertices but 4 have been marked to indicate that they lie on a\n");
  printf("  boundary.\n\n");
  printf("    12  2  0  1\n");
  printf("       1    0   0      5\n");
  printf("       2    0   3      5\n");
  printf("       3    3   0      1\n");
  printf("       4    3   3     33\n");
  printf("       5    1   1      1\n");
  printf("       6    1   2     10\n");
  printf("       7    2   1      1\n");
  printf("       8    2   2     10\n");
  printf("       9    0   1.5    5\n");
  printf("      10    1.5   0    1\n");
  printf("      11    3   1.5    1\n");
  printf("      12    1.5   3    1\n");
  printf("    # Generated by triangle -pqc box.poly\n");
  printf("\n");
  printf("  Here is the output file `box.1.ele', with twelve triangles.\n");
  printf("\n");
  printf("    12  3  0\n");
  printf("       1     5   6   9\n");
  printf("       2    10   3   7\n");
  printf("       3     6   8  12\n");
  printf("       4     9   1   5\n");
  printf("       5     6   2   9\n");
  printf("       6     7   3  11\n");
  printf("       7    11   4   8\n");
  printf("       8     7   5  10\n");
  printf("       9    12   2   6\n");
  printf("      10     8   7  11\n");
  printf("      11     5   1  10\n");
  printf("      12     8   4  12\n");
  printf("    # Generated by triangle -pqc box.poly\n\n");
  printf(
"  Here is the output file `box.1.poly'.  Note that segments have been added\n"
);
  printf(
"  to represent the convex hull, and some segments have been subdivided by\n");
  printf(
"  newly added vertices.  Note also that <# of vertices> is set to zero to\n");
  printf("  indicate that the vertices should be read from the .node file.\n");
  printf("\n");
  printf("    0  2  0  1\n");
  printf("    12  1\n");
  printf("       1     1   9     5\n");
  printf("       2     5   7     1\n");
  printf("       3     8   7     1\n");
  printf("       4     6   8    10\n");
  printf("       5     5   6     1\n");
  printf("       6     3  10     1\n");
  printf("       7     4  11     1\n");
  printf("       8     2  12     1\n");
  printf("       9     9   2     5\n");
  printf("      10    10   1     1\n");
  printf("      11    11   3     1\n");
  printf("      12    12   4     1\n");
  printf("    1\n");
  printf("       1   1.5 1.5\n");
  printf("    # Generated by triangle -pqc box.poly\n");
  printf("\n");
  printf("Refinement and Area Constraints:\n");
  printf("\n");
  printf(
"  The -r switch causes a mesh (.node and .ele files) to be read and\n");
  printf(
"  refined.  If the -p switch is also used, a .poly file is read and used to\n"
);
  printf(
"  specify edges that are constrained and cannot be eliminated (although\n");
  printf(
"  they can be subdivided into smaller edges) by the refinement process.\n");
  printf("\n");
  printf(
"  When you refine a mesh, you generally want to impose tighter constraints.\n"
);
  printf(
"  One way to accomplish this is to use -q with a larger angle, or -a\n");
  printf(
"  followed by a smaller area than you used to generate the mesh you are\n");
  printf(
"  refining.  Another way to do this is to create an .area file, which\n");
  printf(
"  specifies a maximum area for each triangle, and use the -a switch\n");
  printf(
"  (without a number following).  Each triangle's area constraint is applied\n"
);
  printf(
"  to that triangle.  Area constraints tend to diffuse as the mesh is\n");
  printf(
"  refined, so if there are large variations in area constraint between\n");
  printf(
"  adjacent triangles, you may not get the results you want.  In that case,\n"
);
  printf(
"  consider instead using the -u switch and writing a C procedure that\n");
  printf("  determines which triangles are too large.\n\n");
  printf(
"  If you are refining a mesh composed of linear (three-node) elements, the\n"
);
  printf(
"  output mesh contains all the nodes present in the input mesh, in the same\n"
);
  printf(
"  order, with new nodes added at the end of the .node file.  However, the\n");
  printf(
"  refinement is not hierarchical: there is no guarantee that each output\n");
  printf(
"  element is contained in a single input element.  Often, an output element\n"
);
  printf(
"  can overlap two or three input elements, and some input edges are not\n");
  printf(
"  present in the output mesh.  Hence, a sequence of refined meshes forms a\n"
);
  printf(
"  hierarchy of nodes, but not a hierarchy of elements.  If you refine a\n");
  printf(
"  mesh of higher-order elements, the hierarchical property applies only to\n"
);
  printf(
"  the nodes at the corners of an element; the midpoint nodes on each edge\n");
  printf("  are discarded before the mesh is refined.\n\n");
  printf(
"  Maximum area constraints in .poly files operate differently from those in\n"
);
  printf(
"  .area files.  A maximum area in a .poly file applies to the whole\n");
  printf(
"  (segment-bounded) region in which a point falls, whereas a maximum area\n");
  printf(
"  in an .area file applies to only one triangle.  Area constraints in .poly\n"
);
  printf(
"  files are used only when a mesh is first generated, whereas area\n");
  printf(
"  constraints in .area files are used only to refine an existing mesh, and\n"
);
  printf(
"  are typically based on a posteriori error estimates resulting from a\n");
  printf("  finite element simulation on that mesh.\n\n");
  printf(
"  `triangle -rq25 object.1' reads object.1.node and object.1.ele, then\n");
  printf(
"  refines the triangulation to enforce a 25 degree minimum angle, and then\n"
);
  printf(
"  writes the refined triangulation to object.2.node and object.2.ele.\n");
  printf("\n");
  printf(
"  `triangle -rpaa6.2 z.3' reads z.3.node, z.3.ele, z.3.poly, and z.3.area.\n"
);
  printf(
"  After reconstructing the mesh and its subsegments, Triangle refines the\n");
  printf(
"  mesh so that no triangle has area greater than 6.2, and furthermore the\n");
  printf(
"  triangles satisfy the maximum area constraints in z.3.area.  No angle\n");
  printf(
"  bound is imposed at all.  The output is written to z.4.node, z.4.ele, and\n"
);
  printf("  z.4.poly.\n\n");
  printf(
"  The sequence `triangle -qa1 x', `triangle -rqa.3 x.1', `triangle -rqa.1\n");
  printf(
"  x.2' creates a sequence of successively finer meshes x.1, x.2, and x.3,\n");
  printf("  suitable for multigrid.\n\n");
  printf("Convex Hulls and Mesh Boundaries:\n\n");
  printf(
"  If the input is a vertex set (not a PSLG), Triangle produces its convex\n");
  printf(
"  hull as a by-product in the output .poly file if you use the -c switch.\n");
  printf(
"  There are faster algorithms for finding a two-dimensional convex hull\n");
  printf("  than triangulation, of course, but this one comes for free.\n\n");
  printf(
"  If the input is an unconstrained mesh (you are using the -r switch but\n");
  printf(
"  not the -p switch), Triangle produces a list of its boundary edges\n");
  printf(
"  (including hole boundaries) as a by-product when you use the -c switch.\n");
  printf(
"  If you also use the -p switch, the output .poly file contains all the\n");
  printf("  segments from the input .poly file as well.\n\n");
  printf("Voronoi Diagrams:\n\n");
  printf(
"  The -v switch produces a Voronoi diagram, in files suffixed .v.node and\n");
  printf(
"  .v.edge.  For example, `triangle -v points' reads points.node, produces\n");
  printf(
"  its Delaunay triangulation in points.1.node and points.1.ele, and\n");
  printf(
"  produces its Voronoi diagram in points.1.v.node and points.1.v.edge.  The\n"
);
  printf(
"  .v.node file contains a list of all Voronoi vertices, and the .v.edge\n");
  printf(
"  file contains a list of all Voronoi edges, some of which may be infinite\n"
);
  printf(
"  rays.  (The choice of filenames makes it easy to run the set of Voronoi\n");
  printf("  vertices through Triangle, if so desired.)\n\n");
  printf(
"  This implementation does not use exact arithmetic to compute the Voronoi\n"
);
  printf(
"  vertices, and does not check whether neighboring vertices are identical.\n"
);
  printf(
"  Be forewarned that if the Delaunay triangulation is degenerate or\n");
  printf(
"  near-degenerate, the Voronoi diagram may have duplicate vertices or\n");
  printf("  crossing edges.\n\n");
  printf(
"  The result is a valid Voronoi diagram only if Triangle's output is a true\n"
);
  printf(
"  Delaunay triangulation.  The Voronoi output is usually meaningless (and\n");
  printf(
"  may contain crossing edges and other pathology) if the output is a CDT or\n"
);
  printf(
"  CCDT, or if it has holes or concavities.  If the triangulated domain is\n");
  printf(
"  convex and has no holes, you can use -D switch to force Triangle to\n");
  printf(
"  construct a conforming Delaunay triangulation instead of a CCDT, so the\n");
  printf("  Voronoi diagram will be valid.\n\n");
  printf("Mesh Topology:\n\n");
  printf(
"  You may wish to know which triangles are adjacent to a certain Delaunay\n");
  printf(
"  edge in an .edge file, which Voronoi cells are adjacent to a certain\n");
  printf(
"  Voronoi edge in a .v.edge file, or which Voronoi cells are adjacent to\n");
  printf(
"  each other.  All of this information can be found by cross-referencing\n");
  printf(
"  output files with the recollection that the Delaunay triangulation and\n");
  printf("  the Voronoi diagram are planar duals.\n\n");
  printf(
"  Specifically, edge i of an .edge file is the dual of Voronoi edge i of\n");
  printf(
"  the corresponding .v.edge file, and is rotated 90 degrees counterclock-\n");
  printf(
"  wise from the Voronoi edge.  Triangle j of an .ele file is the dual of\n");
  printf(
"  vertex j of the corresponding .v.node file.  Voronoi cell k is the dual\n");
  printf("  of vertex k of the corresponding .node file.\n\n");
  printf(
"  Hence, to find the triangles adjacent to a Delaunay edge, look at the\n");
  printf(
"  vertices of the corresponding Voronoi edge.  If the endpoints of a\n");
  printf(
"  Voronoi edge are Voronoi vertices 2 and 6 respectively, then triangles 2\n"
);
  printf(
"  and 6 adjoin the left and right sides of the corresponding Delaunay edge,\n"
);
  printf(
"  respectively.  To find the Voronoi cells adjacent to a Voronoi edge, look\n"
);
  printf(
"  at the endpoints of the corresponding Delaunay edge.  If the endpoints of\n"
);
  printf(
"  a Delaunay edge are input vertices 7 and 12, then Voronoi cells 7 and 12\n"
);
  printf(
"  adjoin the right and left sides of the corresponding Voronoi edge,\n");
  printf(
"  respectively.  To find which Voronoi cells are adjacent to each other,\n");
  printf("  just read the list of Delaunay edges.\n\n");
  printf(
"  Triangle does not write a list of the edges adjoining each Voronoi cell,\n"
);
  printf(
"  but you can reconstructed it straightforwardly.  For instance, to find\n");
  printf(
"  all the edges of Voronoi cell 1, search the output .edge file for every\n");
  printf(
"  edge that has input vertex 1 as an endpoint.  The corresponding dual\n");
  printf(
"  edges in the output .v.edge file form the boundary of Voronoi cell 1.\n");
  printf("\n");
  printf(
"  For each Voronoi vertex, the .neigh file gives a list of the three\n");
  printf(
"  Voronoi vertices attached to it.  You might find this more convenient\n");
  printf("  than the .v.edge file.\n\n");
  printf("Quadratic Elements:\n\n");
  printf(
"  Triangle generates meshes with subparametric quadratic elements if the\n");
  printf(
"  -o2 switch is specified.  Quadratic elements have six nodes per element,\n"
);
  printf(
"  rather than three.  `Subparametric' means that the edges of the triangles\n"
);
  printf(
"  are always straight, so that subparametric quadratic elements are\n");
  printf(
"  geometrically identical to linear elements, even though they can be used\n"
);
  printf(
"  with quadratic interpolating functions.  The three extra nodes of an\n");
  printf(
"  element fall at the midpoints of the three edges, with the fourth, fifth,\n"
);
  printf(
"  and sixth nodes appearing opposite the first, second, and third corners\n");
  printf("  respectively.\n\n");
  printf("Domains with Small Angles:\n\n");
  printf(
"  If two input segments adjoin each other at a small angle, clearly the -q\n"
);
  printf(
"  switch cannot remove the small angle.  Moreover, Triangle may have no\n");
  printf(
"  choice but to generate additional triangles whose smallest angles are\n");
  printf(
"  smaller than the specified bound.  However, these triangles only appear\n");
  printf(
"  between input segments separated by small angles.  Moreover, if you\n");
  printf(
"  request a minimum angle of theta degrees, Triangle will generally produce\n"
);
  printf(
"  no angle larger than 180 - 2 theta, even if it is forced to compromise on\n"
);
  printf("  the minimum angle.\n\n");
  printf("Statistics:\n\n");
  printf(
"  After generating a mesh, Triangle prints a count of entities in the\n");
  printf(
"  output mesh, including the number of vertices, triangles, edges, exterior\n"
);
  printf(
"  boundary edges (i.e. subsegments on the boundary of the triangulation,\n");
  printf(
"  including hole boundaries), interior boundary edges (i.e. subsegments of\n"
);
  printf(
"  input segments not on the boundary), and total subsegments.  If you've\n");
  printf(
"  forgotten the statistics for an existing mesh, run Triangle on that mesh\n"
);
  printf(
"  with the -rNEP switches to read the mesh and print the statistics without\n"
);
  printf(
"  writing any files.  Use -rpNEP if you've got a .poly file for the mesh.\n");
  printf("\n");
  printf(
"  The -V switch produces extended statistics, including a rough estimate\n");
  printf(
"  of memory use, the number of calls to geometric predicates, and\n");
  printf(
"  histograms of the angles and the aspect ratios of the triangles in the\n");
  printf("  mesh.\n\n");
  printf("Exact Arithmetic:\n\n");
  printf(
"  Triangle uses adaptive exact arithmetic to perform what computational\n");
  printf(
"  geometers call the `orientation' and `incircle' tests.  If the floating-\n"
);
  printf(
"  point arithmetic of your machine conforms to the IEEE 754 standard (as\n");
  printf(
"  most workstations do), and does not use extended precision internal\n");
  printf(
"  floating-point registers, then your output is guaranteed to be an\n");
  printf(
"  absolutely true Delaunay or constrained Delaunay triangulation, roundoff\n"
);
  printf(
"  error notwithstanding.  The word `adaptive' implies that these arithmetic\n"
);
  printf(
"  routines compute the result only to the precision necessary to guarantee\n"
);
  printf(
"  correctness, so they are usually nearly as fast as their approximate\n");
  printf("  counterparts.\n\n");
  printf(
"  May CPUs, including Intel x86 processors, have extended precision\n");
  printf(
"  floating-point registers.  These must be reconfigured so their precision\n"
);
  printf(
"  is reduced to memory precision.  Triangle does this if it is compiled\n");
  printf("  correctly.  See the makefile for details.\n\n");
  printf(
"  The exact tests can be disabled with the -X switch.  On most inputs, this\n"
);
  printf(
"  switch reduces the computation time by about eight percent--it's not\n");
  printf(
"  worth the risk.  There are rare difficult inputs (having many collinear\n");
  printf(
"  and cocircular vertices), however, for which the difference in speed\n");
  printf(
"  could be a factor of two.  Be forewarned that these are precisely the\n");
  printf(
"  inputs most likely to cause errors if you use the -X switch.  Hence, the\n"
);
  printf("  -X switch is not recommended.\n\n");
  printf(
"  Unfortunately, the exact tests don't solve every numerical problem.\n");
  printf(
"  Exact arithmetic is not used to compute the positions of new vertices,\n");
  printf(
"  because the bit complexity of vertex coordinates would grow without\n");
  printf(
"  bound.  Hence, segment intersections aren't computed exactly; in very\n");
  printf(
"  unusual cases, roundoff error in computing an intersection point might\n");
  printf(
"  actually lead to an inverted triangle and an invalid triangulation.\n");
  printf(
"  (This is one reason to specify your own intersection points in your .poly\n"
);
  printf(
"  files.)  Similarly, exact arithmetic is not used to compute the vertices\n"
);
  printf("  of the Voronoi diagram.\n\n");
  printf(
"  Another pair of problems not solved by the exact arithmetic routines is\n");
  printf(
"  underflow and overflow.  If Triangle is compiled for double precision\n");
  printf(
"  arithmetic, I believe that Triangle's geometric predicates work correctly\n"
);
  printf(
"  if the exponent of every input coordinate falls in the range [-148, 201].\n"
);
  printf(
"  Underflow can silently prevent the orientation and incircle tests from\n");
  printf(
"  being performed exactly, while overflow typically causes a floating\n");
  printf("  exception.\n\n");
  printf("Calling Triangle from Another Program:\n\n");
  printf("  Read the file triangle.h for details.\n\n");
  printf("Troubleshooting:\n\n");
  printf("  Please read this section before mailing me bugs.\n\n");
  printf("  `My output mesh has no triangles!'\n\n");
  printf(
"    If you're using a PSLG, you've probably failed to specify a proper set\n"
);
  printf(
"    of bounding segments, or forgotten to use the -c switch.  Or you may\n");
  printf(
"    have placed a hole badly, thereby eating all your triangles.  To test\n");
  printf("    these possibilities, try again with the -c and -O switches.\n");
  printf(
"    Alternatively, all your input vertices may be collinear, in which case\n"
);
  printf("    you can hardly expect to triangulate them.\n\n");
  printf("  `Triangle doesn't terminate, or just crashes.'\n\n");
  printf(
"    Bad things can happen when triangles get so small that the distance\n");
  printf(
"    between their vertices isn't much larger than the precision of your\n");
  printf(
"    machine's arithmetic.  If you've compiled Triangle for single-precision\n"
);
  printf(
"    arithmetic, you might do better by recompiling it for double-precision.\n"
);
  printf(
"    Then again, you might just have to settle for more lenient constraints\n"
);
  printf(
"    on the minimum angle and the maximum area than you had planned.\n");
  printf("\n");
  printf(
"    You can minimize precision problems by ensuring that the origin lies\n");
  printf(
"    inside your vertex set, or even inside the densest part of your\n");
  printf(
"    mesh.  If you're triangulating an object whose x-coordinates all fall\n");
  printf(
"    between 6247133 and 6247134, you're not leaving much floating-point\n");
  printf("    precision for Triangle to work with.\n\n");
  printf(
"    Precision problems can occur covertly if the input PSLG contains two\n");
  printf(
"    segments that meet (or intersect) at an extremely small angle, or if\n");
  printf(
"    such an angle is introduced by the -c switch.  If you don't realize\n");
  printf(
"    that a tiny angle is being formed, you might never discover why\n");
  printf(
"    Triangle is crashing.  To check for this possibility, use the -S switch\n"
);
  printf(
"    (with an appropriate limit on the number of Steiner points, found by\n");
  printf(
"    trial-and-error) to stop Triangle early, and view the output .poly file\n"
);
  printf(
"    with Show Me (described below).  Look carefully for regions where dense\n"
);
  printf(
"    clusters of vertices are forming and for small angles between segments.\n"
);
  printf(
"    Zoom in closely, as such segments might look like a single segment from\n"
);
  printf("    a distance.\n\n");
  printf(
"    If some of the input values are too large, Triangle may suffer a\n");
  printf(
"    floating exception due to overflow when attempting to perform an\n");
  printf(
"    orientation or incircle test.  (Read the section on exact arithmetic\n");
  printf(
"    above.)  Again, I recommend compiling Triangle for double (rather\n");
  printf("    than single) precision arithmetic.\n\n");
  printf(
"    Unexpected problems can arise if you use quality meshing (-q, -a, or\n");
  printf(
"    -u) with an input that is not segment-bounded--that is, if your input\n");
  printf(
"    is a vertex set, or you're using the -c switch.  If the convex hull of\n"
);
  printf(
"    your input vertices has collinear vertices on its boundary, an input\n");
  printf(
"    vertex that you think lies on the convex hull might actually lie just\n");
  printf(
"    inside the convex hull.  If so, the vertex and the nearby convex hull\n");
  printf(
"    edge form an extremely thin triangle.  When Triangle tries to refine\n");
  printf(
"    the mesh to enforce angle and area constraints, Triangle might generate\n"
);
  printf(
"    extremely tiny triangles, or it might fail because of insufficient\n");
  printf("    floating-point precision.\n\n");
  printf(
"  `The numbering of the output vertices doesn't match the input vertices.'\n"
);
  printf("\n");
  printf(
"    You may have had duplicate input vertices, or you may have eaten some\n");
  printf(
"    of your input vertices with a hole, or by placing them outside the area\n"
);
  printf(
"    enclosed by segments.  In any case, you can solve the problem by not\n");
  printf("    using the -j switch.\n\n");
  printf(
"  `Triangle executes without incident, but when I look at the resulting\n");
  printf(
"  mesh, it has overlapping triangles or other geometric inconsistencies.'\n");
  printf("\n");
  printf(
"    If you select the -X switch, Triangle occasionally makes mistakes due\n");
  printf(
"    to floating-point roundoff error.  Although these errors are rare,\n");
  printf(
"    don't use the -X switch.  If you still have problems, please report the\n"
);
  printf("    bug.\n\n");
  printf(
"  `Triangle executes without incident, but when I look at the resulting\n");
  printf("  Voronoi diagram, it has overlapping edges or other geometric\n");
  printf("  inconsistencies.'\n");
  printf("\n");
  printf(
"    If your input is a PSLG (-p), you can only expect a meaningful Voronoi\n"
);
  printf(
"    diagram if the domain you are triangulating is convex and free of\n");
  printf(
"    holes, and you use the -D switch to construct a conforming Delaunay\n");
  printf("    triangulation (instead of a CDT or CCDT).\n\n");
  printf(
"  Strange things can happen if you've taken liberties with your PSLG.  Do\n");
  printf(
"  you have a vertex lying in the middle of a segment?  Triangle sometimes\n");
  printf(
"  copes poorly with that sort of thing.  Do you want to lay out a collinear\n"
);
  printf(
"  row of evenly spaced, segment-connected vertices?  Have you simply\n");
  printf(
"  defined one long segment connecting the leftmost vertex to the rightmost\n"
);
  printf(
"  vertex, and a bunch of vertices lying along it?  This method occasionally\n"
);
  printf(
"  works, especially with horizontal and vertical lines, but often it\n");
  printf(
"  doesn't, and you'll have to connect each adjacent pair of vertices with a\n"
);
  printf("  separate segment.  If you don't like it, tough.\n\n");
  printf(
"  Furthermore, if you have segments that intersect other than at their\n");
  printf(
"  endpoints, try not to let the intersections fall extremely close to PSLG\n"
);
  printf("  vertices or each other.\n\n");
  printf(
"  If you have problems refining a triangulation not produced by Triangle:\n");
  printf(
"  Are you sure the triangulation is geometrically valid?  Is it formatted\n");
  printf(
"  correctly for Triangle?  Are the triangles all listed so the first three\n"
);
  printf(
"  vertices are their corners in counterclockwise order?  Are all of the\n");
  printf(
"  triangles constrained Delaunay?  Triangle's Delaunay refinement algorithm\n"
);
  printf("  assumes that it starts with a CDT.\n\n");
  printf("Show Me:\n\n");
  printf(
"  Triangle comes with a separate program named `Show Me', whose primary\n");
  printf(
"  purpose is to draw meshes on your screen or in PostScript.  Its secondary\n"
);
  printf(
"  purpose is to check the validity of your input files, and do so more\n");
  printf(
"  thoroughly than Triangle does.  Unlike Triangle, Show Me requires that\n");
  printf(
"  you have the X Windows system.  Sorry, Microsoft Windows users.\n");
  printf("\n");
  printf("Triangle on the Web:\n");
  printf("\n");
  printf("  To see an illustrated version of these instructions, check out\n");
  printf("\n");
  printf("    http://www.cs.cmu.edu/~quake/triangle.html\n");
  printf("\n");
  printf("A Brief Plea:\n");
  printf("\n");
  printf(
"  If you use Triangle, and especially if you use it to accomplish real\n");
  printf(
"  work, I would like very much to hear from you.  A short letter or email\n");
  printf(
"  (to jrs@cs.berkeley.edu) describing how you use Triangle will mean a lot\n"
);
  printf(
"  to me.  The more people I know are using this program, the more easily I\n"
);
  printf(
"  can justify spending time on improvements, which in turn will benefit\n");
  printf(
"  you.  Also, I can put you on a list to receive email whenever a new\n");
  printf("  version of Triangle is available.\n\n");
  printf(
"  If you use a mesh generated by Triangle in a publication, please include\n"
);
  printf(
"  an acknowledgment as well.  And please spell Triangle with a capital `T'!\n"
);
  printf(
"  If you want to include a citation, use `Jonathan Richard Shewchuk,\n");
  printf(
"  ``Triangle: Engineering a 2D Quality Mesh Generator and Delaunay\n");
  printf(
"  Triangulator,'' in Applied Computational Geometry:  Towards Geometric\n");
  printf(
"  Engineering (Ming C. Lin and Dinesh Manocha, editors), volume 1148 of\n");
  printf(
"  Lecture Notes in Computer Science, pages 203-222, Springer-Verlag,\n");
  printf(
"  Berlin, May 1996.  (From the First ACM Workshop on Applied Computational\n"
);
  printf("  Geometry.)'\n\n");
  printf("Research credit:\n\n");
  printf(
"  Of course, I can take credit for only a fraction of the ideas that made\n");
  printf(
"  this mesh generator possible.  Triangle owes its existence to the efforts\n"
);
  printf(
"  of many fine computational geometers and other researchers, including\n");
  printf(
"  Marshall Bern, L. Paul Chew, Kenneth L. Clarkson, Boris Delaunay, Rex A.\n"
);
  printf(
"  Dwyer, David Eppstein, Steven Fortune, Leonidas J. Guibas, Donald E.\n");
  printf(
"  Knuth, Charles L. Lawson, Der-Tsai Lee, Gary L. Miller, Ernst P. Mucke,\n");
  printf(
"  Steven E. Pav, Douglas M. Priest, Jim Ruppert, Isaac Saias, Bruce J.\n");
  printf(
"  Schachter, Micha Sharir, Peter W. Shor, Daniel D. Sleator, Jorge Stolfi,\n"
);
  printf("  Robert E. Tarjan, Alper Ungor, Christopher J. Van Wyk, Noel J.\n");
  printf(
"  Walkington, and Binhai Zhu.  See the comments at the beginning of the\n");
  printf("  source code for references.\n\n");
  triexit(0);
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  internalerror()   Ask the user to send me the defective product.  Exit.  */
/*                                                                           */
/*****************************************************************************/

void internalerror()
{
  printf("  Please report this bug to jrs@cs.berkeley.edu\n");
  printf("  Include the message above, your input data set, and the exact\n");
  printf("    command line you used to run Triangle.\n");
  triexit(1);
}

/*****************************************************************************/
/*                                                                           */
/*  parsecommandline()   Read the command line, identify switches, and set   */
/*                       up options and file names.                          */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void parsecommandline(int argc, char **argv, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void parsecommandline(argc, argv, b)
int argc;
char **argv;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
#ifdef TRILIBRARY
#define STARTINDEX 0
#else /* not TRILIBRARY */
#define STARTINDEX 1
  int increment;
  int meshnumber;
#endif /* not TRILIBRARY */
  int i, j, k;
  char workstring[FILENAMESIZE];

  b->poly = b->refine = b->quality = 0;
  b->vararea = b->fixedarea = b->usertest = 0;
  b->regionattrib = b->convex = b->weighted = b->jettison = 0;
  b->firstnumber = 1;
  b->edgesout = b->voronoi = b->neighbors = b->geomview = 0;
  b->nobound = b->nopolywritten = b->nonodewritten = b->noelewritten = 0;
  b->noiterationnum = 0;
  b->noholes = b->noexact = 0;
  b->incremental = b->sweepline = 0;
  b->dwyer = 1;
  b->splitseg = 0;
  b->docheck = 0;
  b->nobisect = 0;
  b->conformdel = 0;
  b->steiner = -1;
  b->order = 1;
  b->minangle = 0.0;
  b->maxarea = -1.0;
  b->quiet = b->verbose = 0;
#ifndef TRILIBRARY
  b->innodefilename[0] = '\0';
#endif /* not TRILIBRARY */

  for (i = STARTINDEX; i < argc; i++) {
#ifndef TRILIBRARY
    if (argv[i][0] == '-') {
#endif /* not TRILIBRARY */
      for (j = STARTINDEX; argv[i][j] != '\0'; j++) {
        if (argv[i][j] == 'p') {
          b->poly = 1;
	}
#ifndef CDT_ONLY
        if (argv[i][j] == 'r') {
          b->refine = 1;
	}
        if (argv[i][j] == 'q') {
          b->quality = 1;
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
            b->minangle = (REAL) strtod(workstring, (char **) NULL);
	  } else {
            b->minangle = 20.0;
	  }
	}
        if (argv[i][j] == 'a') {
          b->quality = 1;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            b->fixedarea = 1;
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            b->maxarea = (REAL) strtod(workstring, (char **) NULL);
            if (b->maxarea <= 0.0) {
              printf("Error:  Maximum area must be greater than zero.\n");
              triexit(1);
	    }
	  } else {
            b->vararea = 1;
	  }
	}
        if (argv[i][j] == 'u') {
          b->quality = 1;
          b->usertest = 1;
        }
#endif /* not CDT_ONLY */
        if (argv[i][j] == 'A') {
          b->regionattrib = 1;
        }
        if (argv[i][j] == 'c') {
          b->convex = 1;
        }
        if (argv[i][j] == 'w') {
          b->weighted = 1;
        }
        if (argv[i][j] == 'W') {
          b->weighted = 2;
        }
        if (argv[i][j] == 'j') {
          b->jettison = 1;
        }
        if (argv[i][j] == 'z') {
          b->firstnumber = 0;
        }
        if (argv[i][j] == 'e') {
          b->edgesout = 1;
	}
        if (argv[i][j] == 'v') {
          b->voronoi = 1;
	}
        if (argv[i][j] == 'n') {
          b->neighbors = 1;
	}
        if (argv[i][j] == 'g') {
          b->geomview = 1;
	}
        if (argv[i][j] == 'B') {
          b->nobound = 1;
	}
        if (argv[i][j] == 'P') {
          b->nopolywritten = 1;
	}
        if (argv[i][j] == 'N') {
          b->nonodewritten = 1;
	}
        if (argv[i][j] == 'E') {
          b->noelewritten = 1;
	}
#ifndef TRILIBRARY
        if (argv[i][j] == 'I') {
          b->noiterationnum = 1;
	}
#endif /* not TRILIBRARY */
        if (argv[i][j] == 'O') {
          b->noholes = 1;
	}
        if (argv[i][j] == 'X') {
          b->noexact = 1;
	}
        if (argv[i][j] == 'o') {
          if (argv[i][j + 1] == '2') {
            j++;
            b->order = 2;
          }
	}
#ifndef CDT_ONLY
        if (argv[i][j] == 'Y') {
          b->nobisect++;
	}
        if (argv[i][j] == 'S') {
          b->steiner = 0;
          while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
            j++;
            b->steiner = b->steiner * 10 + (int) (argv[i][j] - '0');
          }
        }
#endif /* not CDT_ONLY */
#ifndef REDUCED
        if (argv[i][j] == 'i') {
          b->incremental = 1;
        }
        if (argv[i][j] == 'F') {
          b->sweepline = 1;
        }
#endif /* not REDUCED */
        if (argv[i][j] == 'l') {
          b->dwyer = 0;
        }
#ifndef REDUCED
#ifndef CDT_ONLY
        if (argv[i][j] == 's') {
          b->splitseg = 1;
        }
        if ((argv[i][j] == 'D') || (argv[i][j] == 'L')) {
          b->quality = 1;
          b->conformdel = 1;
        }
#endif /* not CDT_ONLY */
        if (argv[i][j] == 'C') {
          b->docheck = 1;
        }
#endif /* not REDUCED */
        if (argv[i][j] == 'Q') {
          b->quiet = 1;
        }
        if (argv[i][j] == 'V') {
          b->verbose++;
        }
#ifndef TRILIBRARY
        if ((argv[i][j] == 'h') || (argv[i][j] == 'H') ||
            (argv[i][j] == '?')) {
          info();
	}
#endif /* not TRILIBRARY */
      }
#ifndef TRILIBRARY
    } else {
      strncpy(b->innodefilename, argv[i], FILENAMESIZE - 1);
      b->innodefilename[FILENAMESIZE - 1] = '\0';
    }
#endif /* not TRILIBRARY */
  }
#ifndef TRILIBRARY
  if (b->innodefilename[0] == '\0') {
    syntax();
  }
  if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".node")) {
    b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
  }
  if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".poly")) {
    b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
    b->poly = 1;
  }
#ifndef CDT_ONLY
  if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 4], ".ele")) {
    b->innodefilename[strlen(b->innodefilename) - 4] = '\0';
    b->refine = 1;
  }
  if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".area")) {
    b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
    b->refine = 1;
    b->quality = 1;
    b->vararea = 1;
  }
#endif /* not CDT_ONLY */
#endif /* not TRILIBRARY */
  b->usesegments = b->poly || b->refine || b->quality || b->convex;
  b->goodangle = cos(b->minangle * PI / 180.0);
  if (b->goodangle == 1.0) {
    b->offconstant = 0.0;
  } else {
    b->offconstant = 0.475 * sqrt((1.0 + b->goodangle) / (1.0 - b->goodangle));
  }
  b->goodangle *= b->goodangle;
  if (b->refine && b->noiterationnum) {
    printf(
      "Error:  You cannot use the -I switch when refining a triangulation.\n");
    triexit(1);
  }
  /* Be careful not to allocate space for element area constraints that */
  /*   will never be assigned any value (other than the default -1.0).  */
  if (!b->refine && !b->poly) {
    b->vararea = 0;
  }
  /* Be careful not to add an extra attribute to each element unless the */
  /*   input supports it (PSLG in, but not refining a preexisting mesh). */
  if (b->refine || !b->poly) {
    b->regionattrib = 0;
  }
  /* Regular/weighted triangulations are incompatible with PSLGs */
  /*   and meshing.                                              */
  if (b->weighted && (b->poly || b->quality)) {
    b->weighted = 0;
    if (!b->quiet) {
      printf("Warning:  weighted triangulations (-w, -W) are incompatible\n");
      printf("  with PSLGs (-p) and meshing (-q, -a, -u).  Weights ignored.\n"
             );
    }
  }
  if (b->jettison && b->nonodewritten && !b->quiet) {
    printf("Warning:  -j and -N switches are somewhat incompatible.\n");
    printf("  If any vertices are jettisoned, you will need the output\n");
    printf("  .node file to reconstruct the new node indices.");
  }

#ifndef TRILIBRARY
  strcpy(b->inpolyfilename, b->innodefilename);
  strcpy(b->inelefilename, b->innodefilename);
  strcpy(b->areafilename, b->innodefilename);
  increment = 0;
  strcpy(workstring, b->innodefilename);
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
  if (b->noiterationnum) {
    strcpy(b->outnodefilename, b->innodefilename);
    strcpy(b->outelefilename, b->innodefilename);
    strcpy(b->edgefilename, b->innodefilename);
    strcpy(b->vnodefilename, b->innodefilename);
    strcpy(b->vedgefilename, b->innodefilename);
    strcpy(b->neighborfilename, b->innodefilename);
    strcpy(b->offfilename, b->innodefilename);
    strcat(b->outnodefilename, ".node");
    strcat(b->outelefilename, ".ele");
    strcat(b->edgefilename, ".edge");
    strcat(b->vnodefilename, ".v.node");
    strcat(b->vedgefilename, ".v.edge");
    strcat(b->neighborfilename, ".neigh");
    strcat(b->offfilename, ".off");
  } else if (increment == 0) {
    strcpy(b->outnodefilename, b->innodefilename);
    strcpy(b->outpolyfilename, b->innodefilename);
    strcpy(b->outelefilename, b->innodefilename);
    strcpy(b->edgefilename, b->innodefilename);
    strcpy(b->vnodefilename, b->innodefilename);
    strcpy(b->vedgefilename, b->innodefilename);
    strcpy(b->neighborfilename, b->innodefilename);
    strcpy(b->offfilename, b->innodefilename);
    strcat(b->outnodefilename, ".1.node");
    strcat(b->outpolyfilename, ".1.poly");
    strcat(b->outelefilename, ".1.ele");
    strcat(b->edgefilename, ".1.edge");
    strcat(b->vnodefilename, ".1.v.node");
    strcat(b->vedgefilename, ".1.v.edge");
    strcat(b->neighborfilename, ".1.neigh");
    strcat(b->offfilename, ".1.off");
  } else {
    workstring[increment] = '%';
    workstring[increment + 1] = 'd';
    workstring[increment + 2] = '\0';
    sprintf(b->outnodefilename, workstring, meshnumber + 1);
    strcpy(b->outpolyfilename, b->outnodefilename);
    strcpy(b->outelefilename, b->outnodefilename);
    strcpy(b->edgefilename, b->outnodefilename);
    strcpy(b->vnodefilename, b->outnodefilename);
    strcpy(b->vedgefilename, b->outnodefilename);
    strcpy(b->neighborfilename, b->outnodefilename);
    strcpy(b->offfilename, b->outnodefilename);
    strcat(b->outnodefilename, ".node");
    strcat(b->outpolyfilename, ".poly");
    strcat(b->outelefilename, ".ele");
    strcat(b->edgefilename, ".edge");
    strcat(b->vnodefilename, ".v.node");
    strcat(b->vedgefilename, ".v.edge");
    strcat(b->neighborfilename, ".neigh");
    strcat(b->offfilename, ".off");
  }
  strcat(b->innodefilename, ".node");
  strcat(b->inpolyfilename, ".poly");
  strcat(b->inelefilename, ".ele");
  strcat(b->areafilename, ".area");
#endif /* not TRILIBRARY */
}

/**                                                                         **/
/**                                                                         **/
/********* User interaction routines begin here                      *********/

/********* Debugging routines begin here                             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  printtriangle()   Print out the details of an oriented triangle.         */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about an     */
/*  oriented triangle in digestible form.  It's also used when the           */
/*  highest level of verbosity (`-VVV') is specified.                        */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void printtriangle(struct mesh *m, struct behavior *b, struct otri *t)
#else /* not ANSI_DECLARATORS */
void printtriangle(m, b, t)
struct mesh *m;
struct behavior *b;
struct otri *t;
#endif /* not ANSI_DECLARATORS */

{
  struct otri printtri;
  struct osub printsh;
  vertex printvertex;

  printf("triangle x%lx with orientation %d:\n", (unsigned long) t->tri,
         t->orient);
  decode(t->tri[0], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [0] = Outer space\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[1], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [1] = Outer space\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(t->tri[2], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [2] = Outer space\n");
  } else {
    printf("    [2] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }

  org(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Origin[%d] = NULL\n", (t->orient + 1) % 3 + 3);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 1) % 3 + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  dest(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Dest  [%d] = NULL\n", (t->orient + 2) % 3 + 3);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           (t->orient + 2) % 3 + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  apex(*t, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Apex  [%d] = NULL\n", t->orient + 3);
  else
    printf("    Apex  [%d] = x%lx  (%.12g, %.12g)\n",
           t->orient + 3, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);

  if (b->usesegments) {
    sdecode(t->tri[6], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [6] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
    sdecode(t->tri[7], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [7] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
    sdecode(t->tri[8], printsh);
    if (printsh.ss != m->dummysub) {
      printf("    [8] = x%lx  %d\n", (unsigned long) printsh.ss,
             printsh.ssorient);
    }
  }

  if (b->vararea) {
    printf("    Area constraint:  %.4g\n", areabound(*t));
  }
}

/*****************************************************************************/
/*                                                                           */
/*  printsubseg()   Print out the details of an oriented subsegment.         */
/*                                                                           */
/*  I originally wrote this procedure to simplify debugging; it can be       */
/*  called directly from the debugger, and presents information about an     */
/*  oriented subsegment in digestible form.  It's also used when the highest */
/*  level of verbosity (`-VVV') is specified.                                */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void printsubseg(struct mesh *m, struct behavior *b, struct osub *s)
#else /* not ANSI_DECLARATORS */
void printsubseg(m, b, s)
struct mesh *m;
struct behavior *b;
struct osub *s;
#endif /* not ANSI_DECLARATORS */

{
  struct osub printsh;
  struct otri printtri;
  vertex printvertex;

  printf("subsegment x%lx with orientation %d and mark %d:\n",
         (unsigned long) s->ss, s->ssorient, mark(*s));
  sdecode(s->ss[0], printsh);
  if (printsh.ss == m->dummysub) {
    printf("    [0] = No subsegment\n");
  } else {
    printf("    [0] = x%lx  %d\n", (unsigned long) printsh.ss,
           printsh.ssorient);
  }
  sdecode(s->ss[1], printsh);
  if (printsh.ss == m->dummysub) {
    printf("    [1] = No subsegment\n");
  } else {
    printf("    [1] = x%lx  %d\n", (unsigned long) printsh.ss,
           printsh.ssorient);
  }

  sorg(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Origin[%d] = NULL\n", 2 + s->ssorient);
  else
    printf("    Origin[%d] = x%lx  (%.12g, %.12g)\n",
           2 + s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  sdest(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Dest  [%d] = NULL\n", 3 - s->ssorient);
  else
    printf("    Dest  [%d] = x%lx  (%.12g, %.12g)\n",
           3 - s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);

  decode(s->ss[6], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [6] = Outer space\n");
  } else {
    printf("    [6] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }
  decode(s->ss[7], printtri);
  if (printtri.tri == m->dummytri) {
    printf("    [7] = Outer space\n");
  } else {
    printf("    [7] = x%lx  %d\n", (unsigned long) printtri.tri,
           printtri.orient);
  }

  segorg(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Segment origin[%d] = NULL\n", 4 + s->ssorient);
  else
    printf("    Segment origin[%d] = x%lx  (%.12g, %.12g)\n",
           4 + s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
  segdest(*s, printvertex);
  if (printvertex == (vertex) NULL)
    printf("    Segment dest  [%d] = NULL\n", 5 - s->ssorient);
  else
    printf("    Segment dest  [%d] = x%lx  (%.12g, %.12g)\n",
           5 - s->ssorient, (unsigned long) printvertex,
           printvertex[0], printvertex[1]);
}

/**                                                                         **/
/**                                                                         **/
/********* Debugging routines end here                               *********/

/********* Memory management routines begin here                     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  poolzero()   Set all of a pool's fields to zero.                         */
/*                                                                           */
/*  This procedure should never be called on a pool that has any memory      */
/*  allocated to it, as that memory would leak.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void poolzero(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
void poolzero(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  pool->firstblock = (VOID **) NULL;
  pool->nowblock = (VOID **) NULL;
  pool->nextitem = (VOID *) NULL;
  pool->deaditemstack = (VOID *) NULL;
  pool->pathblock = (VOID **) NULL;
  pool->pathitem = (VOID *) NULL;
  pool->alignbytes = 0;
  pool->itembytes = 0;
  pool->itemsperblock = 0;
  pool->itemsfirstblock = 0;
  pool->items = 0;
  pool->maxitems = 0;
  pool->unallocateditems = 0;
  pool->pathitemsleft = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  poolrestart()   Deallocate all items in a pool.                          */
/*                                                                           */
/*  The pool is returned to its starting state, except that no memory is     */
/*  freed to the operating system.  Rather, the previously allocated blocks  */
/*  are ready to be reused.                                                  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void poolrestart(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
void poolrestart(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  unsigned long alignptr;

  pool->items = 0;
  pool->maxitems = 0;

  /* Set the currently active block. */
  pool->nowblock = pool->firstblock;
  /* Find the first item in the pool.  Increment by the size of (VOID *). */
  alignptr = (unsigned long) (pool->nowblock + 1);
  /* Align the item on an `alignbytes'-byte boundary. */
  pool->nextitem = (VOID *)
    (alignptr + (unsigned long) pool->alignbytes -
     (alignptr % (unsigned long) pool->alignbytes));
  /* There are lots of unallocated items left in this block. */
  pool->unallocateditems = pool->itemsfirstblock;
  /* The stack of deallocated items is empty. */
  pool->deaditemstack = (VOID *) NULL;
}

/*****************************************************************************/
/*                                                                           */
/*  poolinit()   Initialize a pool of memory for allocation of items.        */
/*                                                                           */
/*  This routine initializes the machinery for allocating items.  A `pool'   */
/*  is created whose records have size at least `bytecount'.  Items will be  */
/*  allocated in `itemcount'-item blocks.  Each item is assumed to be a      */
/*  collection of words, and either pointers or floating-point values are    */
/*  assumed to be the "primary" word type.  (The "primary" word type is used */
/*  to determine alignment of items.)  If `alignment' isn't zero, all items  */
/*  will be `alignment'-byte aligned in memory.  `alignment' must be either  */
/*  a multiple or a factor of the primary word size; powers of two are safe. */
/*  `alignment' is normally used to create a few unused bits at the bottom   */
/*  of each item's pointer, in which information may be stored.              */
/*                                                                           */
/*  Don't change this routine unless you understand it.                      */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void poolinit(struct memorypool *pool, int bytecount, int itemcount,
              int firstitemcount, int alignment)
#else /* not ANSI_DECLARATORS */
void poolinit(pool, bytecount, itemcount, firstitemcount, alignment)
struct memorypool *pool;
int bytecount;
int itemcount;
int firstitemcount;
int alignment;
#endif /* not ANSI_DECLARATORS */

{
  /* Find the proper alignment, which must be at least as large as:   */
  /*   - The parameter `alignment'.                                   */
  /*   - sizeof(VOID *), so the stack of dead items can be maintained */
  /*       without unaligned accesses.                                */
  if (alignment > sizeof(VOID *)) {
    pool->alignbytes = alignment;
  } else {
    pool->alignbytes = sizeof(VOID *);
  }
  pool->itembytes = ((bytecount - 1) / pool->alignbytes + 1) *
                    pool->alignbytes;
  pool->itemsperblock = itemcount;
  if (firstitemcount == 0) {
    pool->itemsfirstblock = itemcount;
  } else {
    pool->itemsfirstblock = firstitemcount;
  }

  /* Allocate a block of items.  Space for `itemsfirstblock' items and one  */
  /*   pointer (to point to the next block) are allocated, as well as space */
  /*   to ensure alignment of the items.                                    */
  pool->firstblock = (VOID **)
    trimalloc(pool->itemsfirstblock * pool->itembytes + (int) sizeof(VOID *) +
              pool->alignbytes);
  /* Set the next block pointer to NULL. */
  *(pool->firstblock) = (VOID *) NULL;
  poolrestart(pool);
}

/*****************************************************************************/
/*                                                                           */
/*  pooldeinit()   Free to the operating system all memory taken by a pool.  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void pooldeinit(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
void pooldeinit(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  while (pool->firstblock != (VOID **) NULL) {
    pool->nowblock = (VOID **) *(pool->firstblock);
    trifree((VOID *) pool->firstblock);
    pool->firstblock = pool->nowblock;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  poolalloc()   Allocate space for an item.                                */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
VOID *poolalloc(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
VOID *poolalloc(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  VOID *newitem;
  VOID **newblock;
  unsigned long alignptr;

  /* First check the linked list of dead items.  If the list is not   */
  /*   empty, allocate an item from the list rather than a fresh one. */
  if (pool->deaditemstack != (VOID *) NULL) {
    newitem = pool->deaditemstack;               /* Take first item in list. */
    pool->deaditemstack = * (VOID **) pool->deaditemstack;
  } else {
    /* Check if there are any free items left in the current block. */
    if (pool->unallocateditems == 0) {
      /* Check if another block must be allocated. */
      if (*(pool->nowblock) == (VOID *) NULL) {
        /* Allocate a new block of items, pointed to by the previous block. */
        newblock = (VOID **) trimalloc(pool->itemsperblock * pool->itembytes +
                                       (int) sizeof(VOID *) +
                                       pool->alignbytes);
        *(pool->nowblock) = (VOID *) newblock;
        /* The next block pointer is NULL. */
        *newblock = (VOID *) NULL;
      }

      /* Move to the new block. */
      pool->nowblock = (VOID **) *(pool->nowblock);
      /* Find the first item in the block.    */
      /*   Increment by the size of (VOID *). */
      alignptr = (unsigned long) (pool->nowblock + 1);
      /* Align the item on an `alignbytes'-byte boundary. */
      pool->nextitem = (VOID *)
        (alignptr + (unsigned long) pool->alignbytes -
         (alignptr % (unsigned long) pool->alignbytes));
      /* There are lots of unallocated items left in this block. */
      pool->unallocateditems = pool->itemsperblock;
    }

    /* Allocate a new item. */
    newitem = pool->nextitem;
    /* Advance `nextitem' pointer to next free item in block. */
    pool->nextitem = (VOID *) ((char *) pool->nextitem + pool->itembytes);
    pool->unallocateditems--;
    pool->maxitems++;
  }
  pool->items++;
  return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  pooldealloc()   Deallocate space for an item.                            */
/*                                                                           */
/*  The deallocated space is stored in a queue for later reuse.              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void pooldealloc(struct memorypool *pool, VOID *dyingitem)
#else /* not ANSI_DECLARATORS */
void pooldealloc(pool, dyingitem)
struct memorypool *pool;
VOID *dyingitem;
#endif /* not ANSI_DECLARATORS */

{
  /* Push freshly killed item onto stack. */
  *((VOID **) dyingitem) = pool->deaditemstack;
  pool->deaditemstack = dyingitem;
  pool->items--;
}

/*****************************************************************************/
/*                                                                           */
/*  traversalinit()   Prepare to traverse the entire list of items.          */
/*                                                                           */
/*  This routine is used in conjunction with traverse().                     */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void traversalinit(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
void traversalinit(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  unsigned long alignptr;

  /* Begin the traversal in the first block. */
  pool->pathblock = pool->firstblock;
  /* Find the first item in the block.  Increment by the size of (VOID *). */
  alignptr = (unsigned long) (pool->pathblock + 1);
  /* Align with item on an `alignbytes'-byte boundary. */
  pool->pathitem = (VOID *)
    (alignptr + (unsigned long) pool->alignbytes -
     (alignptr % (unsigned long) pool->alignbytes));
  /* Set the number of items left in the current block. */
  pool->pathitemsleft = pool->itemsfirstblock;
}

/*****************************************************************************/
/*                                                                           */
/*  traverse()   Find the next item in the list.                             */
/*                                                                           */
/*  This routine is used in conjunction with traversalinit().  Be forewarned */
/*  that this routine successively returns all items in the list, including  */
/*  deallocated ones on the deaditemqueue.  It's up to you to figure out     */
/*  which ones are actually dead.  Why?  I don't want to allocate extra      */
/*  space just to demarcate dead items.  It can usually be done more         */
/*  space-efficiently by a routine that knows something about the structure  */
/*  of the item.                                                             */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
VOID *traverse(struct memorypool *pool)
#else /* not ANSI_DECLARATORS */
VOID *traverse(pool)
struct memorypool *pool;
#endif /* not ANSI_DECLARATORS */

{
  VOID *newitem;
  unsigned long alignptr;

  /* Stop upon exhausting the list of items. */
  if (pool->pathitem == pool->nextitem) {
    return (VOID *) NULL;
  }

  /* Check whether any untraversed items remain in the current block. */
  if (pool->pathitemsleft == 0) {
    /* Find the next block. */
    pool->pathblock = (VOID **) *(pool->pathblock);
    /* Find the first item in the block.  Increment by the size of (VOID *). */
    alignptr = (unsigned long) (pool->pathblock + 1);
    /* Align with item on an `alignbytes'-byte boundary. */
    pool->pathitem = (VOID *)
      (alignptr + (unsigned long) pool->alignbytes -
       (alignptr % (unsigned long) pool->alignbytes));
    /* Set the number of items left in the current block. */
    pool->pathitemsleft = pool->itemsperblock;
  }

  newitem = pool->pathitem;
  /* Find the next item in the block. */
  pool->pathitem = (VOID *) ((char *) pool->pathitem + pool->itembytes);
  pool->pathitemsleft--;
  return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  dummyinit()   Initialize the triangle that fills "outer space" and the   */
/*                omnipresent subsegment.                                    */
/*                                                                           */
/*  The triangle that fills "outer space," called `dummytri', is pointed to  */
/*  by every triangle and subsegment on a boundary (be it outer or inner) of */
/*  the triangulation.  Also, `dummytri' points to one of the triangles on   */
/*  the convex hull (until the holes and concavities are carved), making it  */
/*  possible to find a starting triangle for point location.                 */
/*                                                                           */
/*  The omnipresent subsegment, `dummysub', is pointed to by every triangle  */
/*  or subsegment that doesn't have a full complement of real subsegments    */
/*  to point to.                                                             */
/*                                                                           */
/*  `dummytri' and `dummysub' are generally required to fulfill only a few   */
/*  invariants:  their vertices must remain NULL and `dummytri' must always  */
/*  be bonded (at offset zero) to some triangle on the convex hull of the    */
/*  mesh, via a boundary edge.  Otherwise, the connections of `dummytri' and */
/*  `dummysub' may change willy-nilly.  This makes it possible to avoid      */
/*  writing a good deal of special-case code (in the edge flip, for example) */
/*  for dealing with the boundary of the mesh, places where no subsegment is */
/*  present, and so forth.  Other entities are frequently bonded to          */
/*  `dummytri' and `dummysub' as if they were real mesh entities, with no    */
/*  harm done.                                                               */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void dummyinit(struct mesh *m, struct behavior *b, int trianglebytes,
               int subsegbytes)
#else /* not ANSI_DECLARATORS */
void dummyinit(m, b, trianglebytes, subsegbytes)
struct mesh *m;
struct behavior *b;
int trianglebytes;
int subsegbytes;
#endif /* not ANSI_DECLARATORS */

{
  unsigned long alignptr;

  /* Set up `dummytri', the `triangle' that occupies "outer space." */
  m->dummytribase = (triangle *) trimalloc(trianglebytes +
                                           m->triangles.alignbytes);
  /* Align `dummytri' on a `triangles.alignbytes'-byte boundary. */
  alignptr = (unsigned long) m->dummytribase;
  m->dummytri = (triangle *)
    (alignptr + (unsigned long) m->triangles.alignbytes -
     (alignptr % (unsigned long) m->triangles.alignbytes));
  /* Initialize the three adjoining triangles to be "outer space."  These  */
  /*   will eventually be changed by various bonding operations, but their */
  /*   values don't really matter, as long as they can legally be          */
  /*   dereferenced.                                                       */
  m->dummytri[0] = (triangle) m->dummytri;
  m->dummytri[1] = (triangle) m->dummytri;
  m->dummytri[2] = (triangle) m->dummytri;
  /* Three NULL vertices. */
  m->dummytri[3] = (triangle) NULL;
  m->dummytri[4] = (triangle) NULL;
  m->dummytri[5] = (triangle) NULL;

  if (b->usesegments) {
    /* Set up `dummysub', the omnipresent subsegment pointed to by any */
    /*   triangle side or subsegment end that isn't attached to a real */
    /*   subsegment.                                                   */
    m->dummysubbase = (subseg *) trimalloc(subsegbytes +
                                           m->subsegs.alignbytes);
    /* Align `dummysub' on a `subsegs.alignbytes'-byte boundary. */
    alignptr = (unsigned long) m->dummysubbase;
    m->dummysub = (subseg *)
      (alignptr + (unsigned long) m->subsegs.alignbytes -
       (alignptr % (unsigned long) m->subsegs.alignbytes));
    /* Initialize the two adjoining subsegments to be the omnipresent      */
    /*   subsegment.  These will eventually be changed by various bonding  */
    /*   operations, but their values don't really matter, as long as they */
    /*   can legally be dereferenced.                                      */
    m->dummysub[0] = (subseg) m->dummysub;
    m->dummysub[1] = (subseg) m->dummysub;
    /* Four NULL vertices. */
    m->dummysub[2] = (subseg) NULL;
    m->dummysub[3] = (subseg) NULL;
    m->dummysub[4] = (subseg) NULL;
    m->dummysub[5] = (subseg) NULL;
    /* Initialize the two adjoining triangles to be "outer space." */
    m->dummysub[6] = (subseg) m->dummytri;
    m->dummysub[7] = (subseg) m->dummytri;
    /* Set the boundary marker to zero. */
    * (int *) (m->dummysub + 8) = 0;

    /* Initialize the three adjoining subsegments of `dummytri' to be */
    /*   the omnipresent subsegment.                                  */
    m->dummytri[6] = (triangle) m->dummysub;
    m->dummytri[7] = (triangle) m->dummysub;
    m->dummytri[8] = (triangle) m->dummysub;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  initializevertexpool()   Calculate the size of the vertex data structure */
/*                           and initialize its memory pool.                 */
/*                                                                           */
/*  This routine also computes the `vertexmarkindex' and `vertex2triindex'   */
/*  indices used to find values within each vertex.                          */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void initializevertexpool(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void initializevertexpool(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  int vertexsize;

  /* The index within each vertex at which the boundary marker is found,    */
  /*   followed by the vertex type.  Ensure the vertex marker is aligned to */
  /*   a sizeof(int)-byte address.                                          */
  m->vertexmarkindex = ((m->mesh_dim + m->nextras) * sizeof(REAL) +
                        sizeof(int) - 1) /
                       sizeof(int);
  vertexsize = (m->vertexmarkindex + 2) * sizeof(int);
  if (b->poly) {
    /* The index within each vertex at which a triangle pointer is found.  */
    /*   Ensure the pointer is aligned to a sizeof(triangle)-byte address. */
    m->vertex2triindex = (vertexsize + sizeof(triangle) - 1) /
                         sizeof(triangle);
    vertexsize = (m->vertex2triindex + 1) * sizeof(triangle);
  }

  /* Initialize the pool of vertices. */
  poolinit(&m->vertices, vertexsize, VERTEXPERBLOCK,
           m->invertices > VERTEXPERBLOCK ? m->invertices : VERTEXPERBLOCK,
           sizeof(REAL));
}

/*****************************************************************************/
/*                                                                           */
/*  initializetrisubpools()   Calculate the sizes of the triangle and        */
/*                            subsegment data structures and initialize      */
/*                            their memory pools.                            */
/*                                                                           */
/*  This routine also computes the `highorderindex', `elemattribindex', and  */
/*  `areaboundindex' indices used to find values within each triangle.       */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void initializetrisubpools(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void initializetrisubpools(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  int trisize;

  /* The index within each triangle at which the extra nodes (above three)  */
  /*   associated with high order elements are found.  There are three      */
  /*   pointers to other triangles, three pointers to corners, and possibly */
  /*   three pointers to subsegments before the extra nodes.                */
  m->highorderindex = 6 + (b->usesegments * 3);
  /* The number of bytes occupied by a triangle. */
  trisize = ((b->order + 1) * (b->order + 2) / 2 + (m->highorderindex - 3)) *
            sizeof(triangle);
  /* The index within each triangle at which its attributes are found, */
  /*   where the index is measured in REALs.                           */
  m->elemattribindex = (trisize + sizeof(REAL) - 1) / sizeof(REAL);
  /* The index within each triangle at which the maximum area constraint  */
  /*   is found, where the index is measured in REALs.  Note that if the  */
  /*   `regionattrib' flag is set, an additional attribute will be added. */
  m->areaboundindex = m->elemattribindex + m->eextras + b->regionattrib;
  /* If triangle attributes or an area bound are needed, increase the number */
  /*   of bytes occupied by a triangle.                                      */
  if (b->vararea) {
    trisize = (m->areaboundindex + 1) * sizeof(REAL);
  } else if (m->eextras + b->regionattrib > 0) {
    trisize = m->areaboundindex * sizeof(REAL);
  }
  /* If a Voronoi diagram or triangle neighbor graph is requested, make    */
  /*   sure there's room to store an integer index in each triangle.  This */
  /*   integer index can occupy the same space as the subsegment pointers  */
  /*   or attributes or area constraint or extra nodes.                    */
  if ((b->voronoi || b->neighbors) &&
      (trisize < 6 * sizeof(triangle) + sizeof(int))) {
    trisize = 6 * sizeof(triangle) + sizeof(int);
  }

  /* Having determined the memory size of a triangle, initialize the pool. */
  poolinit(&m->triangles, trisize, TRIPERBLOCK,
           (2 * m->invertices - 2) > TRIPERBLOCK ? (2 * m->invertices - 2) :
           TRIPERBLOCK, 4);

  if (b->usesegments) {
    /* Initialize the pool of subsegments.  Take into account all eight */
    /*   pointers and one boundary marker.                              */
    poolinit(&m->subsegs, 8 * sizeof(triangle) + sizeof(int),
             SUBSEGPERBLOCK, SUBSEGPERBLOCK, 4);

    /* Initialize the "outer space" triangle and omnipresent subsegment. */
    dummyinit(m, b, m->triangles.itembytes, m->subsegs.itembytes);
  } else {
    /* Initialize the "outer space" triangle. */
    dummyinit(m, b, m->triangles.itembytes, 0);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  triangledealloc()   Deallocate space for a triangle, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void triangledealloc(struct mesh *m, triangle *dyingtriangle)
#else /* not ANSI_DECLARATORS */
void triangledealloc(m, dyingtriangle)
struct mesh *m;
triangle *dyingtriangle;
#endif /* not ANSI_DECLARATORS */

{
  /* Mark the triangle as dead.  This makes it possible to detect dead */
  /*   triangles when traversing the list of all triangles.            */
  killtri(dyingtriangle);
  pooldealloc(&m->triangles, (VOID *) dyingtriangle);
}

/*****************************************************************************/
/*                                                                           */
/*  triangletraverse()   Traverse the triangles, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
triangle *triangletraverse(struct mesh *m)
#else /* not ANSI_DECLARATORS */
triangle *triangletraverse(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  triangle *newtriangle;

  do {
    newtriangle = (triangle *) traverse(&m->triangles);
    if (newtriangle == (triangle *) NULL) {
      return (triangle *) NULL;
    }
  } while (deadtri(newtriangle));                         /* Skip dead ones. */
  return newtriangle;
}

/*****************************************************************************/
/*                                                                           */
/*  subsegdealloc()   Deallocate space for a subsegment, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void subsegdealloc(struct mesh *m, subseg *dyingsubseg)
#else /* not ANSI_DECLARATORS */
void subsegdealloc(m, dyingsubseg)
struct mesh *m;
subseg *dyingsubseg;
#endif /* not ANSI_DECLARATORS */

{
  /* Mark the subsegment as dead.  This makes it possible to detect dead */
  /*   subsegments when traversing the list of all subsegments.          */
  killsubseg(dyingsubseg);
  pooldealloc(&m->subsegs, (VOID *) dyingsubseg);
}

/*****************************************************************************/
/*                                                                           */
/*  subsegtraverse()   Traverse the subsegments, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
subseg *subsegtraverse(struct mesh *m)
#else /* not ANSI_DECLARATORS */
subseg *subsegtraverse(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  subseg *newsubseg;

  do {
    newsubseg = (subseg *) traverse(&m->subsegs);
    if (newsubseg == (subseg *) NULL) {
      return (subseg *) NULL;
    }
  } while (deadsubseg(newsubseg));                        /* Skip dead ones. */
  return newsubseg;
}

/*****************************************************************************/
/*                                                                           */
/*  vertexdealloc()   Deallocate space for a vertex, marking it dead.        */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void vertexdealloc(struct mesh *m, vertex dyingvertex)
#else /* not ANSI_DECLARATORS */
void vertexdealloc(m, dyingvertex)
struct mesh *m;
vertex dyingvertex;
#endif /* not ANSI_DECLARATORS */

{
  /* Mark the vertex as dead.  This makes it possible to detect dead */
  /*   vertices when traversing the list of all vertices.            */
  setvertextype(dyingvertex, DEADVERTEX);
  pooldealloc(&m->vertices, (VOID *) dyingvertex);
}

/*****************************************************************************/
/*                                                                           */
/*  vertextraverse()   Traverse the vertices, skipping dead ones.            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
vertex vertextraverse(struct mesh *m)
#else /* not ANSI_DECLARATORS */
vertex vertextraverse(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  vertex newvertex;

  do {
    newvertex = (vertex) traverse(&m->vertices);
    if (newvertex == (vertex) NULL) {
      return (vertex) NULL;
    }
  } while (vertextype(newvertex) == DEADVERTEX);          /* Skip dead ones. */
  return newvertex;
}

/*****************************************************************************/
/*                                                                           */
/*  badsubsegdealloc()   Deallocate space for a bad subsegment, marking it   */
/*                       dead.                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void badsubsegdealloc(struct mesh *m, struct badsubseg *dyingseg)
#else /* not ANSI_DECLARATORS */
void badsubsegdealloc(m, dyingseg)
struct mesh *m;
struct badsubseg *dyingseg;
#endif /* not ANSI_DECLARATORS */

{
  /* Set subsegment's origin to NULL.  This makes it possible to detect dead */
  /*   badsubsegs when traversing the list of all badsubsegs             .   */
  dyingseg->subsegorg = (vertex) NULL;
  pooldealloc(&m->badsubsegs, (VOID *) dyingseg);
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  badsubsegtraverse()   Traverse the bad subsegments, skipping dead ones.  */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
struct badsubseg *badsubsegtraverse(struct mesh *m)
#else /* not ANSI_DECLARATORS */
struct badsubseg *badsubsegtraverse(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  struct badsubseg *newseg;

  do {
    newseg = (struct badsubseg *) traverse(&m->badsubsegs);
    if (newseg == (struct badsubseg *) NULL) {
      return (struct badsubseg *) NULL;
    }
  } while (newseg->subsegorg == (vertex) NULL);           /* Skip dead ones. */
  return newseg;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  getvertex()   Get a specific vertex, by number, from the list.           */
/*                                                                           */
/*  The first vertex is number 'firstnumber'.                                */
/*                                                                           */
/*  Note that this takes O(n) time (with a small constant, if VERTEXPERBLOCK */
/*  is large).  I don't care to take the trouble to make it work in constant */
/*  time.                                                                    */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
vertex getvertex(struct mesh *m, struct behavior *b, int number)
#else /* not ANSI_DECLARATORS */
vertex getvertex(m, b, number)
struct mesh *m;
struct behavior *b;
int number;
#endif /* not ANSI_DECLARATORS */

{
  VOID **getblock;
  char *foundvertex;
  unsigned long alignptr;
  int current;

  getblock = m->vertices.firstblock;
  current = b->firstnumber;

  /* Find the right block. */
  if (current + m->vertices.itemsfirstblock <= number) {
    getblock = (VOID **) *getblock;
    current += m->vertices.itemsfirstblock;
    while (current + m->vertices.itemsperblock <= number) {
      getblock = (VOID **) *getblock;
      current += m->vertices.itemsperblock;
    }
  }

  /* Now find the right vertex. */
  alignptr = (unsigned long) (getblock + 1);
  foundvertex = (char *) (alignptr + (unsigned long) m->vertices.alignbytes -
                          (alignptr % (unsigned long) m->vertices.alignbytes));
  return (vertex) (foundvertex + m->vertices.itembytes * (number - current));
}

/*****************************************************************************/
/*                                                                           */
/*  triangledeinit()   Free all remaining allocated memory.                  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void triangledeinit(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void triangledeinit(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  pooldeinit(&m->triangles);
  trifree((VOID *) m->dummytribase);
  if (b->usesegments) {
    pooldeinit(&m->subsegs);
    trifree((VOID *) m->dummysubbase);
  }
  pooldeinit(&m->vertices);
#ifndef CDT_ONLY
  if (b->quality) {
    pooldeinit(&m->badsubsegs);
    if ((b->minangle > 0.0) || b->vararea || b->fixedarea || b->usertest) {
      pooldeinit(&m->badtriangles);
      pooldeinit(&m->flipstackers);
    }
  }
#endif /* not CDT_ONLY */
}

/**                                                                         **/
/**                                                                         **/
/********* Memory management routines end here                       *********/

/********* Constructors begin here                                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  maketriangle()   Create a new triangle with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void maketriangle(struct mesh *m, struct behavior *b, struct otri *newotri)
#else /* not ANSI_DECLARATORS */
void maketriangle(m, b, newotri)
struct mesh *m;
struct behavior *b;
struct otri *newotri;
#endif /* not ANSI_DECLARATORS */

{
  int i;

  newotri->tri = (triangle *) poolalloc(&m->triangles);
  /* Initialize the three adjoining triangles to be "outer space". */
  newotri->tri[0] = (triangle) m->dummytri;
  newotri->tri[1] = (triangle) m->dummytri;
  newotri->tri[2] = (triangle) m->dummytri;
  /* Three NULL vertices. */
  newotri->tri[3] = (triangle) NULL;
  newotri->tri[4] = (triangle) NULL;
  newotri->tri[5] = (triangle) NULL;
  if (b->usesegments) {
    /* Initialize the three adjoining subsegments to be the omnipresent */
    /*   subsegment.                                                    */
    newotri->tri[6] = (triangle) m->dummysub;
    newotri->tri[7] = (triangle) m->dummysub;
    newotri->tri[8] = (triangle) m->dummysub;
  }
  for (i = 0; i < m->eextras; i++) {
    setelemattribute(*newotri, i, 0.0);
  }
  if (b->vararea) {
    setareabound(*newotri, -1.0);
  }

  newotri->orient = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  makesubseg()   Create a new subsegment with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void makesubseg(struct mesh *m, struct osub *newsubseg)
#else /* not ANSI_DECLARATORS */
void makesubseg(m, newsubseg)
struct mesh *m;
struct osub *newsubseg;
#endif /* not ANSI_DECLARATORS */

{
  newsubseg->ss = (subseg *) poolalloc(&m->subsegs);
  /* Initialize the two adjoining subsegments to be the omnipresent */
  /*   subsegment.                                                  */
  newsubseg->ss[0] = (subseg) m->dummysub;
  newsubseg->ss[1] = (subseg) m->dummysub;
  /* Four NULL vertices. */
  newsubseg->ss[2] = (subseg) NULL;
  newsubseg->ss[3] = (subseg) NULL;
  newsubseg->ss[4] = (subseg) NULL;
  newsubseg->ss[5] = (subseg) NULL;
  /* Initialize the two adjoining triangles to be "outer space." */
  newsubseg->ss[6] = (subseg) m->dummytri;
  newsubseg->ss[7] = (subseg) m->dummytri;
  /* Set the boundary marker to zero. */
  setmark(*newsubseg, 0);

  newsubseg->ssorient = 0;
}

/**                                                                         **/
/**                                                                         **/
/********* Constructors end here                                     *********/

/********* Geometric primitives begin here                           *********/
/**                                                                         **/
/**                                                                         **/

/* The adaptive exact arithmetic geometric predicates implemented herein are */
/*   described in detail in my paper, "Adaptive Precision Floating-Point     */
/*   Arithmetic and Fast Robust Geometric Predicates."  See the header for a */
/*   full citation.                                                          */

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C without       */
/*   forcing the value to be stored to memory (rather than be kept in the    */
/*   register to which the optimizer assigned it).                           */

#define Absolute(a)  ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a; \
  y = b - bvirt

#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL) (x - a); \
  avirt = x - bvirt; \
  bround = b - bvirt; \
  around = a - avirt; \
  y = around + bround

#define Two_Sum(a, b, x, y) \
  x = (REAL) (a + b); \
  Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL) (a - x); \
  avirt = x + bvirt; \
  bround = bvirt - b; \
  around = a - avirt; \
  y = around + bround

#define Two_Diff(a, b, x, y) \
  x = (REAL) (a - b); \
  Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo) \
  c = (REAL) (splitter * a); \
  abig = (REAL) (c - a); \
  ahi = c - abig; \
  alo = a - ahi

#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo); \
  Split(b, bhi, blo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

#define Two_Product(a, b, x, y) \
  x = (REAL) (a * b); \
  Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL) (a * b); \
  Split(a, ahi, alo); \
  err1 = x - (ahi * bhi); \
  err2 = err1 - (alo * bhi); \
  err3 = err2 - (ahi * blo); \
  y = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y) \
  Split(a, ahi, alo); \
  err1 = x - (ahi * ahi); \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

#define Square(a, x, y) \
  x = (REAL) (a * a); \
  Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b , _i, x0); \
  Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b , _i, x0); \
  Two_Sum( a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0); \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0); \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

/* Macro for multiplying a two-component expansion by a single component.    */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo); \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
  Two_Sum(_i, _0, _k, x1); \
  Fast_Two_Sum(_j, _k, x3, x2)

/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void exactinit()
{
  REAL half;
  REAL check, lastcheck;
  int every_other;
#ifdef LINUX
  int cword;
#endif /* LINUX */

#ifdef CPU86
#ifdef SINGLE
  _control87(_PC_24, _MCW_PC); /* Set FPU control word for single precision. */
#else /* not SINGLE */
  _control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. */
#endif /* not SINGLE */
#endif /* CPU86 */
#ifdef LINUX
#ifdef SINGLE
  /*  cword = 4223; */
  cword = 4210;                 /* set FPU control word for single precision */
#else /* not SINGLE */
  /*  cword = 4735; */
  cword = 4722;                 /* set FPU control word for double precision */
#endif /* not SINGLE */
  _FPU_SETCW(cword);
#endif /* LINUX */

  every_other = 1;
  half = 0.5;
  epsilon = 1.0;
  splitter = 1.0;
  check = 1.0;
  /* Repeatedly divide `epsilon' by two until it is too small to add to      */
  /*   one without causing roundoff.  (Also check if the sum is equal to     */
  /*   the previous sum, for machines that round up instead of using exact   */
  /*   rounding.  Not that these routines will work on such machines.)       */
  do {
    lastcheck = check;
    epsilon *= half;
    if (every_other) {
      splitter *= 2.0;
    }
    every_other = !every_other;
    check = 1.0 + epsilon;
  } while ((check != 1.0) && (check != lastcheck));
  splitter += 1.0;
  /* Error bounds for orientation and incircle tests. */
  resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
  ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
  ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
  ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
  iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
  iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
  iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
  o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
  o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
  o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See my Robust Predicates paper for details.             */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
int fast_expansion_sum_zeroelim(int elen, REAL *e, int flen, REAL *f, REAL *h)
#else /* not ANSI_DECLARATORS */
int fast_expansion_sum_zeroelim(elen, e, flen, f, h)  /* h cannot be e or f. */
int elen;
REAL *e;
int flen;
REAL *f;
REAL *h;
#endif /* not ANSI_DECLARATORS */

{
  REAL Q;
  INEXACT REAL Qnew;
  INEXACT REAL hh;
  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  int eindex, findex, hindex;
  REAL enow, fnow;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      Fast_Two_Sum(enow, Q, Qnew, hh);
      enow = e[++eindex];
    } else {
      Fast_Two_Sum(fnow, Q, Qnew, hh);
      fnow = f[++findex];
    }
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        Two_Sum(Q, enow, Qnew, hh);
        enow = e[++eindex];
      } else {
        Two_Sum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
      }
      Q = Qnew;
      if (hh != 0.0) {
        h[hindex++] = hh;
      }
    }
  }
  while (eindex < elen) {
    Two_Sum(Q, enow, Qnew, hh);
    enow = e[++eindex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  while (findex < flen) {
    Two_Sum(Q, fnow, Qnew, hh);
    fnow = f[++findex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See my Robust Predicates paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
int scale_expansion_zeroelim(int elen, REAL *e, REAL b, REAL *h)
#else /* not ANSI_DECLARATORS */
int scale_expansion_zeroelim(elen, e, b, h)   /* e and h cannot be the same. */
int elen;
REAL *e;
REAL b;
REAL *h;
#endif /* not ANSI_DECLARATORS */

{
  INEXACT REAL Q, sum;
  REAL hh;
  INEXACT REAL product1;
  REAL product0;
  int eindex, hindex;
  REAL enow;
  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;

  Split(b, bhi, blo);
  Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
  hindex = 0;
  if (hh != 0) {
    h[hindex++] = hh;
  }
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
    Two_Sum(Q, product0, sum, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
    Fast_Two_Sum(product1, sum, Q, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL estimate(int elen, REAL *e)
#else /* not ANSI_DECLARATORS */
REAL estimate(elen, e)
int elen;
REAL *e;
#endif /* not ANSI_DECLARATORS */

{
  REAL Q;
  int eindex;

  Q = e[0];
  for (eindex = 1; eindex < elen; eindex++) {
    Q += e[eindex];
  }
  return Q;
}

/*****************************************************************************/
/*                                                                           */
/*  counterclockwise()   Return a positive value if the points pa, pb, and   */
/*                       pc occur in counterclockwise order; a negative      */
/*                       value if they occur in clockwise order; and zero    */
/*                       if they are collinear.  The result is also a rough  */
/*                       approximation of twice the signed area of the       */
/*                       triangle defined by the three points.               */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are collinear or nearly so.            */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL counterclockwiseadapt(vertex pa, vertex pb, vertex pc, REAL detsum)
#else /* not ANSI_DECLARATORS */
REAL counterclockwiseadapt(pa, pb, pc, detsum)
vertex pa;
vertex pb;
vertex pc;
REAL detsum;
#endif /* not ANSI_DECLARATORS */

{
  INEXACT REAL acx, acy, bcx, bcy;
  REAL acxtail, acytail, bcxtail, bcytail;
  INEXACT REAL detleft, detright;
  REAL detlefttail, detrighttail;
  REAL det, errbound;
  REAL B[4], C1[8], C2[12], D[16];
  INEXACT REAL B3;
  int C1length, C2length, Dlength;
  REAL u[4];
  INEXACT REAL u3;
  INEXACT REAL s1, t1;
  REAL s0, t0;

  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;
  INEXACT REAL _i, _j;
  REAL _0;

  acx = (REAL) (pa[0] - pc[0]);
  bcx = (REAL) (pb[0] - pc[0]);
  acy = (REAL) (pa[1] - pc[1]);
  bcy = (REAL) (pb[1] - pc[1]);

  Two_Product(acx, bcy, detleft, detlefttail);
  Two_Product(acy, bcx, detright, detrighttail);

  Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
               B3, B[2], B[1], B[0]);
  B[3] = B3;

  det = estimate(4, B);
  errbound = ccwerrboundB * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
  Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
  Two_Diff_Tail(pa[1], pc[1], acy, acytail);
  Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det;
  }

  errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
  det += (acx * bcytail + bcy * acxtail)
       - (acy * bcxtail + bcx * acytail);
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Product(acxtail, bcy, s1, s0);
  Two_Product(acytail, bcx, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

  Two_Product(acx, bcytail, s1, s0);
  Two_Product(acy, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

  Two_Product(acxtail, bcytail, s1, s0);
  Two_Product(acytail, bcxtail, t1, t0);
  Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

  return(D[Dlength - 1]);
}

#ifdef ANSI_DECLARATORS
REAL counterclockwise(struct mesh *m, struct behavior *b,
                      vertex pa, vertex pb, vertex pc)
#else /* not ANSI_DECLARATORS */
REAL counterclockwise(m, b, pa, pb, pc)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
#endif /* not ANSI_DECLARATORS */

{
  REAL detleft, detright, det;
  REAL detsum, errbound;

  m->counterclockcount++;

  detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
  detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
  det = detleft - detright;

  if (b->noexact) {
    return det;
  }

  if (detleft > 0.0) {
    if (detright <= 0.0) {
      return det;
    } else {
      detsum = detleft + detright;
    }
  } else if (detleft < 0.0) {
    if (detright >= 0.0) {
      return det;
    } else {
      detsum = -detleft - detright;
    }
  } else {
    return det;
  }

  errbound = ccwerrboundA * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  return counterclockwiseadapt(pa, pb, pc, detsum);
}

/*****************************************************************************/
/*                                                                           */
/*  incircle()   Return a positive value if the point pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are cocircular or nearly so.           */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL incircleadapt(vertex pa, vertex pb, vertex pc, vertex pd, REAL permanent)
#else /* not ANSI_DECLARATORS */
REAL incircleadapt(pa, pb, pc, pd, permanent)
vertex pa;
vertex pb;
vertex pc;
vertex pd;
REAL permanent;
#endif /* not ANSI_DECLARATORS */

{
  INEXACT REAL adx, bdx, cdx, ady, bdy, cdy;
  REAL det, errbound;

  INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  REAL bc[4], ca[4], ab[4];
  INEXACT REAL bc3, ca3, ab3;
  REAL axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
  int axbclen, axxbclen, aybclen, ayybclen, alen;
  REAL bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
  int bxcalen, bxxcalen, bycalen, byycalen, blen;
  REAL cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
  int cxablen, cxxablen, cyablen, cyyablen, clen;
  REAL abdet[64];
  int ablen;
  REAL fin1[1152], fin2[1152];
  REAL *finnow, *finother, *finswap;
  int finlength;

  REAL adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
  INEXACT REAL adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
  REAL adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
  REAL aa[4], bb[4], cc[4];
  INEXACT REAL aa3, bb3, cc3;
  INEXACT REAL ti1, tj1;
  REAL ti0, tj0;
  REAL u[4], v[4];
  INEXACT REAL u3, v3;
  REAL temp8[8], temp16a[16], temp16b[16], temp16c[16];
  REAL temp32a[32], temp32b[32], temp48[48], temp64[64];
  int temp8len, temp16alen, temp16blen, temp16clen;
  int temp32alen, temp32blen, temp48len, temp64len;
  REAL axtbb[8], axtcc[8], aytbb[8], aytcc[8];
  int axtbblen, axtcclen, aytbblen, aytcclen;
  REAL bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
  int bxtaalen, bxtcclen, bytaalen, bytcclen;
  REAL cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
  int cxtaalen, cxtbblen, cytaalen, cytbblen;
  REAL axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
  int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
  REAL axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16], cytabt[16];
  int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
  REAL axtbctt[8], aytbctt[8], bxtcatt[8];
  REAL bytcatt[8], cxtabtt[8], cytabtt[8];
  int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
  REAL abt[8], bct[8], cat[8];
  int abtlen, bctlen, catlen;
  REAL abtt[4], bctt[4], catt[4];
  int abttlen, bcttlen, cattlen;
  INEXACT REAL abtt3, bctt3, catt3;
  REAL negate;

  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;
  INEXACT REAL _i, _j;
  REAL _0;

  adx = (REAL) (pa[0] - pd[0]);
  bdx = (REAL) (pb[0] - pd[0]);
  cdx = (REAL) (pc[0] - pd[0]);
  ady = (REAL) (pa[1] - pd[1]);
  bdy = (REAL) (pb[1] - pd[1]);
  cdy = (REAL) (pc[1] - pd[1]);

  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  axbclen = scale_expansion_zeroelim(4, bc, adx, axbc);
  axxbclen = scale_expansion_zeroelim(axbclen, axbc, adx, axxbc);
  aybclen = scale_expansion_zeroelim(4, bc, ady, aybc);
  ayybclen = scale_expansion_zeroelim(aybclen, aybc, ady, ayybc);
  alen = fast_expansion_sum_zeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

  Two_Product(cdx, ady, cdxady1, cdxady0);
  Two_Product(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  bxcalen = scale_expansion_zeroelim(4, ca, bdx, bxca);
  bxxcalen = scale_expansion_zeroelim(bxcalen, bxca, bdx, bxxca);
  bycalen = scale_expansion_zeroelim(4, ca, bdy, byca);
  byycalen = scale_expansion_zeroelim(bycalen, byca, bdy, byyca);
  blen = fast_expansion_sum_zeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);

  Two_Product(adx, bdy, adxbdy1, adxbdy0);
  Two_Product(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  cxablen = scale_expansion_zeroelim(4, ab, cdx, cxab);
  cxxablen = scale_expansion_zeroelim(cxablen, cxab, cdx, cxxab);
  cyablen = scale_expansion_zeroelim(4, ab, cdy, cyab);
  cyyablen = scale_expansion_zeroelim(cyablen, cyab, cdy, cyyab);
  clen = fast_expansion_sum_zeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = iccerrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
  Two_Diff_Tail(pa[1], pd[1], ady, adytail);
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)) {
    return det;
  }

  errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
  det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
                                     - (bdy * cdxtail + cdx * bdytail))
          + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
       + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
                                     - (cdy * adxtail + adx * cdytail))
          + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
       + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
                                     - (ady * bdxtail + bdx * adytail))
          + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if ((bdxtail != 0.0) || (bdytail != 0.0)
      || (cdxtail != 0.0) || (cdytail != 0.0)) {
    Square(adx, adxadx1, adxadx0);
    Square(ady, adyady1, adyady0);
    Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
    aa[3] = aa3;
  }
  if ((cdxtail != 0.0) || (cdytail != 0.0)
      || (adxtail != 0.0) || (adytail != 0.0)) {
    Square(bdx, bdxbdx1, bdxbdx0);
    Square(bdy, bdybdy1, bdybdy0);
    Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
    bb[3] = bb3;
  }
  if ((adxtail != 0.0) || (adytail != 0.0)
      || (bdxtail != 0.0) || (bdytail != 0.0)) {
    Square(cdx, cdxcdx1, cdxcdx0);
    Square(cdy, cdycdy1, cdycdy0);
    Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
    cc[3] = cc3;
  }

  if (adxtail != 0.0) {
    axtbclen = scale_expansion_zeroelim(4, bc, adxtail, axtbc);
    temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, 2.0 * adx,
                                          temp16a);

    axtcclen = scale_expansion_zeroelim(4, cc, adxtail, axtcc);
    temp16blen = scale_expansion_zeroelim(axtcclen, axtcc, bdy, temp16b);

    axtbblen = scale_expansion_zeroelim(4, bb, adxtail, axtbb);
    temp16clen = scale_expansion_zeroelim(axtbblen, axtbb, -cdy, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (adytail != 0.0) {
    aytbclen = scale_expansion_zeroelim(4, bc, adytail, aytbc);
    temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, 2.0 * ady,
                                          temp16a);

    aytbblen = scale_expansion_zeroelim(4, bb, adytail, aytbb);
    temp16blen = scale_expansion_zeroelim(aytbblen, aytbb, cdx, temp16b);

    aytcclen = scale_expansion_zeroelim(4, cc, adytail, aytcc);
    temp16clen = scale_expansion_zeroelim(aytcclen, aytcc, -bdx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdxtail != 0.0) {
    bxtcalen = scale_expansion_zeroelim(4, ca, bdxtail, bxtca);
    temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, 2.0 * bdx,
                                          temp16a);

    bxtaalen = scale_expansion_zeroelim(4, aa, bdxtail, bxtaa);
    temp16blen = scale_expansion_zeroelim(bxtaalen, bxtaa, cdy, temp16b);

    bxtcclen = scale_expansion_zeroelim(4, cc, bdxtail, bxtcc);
    temp16clen = scale_expansion_zeroelim(bxtcclen, bxtcc, -ady, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdytail != 0.0) {
    bytcalen = scale_expansion_zeroelim(4, ca, bdytail, bytca);
    temp16alen = scale_expansion_zeroelim(bytcalen, bytca, 2.0 * bdy,
                                          temp16a);

    bytcclen = scale_expansion_zeroelim(4, cc, bdytail, bytcc);
    temp16blen = scale_expansion_zeroelim(bytcclen, bytcc, adx, temp16b);

    bytaalen = scale_expansion_zeroelim(4, aa, bdytail, bytaa);
    temp16clen = scale_expansion_zeroelim(bytaalen, bytaa, -cdx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdxtail != 0.0) {
    cxtablen = scale_expansion_zeroelim(4, ab, cdxtail, cxtab);
    temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, 2.0 * cdx,
                                          temp16a);

    cxtbblen = scale_expansion_zeroelim(4, bb, cdxtail, cxtbb);
    temp16blen = scale_expansion_zeroelim(cxtbblen, cxtbb, ady, temp16b);

    cxtaalen = scale_expansion_zeroelim(4, aa, cdxtail, cxtaa);
    temp16clen = scale_expansion_zeroelim(cxtaalen, cxtaa, -bdy, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdytail != 0.0) {
    cytablen = scale_expansion_zeroelim(4, ab, cdytail, cytab);
    temp16alen = scale_expansion_zeroelim(cytablen, cytab, 2.0 * cdy,
                                          temp16a);

    cytaalen = scale_expansion_zeroelim(4, aa, cdytail, cytaa);
    temp16blen = scale_expansion_zeroelim(cytaalen, cytaa, bdx, temp16b);

    cytbblen = scale_expansion_zeroelim(4, bb, cdytail, cytbb);
    temp16clen = scale_expansion_zeroelim(cytbblen, cytbb, -adx, temp16c);

    temp32alen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                            temp16blen, temp16b, temp32a);
    temp48len = fast_expansion_sum_zeroelim(temp16clen, temp16c,
                                            temp32alen, temp32a, temp48);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48, finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if ((adxtail != 0.0) || (adytail != 0.0)) {
    if ((bdxtail != 0.0) || (bdytail != 0.0)
        || (cdxtail != 0.0) || (cdytail != 0.0)) {
      Two_Product(bdxtail, cdy, ti1, ti0);
      Two_Product(bdx, cdytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -bdy;
      Two_Product(cdxtail, negate, ti1, ti0);
      negate = -bdytail;
      Two_Product(cdx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      bctlen = fast_expansion_sum_zeroelim(4, u, 4, v, bct);

      Two_Product(bdxtail, cdytail, ti1, ti0);
      Two_Product(cdxtail, bdytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
      bctt[3] = bctt3;
      bcttlen = 4;
    } else {
      bct[0] = 0.0;
      bctlen = 1;
      bctt[0] = 0.0;
      bcttlen = 1;
    }

    if (adxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(axtbclen, axtbc, adxtail, temp16a);
      axtbctlen = scale_expansion_zeroelim(bctlen, bct, adxtail, axtbct);
      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, 2.0 * adx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, cc, adxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (cdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, bb, -adxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(axtbctlen, axtbct, adxtail,
                                            temp32a);
      axtbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adxtail, axtbctt);
      temp16alen = scale_expansion_zeroelim(axtbcttlen, axtbctt, 2.0 * adx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(axtbcttlen, axtbctt, adxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (adytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(aytbclen, aytbc, adytail, temp16a);
      aytbctlen = scale_expansion_zeroelim(bctlen, bct, adytail, aytbct);
      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, 2.0 * ady,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(aytbctlen, aytbct, adytail,
                                            temp32a);
      aytbcttlen = scale_expansion_zeroelim(bcttlen, bctt, adytail, aytbctt);
      temp16alen = scale_expansion_zeroelim(aytbcttlen, aytbctt, 2.0 * ady,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(aytbcttlen, aytbctt, adytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }
  if ((bdxtail != 0.0) || (bdytail != 0.0)) {
    if ((cdxtail != 0.0) || (cdytail != 0.0)
        || (adxtail != 0.0) || (adytail != 0.0)) {
      Two_Product(cdxtail, ady, ti1, ti0);
      Two_Product(cdx, adytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -cdy;
      Two_Product(adxtail, negate, ti1, ti0);
      negate = -cdytail;
      Two_Product(adx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      catlen = fast_expansion_sum_zeroelim(4, u, 4, v, cat);

      Two_Product(cdxtail, adytail, ti1, ti0);
      Two_Product(adxtail, cdytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
      catt[3] = catt3;
      cattlen = 4;
    } else {
      cat[0] = 0.0;
      catlen = 1;
      catt[0] = 0.0;
      cattlen = 1;
    }

    if (bdxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(bxtcalen, bxtca, bdxtail, temp16a);
      bxtcatlen = scale_expansion_zeroelim(catlen, cat, bdxtail, bxtcat);
      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, aa, bdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, cdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (adytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, cc, -bdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(bxtcatlen, bxtcat, bdxtail,
                                            temp32a);
      bxtcattlen = scale_expansion_zeroelim(cattlen, catt, bdxtail, bxtcatt);
      temp16alen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(bxtcattlen, bxtcatt, bdxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (bdytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(bytcalen, bytca, bdytail, temp16a);
      bytcatlen = scale_expansion_zeroelim(catlen, cat, bdytail, bytcat);
      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, 2.0 * bdy,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(bytcatlen, bytcat, bdytail,
                                            temp32a);
      bytcattlen = scale_expansion_zeroelim(cattlen, catt, bdytail, bytcatt);
      temp16alen = scale_expansion_zeroelim(bytcattlen, bytcatt, 2.0 * bdy,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(bytcattlen, bytcatt, bdytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }
  if ((cdxtail != 0.0) || (cdytail != 0.0)) {
    if ((adxtail != 0.0) || (adytail != 0.0)
        || (bdxtail != 0.0) || (bdytail != 0.0)) {
      Two_Product(adxtail, bdy, ti1, ti0);
      Two_Product(adx, bdytail, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
      u[3] = u3;
      negate = -ady;
      Two_Product(bdxtail, negate, ti1, ti0);
      negate = -adytail;
      Two_Product(bdx, negate, tj1, tj0);
      Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
      v[3] = v3;
      abtlen = fast_expansion_sum_zeroelim(4, u, 4, v, abt);

      Two_Product(adxtail, bdytail, ti1, ti0);
      Two_Product(bdxtail, adytail, tj1, tj0);
      Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
      abtt[3] = abtt3;
      abttlen = 4;
    } else {
      abt[0] = 0.0;
      abtlen = 1;
      abtt[0] = 0.0;
      abttlen = 1;
    }

    if (cdxtail != 0.0) {
      temp16alen = scale_expansion_zeroelim(cxtablen, cxtab, cdxtail, temp16a);
      cxtabtlen = scale_expansion_zeroelim(abtlen, abt, cdxtail, cxtabt);
      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, bb, cdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, adytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
      if (bdytail != 0.0) {
        temp8len = scale_expansion_zeroelim(4, aa, -cdxtail, temp8);
        temp16alen = scale_expansion_zeroelim(temp8len, temp8, bdytail,
                                              temp16a);
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a, finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }

      temp32alen = scale_expansion_zeroelim(cxtabtlen, cxtabt, cdxtail,
                                            temp32a);
      cxtabttlen = scale_expansion_zeroelim(abttlen, abtt, cdxtail, cxtabtt);
      temp16alen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(cxtabttlen, cxtabtt, cdxtail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
    if (cdytail != 0.0) {
      temp16alen = scale_expansion_zeroelim(cytablen, cytab, cdytail, temp16a);
      cytabtlen = scale_expansion_zeroelim(abtlen, abt, cdytail, cytabt);
      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, 2.0 * cdy,
                                            temp32a);
      temp48len = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp32alen, temp32a, temp48);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                              temp48, finother);
      finswap = finnow; finnow = finother; finother = finswap;


      temp32alen = scale_expansion_zeroelim(cytabtlen, cytabt, cdytail,
                                            temp32a);
      cytabttlen = scale_expansion_zeroelim(abttlen, abtt, cdytail, cytabtt);
      temp16alen = scale_expansion_zeroelim(cytabttlen, cytabtt, 2.0 * cdy,
                                            temp16a);
      temp16blen = scale_expansion_zeroelim(cytabttlen, cytabtt, cdytail,
                                            temp16b);
      temp32blen = fast_expansion_sum_zeroelim(temp16alen, temp16a,
                                              temp16blen, temp16b, temp32b);
      temp64len = fast_expansion_sum_zeroelim(temp32alen, temp32a,
                                              temp32blen, temp32b, temp64);
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                              temp64, finother);
      finswap = finnow; finnow = finother; finother = finswap;
    }
  }

  return finnow[finlength - 1];
}

#ifdef ANSI_DECLARATORS
REAL incircle(struct mesh *m, struct behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd)
#else /* not ANSI_DECLARATORS */
REAL incircle(m, b, pa, pb, pc, pd)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
vertex pd;
#endif /* not ANSI_DECLARATORS */

{
  REAL adx, bdx, cdx, ady, bdy, cdy;
  REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  REAL alift, blift, clift;
  REAL det;
  REAL permanent, errbound;

  m->incirclecount++;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;
  alift = adx * adx + ady * ady;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;
  blift = bdx * bdx + bdy * bdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;
  clift = cdx * cdx + cdy * cdy;

  det = alift * (bdxcdy - cdxbdy)
      + blift * (cdxady - adxcdy)
      + clift * (adxbdy - bdxady);

  if (b->noexact) {
    return det;
  }

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift
            + (Absolute(cdxady) + Absolute(adxcdy)) * blift
            + (Absolute(adxbdy) + Absolute(bdxady)) * clift;
  errbound = iccerrboundA * permanent;
  if ((det > errbound) || (-det > errbound)) {
    return det;
  }

  return incircleadapt(pa, pb, pc, pd, permanent);
}

/*****************************************************************************/
/*                                                                           */
/*  orient3d()   Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are coplanar or nearly so.             */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL orient3dadapt(vertex pa, vertex pb, vertex pc, vertex pd,
                   REAL aheight, REAL bheight, REAL cheight, REAL dheight,
                   REAL permanent)
#else /* not ANSI_DECLARATORS */
REAL orient3dadapt(pa, pb, pc, pd,
                   aheight, bheight, cheight, dheight, permanent)
vertex pa;
vertex pb;
vertex pc;
vertex pd;
REAL aheight;
REAL bheight;
REAL cheight;
REAL dheight;
REAL permanent;
#endif /* not ANSI_DECLARATORS */

{
  INEXACT REAL adx, bdx, cdx, ady, bdy, cdy, adheight, bdheight, cdheight;
  REAL det, errbound;

  INEXACT REAL bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  REAL bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  REAL bc[4], ca[4], ab[4];
  INEXACT REAL bc3, ca3, ab3;
  REAL adet[8], bdet[8], cdet[8];
  int alen, blen, clen;
  REAL abdet[16];
  int ablen;
  REAL *finnow, *finother, *finswap;
  REAL fin1[192], fin2[192];
  int finlength;

  REAL adxtail, bdxtail, cdxtail;
  REAL adytail, bdytail, cdytail;
  REAL adheighttail, bdheighttail, cdheighttail;
  INEXACT REAL at_blarge, at_clarge;
  INEXACT REAL bt_clarge, bt_alarge;
  INEXACT REAL ct_alarge, ct_blarge;
  REAL at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
  int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
  INEXACT REAL bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
  INEXACT REAL adxt_cdy1, adxt_bdy1, bdxt_ady1;
  REAL bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
  REAL adxt_cdy0, adxt_bdy0, bdxt_ady0;
  INEXACT REAL bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
  INEXACT REAL adyt_cdx1, adyt_bdx1, bdyt_adx1;
  REAL bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
  REAL adyt_cdx0, adyt_bdx0, bdyt_adx0;
  REAL bct[8], cat[8], abt[8];
  int bctlen, catlen, abtlen;
  INEXACT REAL bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
  INEXACT REAL adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
  REAL bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
  REAL adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
  REAL u[4], v[12], w[16];
  INEXACT REAL u3;
  int vlength, wlength;
  REAL negate;

  INEXACT REAL bvirt;
  REAL avirt, bround, around;
  INEXACT REAL c;
  INEXACT REAL abig;
  REAL ahi, alo, bhi, blo;
  REAL err1, err2, err3;
  INEXACT REAL _i, _j, _k;
  REAL _0;

  adx = (REAL) (pa[0] - pd[0]);
  bdx = (REAL) (pb[0] - pd[0]);
  cdx = (REAL) (pc[0] - pd[0]);
  ady = (REAL) (pa[1] - pd[1]);
  bdy = (REAL) (pb[1] - pd[1]);
  cdy = (REAL) (pc[1] - pd[1]);
  adheight = (REAL) (aheight - dheight);
  bdheight = (REAL) (bheight - dheight);
  cdheight = (REAL) (cheight - dheight);

  Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
  Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  alen = scale_expansion_zeroelim(4, bc, adheight, adet);

  Two_Product(cdx, ady, cdxady1, cdxady0);
  Two_Product(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  blen = scale_expansion_zeroelim(4, ca, bdheight, bdet);

  Two_Product(adx, bdy, adxbdy1, adxbdy0);
  Two_Product(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  clen = scale_expansion_zeroelim(4, ab, cdheight, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = o3derrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
  Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
  Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
  Two_Diff_Tail(pa[1], pd[1], ady, adytail);
  Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
  Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
  Two_Diff_Tail(aheight, dheight, adheight, adheighttail);
  Two_Diff_Tail(bheight, dheight, bdheight, bdheighttail);
  Two_Diff_Tail(cheight, dheight, cdheight, cdheighttail);

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) &&
      (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0) &&
      (adheighttail == 0.0) &&
      (bdheighttail == 0.0) &&
      (cdheighttail == 0.0)) {
    return det;
  }

  errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
  det += (adheight * ((bdx * cdytail + cdy * bdxtail) -
                      (bdy * cdxtail + cdx * bdytail)) +
          adheighttail * (bdx * cdy - bdy * cdx)) +
         (bdheight * ((cdx * adytail + ady * cdxtail) -
                      (cdy * adxtail + adx * cdytail)) +
          bdheighttail * (cdx * ady - cdy * adx)) +
         (cdheight * ((adx * bdytail + bdy * adxtail) -
                      (ady * bdxtail + bdx * adytail)) +
          cdheighttail * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b[0] = 0.0;
      at_blen = 1;
      at_c[0] = 0.0;
      at_clen = 1;
    } else {
      negate = -adytail;
      Two_Product(negate, bdx, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      Two_Product(adytail, cdx, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    }
  } else {
    if (adytail == 0.0) {
      Two_Product(adxtail, bdy, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      negate = -adxtail;
      Two_Product(negate, cdy, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    } else {
      Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
      Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
      Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   at_blarge, at_b[2], at_b[1], at_b[0]);
      at_b[3] = at_blarge;
      at_blen = 4;
      Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
      Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
      Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   at_clarge, at_c[2], at_c[1], at_c[0]);
      at_c[3] = at_clarge;
      at_clen = 4;
    }
  }
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c[0] = 0.0;
      bt_clen = 1;
      bt_a[0] = 0.0;
      bt_alen = 1;
    } else {
      negate = -bdytail;
      Two_Product(negate, cdx, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    }
  } else {
    if (bdytail == 0.0) {
      Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      negate = -bdxtail;
      Two_Product(negate, ady, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    } else {
      Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
      Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
      Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                   bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
      bt_c[3] = bt_clarge;
      bt_clen = 4;
      Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
      Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
      Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                  bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
      bt_a[3] = bt_alarge;
      bt_alen = 4;
    }
  }
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a[0] = 0.0;
      ct_alen = 1;
      ct_b[0] = 0.0;
      ct_blen = 1;
    } else {
      negate = -cdytail;
      Two_Product(negate, adx, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    }
  } else {
    if (cdytail == 0.0) {
      Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      negate = -cdxtail;
      Two_Product(negate, bdy, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    } else {
      Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
      Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
      Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
      ct_a[3] = ct_alarge;
      ct_alen = 4;
      Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
      Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
      Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
      ct_b[3] = ct_blarge;
      ct_blen = 4;
    }
  }

  bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
  wlength = scale_expansion_zeroelim(bctlen, bct, adheight, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
  wlength = scale_expansion_zeroelim(catlen, cat, bdheight, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
  wlength = scale_expansion_zeroelim(abtlen, abt, cdheight, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  if (adheighttail != 0.0) {
    vlength = scale_expansion_zeroelim(4, bc, adheighttail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdheighttail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ca, bdheighttail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdheighttail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ab, cdheighttail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
      Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdheighttail != 0.0) {
        Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (cdytail != 0.0) {
      negate = -adxtail;
      Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
      Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdheighttail != 0.0) {
        Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
      Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adheighttail != 0.0) {
        Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (adytail != 0.0) {
      negate = -bdxtail;
      Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
      Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdheighttail != 0.0) {
        Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
      Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdheighttail != 0.0) {
        Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (bdytail != 0.0) {
      negate = -cdxtail;
      Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
      Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adheight, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adheighttail != 0.0) {
        Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adheighttail,
                        u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }

  if (adheighttail != 0.0) {
    wlength = scale_expansion_zeroelim(bctlen, bct, adheighttail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdheighttail != 0.0) {
    wlength = scale_expansion_zeroelim(catlen, cat, bdheighttail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdheighttail != 0.0) {
    wlength = scale_expansion_zeroelim(abtlen, abt, cdheighttail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  return finnow[finlength - 1];
}

#ifdef ANSI_DECLARATORS
REAL orient3d(struct mesh *m, struct behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd,
              REAL aheight, REAL bheight, REAL cheight, REAL dheight)
#else /* not ANSI_DECLARATORS */
REAL orient3d(m, b, pa, pb, pc, pd, aheight, bheight, cheight, dheight)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
vertex pd;
REAL aheight;
REAL bheight;
REAL cheight;
REAL dheight;
#endif /* not ANSI_DECLARATORS */

{
  REAL adx, bdx, cdx, ady, bdy, cdy, adheight, bdheight, cdheight;
  REAL bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
  REAL det;
  REAL permanent, errbound;

  m->orient3dcount++;

  adx = pa[0] - pd[0];
  bdx = pb[0] - pd[0];
  cdx = pc[0] - pd[0];
  ady = pa[1] - pd[1];
  bdy = pb[1] - pd[1];
  cdy = pc[1] - pd[1];
  adheight = aheight - dheight;
  bdheight = bheight - dheight;
  cdheight = cheight - dheight;

  bdxcdy = bdx * cdy;
  cdxbdy = cdx * bdy;

  cdxady = cdx * ady;
  adxcdy = adx * cdy;

  adxbdy = adx * bdy;
  bdxady = bdx * ady;

  det = adheight * (bdxcdy - cdxbdy) 
      + bdheight * (cdxady - adxcdy)
      + cdheight * (adxbdy - bdxady);

  if (b->noexact) {
    return det;
  }

  permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adheight)
            + (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdheight)
            + (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdheight);
  errbound = o3derrboundA * permanent;
  if ((det > errbound) || (-det > errbound)) {
    return det;
  }

  return orient3dadapt(pa, pb, pc, pd, aheight, bheight, cheight, dheight,
                       permanent);
}

/*****************************************************************************/
/*                                                                           */
/*  nonregular()   Return a positive value if the point pd is incompatible   */
/*                 with the circle or plane passing through pa, pb, and pc   */
/*                 (meaning that pd is inside the circle or below the        */
/*                 plane); a negative value if it is compatible; and zero if */
/*                 the four points are cocircular/coplanar.  The points pa,  */
/*                 pb, and pc must be in counterclockwise order, or the sign */
/*                 of the result will be reversed.                           */
/*                                                                           */
/*  If the -w switch is used, the points are lifted onto the parabolic       */
/*  lifting map, then they are dropped according to their weights, then the  */
/*  3D orientation test is applied.  If the -W switch is used, the points'   */
/*  heights are already provided, so the 3D orientation test is applied      */
/*  directly.  If neither switch is used, the incircle test is applied.      */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
REAL nonregular(struct mesh *m, struct behavior *b,
                vertex pa, vertex pb, vertex pc, vertex pd)
#else /* not ANSI_DECLARATORS */
REAL nonregular(m, b, pa, pb, pc, pd)
struct mesh *m;
struct behavior *b;
vertex pa;
vertex pb;
vertex pc;
vertex pd;
#endif /* not ANSI_DECLARATORS */

{
  if (b->weighted == 0) {
    return incircle(m, b, pa, pb, pc, pd);
  } else if (b->weighted == 1) {
    return orient3d(m, b, pa, pb, pc, pd,
                    pa[0] * pa[0] + pa[1] * pa[1] - pa[2],
                    pb[0] * pb[0] + pb[1] * pb[1] - pb[2],
                    pc[0] * pc[0] + pc[1] * pc[1] - pc[2],
                    pd[0] * pd[0] + pd[1] * pd[1] - pd[2]);
  } else {
    return orient3d(m, b, pa, pb, pc, pd, pa[2], pb[2], pc[2], pd[2]);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  findcircumcenter()   Find the circumcenter of a triangle.                */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  (barycentric) coordinates.  The xi-eta coordinate system is defined in   */
/*  terms of the triangle:  the origin of the triangle is the origin of the  */
/*  coordinate system; the destination of the triangle is one unit along the */
/*  xi axis; and the apex of the triangle is one unit along the eta axis.    */
/*  This procedure also returns the square of the length of the triangle's   */
/*  shortest edge.                                                           */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void findcircumcenter(struct mesh *m, struct behavior *b,
                      vertex torg, vertex tdest, vertex tapex,
                      vertex circumcenter, REAL *xi, REAL *eta, int offcenter)
#else /* not ANSI_DECLARATORS */
void findcircumcenter(m, b, torg, tdest, tapex, circumcenter, xi, eta,
                      offcenter)
struct mesh *m;
struct behavior *b;
vertex torg;
vertex tdest;
vertex tapex;
vertex circumcenter;
REAL *xi;
REAL *eta;
int offcenter;
#endif /* not ANSI_DECLARATORS */

{
  REAL xdo, ydo, xao, yao;
  REAL dodist, aodist, dadist;
  REAL denominator;
  REAL dx, dy, dxoff, dyoff;

  m->circumcentercount++;

  /* Compute the circumcenter of the triangle. */
  xdo = tdest[0] - torg[0];
  ydo = tdest[1] - torg[1];
  xao = tapex[0] - torg[0];
  yao = tapex[1] - torg[1];
  dodist = xdo * xdo + ydo * ydo;
  aodist = xao * xao + yao * yao;
  dadist = (tdest[0] - tapex[0]) * (tdest[0] - tapex[0]) +
           (tdest[1] - tapex[1]) * (tdest[1] - tapex[1]);
  if (b->noexact) {
    denominator = 0.5 / (xdo * yao - xao * ydo);
  } else {
    /* Use the counterclockwise() routine to ensure a positive (and */
    /*   reasonably accurate) result, avoiding any possibility of   */
    /*   division by zero.                                          */
    denominator = 0.5 / counterclockwise(m, b, tdest, tapex, torg);
    /* Don't count the above as an orientation test. */
    m->counterclockcount--;
  }
  dx = (yao * dodist - ydo * aodist) * denominator;
  dy = (xdo * aodist - xao * dodist) * denominator;

  /* Find the (squared) length of the triangle's shortest edge.  This   */
  /*   serves as a conservative estimate of the insertion radius of the */
  /*   circumcenter's parent.  The estimate is used to ensure that      */
  /*   the algorithm terminates even if very small angles appear in     */
  /*   the input PSLG.                                                  */
  if ((dodist < aodist) && (dodist < dadist)) {
    if (offcenter && (b->offconstant > 0.0)) {
      /* Find the position of the off-center, as described by Alper Ungor. */
      dxoff = 0.5 * xdo - b->offconstant * ydo;
      dyoff = 0.5 * ydo + b->offconstant * xdo;
      /* If the off-center is closer to the origin than the */
      /*   circumcenter, use the off-center instead.        */
      if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) {
        dx = dxoff;
        dy = dyoff;
      }
    }
  } else if (aodist < dadist) {
    if (offcenter && (b->offconstant > 0.0)) {
      dxoff = 0.5 * xao + b->offconstant * yao;
      dyoff = 0.5 * yao - b->offconstant * xao;
      /* If the off-center is closer to the origin than the */
      /*   circumcenter, use the off-center instead.        */
      if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) {
        dx = dxoff;
        dy = dyoff;
      }
    }
  } else {
    if (offcenter && (b->offconstant > 0.0)) {
      dxoff = 0.5 * (tapex[0] - tdest[0]) -
              b->offconstant * (tapex[1] - tdest[1]);
      dyoff = 0.5 * (tapex[1] - tdest[1]) +
              b->offconstant * (tapex[0] - tdest[0]);
      /* If the off-center is closer to the destination than the */
      /*   circumcenter, use the off-center instead.             */
      if (dxoff * dxoff + dyoff * dyoff <
          (dx - xdo) * (dx - xdo) + (dy - ydo) * (dy - ydo)) {
        dx = xdo + dxoff;
        dy = ydo + dyoff;
      }
    }
  }

  circumcenter[0] = torg[0] + dx;
  circumcenter[1] = torg[1] + dy;

  /* To interpolate vertex attributes for the new vertex inserted at */
  /*   the circumcenter, define a coordinate system with a xi-axis,  */
  /*   directed from the triangle's origin to its destination, and   */
  /*   an eta-axis, directed from its origin to its apex.            */
  /*   Calculate the xi and eta coordinates of the circumcenter.     */
  *xi = (yao * dx - xao * dy) * (2.0 * denominator);
  *eta = (xdo * dy - ydo * dx) * (2.0 * denominator);
}

/**                                                                         **/
/**                                                                         **/
/********* Geometric primitives end here                             *********/

/*****************************************************************************/
/*                                                                           */
/*  triangleinit()   Initialize some variables.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void triangleinit(struct mesh *m)
#else /* not ANSI_DECLARATORS */
void triangleinit(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  poolzero(&m->vertices);
  poolzero(&m->triangles);
  poolzero(&m->subsegs);
  poolzero(&m->viri);
  poolzero(&m->badsubsegs);
  poolzero(&m->badtriangles);
  poolzero(&m->flipstackers);
  poolzero(&m->splaynodes);

  m->recenttri.tri = (triangle *) NULL; /* No triangle has been visited yet. */
  m->undeads = 0;                       /* No eliminated input vertices yet. */
  m->samples = 1;         /* Point location should take at least one sample. */
  m->checksegments = 0;   /* There are no segments in the triangulation yet. */
  m->checkquality = 0;     /* The quality triangulation stage has not begun. */
  m->incirclecount = m->counterclockcount = m->orient3dcount = 0;
  m->hyperbolacount = m->circletopcount = m->circumcentercount = 0;
  randomseed = 1;

  exactinit();                     /* Initialize exact arithmetic constants. */
}

/*****************************************************************************/
/*                                                                           */
/*  randomnation()   Generate a random number between 0 and `choices' - 1.   */
/*                                                                           */
/*  This is a simple linear congruential random number generator.  Hence, it */
/*  is a bad random number generator, but good enough for most randomized    */
/*  geometric algorithms.                                                    */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
unsigned long randomnation(unsigned int choices)
#else /* not ANSI_DECLARATORS */
unsigned long randomnation(choices)
unsigned int choices;
#endif /* not ANSI_DECLARATORS */

{
  randomseed = (randomseed * 1366l + 150889l) % 714025l;
  return randomseed / (714025l / choices + 1);
}

/********* Mesh quality testing routines begin here                  *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  checkmesh()   Test the mesh for topological consistency.                 */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void checkmesh(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void checkmesh(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop;
  struct otri oppotri, oppooppotri;
  vertex triorg, tridest, triapex;
  vertex oppoorg, oppodest;
  int horrors;
  int saveexact;
  triangle ptr;                         /* Temporary variable used by sym(). */

  /* Temporarily turn on exact arithmetic if it's off. */
  saveexact = b->noexact;
  b->noexact = 0;
  if (!b->quiet) {
    printf("  Checking consistency of mesh...\n");
  }
  horrors = 0;
  /* Run through the list of triangles, checking each one. */
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  while (triangleloop.tri != (triangle *) NULL) {
    /* Check all three edges of the triangle. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      org(triangleloop, triorg);
      dest(triangleloop, tridest);
      if (triangleloop.orient == 0) {       /* Only test for inversion once. */
        /* Test if the triangle is flat or inverted. */
        apex(triangleloop, triapex);
        if (counterclockwise(m, b, triorg, tridest, triapex) <= 0.0) {
          printf("  !! !! Inverted ");
          printtriangle(m, b, &triangleloop);
          horrors++;
        }
      }
      /* Find the neighboring triangle on this edge. */
      sym(triangleloop, oppotri);
      if (oppotri.tri != m->dummytri) {
        /* Check that the triangle's neighbor knows it's a neighbor. */
        sym(oppotri, oppooppotri);
        if ((triangleloop.tri != oppooppotri.tri)
            || (triangleloop.orient != oppooppotri.orient)) {
          printf("  !! !! Asymmetric triangle-triangle bond:\n");
          if (triangleloop.tri == oppooppotri.tri) {
            printf("   (Right triangle, wrong orientation)\n");
          }
          printf("    First ");
          printtriangle(m, b, &triangleloop);
          printf("    Second (nonreciprocating) ");
          printtriangle(m, b, &oppotri);
          horrors++;
        }
        /* Check that both triangles agree on the identities */
        /*   of their shared vertices.                       */
        org(oppotri, oppoorg);
        dest(oppotri, oppodest);
        if ((triorg != oppodest) || (tridest != oppoorg)) {
          printf("  !! !! Mismatched edge coordinates between two triangles:\n"
                 );
          printf("    First mismatched ");
          printtriangle(m, b, &triangleloop);
          printf("    Second mismatched ");
          printtriangle(m, b, &oppotri);
          horrors++;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }
  if (horrors == 0) {
    if (!b->quiet) {
      printf("  In my studied opinion, the mesh appears to be consistent.\n");
    }
  } else if (horrors == 1) {
    printf("  !! !! !! !! Precisely one festering wound discovered.\n");
  } else {
    printf("  !! !! !! !! %d abominations witnessed.\n", horrors);
  }
  /* Restore the status of exact arithmetic. */
  b->noexact = saveexact;
}

#endif /* not REDUCED */

/*****************************************************************************/
/*                                                                           */
/*  checkdelaunay()   Ensure that the mesh is (constrained) Delaunay.        */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void checkdelaunay(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void checkdelaunay(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop;
  struct otri oppotri;
  struct osub opposubseg;
  vertex triorg, tridest, triapex;
  vertex oppoapex;
  int shouldbedelaunay;
  int horrors;
  int saveexact;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  /* Temporarily turn on exact arithmetic if it's off. */
  saveexact = b->noexact;
  b->noexact = 0;
  if (!b->quiet) {
    printf("  Checking Delaunay property of mesh...\n");
  }
  horrors = 0;
  /* Run through the list of triangles, checking each one. */
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  while (triangleloop.tri != (triangle *) NULL) {
    /* Check all three edges of the triangle. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      org(triangleloop, triorg);
      dest(triangleloop, tridest);
      apex(triangleloop, triapex);
      sym(triangleloop, oppotri);
      apex(oppotri, oppoapex);
      /* Only test that the edge is locally Delaunay if there is an   */
      /*   adjoining triangle whose pointer is larger (to ensure that */
      /*   each pair isn't tested twice).                             */
      shouldbedelaunay = (oppotri.tri != m->dummytri) &&
            !deadtri(oppotri.tri) && (triangleloop.tri < oppotri.tri) &&
            (triorg != m->infvertex1) && (triorg != m->infvertex2) &&
            (triorg != m->infvertex3) &&
            (tridest != m->infvertex1) && (tridest != m->infvertex2) &&
            (tridest != m->infvertex3) &&
            (triapex != m->infvertex1) && (triapex != m->infvertex2) &&
            (triapex != m->infvertex3) &&
            (oppoapex != m->infvertex1) && (oppoapex != m->infvertex2) &&
            (oppoapex != m->infvertex3);
      if (m->checksegments && shouldbedelaunay) {
        /* If a subsegment separates the triangles, then the edge is */
        /*   constrained, so no local Delaunay test should be done.  */
        tspivot(triangleloop, opposubseg);
        if (opposubseg.ss != m->dummysub){
          shouldbedelaunay = 0;
        }
      }
      if (shouldbedelaunay) {
        if (nonregular(m, b, triorg, tridest, triapex, oppoapex) > 0.0) {
          if (!b->weighted) {
            printf("  !! !! Non-Delaunay pair of triangles:\n");
            printf("    First non-Delaunay ");
            printtriangle(m, b, &triangleloop);
            printf("    Second non-Delaunay ");
          } else {
            printf("  !! !! Non-regular pair of triangles:\n");
            printf("    First non-regular ");
            printtriangle(m, b, &triangleloop);
            printf("    Second non-regular ");
          }
          printtriangle(m, b, &oppotri);
          horrors++;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }
  if (horrors == 0) {
    if (!b->quiet) {
      printf(
  "  By virtue of my perceptive intelligence, I declare the mesh Delaunay.\n");
    }
  } else if (horrors == 1) {
    printf(
         "  !! !! !! !! Precisely one terrifying transgression identified.\n");
  } else {
    printf("  !! !! !! !! %d obscenities viewed with horror.\n", horrors);
  }
  /* Restore the status of exact arithmetic. */
  b->noexact = saveexact;
}

#endif /* not REDUCED */

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtriang()   Add a bad triangle data structure to the end of a   */
/*                       queue.                                              */
/*                                                                           */
/*  The queue is actually a set of 4096 queues.  I use multiple queues to    */
/*  give priority to smaller angles.  I originally implemented a heap, but   */
/*  the queues are faster by a larger margin than I'd suspected.             */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void enqueuebadtriang(struct mesh *m, struct behavior *b,
                      struct badtriang *badtri)
#else /* not ANSI_DECLARATORS */
void enqueuebadtriang(m, b, badtri)
struct mesh *m;
struct behavior *b;
struct badtriang *badtri;
#endif /* not ANSI_DECLARATORS */

{
  REAL length, multiplier;
  int exponent, expincrement;
  int queuenumber;
  int posexponent;
  int i;

  if (b->verbose > 2) {
    printf("  Queueing bad triangle:\n");
    printf("    (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
           badtri->triangorg[0], badtri->triangorg[1],
           badtri->triangdest[0], badtri->triangdest[1],
           badtri->triangapex[0], badtri->triangapex[1]);
  }

  /* Determine the appropriate queue to put the bad triangle into.    */
  /*   Recall that the key is the square of its shortest edge length. */
  if (badtri->key >= 1.0) {
    length = badtri->key;
    posexponent = 1;
  } else {
    /* `badtri->key' is 2.0 to a negative exponent, so we'll record that */
    /*   fact and use the reciprocal of `badtri->key', which is > 1.0.   */
    length = 1.0 / badtri->key;
    posexponent = 0;
  }
  /* `length' is approximately 2.0 to what exponent?  The following code */
  /*   determines the answer in time logarithmic in the exponent.        */
  exponent = 0;
  while (length > 2.0) {
    /* Find an approximation by repeated squaring of two. */
    expincrement = 1;
    multiplier = 0.5;
    while (length * multiplier * multiplier > 1.0) {
      expincrement *= 2;
      multiplier *= multiplier;
    }
    /* Reduce the value of `length', then iterate if necessary. */
    exponent += expincrement;
    length *= multiplier;
  }
  /* `length' is approximately squareroot(2.0) to what exponent? */
  exponent = 2.0 * exponent + (length > SQUAREROOTTWO);
  /* `exponent' is now in the range 0...2047 for IEEE double precision.   */
  /*   Choose a queue in the range 0...4095.  The shortest edges have the */
  /*   highest priority (queue 4095).                                     */
  if (posexponent) {
    queuenumber = 2047 - exponent;
  } else {
    queuenumber = 2048 + exponent;
  }

  /* Are we inserting into an empty queue? */
  if (m->queuefront[queuenumber] == (struct badtriang *) NULL) {
    /* Yes, we are inserting into an empty queue.     */
    /*   Will this become the highest-priority queue? */
    if (queuenumber > m->firstnonemptyq) {
      /* Yes, this is the highest-priority queue. */
      m->nextnonemptyq[queuenumber] = m->firstnonemptyq;
      m->firstnonemptyq = queuenumber;
    } else {
      /* No, this is not the highest-priority queue. */
      /*   Find the queue with next higher priority. */
      i = queuenumber + 1;
      while (m->queuefront[i] == (struct badtriang *) NULL) {
        i++;
      }
      /* Mark the newly nonempty queue as following a higher-priority queue. */
      m->nextnonemptyq[queuenumber] = m->nextnonemptyq[i];
      m->nextnonemptyq[i] = queuenumber;
    }
    /* Put the bad triangle at the beginning of the (empty) queue. */
    m->queuefront[queuenumber] = badtri;
  } else {
    /* Add the bad triangle to the end of an already nonempty queue. */
    m->queuetail[queuenumber]->nexttriang = badtri;
  }
  /* Maintain a pointer to the last triangle of the queue. */
  m->queuetail[queuenumber] = badtri;
  /* Newly enqueued bad triangle has no successor in the queue. */
  badtri->nexttriang = (struct badtriang *) NULL;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtri()   Add a bad triangle to the end of a queue.              */
/*                                                                           */
/*  Allocates a badtriang data structure for the triangle, then passes it to */
/*  enqueuebadtriang().                                                      */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void enqueuebadtri(struct mesh *m, struct behavior *b, struct otri *enqtri,
                   REAL minedge, vertex enqapex, vertex enqorg, vertex enqdest)
#else /* not ANSI_DECLARATORS */
void enqueuebadtri(m, b, enqtri, minedge, enqapex, enqorg, enqdest)
struct mesh *m;
struct behavior *b;
struct otri *enqtri;
REAL minedge;
vertex enqapex;
vertex enqorg;
vertex enqdest;
#endif /* not ANSI_DECLARATORS */

{
  struct badtriang *newbad;

  /* Allocate space for the bad triangle. */
  newbad = (struct badtriang *) poolalloc(&m->badtriangles);
  newbad->poortri = encode(*enqtri);
  newbad->key = minedge;
  newbad->triangapex = enqapex;
  newbad->triangorg = enqorg;
  newbad->triangdest = enqdest;
  enqueuebadtriang(m, b, newbad);
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  dequeuebadtriang()   Remove a triangle from the front of the queue.      */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
struct badtriang *dequeuebadtriang(struct mesh *m)
#else /* not ANSI_DECLARATORS */
struct badtriang *dequeuebadtriang(m)
struct mesh *m;
#endif /* not ANSI_DECLARATORS */

{
  struct badtriang *result;

  /* If no queues are nonempty, return NULL. */
  if (m->firstnonemptyq < 0) {
    return (struct badtriang *) NULL;
  }
  /* Find the first triangle of the highest-priority queue. */
  result = m->queuefront[m->firstnonemptyq];
  /* Remove the triangle from the queue. */
  m->queuefront[m->firstnonemptyq] = result->nexttriang;
  /* If this queue is now empty, note the new highest-priority */
  /*   nonempty queue.                                         */
  if (result == m->queuetail[m->firstnonemptyq]) {
    m->firstnonemptyq = m->nextnonemptyq[m->firstnonemptyq];
  }
  return result;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  checkseg4encroach()   Check a subsegment to see if it is encroached; add */
/*                        it to the list if it is.                           */
/*                                                                           */
/*  A subsegment is encroached if there is a vertex in its diametral lens.   */
/*  For Ruppert's algorithm (-D switch), the "diametral lens" is the         */
/*  diametral circle.  For Chew's algorithm (default), the diametral lens is */
/*  just big enough to enclose two isosceles triangles whose bases are the   */
/*  subsegment.  Each of the two isosceles triangles has two angles equal    */
/*  to `b->minangle'.                                                        */
/*                                                                           */
/*  Chew's algorithm does not require diametral lenses at all--but they save */
/*  time.  Any vertex inside a subsegment's diametral lens implies that the  */
/*  triangle adjoining the subsegment will be too skinny, so it's only a     */
/*  matter of time before the encroaching vertex is deleted by Chew's        */
/*  algorithm.  It's faster to simply not insert the doomed vertex in the    */
/*  first place, which is why I use diametral lenses with Chew's algorithm.  */
/*                                                                           */
/*  Returns a nonzero value if the subsegment is encroached.                 */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
int checkseg4encroach(struct mesh *m, struct behavior *b,
                      struct osub *testsubseg)
#else /* not ANSI_DECLARATORS */
int checkseg4encroach(m, b, testsubseg)
struct mesh *m;
struct behavior *b;
struct osub *testsubseg;
#endif /* not ANSI_DECLARATORS */

{
  struct otri neighbortri;
  struct osub testsym;
  struct badsubseg *encroachedseg;
  REAL dotproduct;
  int encroached;
  int sides;
  vertex eorg, edest, eapex;
  triangle ptr;                     /* Temporary variable used by stpivot(). */

  encroached = 0;
  sides = 0;

  sorg(*testsubseg, eorg);
  sdest(*testsubseg, edest);
  /* Check one neighbor of the subsegment. */
  stpivot(*testsubseg, neighbortri);
  /* Does the neighbor exist, or is this a boundary edge? */
  if (neighbortri.tri != m->dummytri) {
    sides++;
    /* Find a vertex opposite this subsegment. */
    apex(neighbortri, eapex);
    /* Check whether the apex is in the diametral lens of the subsegment */
    /*   (the diametral circle if `conformdel' is set).  A dot product   */
    /*   of two sides of the triangle is used to check whether the angle */
    /*   at the apex is greater than (180 - 2 `minangle') degrees (for   */
    /*   lenses; 90 degrees for diametral circles).                      */
    dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                 (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
    if (dotproduct < 0.0) {
      if (b->conformdel ||
          (dotproduct * dotproduct >=
           (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
           ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
            (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
           ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
            (edest[1] - eapex[1]) * (edest[1] - eapex[1])))) {
        encroached = 1;
      }
    }
  }
  /* Check the other neighbor of the subsegment. */
  ssym(*testsubseg, testsym);
  stpivot(testsym, neighbortri);
  /* Does the neighbor exist, or is this a boundary edge? */
  if (neighbortri.tri != m->dummytri) {
    sides++;
    /* Find the other vertex opposite this subsegment. */
    apex(neighbortri, eapex);
    /* Check whether the apex is in the diametral lens of the subsegment */
    /*   (or the diametral circle, if `conformdel' is set).              */
    dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                 (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
    if (dotproduct < 0.0) {
      if (b->conformdel ||
          (dotproduct * dotproduct >=
           (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
           ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
            (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
           ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
            (edest[1] - eapex[1]) * (edest[1] - eapex[1])))) {
        encroached += 2;
      }
    }
  }

  if (encroached && (!b->nobisect || ((b->nobisect == 1) && (sides == 2)))) {
    if (b->verbose > 2) {
      printf(
        "  Queueing encroached subsegment (%.12g, %.12g) (%.12g, %.12g).\n",
        eorg[0], eorg[1], edest[0], edest[1]);
    }
    /* Add the subsegment to the list of encroached subsegments. */
    /*   Be sure to get the orientation right.                   */
    encroachedseg = (struct badsubseg *) poolalloc(&m->badsubsegs);
    if (encroached == 1) {
      encroachedseg->encsubseg = sencode(*testsubseg);
      encroachedseg->subsegorg = eorg;
      encroachedseg->subsegdest = edest;
    } else {
      encroachedseg->encsubseg = sencode(testsym);
      encroachedseg->subsegorg = edest;
      encroachedseg->subsegdest = eorg;
    }
  }

  return encroached;
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  testtriangle()   Test a triangle for quality and size.                   */
/*                                                                           */
/*  Tests a triangle to see if it satisfies the minimum angle condition and  */
/*  the maximum area condition.  Triangles that aren't up to spec are added  */
/*  to the bad triangle queue.                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void testtriangle(struct mesh *m, struct behavior *b, struct otri *testtri)
#else /* not ANSI_DECLARATORS */
void testtriangle(m, b, testtri)
struct mesh *m;
struct behavior *b;
struct otri *testtri;
#endif /* not ANSI_DECLARATORS */

{
  struct otri tri1, tri2;
  struct osub testsub;
  vertex torg, tdest, tapex;
  vertex base1, base2;
  vertex org1, dest1, org2, dest2;
  vertex joinvertex;
  REAL dxod, dyod, dxda, dyda, dxao, dyao;
  REAL dxod2, dyod2, dxda2, dyda2, dxao2, dyao2;
  REAL apexlen, orglen, destlen, minedge;
  REAL angle;
  REAL area;
  REAL dist1, dist2;
  subseg sptr;                      /* Temporary variable used by tspivot(). */
  triangle ptr;           /* Temporary variable used by oprev() and dnext(). */

  org(*testtri, torg);
  dest(*testtri, tdest);
  apex(*testtri, tapex);
  dxod = torg[0] - tdest[0];
  dyod = torg[1] - tdest[1];
  dxda = tdest[0] - tapex[0];
  dyda = tdest[1] - tapex[1];
  dxao = tapex[0] - torg[0];
  dyao = tapex[1] - torg[1];
  dxod2 = dxod * dxod;
  dyod2 = dyod * dyod;
  dxda2 = dxda * dxda;
  dyda2 = dyda * dyda;
  dxao2 = dxao * dxao;
  dyao2 = dyao * dyao;
  /* Find the lengths of the triangle's three edges. */
  apexlen = dxod2 + dyod2;
  orglen = dxda2 + dyda2;
  destlen = dxao2 + dyao2;

  if ((apexlen < orglen) && (apexlen < destlen)) {
    /* The edge opposite the apex is shortest. */
    minedge = apexlen;
    /* Find the square of the cosine of the angle at the apex. */
    angle = dxda * dxao + dyda * dyao;
    angle = angle * angle / (orglen * destlen);
    base1 = torg;
    base2 = tdest;
    otricopy(*testtri, tri1);
  } else if (orglen < destlen) {
    /* The edge opposite the origin is shortest. */
    minedge = orglen;
    /* Find the square of the cosine of the angle at the origin. */
    angle = dxod * dxao + dyod * dyao;
    angle = angle * angle / (apexlen * destlen);
    base1 = tdest;
    base2 = tapex;
    lnext(*testtri, tri1);
  } else {
    /* The edge opposite the destination is shortest. */
    minedge = destlen;
    /* Find the square of the cosine of the angle at the destination. */
    angle = dxod * dxda + dyod * dyda;
    angle = angle * angle / (apexlen * orglen);
    base1 = tapex;
    base2 = torg;
    lprev(*testtri, tri1);
  }

  if (b->vararea || b->fixedarea || b->usertest) {
    /* Check whether the area is larger than permitted. */
    area = 0.5 * (dxod * dyda - dyod * dxda);
    if (b->fixedarea && (area > b->maxarea)) {
      /* Add this triangle to the list of bad triangles. */
      enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
      return;
    }

    /* Nonpositive area constraints are treated as unconstrained. */
    if ((b->vararea) && (area > areabound(*testtri)) &&
        (areabound(*testtri) > 0.0)) {
      /* Add this triangle to the list of bad triangles. */
      enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
      return;
    }

    if (b->usertest) {
      /* Check whether the user thinks this triangle is too large. */
      if (triunsuitable(torg, tdest, tapex, area)) {
        enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
        return;
      }
    }
  }

  /* Check whether the angle is smaller than permitted. */
  if (angle > b->goodangle) {
    /* Use the rules of Miller, Pav, and Walkington to decide that certain */
    /*   triangles should not be split, even if they have bad angles.      */
    /*   A skinny triangle is not split if its shortest edge subtends a    */
    /*   small input angle, and both endpoints of the edge lie on a        */
    /*   concentric circular shell.  For convenience, I make a small       */
    /*   adjustment to that rule:  I check if the endpoints of the edge    */
    /*   both lie in segment interiors, equidistant from the apex where    */
    /*   the two segments meet.                                            */
    /* First, check if both points lie in segment interiors.               */
    if ((vertextype(base1) == SEGMENTVERTEX) &&
        (vertextype(base2) == SEGMENTVERTEX)) {
      /* Check if both points lie in a common segment.  If they do, the */
      /*   skinny triangle is enqueued to be split as usual.            */
      tspivot(tri1, testsub);
      if (testsub.ss == m->dummysub) {
        /* No common segment.  Find a subsegment that contains `torg'. */
        otricopy(tri1, tri2);
        do {
          oprevself(tri1);
          tspivot(tri1, testsub);
        } while (testsub.ss == m->dummysub);
        /* Find the endpoints of the containing segment. */
        segorg(testsub, org1);
        segdest(testsub, dest1);
        /* Find a subsegment that contains `tdest'. */
        do {
          dnextself(tri2);
          tspivot(tri2, testsub);
        } while (testsub.ss == m->dummysub);
        /* Find the endpoints of the containing segment. */
        segorg(testsub, org2);
        segdest(testsub, dest2);
        /* Check if the two containing segments have an endpoint in common. */
        joinvertex = (vertex) NULL;
        if ((dest1[0] == org2[0]) && (dest1[1] == org2[1])) {
          joinvertex = dest1;
        } else if ((org1[0] == dest2[0]) && (org1[1] == dest2[1])) {
          joinvertex = org1;
        }
        if (joinvertex != (vertex) NULL) {
          /* Compute the distance from the common endpoint (of the two  */
          /*   segments) to each of the endpoints of the shortest edge. */
          dist1 = ((base1[0] - joinvertex[0]) * (base1[0] - joinvertex[0]) +
                   (base1[1] - joinvertex[1]) * (base1[1] - joinvertex[1]));
          dist2 = ((base2[0] - joinvertex[0]) * (base2[0] - joinvertex[0]) +
                   (base2[1] - joinvertex[1]) * (base2[1] - joinvertex[1]));
          /* If the two distances are equal, don't split the triangle. */
          if ((dist1 < 1.001 * dist2) && (dist1 > 0.999 * dist2)) {
            /* Return now to avoid enqueueing the bad triangle. */
            return;
          }
        }
      }
    }

    /* Add this triangle to the list of bad triangles. */
    enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
  }
}

#endif /* not CDT_ONLY */

/**                                                                         **/
/**                                                                         **/
/********* Mesh quality testing routines end here                    *********/

/********* Point location routines begin here                        *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  makevertexmap()   Construct a mapping from vertices to triangles to      */
/*                    improve the speed of point location for segment        */
/*                    insertion.                                             */
/*                                                                           */
/*  Traverses all the triangles, and provides each corner of each triangle   */
/*  with a pointer to that triangle.  Of course, pointers will be            */
/*  overwritten by other pointers because (almost) each vertex is a corner   */
/*  of several triangles, but in the end every vertex will point to some     */
/*  triangle that contains it.                                               */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void makevertexmap(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void makevertexmap(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop;
  vertex triorg;

  if (b->verbose) {
    printf("    Constructing mapping from vertices to triangles.\n");
  }
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  while (triangleloop.tri != (triangle *) NULL) {
    /* Check all three vertices of the triangle. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      org(triangleloop, triorg);
      setvertex2tri(triorg, encode(triangleloop));
    }
    triangleloop.tri = triangletraverse(m);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  preciselocate()   Find a triangle or edge containing a given point.      */
/*                                                                           */
/*  Begins its search from `searchtri'.  It is important that `searchtri'    */
/*  be a handle with the property that `searchpoint' is strictly to the left */
/*  of the edge denoted by `searchtri', or is collinear with that edge and   */
/*  does not intersect that edge.  (In particular, `searchpoint' should not  */
/*  be the origin or destination of that edge.)                              */
/*                                                                           */
/*  These conditions are imposed because preciselocate() is normally used in */
/*  one of two situations:                                                   */
/*                                                                           */
/*  (1)  To try to find the location to insert a new point.  Normally, we    */
/*       know an edge that the point is strictly to the left of.  In the     */
/*       incremental Delaunay algorithm, that edge is a bounding box edge.   */
/*       In Ruppert's Delaunay refinement algorithm for quality meshing,     */
/*       that edge is the shortest edge of the triangle whose circumcenter   */
/*       is being inserted.                                                  */
/*                                                                           */
/*  (2)  To try to find an existing point.  In this case, any edge on the    */
/*       convex hull is a good starting edge.  You must screen out the       */
/*       possibility that the vertex sought is an endpoint of the starting   */
/*       edge before you call preciselocate().                               */
/*                                                                           */
/*  On completion, `searchtri' is a triangle that contains `searchpoint'.    */
/*                                                                           */
/*  This implementation differs from that given by Guibas and Stolfi.  It    */
/*  walks from triangle to triangle, crossing an edge only if `searchpoint'  */
/*  is on the other side of the line containing that edge.  After entering   */
/*  a triangle, there are two edges by which one can leave that triangle.    */
/*  If both edges are valid (`searchpoint' is on the other side of both      */
/*  edges), one of the two is chosen by drawing a line perpendicular to      */
/*  the entry edge (whose endpoints are `forg' and `fdest') passing through  */
/*  `fapex'.  Depending on which side of this perpendicular `searchpoint'    */
/*  falls on, an exit edge is chosen.                                        */
/*                                                                           */
/*  This implementation is empirically faster than the Guibas and Stolfi     */
/*  point location routine (which I originally used), which tends to spiral  */
/*  in toward its target.                                                    */
/*                                                                           */
/*  Returns ONVERTEX if the point lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the point lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the point lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the point lies strictly within a triangle.         */
/*  `searchtri' is a handle on the triangle that contains the point.         */
/*                                                                           */
/*  Returns OUTSIDE if the point lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the point is to the right of.  This might      */
/*  occur when the circumcenter of a triangle falls just slightly outside    */
/*  the mesh due to floating-point roundoff error.  It also occurs when      */
/*  seeking a hole or region point that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  If `stopatsubsegment' is nonzero, the search will stop if it tries to    */
/*  walk through a subsegment, and will return OUTSIDE.                      */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*  However, it can still be used to find the circumcenter of a triangle, as */
/*  long as the search is begun from the triangle in question.               */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
enum locateresult preciselocate(struct mesh *m, struct behavior *b,
                                vertex searchpoint, struct otri *searchtri,
                                int stopatsubsegment)
#else /* not ANSI_DECLARATORS */
enum locateresult preciselocate(m, b, searchpoint, searchtri, stopatsubsegment)
struct mesh *m;
struct behavior *b;
vertex searchpoint;
struct otri *searchtri;
int stopatsubsegment;
#endif /* not ANSI_DECLARATORS */

{
  struct otri backtracktri;
  struct osub checkedge;
  vertex forg, fdest, fapex;
  REAL orgorient, destorient;
  int moveleft;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (b->verbose > 2) {
    printf("  Searching for point (%.12g, %.12g).\n",
           searchpoint[0], searchpoint[1]);
  }
  /* Where are we? */
  org(*searchtri, forg);
  dest(*searchtri, fdest);
  apex(*searchtri, fapex);
  while (1) {
    if (b->verbose > 2) {
      printf("    At (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             forg[0], forg[1], fdest[0], fdest[1], fapex[0], fapex[1]);
    }
    /* Check whether the apex is the point we seek. */
    if ((fapex[0] == searchpoint[0]) && (fapex[1] == searchpoint[1])) {
      lprevself(*searchtri);
      return ONVERTEX;
    }
    /* Does the point lie on the other side of the line defined by the */
    /*   triangle edge opposite the triangle's destination?            */
    destorient = counterclockwise(m, b, forg, fapex, searchpoint);
    /* Does the point lie on the other side of the line defined by the */
    /*   triangle edge opposite the triangle's origin?                 */
    orgorient = counterclockwise(m, b, fapex, fdest, searchpoint);
    if (destorient > 0.0) {
      if (orgorient > 0.0) {
        /* Move left if the inner product of (fapex - searchpoint) and  */
        /*   (fdest - forg) is positive.  This is equivalent to drawing */
        /*   a line perpendicular to the line (forg, fdest) and passing */
        /*   through `fapex', and determining which side of this line   */
        /*   `searchpoint' falls on.                                    */
        moveleft = (fapex[0] - searchpoint[0]) * (fdest[0] - forg[0]) +
                   (fapex[1] - searchpoint[1]) * (fdest[1] - forg[1]) > 0.0;
      } else {
        moveleft = 1;
      }
    } else {
      if (orgorient > 0.0) {
        moveleft = 0;
      } else {
        /* The point we seek must be on the boundary of or inside this */
        /*   triangle.                                                 */
        if (destorient == 0.0) {
          lprevself(*searchtri);
          return ONEDGE;
        }
        if (orgorient == 0.0) {
          lnextself(*searchtri);
          return ONEDGE;
        }
        return INTRIANGLE;
      }
    }

    /* Move to another triangle.  Leave a trace `backtracktri' in case */
    /*   floating-point roundoff or some such bogey causes us to walk  */
    /*   off a boundary of the triangulation.                          */
    if (moveleft) {
      lprev(*searchtri, backtracktri);
      fdest = fapex;
    } else {
      lnext(*searchtri, backtracktri);
      forg = fapex;
    }
    sym(backtracktri, *searchtri);

    if (m->checksegments && stopatsubsegment) {
      /* Check for walking through a subsegment. */
      tspivot(backtracktri, checkedge);
      if (checkedge.ss != m->dummysub) {
        /* Go back to the last triangle. */
        otricopy(backtracktri, *searchtri);
        return OUTSIDE;
      }
    }
    /* Check for walking right out of the triangulation. */
    if (searchtri->tri == m->dummytri) {
      /* Go back to the last triangle. */
      otricopy(backtracktri, *searchtri);
      return OUTSIDE;
    }

    apex(*searchtri, fapex);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  locate()   Find a triangle or edge containing a given point.             */
/*                                                                           */
/*  Searching begins from one of:  the input `searchtri', a recently         */
/*  encountered triangle `recenttri', or from a triangle chosen from a       */
/*  random sample.  The choice is made by determining which triangle's       */
/*  origin is closest to the point we are searching for.  Normally,          */
/*  `searchtri' should be a handle on the convex hull of the triangulation.  */
/*                                                                           */
/*  Details on the random sampling method can be found in the Mucke, Saias,  */
/*  and Zhu paper cited in the header of this code.                          */
/*                                                                           */
/*  On completion, `searchtri' is a triangle that contains `searchpoint'.    */
/*                                                                           */
/*  Returns ONVERTEX if the point lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the point lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the point lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the point lies strictly within a triangle.         */
/*  `searchtri' is a handle on the triangle that contains the point.         */
/*                                                                           */
/*  Returns OUTSIDE if the point lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the point is to the right of.  This might      */
/*  occur when the circumcenter of a triangle falls just slightly outside    */
/*  the mesh due to floating-point roundoff error.  It also occurs when      */
/*  seeking a hole or region point that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
enum locateresult locate(struct mesh *m, struct behavior *b,
                         vertex searchpoint, struct otri *searchtri)
#else /* not ANSI_DECLARATORS */
enum locateresult locate(m, b, searchpoint, searchtri)
struct mesh *m;
struct behavior *b;
vertex searchpoint;
struct otri *searchtri;
#endif /* not ANSI_DECLARATORS */

{
  VOID **sampleblock;
  char *firsttri;
  struct otri sampletri;
  vertex torg, tdest;
  unsigned long alignptr;
  REAL searchdist, dist;
  REAL ahead;
  long samplesperblock, totalsamplesleft, samplesleft;
  long population, totalpopulation;
  triangle ptr;                         /* Temporary variable used by sym(). */

  if (b->verbose > 2) {
    printf("  Randomly sampling for a triangle near point (%.12g, %.12g).\n",
           searchpoint[0], searchpoint[1]);
  }
  /* Record the distance from the suggested starting triangle to the */
  /*   point we seek.                                                */
  org(*searchtri, torg);
  searchdist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
               (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
  if (b->verbose > 2) {
    printf("    Boundary triangle has origin (%.12g, %.12g).\n",
           torg[0], torg[1]);
  }

  /* If a recently encountered triangle has been recorded and has not been */
  /*   deallocated, test it as a good starting point.                      */
  if (m->recenttri.tri != (triangle *) NULL) {
    if (!deadtri(m->recenttri.tri)) {
      org(m->recenttri, torg);
      if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1])) {
        otricopy(m->recenttri, *searchtri);
        return ONVERTEX;
      }
      dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
             (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
      if (dist < searchdist) {
        otricopy(m->recenttri, *searchtri);
        searchdist = dist;
        if (b->verbose > 2) {
          printf("    Choosing recent triangle with origin (%.12g, %.12g).\n",
                 torg[0], torg[1]);
        }
      }
    }
  }

  /* The number of random samples taken is proportional to the cube root of */
  /*   the number of triangles in the mesh.  The next bit of code assumes   */
  /*   that the number of triangles increases monotonically (or at least    */
  /*   doesn't decrease enough to matter).                                  */
  while (SAMPLEFACTOR * m->samples * m->samples * m->samples <
         m->triangles.items) {
    m->samples++;
  }

  /* We'll draw ceiling(samples * TRIPERBLOCK / maxitems) random samples  */
  /*   from each block of triangles (except the first)--until we meet the */
  /*   sample quota.  The ceiling means that blocks at the end might be   */
  /*   neglected, but I don't care.                                       */
  samplesperblock = (m->samples * TRIPERBLOCK - 1) / m->triangles.maxitems + 1;
  /* We'll draw ceiling(samples * itemsfirstblock / maxitems) random samples */
  /*   from the first block of triangles.                                    */
  samplesleft = (m->samples * m->triangles.itemsfirstblock - 1) /
                m->triangles.maxitems + 1;
  totalsamplesleft = m->samples;
  population = m->triangles.itemsfirstblock;
  totalpopulation = m->triangles.maxitems;
  sampleblock = m->triangles.firstblock;
  sampletri.orient = 0;
  while (totalsamplesleft > 0) {
    /* If we're in the last block, `population' needs to be corrected. */
    if (population > totalpopulation) {
      population = totalpopulation;
    }
    /* Find a pointer to the first triangle in the block. */
    alignptr = (unsigned long) (sampleblock + 1);
    firsttri = (char *) (alignptr +
                         (unsigned long) m->triangles.alignbytes -
                         (alignptr %
                          (unsigned long) m->triangles.alignbytes));

    /* Choose `samplesleft' randomly sampled triangles in this block. */
    do {
      sampletri.tri = (triangle *) (firsttri +
                                    (randomnation((unsigned int) population) *
                                     m->triangles.itembytes));
      if (!deadtri(sampletri.tri)) {
        org(sampletri, torg);
        dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
               (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
        if (dist < searchdist) {
          otricopy(sampletri, *searchtri);
          searchdist = dist;
          if (b->verbose > 2) {
            printf("    Choosing triangle with origin (%.12g, %.12g).\n",
                   torg[0], torg[1]);
          }
        }
      }

      samplesleft--;
      totalsamplesleft--;
    } while ((samplesleft > 0) && (totalsamplesleft > 0));

    if (totalsamplesleft > 0) {
      sampleblock = (VOID **) *sampleblock;
      samplesleft = samplesperblock;
      totalpopulation -= population;
      population = TRIPERBLOCK;
    }
  }

  /* Where are we? */
  org(*searchtri, torg);
  dest(*searchtri, tdest);
  /* Check the starting triangle's vertices. */
  if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1])) {
    return ONVERTEX;
  }
  if ((tdest[0] == searchpoint[0]) && (tdest[1] == searchpoint[1])) {
    lnextself(*searchtri);
    return ONVERTEX;
  }
  /* Orient `searchtri' to fit the preconditions of calling preciselocate(). */
  ahead = counterclockwise(m, b, torg, tdest, searchpoint);
  if (ahead < 0.0) {
    /* Turn around so that `searchpoint' is to the left of the */
    /*   edge specified by `searchtri'.                        */
    symself(*searchtri);
  } else if (ahead == 0.0) {
    /* Check if `searchpoint' is between `torg' and `tdest'. */
    if (((torg[0] < searchpoint[0]) == (searchpoint[0] < tdest[0])) &&
        ((torg[1] < searchpoint[1]) == (searchpoint[1] < tdest[1]))) {
      return ONEDGE;
    }
  }
  return preciselocate(m, b, searchpoint, searchtri, 0);
}

/**                                                                         **/
/**                                                                         **/
/********* Point location routines end here                          *********/

/********* Mesh transformation routines begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  insertsubseg()   Create a new subsegment and insert it between two       */
/*                   triangles.                                              */
/*                                                                           */
/*  The new subsegment is inserted at the edge described by the handle       */
/*  `tri'.  Its vertices are properly initialized.  The marker `subsegmark'  */
/*  is applied to the subsegment and, if appropriate, its vertices.          */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void insertsubseg(struct mesh *m, struct behavior *b, struct otri *tri,
                  int subsegmark)
#else /* not ANSI_DECLARATORS */
void insertsubseg(m, b, tri, subsegmark)
struct mesh *m;
struct behavior *b;
struct otri *tri;             /* Edge at which to insert the new subsegment. */
int subsegmark;                            /* Marker for the new subsegment. */
#endif /* not ANSI_DECLARATORS */

{
  struct otri oppotri;
  struct osub newsubseg;
  vertex triorg, tridest;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  org(*tri, triorg);
  dest(*tri, tridest);
  /* Mark vertices if possible. */
  if (vertexmark(triorg) == 0) {
    setvertexmark(triorg, subsegmark);
  }
  if (vertexmark(tridest) == 0) {
    setvertexmark(tridest, subsegmark);
  }
  /* Check if there's already a subsegment here. */
  tspivot(*tri, newsubseg);
  if (newsubseg.ss == m->dummysub) {
    /* Make new subsegment and initialize its vertices. */
    makesubseg(m, &newsubseg);
    setsorg(newsubseg, tridest);
    setsdest(newsubseg, triorg);
    setsegorg(newsubseg, tridest);
    setsegdest(newsubseg, triorg);
    /* Bond new subsegment to the two triangles it is sandwiched between. */
    /*   Note that the facing triangle `oppotri' might be equal to        */
    /*   `dummytri' (outer space), but the new subsegment is bonded to it */
    /*   all the same.                                                    */
    tsbond(*tri, newsubseg);
    sym(*tri, oppotri);
    ssymself(newsubseg);
    tsbond(oppotri, newsubseg);
    setmark(newsubseg, subsegmark);
    if (b->verbose > 2) {
      printf("  Inserting new ");
      printsubseg(m, b, &newsubseg);
    }
  } else {
    if (mark(newsubseg) == 0) {
      setmark(newsubseg, subsegmark);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  Terminology                                                              */
/*                                                                           */
/*  A "local transformation" replaces a small set of triangles with another  */
/*  set of triangles.  This may or may not involve inserting or deleting a   */
/*  vertex.                                                                  */
/*                                                                           */
/*  The term "casing" is used to describe the set of triangles that are      */
/*  attached to the triangles being transformed, but are not transformed     */
/*  themselves.  Think of the casing as a fixed hollow structure inside      */
/*  which all the action happens.  A "casing" is only defined relative to    */
/*  a single transformation; each occurrence of a transformation will        */
/*  involve a different casing.                                              */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  flip()   Transform two triangles to two different triangles by flipping  */
/*           an edge counterclockwise within a quadrilateral.                */
/*                                                                           */
/*  Imagine the original triangles, abc and bad, oriented so that the        */
/*  shared edge ab lies in a horizontal plane, with the vertex b on the left */
/*  and the vertex a on the right.  The vertex c lies below the edge, and    */
/*  the vertex d lies above the edge.  The `flipedge' handle holds the edge  */
/*  ab of triangle abc, and is directed left, from vertex a to vertex b.     */
/*                                                                           */
/*  The triangles abc and bad are deleted and replaced by the triangles cdb  */
/*  and dca.  The triangles that represent abc and bad are NOT deallocated;  */
/*  they are reused for dca and cdb, respectively.  Hence, any handles that  */
/*  may have held the original triangles are still valid, although not       */
/*  directed as they were before.                                            */
/*                                                                           */
/*  Upon completion of this routine, the `flipedge' handle holds the edge    */
/*  dc of triangle dca, and is directed down, from vertex d to vertex c.     */
/*  (Hence, the two triangles have rotated counterclockwise.)                */
/*                                                                           */
/*  WARNING:  This transformation is geometrically valid only if the         */
/*  quadrilateral adbc is convex.  Furthermore, this transformation is       */
/*  valid only if there is not a subsegment between the triangles abc and    */
/*  bad.  This routine does not check either of these preconditions, and     */
/*  it is the responsibility of the calling routine to ensure that they are  */
/*  met.  If they are not, the streets shall be filled with wailing and      */
/*  gnashing of teeth.                                                       */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void flip(struct mesh *m, struct behavior *b, struct otri *flipedge)
#else /* not ANSI_DECLARATORS */
void flip(m, b, flipedge)
struct mesh *m;
struct behavior *b;
struct otri *flipedge;                    /* Handle for the triangle abc. */
#endif /* not ANSI_DECLARATORS */

{
  struct otri botleft, botright;
  struct otri topleft, topright;
  struct otri top;
  struct otri botlcasing, botrcasing;
  struct otri toplcasing, toprcasing;
  struct osub botlsubseg, botrsubseg;
  struct osub toplsubseg, toprsubseg;
  vertex leftvertex, rightvertex, botvertex;
  vertex farvertex;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  /* Identify the vertices of the quadrilateral. */
  org(*flipedge, rightvertex);
  dest(*flipedge, leftvertex);
  apex(*flipedge, botvertex);
  sym(*flipedge, top);
#ifdef SELF_CHECK
  if (top.tri == m->dummytri) {
    printf("Internal error in flip():  Attempt to flip on boundary.\n");
    lnextself(*flipedge);
    return;
  }
  if (m->checksegments) {
    tspivot(*flipedge, toplsubseg);
    if (toplsubseg.ss != m->dummysub) {
      printf("Internal error in flip():  Attempt to flip a segment.\n");
      lnextself(*flipedge);
      return;
    }
  }
#endif /* SELF_CHECK */
  apex(top, farvertex);

  /* Identify the casing of the quadrilateral. */
  lprev(top, topleft);
  sym(topleft, toplcasing);
  lnext(top, topright);
  sym(topright, toprcasing);
  lnext(*flipedge, botleft);
  sym(botleft, botlcasing);
  lprev(*flipedge, botright);
  sym(botright, botrcasing);
  /* Rotate the quadrilateral one-quarter turn counterclockwise. */
  bond(topleft, botlcasing);
  bond(botleft, botrcasing);
  bond(botright, toprcasing);
  bond(topright, toplcasing);

  if (m->checksegments) {
    /* Check for subsegments and rebond them to the quadrilateral. */
    tspivot(topleft, toplsubseg);
    tspivot(botleft, botlsubseg);
    tspivot(botright, botrsubseg);
    tspivot(topright, toprsubseg);
    if (toplsubseg.ss == m->dummysub) {
      tsdissolve(topright);
    } else {
      tsbond(topright, toplsubseg);
    }
    if (botlsubseg.ss == m->dummysub) {
      tsdissolve(topleft);
    } else {
      tsbond(topleft, botlsubseg);
    }
    if (botrsubseg.ss == m->dummysub) {
      tsdissolve(botleft);
    } else {
      tsbond(botleft, botrsubseg);
    }
    if (toprsubseg.ss == m->dummysub) {
      tsdissolve(botright);
    } else {
      tsbond(botright, toprsubseg);
    }
  }

  /* New vertex assignments for the rotated quadrilateral. */
  setorg(*flipedge, farvertex);
  setdest(*flipedge, botvertex);
  setapex(*flipedge, rightvertex);
  setorg(top, botvertex);
  setdest(top, farvertex);
  setapex(top, leftvertex);
  if (b->verbose > 2) {
    printf("  Edge flip results in left ");
    printtriangle(m, b, &top);
    printf("  and right ");
    printtriangle(m, b, flipedge);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  unflip()   Transform two triangles to two different triangles by         */
/*             flipping an edge clockwise within a quadrilateral.  Reverses  */
/*             the flip() operation so that the data structures representing */
/*             the triangles are back where they were before the flip().     */
/*                                                                           */
/*  Imagine the original triangles, abc and bad, oriented so that the        */
/*  shared edge ab lies in a horizontal plane, with the vertex b on the left */
/*  and the vertex a on the right.  The vertex c lies below the edge, and    */
/*  the vertex d lies above the edge.  The `flipedge' handle holds the edge  */
/*  ab of triangle abc, and is directed left, from vertex a to vertex b.     */
/*                                                                           */
/*  The triangles abc and bad are deleted and replaced by the triangles cdb  */
/*  and dca.  The triangles that represent abc and bad are NOT deallocated;  */
/*  they are reused for cdb and dca, respectively.  Hence, any handles that  */
/*  may have held the original triangles are still valid, although not       */
/*  directed as they were before.                                            */
/*                                                                           */
/*  Upon completion of this routine, the `flipedge' handle holds the edge    */
/*  cd of triangle cdb, and is directed up, from vertex c to vertex d.       */
/*  (Hence, the two triangles have rotated clockwise.)                       */
/*                                                                           */
/*  WARNING:  This transformation is geometrically valid only if the         */
/*  quadrilateral adbc is convex.  Furthermore, this transformation is       */
/*  valid only if there is not a subsegment between the triangles abc and    */
/*  bad.  This routine does not check either of these preconditions, and     */
/*  it is the responsibility of the calling routine to ensure that they are  */
/*  met.  If they are not, the streets shall be filled with wailing and      */
/*  gnashing of teeth.                                                       */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void unflip(struct mesh *m, struct behavior *b, struct otri *flipedge)
#else /* not ANSI_DECLARATORS */
void unflip(m, b, flipedge)
struct mesh *m;
struct behavior *b;
struct otri *flipedge;                    /* Handle for the triangle abc. */
#endif /* not ANSI_DECLARATORS */

{
  struct otri botleft, botright;
  struct otri topleft, topright;
  struct otri top;
  struct otri botlcasing, botrcasing;
  struct otri toplcasing, toprcasing;
  struct osub botlsubseg, botrsubseg;
  struct osub toplsubseg, toprsubseg;
  vertex leftvertex, rightvertex, botvertex;
  vertex farvertex;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  /* Identify the vertices of the quadrilateral. */
  org(*flipedge, rightvertex);
  dest(*flipedge, leftvertex);
  apex(*flipedge, botvertex);
  sym(*flipedge, top);
#ifdef SELF_CHECK
  if (top.tri == m->dummytri) {
    printf("Internal error in unflip():  Attempt to flip on boundary.\n");
    lnextself(*flipedge);
    return;
  }
  if (m->checksegments) {
    tspivot(*flipedge, toplsubseg);
    if (toplsubseg.ss != m->dummysub) {
      printf("Internal error in unflip():  Attempt to flip a subsegment.\n");
      lnextself(*flipedge);
      return;
    }
  }
#endif /* SELF_CHECK */
  apex(top, farvertex);

  /* Identify the casing of the quadrilateral. */
  lprev(top, topleft);
  sym(topleft, toplcasing);
  lnext(top, topright);
  sym(topright, toprcasing);
  lnext(*flipedge, botleft);
  sym(botleft, botlcasing);
  lprev(*flipedge, botright);
  sym(botright, botrcasing);
  /* Rotate the quadrilateral one-quarter turn clockwise. */
  bond(topleft, toprcasing);
  bond(botleft, toplcasing);
  bond(botright, botlcasing);
  bond(topright, botrcasing);

  if (m->checksegments) {
    /* Check for subsegments and rebond them to the quadrilateral. */
    tspivot(topleft, toplsubseg);
    tspivot(botleft, botlsubseg);
    tspivot(botright, botrsubseg);
    tspivot(topright, toprsubseg);
    if (toplsubseg.ss == m->dummysub) {
      tsdissolve(botleft);
    } else {
      tsbond(botleft, toplsubseg);
    }
    if (botlsubseg.ss == m->dummysub) {
      tsdissolve(botright);
    } else {
      tsbond(botright, botlsubseg);
    }
    if (botrsubseg.ss == m->dummysub) {
      tsdissolve(topright);
    } else {
      tsbond(topright, botrsubseg);
    }
    if (toprsubseg.ss == m->dummysub) {
      tsdissolve(topleft);
    } else {
      tsbond(topleft, toprsubseg);
    }
  }

  /* New vertex assignments for the rotated quadrilateral. */
  setorg(*flipedge, botvertex);
  setdest(*flipedge, farvertex);
  setapex(*flipedge, leftvertex);
  setorg(top, farvertex);
  setdest(top, botvertex);
  setapex(top, rightvertex);
  if (b->verbose > 2) {
    printf("  Edge unflip results in left ");
    printtriangle(m, b, flipedge);
    printf("  and right ");
    printtriangle(m, b, &top);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  insertvertex()   Insert a vertex into a Delaunay triangulation,          */
/*                   performing flips as necessary to maintain the Delaunay  */
/*                   property.                                               */
/*                                                                           */
/*  The point `insertvertex' is located.  If `searchtri.tri' is not NULL,    */
/*  the search for the containing triangle begins from `searchtri'.  If      */
/*  `searchtri.tri' is NULL, a full point location procedure is called.      */
/*  If `insertvertex' is found inside a triangle, the triangle is split into */
/*  three; if `insertvertex' lies on an edge, the edge is split in two,      */
/*  thereby splitting the two adjacent triangles into four.  Edge flips are  */
/*  used to restore the Delaunay property.  If `insertvertex' lies on an     */
/*  existing vertex, no action is taken, and the value DUPLICATEVERTEX is    */
/*  returned.  On return, `searchtri' is set to a handle whose origin is the */
/*  existing vertex.                                                         */
/*                                                                           */
/*  Normally, the parameter `splitseg' is set to NULL, implying that no      */
/*  subsegment should be split.  In this case, if `insertvertex' is found to */
/*  lie on a segment, no action is taken, and the value VIOLATINGVERTEX is   */
/*  returned.  On return, `searchtri' is set to a handle whose primary edge  */
/*  is the violated subsegment.                                              */
/*                                                                           */
/*  If the calling routine wishes to split a subsegment by inserting a       */
/*  vertex in it, the parameter `splitseg' should be that subsegment.  In    */
/*  this case, `searchtri' MUST be the triangle handle reached by pivoting   */
/*  from that subsegment; no point location is done.                         */
/*                                                                           */
/*  `segmentflaws' and `triflaws' are flags that indicate whether or not     */
/*  there should be checks for the creation of encroached subsegments or bad */
/*  quality triangles.  If a newly inserted vertex encroaches upon           */
/*  subsegments, these subsegments are added to the list of subsegments to   */
/*  be split if `segmentflaws' is set.  If bad triangles are created, these  */
/*  are added to the queue if `triflaws' is set.                             */
/*                                                                           */
/*  If a duplicate vertex or violated segment does not prevent the vertex    */
/*  from being inserted, the return value will be ENCROACHINGVERTEX if the   */
/*  vertex encroaches upon a subsegment (and checking is enabled), or        */
/*  SUCCESSFULVERTEX otherwise.  In either case, `searchtri' is set to a     */
/*  handle whose origin is the newly inserted vertex.                        */
/*                                                                           */
/*  insertvertex() does not use flip() for reasons of speed; some            */
/*  information can be reused from edge flip to edge flip, like the          */
/*  locations of subsegments.                                                */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
enum insertvertexresult insertvertex(struct mesh *m, struct behavior *b,
                                     vertex newvertex, struct otri *searchtri,
                                     struct osub *splitseg,
                                     int segmentflaws, int triflaws)
#else /* not ANSI_DECLARATORS */
enum insertvertexresult insertvertex(m, b, newvertex, searchtri, splitseg,
                                     segmentflaws, triflaws)
struct mesh *m;
struct behavior *b;
vertex newvertex;
struct otri *searchtri;
struct osub *splitseg;
int segmentflaws;
int triflaws;
#endif /* not ANSI_DECLARATORS */

{
  struct otri horiz;
  struct otri top;
  struct otri botleft, botright;
  struct otri topleft, topright;
  struct otri newbotleft, newbotright;
  struct otri newtopright;
  struct otri botlcasing, botrcasing;
  struct otri toplcasing, toprcasing;
  struct otri testtri;
  struct osub botlsubseg, botrsubseg;
  struct osub toplsubseg, toprsubseg;
  struct osub brokensubseg;
  struct osub checksubseg;
  struct osub rightsubseg;
  struct osub newsubseg;
  struct badsubseg *encroached;
  struct flipstacker *newflip;
  vertex first;
  vertex leftvertex, rightvertex, botvertex, topvertex, farvertex;
  vertex segmentorg, segmentdest;
  REAL attrib;
  REAL area;
  enum insertvertexresult success;
  enum locateresult intersect;
  int doflip;
  int mirrorflag;
  int enq;
  int i;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;         /* Temporary variable used by spivot() and tspivot(). */

  if (b->verbose > 1) {
    printf("  Inserting (%.12g, %.12g).\n", newvertex[0], newvertex[1]);
  }

  if (splitseg == (struct osub *) NULL) {
    /* Find the location of the vertex to be inserted.  Check if a good */
    /*   starting triangle has already been provided by the caller.     */
    if (searchtri->tri == m->dummytri) {
      /* Find a boundary triangle. */
      horiz.tri = m->dummytri;
      horiz.orient = 0;
      symself(horiz);
      /* Search for a triangle containing `newvertex'. */
      intersect = locate(m, b, newvertex, &horiz);
    } else {
      /* Start searching from the triangle provided by the caller. */
      otricopy(*searchtri, horiz);
      intersect = preciselocate(m, b, newvertex, &horiz, 1);
    }
  } else {
    /* The calling routine provides the subsegment in which */
    /*   the vertex is inserted.                             */
    otricopy(*searchtri, horiz);
    intersect = ONEDGE;
  }

  if (intersect == ONVERTEX) {
    /* There's already a vertex there.  Return in `searchtri' a triangle */
    /*   whose origin is the existing vertex.                            */
    otricopy(horiz, *searchtri);
    otricopy(horiz, m->recenttri);
    return DUPLICATEVERTEX;
  }
  if ((intersect == ONEDGE) || (intersect == OUTSIDE)) {
    /* The vertex falls on an edge or boundary. */
    if (m->checksegments && (splitseg == (struct osub *) NULL)) {
      /* Check whether the vertex falls on a subsegment. */
      tspivot(horiz, brokensubseg);
      if (brokensubseg.ss != m->dummysub) {
        /* The vertex falls on a subsegment, and hence will not be inserted. */
        if (segmentflaws) {
          enq = b->nobisect != 2;
          if (enq && (b->nobisect == 1)) {
            /* This subsegment may be split only if it is an */
            /*   internal boundary.                          */
            sym(horiz, testtri);
            enq = testtri.tri != m->dummytri;
          }
          if (enq) {
            /* Add the subsegment to the list of encroached subsegments. */
            encroached = (struct badsubseg *) poolalloc(&m->badsubsegs);
            encroached->encsubseg = sencode(brokensubseg);
            sorg(brokensubseg, encroached->subsegorg);
            sdest(brokensubseg, encroached->subsegdest);
            if (b->verbose > 2) {
              printf(
          "  Queueing encroached subsegment (%.12g, %.12g) (%.12g, %.12g).\n",
                     encroached->subsegorg[0], encroached->subsegorg[1],
                     encroached->subsegdest[0], encroached->subsegdest[1]);
            }
          }
        }
        /* Return a handle whose primary edge contains the vertex, */
        /*   which has not been inserted.                          */
        otricopy(horiz, *searchtri);
        otricopy(horiz, m->recenttri);
        return VIOLATINGVERTEX;
      }
    }

    /* Insert the vertex on an edge, dividing one triangle into two (if */
    /*   the edge lies on a boundary) or two triangles into four.       */
    lprev(horiz, botright);
    sym(botright, botrcasing);
    sym(horiz, topright);
    /* Is there a second triangle?  (Or does this edge lie on a boundary?) */
    mirrorflag = topright.tri != m->dummytri;
    if (mirrorflag) {
      lnextself(topright);
      sym(topright, toprcasing);
      maketriangle(m, b, &newtopright);
    } else {
      /* Splitting a boundary edge increases the number of boundary edges. */
      m->hullsize++;
    }
    maketriangle(m, b, &newbotright);

    /* Set the vertices of changed and new triangles. */
    org(horiz, rightvertex);
    dest(horiz, leftvertex);
    apex(horiz, botvertex);
    setorg(newbotright, botvertex);
    setdest(newbotright, rightvertex);
    setapex(newbotright, newvertex);
    setorg(horiz, newvertex);
    for (i = 0; i < m->eextras; i++) {
      /* Set the element attributes of a new triangle. */
      setelemattribute(newbotright, i, elemattribute(botright, i));
    }
    if (b->vararea) {
      /* Set the area constraint of a new triangle. */
      setareabound(newbotright, areabound(botright));
    }
    if (mirrorflag) {
      dest(topright, topvertex);
      setorg(newtopright, rightvertex);
      setdest(newtopright, topvertex);
      setapex(newtopright, newvertex);
      setorg(topright, newvertex);
      for (i = 0; i < m->eextras; i++) {
        /* Set the element attributes of another new triangle. */
        setelemattribute(newtopright, i, elemattribute(topright, i));
      }
      if (b->vararea) {
        /* Set the area constraint of another new triangle. */
        setareabound(newtopright, areabound(topright));
      }
    }

    /* There may be subsegments that need to be bonded */
    /*   to the new triangle(s).                       */
    if (m->checksegments) {
      tspivot(botright, botrsubseg);
      if (botrsubseg.ss != m->dummysub) {
        tsdissolve(botright);
        tsbond(newbotright, botrsubseg);
      }
      if (mirrorflag) {
        tspivot(topright, toprsubseg);
        if (toprsubseg.ss != m->dummysub) {
          tsdissolve(topright);
          tsbond(newtopright, toprsubseg);
        }
      }
    }

    /* Bond the new triangle(s) to the surrounding triangles. */
    bond(newbotright, botrcasing);
    lprevself(newbotright);
    bond(newbotright, botright);
    lprevself(newbotright);
    if (mirrorflag) {
      bond(newtopright, toprcasing);
      lnextself(newtopright);
      bond(newtopright, topright);
      lnextself(newtopright);
      bond(newtopright, newbotright);
    }

    if (splitseg != (struct osub *) NULL) {
      /* Split the subsegment into two. */
      setsdest(*splitseg, newvertex);
      segorg(*splitseg, segmentorg);
      segdest(*splitseg, segmentdest);
      ssymself(*splitseg);
      spivot(*splitseg, rightsubseg);
      insertsubseg(m, b, &newbotright, mark(*splitseg));
      tspivot(newbotright, newsubseg);
      setsegorg(newsubseg, segmentorg);
      setsegdest(newsubseg, segmentdest);
      sbond(*splitseg, newsubseg);
      ssymself(newsubseg);
      sbond(newsubseg, rightsubseg);
      ssymself(*splitseg);
      /* Transfer the subsegment's boundary marker to the vertex */
      /*   if required.                                          */
      if (vertexmark(newvertex) == 0) {
        setvertexmark(newvertex, mark(*splitseg));
      }
    }

    if (m->checkquality) {
      poolrestart(&m->flipstackers);
      m->lastflip = (struct flipstacker *) poolalloc(&m->flipstackers);
      m->lastflip->flippedtri = encode(horiz);
      m->lastflip->prevflip = (struct flipstacker *) &insertvertex;
    }

#ifdef SELF_CHECK
    if (counterclockwise(m, b, rightvertex, leftvertex, botvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf(
            "  Clockwise triangle prior to edge vertex insertion (bottom).\n");
    }
    if (mirrorflag) {
      if (counterclockwise(m, b, leftvertex, rightvertex, topvertex) < 0.0) {
        printf("Internal error in insertvertex():\n");
        printf("  Clockwise triangle prior to edge vertex insertion (top).\n");
      }
      if (counterclockwise(m, b, rightvertex, topvertex, newvertex) < 0.0) {
        printf("Internal error in insertvertex():\n");
        printf(
            "  Clockwise triangle after edge vertex insertion (top right).\n");
      }
      if (counterclockwise(m, b, topvertex, leftvertex, newvertex) < 0.0) {
        printf("Internal error in insertvertex():\n");
        printf(
            "  Clockwise triangle after edge vertex insertion (top left).\n");
      }
    }
    if (counterclockwise(m, b, leftvertex, botvertex, newvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf(
          "  Clockwise triangle after edge vertex insertion (bottom left).\n");
    }
    if (counterclockwise(m, b, botvertex, rightvertex, newvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf(
        "  Clockwise triangle after edge vertex insertion (bottom right).\n");
    }
#endif /* SELF_CHECK */
    if (b->verbose > 2) {
      printf("  Updating bottom left ");
      printtriangle(m, b, &botright);
      if (mirrorflag) {
        printf("  Updating top left ");
        printtriangle(m, b, &topright);
        printf("  Creating top right ");
        printtriangle(m, b, &newtopright);
      }
      printf("  Creating bottom right ");
      printtriangle(m, b, &newbotright);
    }

    /* Position `horiz' on the first edge to check for */
    /*   the Delaunay property.                        */
    lnextself(horiz);
  } else {
    /* Insert the vertex in a triangle, splitting it into three. */
    lnext(horiz, botleft);
    lprev(horiz, botright);
    sym(botleft, botlcasing);
    sym(botright, botrcasing);
    maketriangle(m, b, &newbotleft);
    maketriangle(m, b, &newbotright);

    /* Set the vertices of changed and new triangles. */
    org(horiz, rightvertex);
    dest(horiz, leftvertex);
    apex(horiz, botvertex);
    setorg(newbotleft, leftvertex);
    setdest(newbotleft, botvertex);
    setapex(newbotleft, newvertex);
    setorg(newbotright, botvertex);
    setdest(newbotright, rightvertex);
    setapex(newbotright, newvertex);
    setapex(horiz, newvertex);
    for (i = 0; i < m->eextras; i++) {
      /* Set the element attributes of the new triangles. */
      attrib = elemattribute(horiz, i);
      setelemattribute(newbotleft, i, attrib);
      setelemattribute(newbotright, i, attrib);
    }
    if (b->vararea) {
      /* Set the area constraint of the new triangles. */
      area = areabound(horiz);
      setareabound(newbotleft, area);
      setareabound(newbotright, area);
    }

    /* There may be subsegments that need to be bonded */
    /*   to the new triangles.                         */
    if (m->checksegments) {
      tspivot(botleft, botlsubseg);
      if (botlsubseg.ss != m->dummysub) {
        tsdissolve(botleft);
        tsbond(newbotleft, botlsubseg);
      }
      tspivot(botright, botrsubseg);
      if (botrsubseg.ss != m->dummysub) {
        tsdissolve(botright);
        tsbond(newbotright, botrsubseg);
      }
    }

    /* Bond the new triangles to the surrounding triangles. */
    bond(newbotleft, botlcasing);
    bond(newbotright, botrcasing);
    lnextself(newbotleft);
    lprevself(newbotright);
    bond(newbotleft, newbotright);
    lnextself(newbotleft);
    bond(botleft, newbotleft);
    lprevself(newbotright);
    bond(botright, newbotright);

    if (m->checkquality) {
      poolrestart(&m->flipstackers);
      m->lastflip = (struct flipstacker *) poolalloc(&m->flipstackers);
      m->lastflip->flippedtri = encode(horiz);
      m->lastflip->prevflip = (struct flipstacker *) NULL;
    }

#ifdef SELF_CHECK
    if (counterclockwise(m, b, rightvertex, leftvertex, botvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf("  Clockwise triangle prior to vertex insertion.\n");
    }
    if (counterclockwise(m, b, rightvertex, leftvertex, newvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf("  Clockwise triangle after vertex insertion (top).\n");
    }
    if (counterclockwise(m, b, leftvertex, botvertex, newvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf("  Clockwise triangle after vertex insertion (left).\n");
    }
    if (counterclockwise(m, b, botvertex, rightvertex, newvertex) < 0.0) {
      printf("Internal error in insertvertex():\n");
      printf("  Clockwise triangle after vertex insertion (right).\n");
    }
#endif /* SELF_CHECK */
    if (b->verbose > 2) {
      printf("  Updating top ");
      printtriangle(m, b, &horiz);
      printf("  Creating left ");
      printtriangle(m, b, &newbotleft);
      printf("  Creating right ");
      printtriangle(m, b, &newbotright);
    }
  }

  /* The insertion is successful by default, unless an encroached */
  /*   subsegment is found.                                       */
  success = SUCCESSFULVERTEX;
  /* Circle around the newly inserted vertex, checking each edge opposite */
  /*   it for the Delaunay property.  Non-Delaunay edges are flipped.     */
  /*   `horiz' is always the edge being checked.  `first' marks where to  */
  /*   stop circling.                                                     */
  org(horiz, first);
  rightvertex = first;
  dest(horiz, leftvertex);
  /* Circle until finished. */
  while (1) {
    /* By default, the edge will be flipped. */
    doflip = 1;

    if (m->checksegments) {
      /* Check for a subsegment, which cannot be flipped. */
      tspivot(horiz, checksubseg);
      if (checksubseg.ss != m->dummysub) {
        /* The edge is a subsegment and cannot be flipped. */
        doflip = 0;
#ifndef CDT_ONLY
        if (segmentflaws) {
          /* Does the new vertex encroach upon this subsegment? */
          if (checkseg4encroach(m, b, &checksubseg)) {
            success = ENCROACHINGVERTEX;
          }
        }
#endif /* not CDT_ONLY */
      }
    }

    if (doflip) {
      /* Check if the edge is a boundary edge. */
      sym(horiz, top);
      if (top.tri == m->dummytri) {
        /* The edge is a boundary edge and cannot be flipped. */
        doflip = 0;
      } else {
        /* Find the vertex on the other side of the edge. */
        apex(top, farvertex);
        /* In the incremental Delaunay triangulation algorithm, any of      */
        /*   `leftvertex', `rightvertex', and `farvertex' could be vertices */
        /*   of the triangular bounding box.  These vertices must be        */
        /*   treated as if they are infinitely distant, even though their   */
        /*   "coordinates" are not.                                         */
        if ((leftvertex == m->infvertex1) || (leftvertex == m->infvertex2) ||
            (leftvertex == m->infvertex3)) {
          /* `leftvertex' is infinitely distant.  Check the convexity of  */
          /*   the boundary of the triangulation.  'farvertex' might be   */
          /*   infinite as well, but trust me, this same condition should */
          /*   be applied.                                                */
          doflip = counterclockwise(m, b, newvertex, rightvertex, farvertex)
                   > 0.0;
        } else if ((rightvertex == m->infvertex1) ||
                   (rightvertex == m->infvertex2) ||
                   (rightvertex == m->infvertex3)) {
          /* `rightvertex' is infinitely distant.  Check the convexity of */
          /*   the boundary of the triangulation.  'farvertex' might be   */
          /*   infinite as well, but trust me, this same condition should */
          /*   be applied.                                                */
          doflip = counterclockwise(m, b, farvertex, leftvertex, newvertex)
                   > 0.0;
        } else if ((farvertex == m->infvertex1) ||
                   (farvertex == m->infvertex2) ||
                   (farvertex == m->infvertex3)) {
          /* `farvertex' is infinitely distant and cannot be inside */
          /*   the circumcircle of the triangle `horiz'.            */
          doflip = 0;
        } else {
          /* Test whether the edge is locally Delaunay. */
          doflip = incircle(m, b, leftvertex, newvertex, rightvertex,
                            farvertex) > 0.0;
        }
        if (doflip) {
          /* We made it!  Flip the edge `horiz' by rotating its containing */
          /*   quadrilateral (the two triangles adjacent to `horiz').      */
          /* Identify the casing of the quadrilateral. */
          lprev(top, topleft);
          sym(topleft, toplcasing);
          lnext(top, topright);
          sym(topright, toprcasing);
          lnext(horiz, botleft);
          sym(botleft, botlcasing);
          lprev(horiz, botright);
          sym(botright, botrcasing);
          /* Rotate the quadrilateral one-quarter turn counterclockwise. */
          bond(topleft, botlcasing);
          bond(botleft, botrcasing);
          bond(botright, toprcasing);
          bond(topright, toplcasing);
          if (m->checksegments) {
            /* Check for subsegments and rebond them to the quadrilateral. */
            tspivot(topleft, toplsubseg);
            tspivot(botleft, botlsubseg);
            tspivot(botright, botrsubseg);
            tspivot(topright, toprsubseg);
            if (toplsubseg.ss == m->dummysub) {
              tsdissolve(topright);
            } else {
              tsbond(topright, toplsubseg);
            }
            if (botlsubseg.ss == m->dummysub) {
              tsdissolve(topleft);
            } else {
              tsbond(topleft, botlsubseg);
            }
            if (botrsubseg.ss == m->dummysub) {
              tsdissolve(botleft);
            } else {
              tsbond(botleft, botrsubseg);
            }
            if (toprsubseg.ss == m->dummysub) {
              tsdissolve(botright);
            } else {
              tsbond(botright, toprsubseg);
            }
          }
          /* New vertex assignments for the rotated quadrilateral. */
          setorg(horiz, farvertex);
          setdest(horiz, newvertex);
          setapex(horiz, rightvertex);
          setorg(top, newvertex);
          setdest(top, farvertex);
          setapex(top, leftvertex);
          for (i = 0; i < m->eextras; i++) {
            /* Take the average of the two triangles' attributes. */
            attrib = 0.5 * (elemattribute(top, i) + elemattribute(horiz, i));
            setelemattribute(top, i, attrib);
            setelemattribute(horiz, i, attrib);
          }
          if (b->vararea) {
            if ((areabound(top) <= 0.0) || (areabound(horiz) <= 0.0)) {
              area = -1.0;
            } else {
              /* Take the average of the two triangles' area constraints.    */
              /*   This prevents small area constraints from migrating a     */
              /*   long, long way from their original location due to flips. */
              area = 0.5 * (areabound(top) + areabound(horiz));
            }
            setareabound(top, area);
            setareabound(horiz, area);
          }

          if (m->checkquality) {
            newflip = (struct flipstacker *) poolalloc(&m->flipstackers);
            newflip->flippedtri = encode(horiz);
            newflip->prevflip = m->lastflip;
            m->lastflip = newflip;
          }

#ifdef SELF_CHECK
          if (newvertex != (vertex) NULL) {
            if (counterclockwise(m, b, leftvertex, newvertex, rightvertex) <
                0.0) {
              printf("Internal error in insertvertex():\n");
              printf("  Clockwise triangle prior to edge flip (bottom).\n");
            }
            /* The following test has been removed because constrainededge() */
            /*   sometimes generates inverted triangles that insertvertex()  */
            /*   removes.                                                    */
/*
            if (counterclockwise(m, b, rightvertex, farvertex, leftvertex) <
                0.0) {
              printf("Internal error in insertvertex():\n");
              printf("  Clockwise triangle prior to edge flip (top).\n");
            }
*/
            if (counterclockwise(m, b, farvertex, leftvertex, newvertex) <
                0.0) {
              printf("Internal error in insertvertex():\n");
              printf("  Clockwise triangle after edge flip (left).\n");
            }
            if (counterclockwise(m, b, newvertex, rightvertex, farvertex) <
                0.0) {
              printf("Internal error in insertvertex():\n");
              printf("  Clockwise triangle after edge flip (right).\n");
            }
          }
#endif /* SELF_CHECK */
          if (b->verbose > 2) {
            printf("  Edge flip results in left ");
            lnextself(topleft);
            printtriangle(m, b, &topleft);
            printf("  and right ");
            printtriangle(m, b, &horiz);
          }
          /* On the next iterations, consider the two edges that were  */
          /*   exposed (this is, are now visible to the newly inserted */
          /*   vertex) by the edge flip.                               */
          lprevself(horiz);
          leftvertex = farvertex;
        }
      }
    }
    if (!doflip) {
      /* The handle `horiz' is accepted as locally Delaunay. */
#ifndef CDT_ONLY
      if (triflaws) {
        /* Check the triangle `horiz' for quality. */
        testtriangle(m, b, &horiz);
      }
#endif /* not CDT_ONLY */
      /* Look for the next edge around the newly inserted vertex. */
      lnextself(horiz);
      sym(horiz, testtri);
      /* Check for finishing a complete revolution about the new vertex, or */
      /*   falling outside  of the triangulation.  The latter will happen   */
      /*   when a vertex is inserted at a boundary.                         */
      if ((leftvertex == first) || (testtri.tri == m->dummytri)) {
        /* We're done.  Return a triangle whose origin is the new vertex. */
        lnext(horiz, *searchtri);
        lnext(horiz, m->recenttri);
        return success;
      }
      /* Finish finding the next edge around the newly inserted vertex. */
      lnext(testtri, horiz);
      rightvertex = leftvertex;
      dest(horiz, leftvertex);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  triangulatepolygon()   Find the Delaunay triangulation of a polygon that */
/*                         has a certain "nice" shape.  This includes the    */
/*                         polygons that result from deletion of a vertex or */
/*                         insertion of a segment.                           */
/*                                                                           */
/*  This is a conceptually difficult routine.  The starting assumption is    */
/*  that we have a polygon with n sides.  n - 1 of these sides are currently */
/*  represented as edges in the mesh.  One side, called the "base", need not */
/*  be.                                                                      */
/*                                                                           */
/*  Inside the polygon is a structure I call a "fan", consisting of n - 1    */
/*  triangles that share a common origin.  For each of these triangles, the  */
/*  edge opposite the origin is one of the sides of the polygon.  The        */
/*  primary edge of each triangle is the edge directed from the origin to    */
/*  the destination; note that this is not the same edge that is a side of   */
/*  the polygon.  `firstedge' is the primary edge of the first triangle.     */
/*  From there, the triangles follow in counterclockwise order about the     */
/*  polygon, until `lastedge', the primary edge of the last triangle.        */
/*  `firstedge' and `lastedge' are probably connected to other triangles     */
/*  beyond the extremes of the fan, but their identity is not important, as  */
/*  long as the fan remains connected to them.                               */
/*                                                                           */
/*  Imagine the polygon oriented so that its base is at the bottom.  This    */
/*  puts `firstedge' on the far right, and `lastedge' on the far left.       */
/*  The right vertex of the base is the destination of `firstedge', and the  */
/*  left vertex of the base is the apex of `lastedge'.                       */
/*                                                                           */
/*  The challenge now is to find the right sequence of edge flips to         */
/*  transform the fan into a Delaunay triangulation of the polygon.  Each    */
/*  edge flip effectively removes one triangle from the fan, committing it   */
/*  to the polygon.  The resulting polygon has one fewer edge.  If `doflip'  */
/*  is set, the final flip will be performed, resulting in a fan of one      */
/*  (useless?) triangle.  If `doflip' is not set, the final flip is not      */
/*  performed, resulting in a fan of two triangles, and an unfinished        */
/*  triangular polygon that is not yet filled out with a single triangle.    */
/*  On completion of the routine, `lastedge' is the last remaining triangle, */
/*  or the leftmost of the last two.                                         */
/*                                                                           */
/*  Although the flips are performed in the order described above, the       */
/*  decisions about what flips to perform are made in precisely the reverse  */
/*  order.  The recursive triangulatepolygon() procedure makes a decision,   */
/*  uses up to two recursive calls to triangulate the "subproblems"          */
/*  (polygons with fewer edges), and then performs an edge flip.             */
/*                                                                           */
/*  The "decision" it makes is which vertex of the polygon should be         */
/*  connected to the base.  This decision is made by testing every possible  */
/*  vertex.  Once the best vertex is found, the two edges that connect this  */
/*  vertex to the base become the bases for two smaller polygons.  These     */
/*  are triangulated recursively.  Unfortunately, this approach can take     */
/*  O(n^2) time not only in the worst case, but in many common cases.  It's  */
/*  rarely a big deal for vertex deletion, where n is rarely larger than     */
/*  ten, but it could be a big deal for segment insertion, especially if     */
/*  there's a lot of long segments that each cut many triangles.  I ought to */
/*  code a faster algorithm some day.                                        */
/*                                                                           */
/*  The `edgecount' parameter is the number of sides of the polygon,         */
/*  including its base.  `triflaws' is a flag that determines whether the    */
/*  new triangles should be tested for quality, and enqueued if they are     */
/*  bad.                                                                     */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void triangulatepolygon(struct mesh *m, struct behavior *b,
                        struct otri *firstedge, struct otri *lastedge,
                        int edgecount, int doflip, int triflaws)
#else /* not ANSI_DECLARATORS */
void triangulatepolygon(m, b, firstedge, lastedge, edgecount, doflip, triflaws)
struct mesh *m;
struct behavior *b;
struct otri *firstedge;
struct otri *lastedge;
int edgecount;
int doflip;
int triflaws;
#endif /* not ANSI_DECLARATORS */

{
  struct otri testtri;
  struct otri besttri;
  struct otri tempedge;
  vertex leftbasevertex, rightbasevertex;
  vertex testvertex;
  vertex bestvertex;
  int bestnumber;
  int i;
  triangle ptr;   /* Temporary variable used by sym(), onext(), and oprev(). */

  /* Identify the base vertices. */
  apex(*lastedge, leftbasevertex);
  dest(*firstedge, rightbasevertex);
  if (b->verbose > 2) {
    printf("  Triangulating interior polygon at edge\n");
    printf("    (%.12g, %.12g) (%.12g, %.12g)\n", leftbasevertex[0],
           leftbasevertex[1], rightbasevertex[0], rightbasevertex[1]);
  }
  /* Find the best vertex to connect the base to. */
  onext(*firstedge, besttri);
  dest(besttri, bestvertex);
  otricopy(besttri, testtri);
  bestnumber = 1;
  for (i = 2; i <= edgecount - 2; i++) {
    onextself(testtri);
    dest(testtri, testvertex);
    /* Is this a better vertex? */
    if (incircle(m, b, leftbasevertex, rightbasevertex, bestvertex,
                 testvertex) > 0.0) {
      otricopy(testtri, besttri);
      bestvertex = testvertex;
      bestnumber = i;
    }
  }
  if (b->verbose > 2) {
    printf("    Connecting edge to (%.12g, %.12g)\n", bestvertex[0],
           bestvertex[1]);
  }
  if (bestnumber > 1) {
    /* Recursively triangulate the smaller polygon on the right. */
    oprev(besttri, tempedge);
    triangulatepolygon(m, b, firstedge, &tempedge, bestnumber + 1, 1,
                       triflaws);
  }
  if (bestnumber < edgecount - 2) {
    /* Recursively triangulate the smaller polygon on the left. */
    sym(besttri, tempedge);
    triangulatepolygon(m, b, &besttri, lastedge, edgecount - bestnumber, 1,
                       triflaws);
    /* Find `besttri' again; it may have been lost to edge flips. */
    sym(tempedge, besttri);
  }
  if (doflip) {
    /* Do one final edge flip. */
    flip(m, b, &besttri);
#ifndef CDT_ONLY
    if (triflaws) {
      /* Check the quality of the newly committed triangle. */
      sym(besttri, testtri);
      testtriangle(m, b, &testtri);
    }
#endif /* not CDT_ONLY */
  }
  /* Return the base triangle. */
  otricopy(besttri, *lastedge);
}

/*****************************************************************************/
/*                                                                           */
/*  deletevertex()   Delete a vertex from a Delaunay triangulation, ensuring */
/*                   that the triangulation remains Delaunay.                */
/*                                                                           */
/*  The origin of `deltri' is deleted.  The union of the triangles adjacent  */
/*  to this vertex is a polygon, for which the Delaunay triangulation is     */
/*  found.  Two triangles are removed from the mesh.                         */
/*                                                                           */
/*  Only interior vertices that do not lie on segments or boundaries may be  */
/*  deleted.                                                                 */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void deletevertex(struct mesh *m, struct behavior *b, struct otri *deltri)
#else /* not ANSI_DECLARATORS */
void deletevertex(m, b, deltri)
struct mesh *m;
struct behavior *b;
struct otri *deltri;
#endif /* not ANSI_DECLARATORS */

{
  struct otri countingtri;
  struct otri firstedge, lastedge;
  struct otri deltriright;
  struct otri lefttri, righttri;
  struct otri leftcasing, rightcasing;
  struct osub leftsubseg, rightsubseg;
  vertex delvertex;
  vertex neworg;
  int edgecount;
  triangle ptr;   /* Temporary variable used by sym(), onext(), and oprev(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  org(*deltri, delvertex);
  if (b->verbose > 1) {
    printf("  Deleting (%.12g, %.12g).\n", delvertex[0], delvertex[1]);
  }
  vertexdealloc(m, delvertex);

  /* Count the degree of the vertex being deleted. */
  onext(*deltri, countingtri);
  edgecount = 1;
  while (!otriequal(*deltri, countingtri)) {
#ifdef SELF_CHECK
    if (countingtri.tri == m->dummytri) {
      printf("Internal error in deletevertex():\n");
      printf("  Attempt to delete boundary vertex.\n");
      internalerror();
    }
#endif /* SELF_CHECK */
    edgecount++;
    onextself(countingtri);
  }

#ifdef SELF_CHECK
  if (edgecount < 3) {
    printf("Internal error in deletevertex():\n  Vertex has degree %d.\n",
           edgecount);
    internalerror();
  }
#endif /* SELF_CHECK */
  if (edgecount > 3) {
    /* Triangulate the polygon defined by the union of all triangles */
    /*   adjacent to the vertex being deleted.  Check the quality of */
    /*   the resulting triangles.                                    */
    onext(*deltri, firstedge);
    oprev(*deltri, lastedge);
    triangulatepolygon(m, b, &firstedge, &lastedge, edgecount, 0,
                       !b->nobisect);
  }
  /* Splice out two triangles. */
  lprev(*deltri, deltriright);
  dnext(*deltri, lefttri);
  sym(lefttri, leftcasing);
  oprev(deltriright, righttri);
  sym(righttri, rightcasing);
  bond(*deltri, leftcasing);
  bond(deltriright, rightcasing);
  tspivot(lefttri, leftsubseg);
  if (leftsubseg.ss != m->dummysub) {
    tsbond(*deltri, leftsubseg);
  }
  tspivot(righttri, rightsubseg);
  if (rightsubseg.ss != m->dummysub) {
    tsbond(deltriright, rightsubseg);
  }

  /* Set the new origin of `deltri' and check its quality. */
  org(lefttri, neworg);
  setorg(*deltri, neworg);
  if (!b->nobisect) {
    testtriangle(m, b, deltri);
  }

  /* Delete the two spliced-out triangles. */
  triangledealloc(m, lefttri.tri);
  triangledealloc(m, righttri.tri);
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  undovertex()   Undo the most recent vertex insertion.                    */
/*                                                                           */
/*  Walks through the list of transformations (flips and a vertex insertion) */
/*  in the reverse of the order in which they were done, and undoes them.    */
/*  The inserted vertex is removed from the triangulation and deallocated.   */
/*  Two triangles (possibly just one) are also deallocated.                  */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void undovertex(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void undovertex(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri fliptri;
  struct otri botleft, botright, topright;
  struct otri botlcasing, botrcasing, toprcasing;
  struct otri gluetri;
  struct osub botlsubseg, botrsubseg, toprsubseg;
  vertex botvertex, rightvertex;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  /* Walk through the list of transformations (flips and a vertex insertion) */
  /*   in the reverse of the order in which they were done, and undo them.   */
  while (m->lastflip != (struct flipstacker *) NULL) {
    /* Find a triangle involved in the last unreversed transformation. */
    decode(m->lastflip->flippedtri, fliptri);

    /* We are reversing one of three transformations:  a trisection of one */
    /*   triangle into three (by inserting a vertex in the triangle), a    */
    /*   bisection of two triangles into four (by inserting a vertex in an */
    /*   edge), or an edge flip.                                           */
    if (m->lastflip->prevflip == (struct flipstacker *) NULL) {
      /* Restore a triangle that was split into three triangles, */
      /*   so it is again one triangle.                          */
      dprev(fliptri, botleft);
      lnextself(botleft);
      onext(fliptri, botright);
      lprevself(botright);
      sym(botleft, botlcasing);
      sym(botright, botrcasing);
      dest(botleft, botvertex);

      setapex(fliptri, botvertex);
      lnextself(fliptri);
      bond(fliptri, botlcasing);
      tspivot(botleft, botlsubseg);
      tsbond(fliptri, botlsubseg);
      lnextself(fliptri);
      bond(fliptri, botrcasing);
      tspivot(botright, botrsubseg);
      tsbond(fliptri, botrsubseg);

      /* Delete the two spliced-out triangles. */
      triangledealloc(m, botleft.tri);
      triangledealloc(m, botright.tri);
    } else if (m->lastflip->prevflip == (struct flipstacker *) &insertvertex) {
      /* Restore two triangles that were split into four triangles, */
      /*   so they are again two triangles.                         */
      lprev(fliptri, gluetri);
      sym(gluetri, botright);
      lnextself(botright);
      sym(botright, botrcasing);
      dest(botright, rightvertex);

      setorg(fliptri, rightvertex);
      bond(gluetri, botrcasing);
      tspivot(botright, botrsubseg);
      tsbond(gluetri, botrsubseg);

      /* Delete the spliced-out triangle. */
      triangledealloc(m, botright.tri);

      sym(fliptri, gluetri);
      if (gluetri.tri != m->dummytri) {
        lnextself(gluetri);
        dnext(gluetri, topright);
        sym(topright, toprcasing);

        setorg(gluetri, rightvertex);
        bond(gluetri, toprcasing);
        tspivot(topright, toprsubseg);
        tsbond(gluetri, toprsubseg);

        /* Delete the spliced-out triangle. */
        triangledealloc(m, topright.tri);
      }

      /* This is the end of the list, sneakily encoded. */
      m->lastflip->prevflip = (struct flipstacker *) NULL;
    } else {
      /* Undo an edge flip. */
      unflip(m, b, &fliptri);
    }

    /* Go on and process the next transformation. */
    m->lastflip = m->lastflip->prevflip;
  }
}

#endif /* not CDT_ONLY */

/**                                                                         **/
/**                                                                         **/
/********* Mesh transformation routines end here                     *********/

/********* Divide-and-conquer Delaunay triangulation begins here     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  The divide-and-conquer bounding box                                      */
/*                                                                           */
/*  I originally implemented the divide-and-conquer and incremental Delaunay */
/*  triangulations using the edge-based data structure presented by Guibas   */
/*  and Stolfi.  Switching to a triangle-based data structure doubled the    */
/*  speed.  However, I had to think of a few extra tricks to maintain the    */
/*  elegance of the original algorithms.                                     */
/*                                                                           */
/*  The "bounding box" used by my variant of the divide-and-conquer          */
/*  algorithm uses one triangle for each edge of the convex hull of the      */
/*  triangulation.  These bounding triangles all share a common apical       */
/*  vertex, which is represented by NULL and which represents nothing.       */
/*  The bounding triangles are linked in a circular fan about this NULL      */
/*  vertex, and the edges on the convex hull of the triangulation appear     */
/*  opposite the NULL vertex.  You might find it easiest to imagine that     */
/*  the NULL vertex is a point in 3D space behind the center of the          */
/*  triangulation, and that the bounding triangles form a sort of cone.      */
/*                                                                           */
/*  This bounding box makes it easy to represent degenerate cases.  For      */
/*  instance, the triangulation of two vertices is a single edge.  This edge */
/*  is represented by two bounding box triangles, one on each "side" of the  */
/*  edge.  These triangles are also linked together in a fan about the NULL  */
/*  vertex.                                                                  */
/*                                                                           */
/*  The bounding box also makes it easy to traverse the convex hull, as the  */
/*  divide-and-conquer algorithm needs to do.                                */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  vertexsort()   Sort an array of vertices by x-coordinate, using the      */
/*                 y-coordinate as a secondary key.                          */
/*                                                                           */
/*  Uses quicksort.  Randomized O(n log n) time.  No, I did not make any of  */
/*  the usual quicksort mistakes.                                            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void vertexsort(vertex *sortarray, int arraysize)
#else /* not ANSI_DECLARATORS */
void vertexsort(sortarray, arraysize)
vertex *sortarray;
int arraysize;
#endif /* not ANSI_DECLARATORS */

{
  int left, right;
  int pivot;
  REAL pivotx, pivoty;
  vertex temp;

  if (arraysize == 2) {
    /* Recursive base case. */
    if ((sortarray[0][0] > sortarray[1][0]) ||
        ((sortarray[0][0] == sortarray[1][0]) &&
         (sortarray[0][1] > sortarray[1][1]))) {
      temp = sortarray[1];
      sortarray[1] = sortarray[0];
      sortarray[0] = temp;
    }
    return;
  }
  /* Choose a random pivot to split the array. */
  pivot = (int) randomnation((unsigned int) arraysize);
  pivotx = sortarray[pivot][0];
  pivoty = sortarray[pivot][1];
  /* Split the array. */
  left = -1;
  right = arraysize;
  while (left < right) {
    /* Search for a vertex whose x-coordinate is too large for the left. */
    do {
      left++;
    } while ((left <= right) && ((sortarray[left][0] < pivotx) ||
                                 ((sortarray[left][0] == pivotx) &&
                                  (sortarray[left][1] < pivoty))));
    /* Search for a vertex whose x-coordinate is too small for the right. */
    do {
      right--;
    } while ((left <= right) && ((sortarray[right][0] > pivotx) ||
                                 ((sortarray[right][0] == pivotx) &&
                                  (sortarray[right][1] > pivoty))));
    if (left < right) {
      /* Swap the left and right vertices. */
      temp = sortarray[left];
      sortarray[left] = sortarray[right];
      sortarray[right] = temp;
    }
  }
  if (left > 1) {
    /* Recursively sort the left subset. */
    vertexsort(sortarray, left);
  }
  if (right < arraysize - 2) {
    /* Recursively sort the right subset. */
    vertexsort(&sortarray[right + 1], arraysize - right - 1);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  vertexmedian()   An order statistic algorithm, almost.  Shuffles an      */
/*                   array of vertices so that the first `median' vertices   */
/*                   occur lexicographically before the remaining vertices.  */
/*                                                                           */
/*  Uses the x-coordinate as the primary key if axis == 0; the y-coordinate  */
/*  if axis == 1.  Very similar to the vertexsort() procedure, but runs in   */
/*  randomized linear time.                                                  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void vertexmedian(vertex *sortarray, int arraysize, int median, int axis)
#else /* not ANSI_DECLARATORS */
void vertexmedian(sortarray, arraysize, median, axis)
vertex *sortarray;
int arraysize;
int median;
int axis;
#endif /* not ANSI_DECLARATORS */

{
  int left, right;
  int pivot;
  REAL pivot1, pivot2;
  vertex temp;

  if (arraysize == 2) {
    /* Recursive base case. */
    if ((sortarray[0][axis] > sortarray[1][axis]) ||
        ((sortarray[0][axis] == sortarray[1][axis]) &&
         (sortarray[0][1 - axis] > sortarray[1][1 - axis]))) {
      temp = sortarray[1];
      sortarray[1] = sortarray[0];
      sortarray[0] = temp;
    }
    return;
  }
  /* Choose a random pivot to split the array. */
  pivot = (int) randomnation((unsigned int) arraysize);
  pivot1 = sortarray[pivot][axis];
  pivot2 = sortarray[pivot][1 - axis];
  /* Split the array. */
  left = -1;
  right = arraysize;
  while (left < right) {
    /* Search for a vertex whose x-coordinate is too large for the left. */
    do {
      left++;
    } while ((left <= right) && ((sortarray[left][axis] < pivot1) ||
                                 ((sortarray[left][axis] == pivot1) &&
                                  (sortarray[left][1 - axis] < pivot2))));
    /* Search for a vertex whose x-coordinate is too small for the right. */
    do {
      right--;
    } while ((left <= right) && ((sortarray[right][axis] > pivot1) ||
                                 ((sortarray[right][axis] == pivot1) &&
                                  (sortarray[right][1 - axis] > pivot2))));
    if (left < right) {
      /* Swap the left and right vertices. */
      temp = sortarray[left];
      sortarray[left] = sortarray[right];
      sortarray[right] = temp;
    }
  }
  /* Unlike in vertexsort(), at most one of the following */
  /*   conditionals is true.                             */
  if (left > median) {
    /* Recursively shuffle the left subset. */
    vertexmedian(sortarray, left, median, axis);
  }
  if (right < median - 1) {
    /* Recursively shuffle the right subset. */
    vertexmedian(&sortarray[right + 1], arraysize - right - 1,
                 median - right - 1, axis);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  alternateaxes()   Sorts the vertices as appropriate for the divide-and-  */
/*                    conquer algorithm with alternating cuts.               */
/*                                                                           */
/*  Partitions by x-coordinate if axis == 0; by y-coordinate if axis == 1.   */
/*  For the base case, subsets containing only two or three vertices are     */
/*  always sorted by x-coordinate.                                           */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void alternateaxes(vertex *sortarray, int arraysize, int axis)
#else /* not ANSI_DECLARATORS */
void alternateaxes(sortarray, arraysize, axis)
vertex *sortarray;
int arraysize;
int axis;
#endif /* not ANSI_DECLARATORS */

{
  int divider;

  divider = arraysize >> 1;
  if (arraysize <= 3) {
    /* Recursive base case:  subsets of two or three vertices will be    */
    /*   handled specially, and should always be sorted by x-coordinate. */
    axis = 0;
  }
  /* Partition with a horizontal or vertical cut. */
  vertexmedian(sortarray, arraysize, divider, axis);
  /* Recursively partition the subsets with a cross cut. */
  if (arraysize - divider >= 2) {
    if (divider >= 2) {
      alternateaxes(sortarray, divider, 1 - axis);
    }
    alternateaxes(&sortarray[divider], arraysize - divider, 1 - axis);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  mergehulls()   Merge two adjacent Delaunay triangulations into a         */
/*                 single Delaunay triangulation.                            */
/*                                                                           */
/*  This is similar to the algorithm given by Guibas and Stolfi, but uses    */
/*  a triangle-based, rather than edge-based, data structure.                */
/*                                                                           */
/*  The algorithm walks up the gap between the two triangulations, knitting  */
/*  them together.  As they are merged, some of their bounding triangles     */
/*  are converted into real triangles of the triangulation.  The procedure   */
/*  pulls each hull's bounding triangles apart, then knits them together     */
/*  like the teeth of two gears.  The Delaunay property determines, at each  */
/*  step, whether the next "tooth" is a bounding triangle of the left hull   */
/*  or the right.  When a bounding triangle becomes real, its apex is        */
/*  changed from NULL to a real vertex.                                      */
/*                                                                           */
/*  Only two new triangles need to be allocated.  These become new bounding  */
/*  triangles at the top and bottom of the seam.  They are used to connect   */
/*  the remaining bounding triangles (those that have not been converted     */
/*  into real triangles) into a single fan.                                  */
/*                                                                           */
/*  On entry, `farleft' and `innerleft' are bounding triangles of the left   */
/*  triangulation.  The origin of `farleft' is the leftmost vertex, and      */
/*  the destination of `innerleft' is the rightmost vertex of the            */
/*  triangulation.  Similarly, `innerright' and `farright' are bounding      */
/*  triangles of the right triangulation.  The origin of `innerright' and    */
/*  destination of `farright' are the leftmost and rightmost vertices.       */
/*                                                                           */
/*  On completion, the origin of `farleft' is the leftmost vertex of the     */
/*  merged triangulation, and the destination of `farright' is the rightmost */
/*  vertex.                                                                  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void mergehulls(struct mesh *m, struct behavior *b, struct otri *farleft,
                struct otri *innerleft, struct otri *innerright,
                struct otri *farright, int axis)
#else /* not ANSI_DECLARATORS */
void mergehulls(m, b, farleft, innerleft, innerright, farright, axis)
struct mesh *m;
struct behavior *b;
struct otri *farleft;
struct otri *innerleft;
struct otri *innerright;
struct otri *farright;
int axis;
#endif /* not ANSI_DECLARATORS */

{
  struct otri leftcand, rightcand;
  struct otri baseedge;
  struct otri nextedge;
  struct otri sidecasing, topcasing, outercasing;
  struct otri checkedge;
  vertex innerleftdest;
  vertex innerrightorg;
  vertex innerleftapex, innerrightapex;
  vertex farleftpt, farrightpt;
  vertex farleftapex, farrightapex;
  vertex lowerleft, lowerright;
  vertex upperleft, upperright;
  vertex nextapex;
  vertex checkvertex;
  int changemade;
  int badedge;
  int leftfinished, rightfinished;
  triangle ptr;                         /* Temporary variable used by sym(). */

  dest(*innerleft, innerleftdest);
  apex(*innerleft, innerleftapex);
  org(*innerright, innerrightorg);
  apex(*innerright, innerrightapex);
  /* Special treatment for horizontal cuts. */
  if (b->dwyer && (axis == 1)) {
    org(*farleft, farleftpt);
    apex(*farleft, farleftapex);
    dest(*farright, farrightpt);
    apex(*farright, farrightapex);
    /* The pointers to the extremal vertices are shifted to point to the */
    /*   topmost and bottommost vertex of each hull, rather than the     */
    /*   leftmost and rightmost vertices.                                */
    while (farleftapex[1] < farleftpt[1]) {
      lnextself(*farleft);
      symself(*farleft);
      farleftpt = farleftapex;
      apex(*farleft, farleftapex);
    }
    sym(*innerleft, checkedge);
    apex(checkedge, checkvertex);
    while (checkvertex[1] > innerleftdest[1]) {
      lnext(checkedge, *innerleft);
      innerleftapex = innerleftdest;
      innerleftdest = checkvertex;
      sym(*innerleft, checkedge);
      apex(checkedge, checkvertex);
    }
    while (innerrightapex[1] < innerrightorg[1]) {
      lnextself(*innerright);
      symself(*innerright);
      innerrightorg = innerrightapex;
      apex(*innerright, innerrightapex);
    }
    sym(*farright, checkedge);
    apex(checkedge, checkvertex);
    while (checkvertex[1] > farrightpt[1]) {
      lnext(checkedge, *farright);
      farrightapex = farrightpt;
      farrightpt = checkvertex;
      sym(*farright, checkedge);
      apex(checkedge, checkvertex);
    }
  }
  /* Find a line tangent to and below both hulls. */
  do {
    changemade = 0;
    /* Make innerleftdest the "bottommost" vertex of the left hull. */
    if (counterclockwise(m, b, innerleftdest, innerleftapex, innerrightorg) >
        0.0) {
      lprevself(*innerleft);
      symself(*innerleft);
      innerleftdest = innerleftapex;
      apex(*innerleft, innerleftapex);
      changemade = 1;
    }
    /* Make innerrightorg the "bottommost" vertex of the right hull. */
    if (counterclockwise(m, b, innerrightapex, innerrightorg, innerleftdest) >
        0.0) {
      lnextself(*innerright);
      symself(*innerright);
      innerrightorg = innerrightapex;
      apex(*innerright, innerrightapex);
      changemade = 1;
    }
  } while (changemade);
  /* Find the two candidates to be the next "gear tooth." */
  sym(*innerleft, leftcand);
  sym(*innerright, rightcand);
  /* Create the bottom new bounding triangle. */
  maketriangle(m, b, &baseedge);
  /* Connect it to the bounding boxes of the left and right triangulations. */
  bond(baseedge, *innerleft);
  lnextself(baseedge);
  bond(baseedge, *innerright);
  lnextself(baseedge);
  setorg(baseedge, innerrightorg);
  setdest(baseedge, innerleftdest);
  /* Apex is intentionally left NULL. */
  if (b->verbose > 2) {
    printf("  Creating base bounding ");
    printtriangle(m, b, &baseedge);
  }
  /* Fix the extreme triangles if necessary. */
  org(*farleft, farleftpt);
  if (innerleftdest == farleftpt) {
    lnext(baseedge, *farleft);
  }
  dest(*farright, farrightpt);
  if (innerrightorg == farrightpt) {
    lprev(baseedge, *farright);
  }
  /* The vertices of the current knitting edge. */
  lowerleft = innerleftdest;
  lowerright = innerrightorg;
  /* The candidate vertices for knitting. */
  apex(leftcand, upperleft);
  apex(rightcand, upperright);
  /* Walk up the gap between the two triangulations, knitting them together. */
  while (1) {
    /* Have we reached the top?  (This isn't quite the right question,       */
    /*   because even though the left triangulation might seem finished now, */
    /*   moving up on the right triangulation might reveal a new vertex of   */
    /*   the left triangulation.  And vice-versa.)                           */
    leftfinished = counterclockwise(m, b, upperleft, lowerleft, lowerright) <=
                   0.0;
    rightfinished = counterclockwise(m, b, upperright, lowerleft, lowerright)
                 <= 0.0;
    if (leftfinished && rightfinished) {
      /* Create the top new bounding triangle. */
      maketriangle(m, b, &nextedge);
      setorg(nextedge, lowerleft);
      setdest(nextedge, lowerright);
      /* Apex is intentionally left NULL. */
      /* Connect it to the bounding boxes of the two triangulations. */
      bond(nextedge, baseedge);
      lnextself(nextedge);
      bond(nextedge, rightcand);
      lnextself(nextedge);
      bond(nextedge, leftcand);
      if (b->verbose > 2) {
        printf("  Creating top bounding ");
        printtriangle(m, b, &nextedge);
      }
      /* Special treatment for horizontal cuts. */
      if (b->dwyer && (axis == 1)) {
        org(*farleft, farleftpt);
        apex(*farleft, farleftapex);
        dest(*farright, farrightpt);
        apex(*farright, farrightapex);
        sym(*farleft, checkedge);
        apex(checkedge, checkvertex);
        /* The pointers to the extremal vertices are restored to the  */
        /*   leftmost and rightmost vertices (rather than topmost and */
        /*   bottommost).                                             */
        while (checkvertex[0] < farleftpt[0]) {
          lprev(checkedge, *farleft);
          farleftapex = farleftpt;
          farleftpt = checkvertex;
          sym(*farleft, checkedge);
          apex(checkedge, checkvertex);
        }
        while (farrightapex[0] > farrightpt[0]) {
          lprevself(*farright);
          symself(*farright);
          farrightpt = farrightapex;
          apex(*farright, farrightapex);
        }
      }
      return;
    }
    /* Consider eliminating edges from the left triangulation. */
    if (!leftfinished) {
      /* What vertex would be exposed if an edge were deleted? */
      lprev(leftcand, nextedge);
      symself(nextedge);
      apex(nextedge, nextapex);
      /* If nextapex is NULL, then no vertex would be exposed; the */
      /*   triangulation would have been eaten right through.      */
      if (nextapex != (vertex) NULL) {
        /* Check whether the edge is Delaunay. */
        badedge = incircle(m, b, lowerleft, lowerright, upperleft, nextapex) >
                  0.0;
        while (badedge) {
          /* Eliminate the edge with an edge flip.  As a result, the    */
          /*   left triangulation will have one more boundary triangle. */
          lnextself(nextedge);
          sym(nextedge, topcasing);
          lnextself(nextedge);
          sym(nextedge, sidecasing);
          bond(nextedge, topcasing);
          bond(leftcand, sidecasing);
          lnextself(leftcand);
          sym(leftcand, outercasing);
          lprevself(nextedge);
          bond(nextedge, outercasing);
          /* Correct the vertices to reflect the edge flip. */
          setorg(leftcand, lowerleft);
          setdest(leftcand, NULL);
          setapex(leftcand, nextapex);
          setorg(nextedge, NULL);
          setdest(nextedge, upperleft);
          setapex(nextedge, nextapex);
          /* Consider the newly exposed vertex. */
          upperleft = nextapex;
          /* What vertex would be exposed if another edge were deleted? */
          otricopy(sidecasing, nextedge);
          apex(nextedge, nextapex);
          if (nextapex != (vertex) NULL) {
            /* Check whether the edge is Delaunay. */
            badedge = incircle(m, b, lowerleft, lowerright, upperleft,
                               nextapex) > 0.0;
          } else {
            /* Avoid eating right through the triangulation. */
            badedge = 0;
          }
        }
      }
    }
    /* Consider eliminating edges from the right triangulation. */
    if (!rightfinished) {
      /* What vertex would be exposed if an edge were deleted? */
      lnext(rightcand, nextedge);
      symself(nextedge);
      apex(nextedge, nextapex);
      /* If nextapex is NULL, then no vertex would be exposed; the */
      /*   triangulation would have been eaten right through.      */
      if (nextapex != (vertex) NULL) {
        /* Check whether the edge is Delaunay. */
        badedge = incircle(m, b, lowerleft, lowerright, upperright, nextapex) >
                  0.0;
        while (badedge) {
          /* Eliminate the edge with an edge flip.  As a result, the     */
          /*   right triangulation will have one more boundary triangle. */
          lprevself(nextedge);
          sym(nextedge, topcasing);
          lprevself(nextedge);
          sym(nextedge, sidecasing);
          bond(nextedge, topcasing);
          bond(rightcand, sidecasing);
          lprevself(rightcand);
          sym(rightcand, outercasing);
          lnextself(nextedge);
          bond(nextedge, outercasing);
          /* Correct the vertices to reflect the edge flip. */
          setorg(rightcand, NULL);
          setdest(rightcand, lowerright);
          setapex(rightcand, nextapex);
          setorg(nextedge, upperright);
          setdest(nextedge, NULL);
          setapex(nextedge, nextapex);
          /* Consider the newly exposed vertex. */
          upperright = nextapex;
          /* What vertex would be exposed if another edge were deleted? */
          otricopy(sidecasing, nextedge);
          apex(nextedge, nextapex);
          if (nextapex != (vertex) NULL) {
            /* Check whether the edge is Delaunay. */
            badedge = incircle(m, b, lowerleft, lowerright, upperright,
                               nextapex) > 0.0;
          } else {
            /* Avoid eating right through the triangulation. */
            badedge = 0;
          }
        }
      }
    }
    if (leftfinished || (!rightfinished &&
           (incircle(m, b, upperleft, lowerleft, lowerright, upperright) >
            0.0))) {
      /* Knit the triangulations, adding an edge from `lowerleft' */
      /*   to `upperright'.                                       */
      bond(baseedge, rightcand);
      lprev(rightcand, baseedge);
      setdest(baseedge, lowerleft);
      lowerright = upperright;
      sym(baseedge, rightcand);
      apex(rightcand, upperright);
    } else {
      /* Knit the triangulations, adding an edge from `upperleft' */
      /*   to `lowerright'.                                       */
      bond(baseedge, leftcand);
      lnext(leftcand, baseedge);
      setorg(baseedge, lowerright);
      lowerleft = upperleft;
      sym(baseedge, leftcand);
      apex(leftcand, upperleft);
    }
    if (b->verbose > 2) {
      printf("  Connecting ");
      printtriangle(m, b, &baseedge);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  divconqrecurse()   Recursively form a Delaunay triangulation by the      */
/*                     divide-and-conquer method.                            */
/*                                                                           */
/*  Recursively breaks down the problem into smaller pieces, which are       */
/*  knitted together by mergehulls().  The base cases (problems of two or    */
/*  three vertices) are handled specially here.                              */
/*                                                                           */
/*  On completion, `farleft' and `farright' are bounding triangles such that */
/*  the origin of `farleft' is the leftmost vertex (breaking ties by         */
/*  choosing the highest leftmost vertex), and the destination of            */
/*  `farright' is the rightmost vertex (breaking ties by choosing the        */
/*  lowest rightmost vertex).                                                */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void divconqrecurse(struct mesh *m, struct behavior *b, vertex *sortarray,
                    int vertices, int axis,
                    struct otri *farleft, struct otri *farright)
#else /* not ANSI_DECLARATORS */
void divconqrecurse(m, b, sortarray, vertices, axis, farleft, farright)
struct mesh *m;
struct behavior *b;
vertex *sortarray;
int vertices;
int axis;
struct otri *farleft;
struct otri *farright;
#endif /* not ANSI_DECLARATORS */

{
  struct otri midtri, tri1, tri2, tri3;
  struct otri innerleft, innerright;
  REAL area;
  int divider;

  if (b->verbose > 2) {
    printf("  Triangulating %d vertices.\n", vertices);
  }
  if (vertices == 2) {
    /* The triangulation of two vertices is an edge.  An edge is */
    /*   represented by two bounding triangles.                  */
    maketriangle(m, b, farleft);
    setorg(*farleft, sortarray[0]);
    setdest(*farleft, sortarray[1]);
    /* The apex is intentionally left NULL. */
    maketriangle(m, b, farright);
    setorg(*farright, sortarray[1]);
    setdest(*farright, sortarray[0]);
    /* The apex is intentionally left NULL. */
    bond(*farleft, *farright);
    lprevself(*farleft);
    lnextself(*farright);
    bond(*farleft, *farright);
    lprevself(*farleft);
    lnextself(*farright);
    bond(*farleft, *farright);
    if (b->verbose > 2) {
      printf("  Creating ");
      printtriangle(m, b, farleft);
      printf("  Creating ");
      printtriangle(m, b, farright);
    }
    /* Ensure that the origin of `farleft' is sortarray[0]. */
    lprev(*farright, *farleft);
    return;
  } else if (vertices == 3) {
    /* The triangulation of three vertices is either a triangle (with */
    /*   three bounding triangles) or two edges (with four bounding   */
    /*   triangles).  In either case, four triangles are created.     */
    maketriangle(m, b, &midtri);
    maketriangle(m, b, &tri1);
    maketriangle(m, b, &tri2);
    maketriangle(m, b, &tri3);
    area = counterclockwise(m, b, sortarray[0], sortarray[1], sortarray[2]);
    if (area == 0.0) {
      /* Three collinear vertices; the triangulation is two edges. */
      setorg(midtri, sortarray[0]);
      setdest(midtri, sortarray[1]);
      setorg(tri1, sortarray[1]);
      setdest(tri1, sortarray[0]);
      setorg(tri2, sortarray[2]);
      setdest(tri2, sortarray[1]);
      setorg(tri3, sortarray[1]);
      setdest(tri3, sortarray[2]);
      /* All apices are intentionally left NULL. */
      bond(midtri, tri1);
      bond(tri2, tri3);
      lnextself(midtri);
      lprevself(tri1);
      lnextself(tri2);
      lprevself(tri3);
      bond(midtri, tri3);
      bond(tri1, tri2);
      lnextself(midtri);
      lprevself(tri1);
      lnextself(tri2);
      lprevself(tri3);
      bond(midtri, tri1);
      bond(tri2, tri3);
      /* Ensure that the origin of `farleft' is sortarray[0]. */
      otricopy(tri1, *farleft);
      /* Ensure that the destination of `farright' is sortarray[2]. */
      otricopy(tri2, *farright);
    } else {
      /* The three vertices are not collinear; the triangulation is one */
      /*   triangle, namely `midtri'.                                   */
      setorg(midtri, sortarray[0]);
      setdest(tri1, sortarray[0]);
      setorg(tri3, sortarray[0]);
      /* Apices of tri1, tri2, and tri3 are left NULL. */
      if (area > 0.0) {
        /* The vertices are in counterclockwise order. */
        setdest(midtri, sortarray[1]);
        setorg(tri1, sortarray[1]);
        setdest(tri2, sortarray[1]);
        setapex(midtri, sortarray[2]);
        setorg(tri2, sortarray[2]);
        setdest(tri3, sortarray[2]);
      } else {
        /* The vertices are in clockwise order. */
        setdest(midtri, sortarray[2]);
        setorg(tri1, sortarray[2]);
        setdest(tri2, sortarray[2]);
        setapex(midtri, sortarray[1]);
        setorg(tri2, sortarray[1]);
        setdest(tri3, sortarray[1]);
      }
      /* The topology does not depend on how the vertices are ordered. */
      bond(midtri, tri1);
      lnextself(midtri);
      bond(midtri, tri2);
      lnextself(midtri);
      bond(midtri, tri3);
      lprevself(tri1);
      lnextself(tri2);
      bond(tri1, tri2);
      lprevself(tri1);
      lprevself(tri3);
      bond(tri1, tri3);
      lnextself(tri2);
      lprevself(tri3);
      bond(tri2, tri3);
      /* Ensure that the origin of `farleft' is sortarray[0]. */
      otricopy(tri1, *farleft);
      /* Ensure that the destination of `farright' is sortarray[2]. */
      if (area > 0.0) {
        otricopy(tri2, *farright);
      } else {
        lnext(*farleft, *farright);
      }
    }
    if (b->verbose > 2) {
      printf("  Creating ");
      printtriangle(m, b, &midtri);
      printf("  Creating ");
      printtriangle(m, b, &tri1);
      printf("  Creating ");
      printtriangle(m, b, &tri2);
      printf("  Creating ");
      printtriangle(m, b, &tri3);
    }
    return;
  } else {
    /* Split the vertices in half. */
    divider = vertices >> 1;
    /* Recursively triangulate each half. */
    divconqrecurse(m, b, sortarray, divider, 1 - axis, farleft, &innerleft);
    divconqrecurse(m, b, &sortarray[divider], vertices - divider, 1 - axis,
                   &innerright, farright);
    if (b->verbose > 1) {
      printf("  Joining triangulations with %d and %d vertices.\n", divider,
             vertices - divider);
    }
    /* Merge the two triangulations into one. */
    mergehulls(m, b, farleft, &innerleft, &innerright, farright, axis);
  }
}

#ifdef ANSI_DECLARATORS
long removeghosts(struct mesh *m, struct behavior *b, struct otri *startghost)
#else /* not ANSI_DECLARATORS */
long removeghosts(m, b, startghost)
struct mesh *m;
struct behavior *b;
struct otri *startghost;
#endif /* not ANSI_DECLARATORS */

{
  struct otri searchedge;
  struct otri dissolveedge;
  struct otri deadtriangle;
  vertex markorg;
  long hullsize;
  triangle ptr;                         /* Temporary variable used by sym(). */

  if (b->verbose) {
    printf("  Removing ghost triangles.\n");
  }
  /* Find an edge on the convex hull to start point location from. */
  lprev(*startghost, searchedge);
  symself(searchedge);
  m->dummytri[0] = encode(searchedge);
  /* Remove the bounding box and count the convex hull edges. */
  otricopy(*startghost, dissolveedge);
  hullsize = 0;
  do {
    hullsize++;
    lnext(dissolveedge, deadtriangle);
    lprevself(dissolveedge);
    symself(dissolveedge);
    /* If no PSLG is involved, set the boundary markers of all the vertices */
    /*   on the convex hull.  If a PSLG is used, this step is done later.   */
    if (!b->poly) {
      /* Watch out for the case where all the input vertices are collinear. */
      if (dissolveedge.tri != m->dummytri) {
        org(dissolveedge, markorg);
        if (vertexmark(markorg) == 0) {
          setvertexmark(markorg, 1);
        }
      }
    }
    /* Remove a bounding triangle from a convex hull triangle. */
    dissolve(dissolveedge);
    /* Find the next bounding triangle. */
    sym(deadtriangle, dissolveedge);
    /* Delete the bounding triangle. */
    triangledealloc(m, deadtriangle.tri);
  } while (!otriequal(dissolveedge, *startghost));
  return hullsize;
}

/*****************************************************************************/
/*                                                                           */
/*  divconqdelaunay()   Form a Delaunay triangulation by the divide-and-     */
/*                      conquer method.                                      */
/*                                                                           */
/*  Sorts the vertices, calls a recursive procedure to triangulate them, and */
/*  removes the bounding box, setting boundary markers as appropriate.       */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
long divconqdelaunay(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
long divconqdelaunay(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  vertex *sortarray;
  struct otri hullleft, hullright;
  int divider;
  int i, j;

  if (b->verbose) {
    printf("  Sorting vertices.\n");
  }

  /* Allocate an array of pointers to vertices for sorting. */
  sortarray = (vertex *) trimalloc(m->invertices * (int) sizeof(vertex));
  traversalinit(&m->vertices);
  for (i = 0; i < m->invertices; i++) {
    sortarray[i] = vertextraverse(m);
  }
  /* Sort the vertices. */
  vertexsort(sortarray, m->invertices);
  /* Discard duplicate vertices, which can really mess up the algorithm. */
  i = 0;
  for (j = 1; j < m->invertices; j++) {
    if ((sortarray[i][0] == sortarray[j][0])
        && (sortarray[i][1] == sortarray[j][1])) {
      if (!b->quiet) {
        printf(
"Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.\n",
               sortarray[j][0], sortarray[j][1]);
      }
      setvertextype(sortarray[j], UNDEADVERTEX);
      m->undeads++;
    } else {
      i++;
      sortarray[i] = sortarray[j];
    }
  }
  i++;
  if (b->dwyer) {
    /* Re-sort the array of vertices to accommodate alternating cuts. */
    divider = i >> 1;
    if (i - divider >= 2) {
      if (divider >= 2) {
        alternateaxes(sortarray, divider, 1);
      }
      alternateaxes(&sortarray[divider], i - divider, 1);
    }
  }

  if (b->verbose) {
    printf("  Forming triangulation.\n");
  }

  /* Form the Delaunay triangulation. */
  divconqrecurse(m, b, sortarray, i, 0, &hullleft, &hullright);
  trifree((VOID *) sortarray);

  return removeghosts(m, b, &hullleft);
}

/**                                                                         **/
/**                                                                         **/
/********* Divide-and-conquer Delaunay triangulation ends here       *********/

/********* Incremental Delaunay triangulation begins here            *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  boundingbox()   Form an "infinite" bounding triangle to insert vertices  */
/*                  into.                                                    */
/*                                                                           */
/*  The vertices at "infinity" are assigned finite coordinates, which are    */
/*  used by the point location routines, but (mostly) ignored by the         */
/*  Delaunay edge flip routines.                                             */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void boundingbox(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void boundingbox(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri inftri;          /* Handle for the triangular bounding box. */
  REAL width;

  if (b->verbose) {
    printf("  Creating triangular bounding box.\n");
  }
  /* Find the width (or height, whichever is larger) of the triangulation. */
  width = m->xmax - m->xmin;
  if (m->ymax - m->ymin > width) {
    width = m->ymax - m->ymin;
  }
  if (width == 0.0) {
    width = 1.0;
  }
  /* Create the vertices of the bounding box. */
  m->infvertex1 = (vertex) trimalloc(m->vertices.itembytes);
  m->infvertex2 = (vertex) trimalloc(m->vertices.itembytes);
  m->infvertex3 = (vertex) trimalloc(m->vertices.itembytes);
  m->infvertex1[0] = m->xmin - 50.0 * width;
  m->infvertex1[1] = m->ymin - 40.0 * width;
  m->infvertex2[0] = m->xmax + 50.0 * width;
  m->infvertex2[1] = m->ymin - 40.0 * width;
  m->infvertex3[0] = 0.5 * (m->xmin + m->xmax);
  m->infvertex3[1] = m->ymax + 60.0 * width;

  /* Create the bounding box. */
  maketriangle(m, b, &inftri);
  setorg(inftri, m->infvertex1);
  setdest(inftri, m->infvertex2);
  setapex(inftri, m->infvertex3);
  /* Link dummytri to the bounding box so we can always find an */
  /*   edge to begin searching (point location) from.           */
  m->dummytri[0] = (triangle) inftri.tri;
  if (b->verbose > 2) {
    printf("  Creating ");
    printtriangle(m, b, &inftri);
  }
}

#endif /* not REDUCED */

/*****************************************************************************/
/*                                                                           */
/*  removebox()   Remove the "infinite" bounding triangle, setting boundary  */
/*                markers as appropriate.                                    */
/*                                                                           */
/*  The triangular bounding box has three boundary triangles (one for each   */
/*  side of the bounding box), and a bunch of triangles fanning out from     */
/*  the three bounding box vertices (one triangle for each edge of the       */
/*  convex hull of the inner mesh).  This routine removes these triangles.   */
/*                                                                           */
/*  Returns the number of edges on the convex hull of the triangulation.     */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
long removebox(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
long removebox(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri deadtriangle;
  struct otri searchedge;
  struct otri checkedge;
  struct otri nextedge, finaledge, dissolveedge;
  vertex markorg;
  long hullsize;
  triangle ptr;                         /* Temporary variable used by sym(). */

  if (b->verbose) {
    printf("  Removing triangular bounding box.\n");
  }
  /* Find a boundary triangle. */
  nextedge.tri = m->dummytri;
  nextedge.orient = 0;
  symself(nextedge);
  /* Mark a place to stop. */
  lprev(nextedge, finaledge);
  lnextself(nextedge);
  symself(nextedge);
  /* Find a triangle (on the boundary of the vertex set) that isn't */
  /*   a bounding box triangle.                                     */
  lprev(nextedge, searchedge);
  symself(searchedge);
  /* Check whether nextedge is another boundary triangle */
  /*   adjacent to the first one.                        */
  lnext(nextedge, checkedge);
  symself(checkedge);
  if (checkedge.tri == m->dummytri) {
    /* Go on to the next triangle.  There are only three boundary   */
    /*   triangles, and this next triangle cannot be the third one, */
    /*   so it's safe to stop here.                                 */
    lprevself(searchedge);
    symself(searchedge);
  }
  /* Find a new boundary edge to search from, as the current search */
  /*   edge lies on a bounding box triangle and will be deleted.    */
  m->dummytri[0] = encode(searchedge);
  hullsize = -2l;
  while (!otriequal(nextedge, finaledge)) {
    hullsize++;
    lprev(nextedge, dissolveedge);
    symself(dissolveedge);
    /* If not using a PSLG, the vertices should be marked now. */
    /*   (If using a PSLG, markhull() will do the job.)        */
    if (!b->poly) {
      /* Be careful!  One must check for the case where all the input     */
      /*   vertices are collinear, and thus all the triangles are part of */
      /*   the bounding box.  Otherwise, the setvertexmark() call below   */
      /*   will cause a bad pointer reference.                            */
      if (dissolveedge.tri != m->dummytri) {
        org(dissolveedge, markorg);
        if (vertexmark(markorg) == 0) {
          setvertexmark(markorg, 1);
        }
      }
    }
    /* Disconnect the bounding box triangle from the mesh triangle. */
    dissolve(dissolveedge);
    lnext(nextedge, deadtriangle);
    sym(deadtriangle, nextedge);
    /* Get rid of the bounding box triangle. */
    triangledealloc(m, deadtriangle.tri);
    /* Do we need to turn the corner? */
    if (nextedge.tri == m->dummytri) {
      /* Turn the corner. */
      otricopy(dissolveedge, nextedge);
    }
  }
  triangledealloc(m, finaledge.tri);

  trifree((VOID *) m->infvertex1);  /* Deallocate the bounding box vertices. */
  trifree((VOID *) m->infvertex2);
  trifree((VOID *) m->infvertex3);

  return hullsize;
}

#endif /* not REDUCED */

/*****************************************************************************/
/*                                                                           */
/*  incrementaldelaunay()   Form a Delaunay triangulation by incrementally   */
/*                          inserting vertices.                              */
/*                                                                           */
/*  Returns the number of edges on the convex hull of the triangulation.     */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
long incrementaldelaunay(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
long incrementaldelaunay(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri starttri;
  vertex vertexloop;

  /* Create a triangular bounding box. */
  boundingbox(m, b);
  if (b->verbose) {
    printf("  Incrementally inserting vertices.\n");
  }
  traversalinit(&m->vertices);
  vertexloop = vertextraverse(m);
  while (vertexloop != (vertex) NULL) {
    starttri.tri = m->dummytri;
    if (insertvertex(m, b, vertexloop, &starttri, (struct osub *) NULL, 0, 0)
        == DUPLICATEVERTEX) {
      if (!b->quiet) {
        printf(
"Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.\n",
               vertexloop[0], vertexloop[1]);
      }
      setvertextype(vertexloop, UNDEADVERTEX);
      m->undeads++;
    }
    vertexloop = vertextraverse(m);
  }
  /* Remove the bounding box. */
  return removebox(m, b);
}

#endif /* not REDUCED */

/**                                                                         **/
/**                                                                         **/
/********* Incremental Delaunay triangulation ends here              *********/

/********* Sweepline Delaunay triangulation begins here              *********/
/**                                                                         **/
/**                                                                         **/

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void eventheapinsert(struct event **heap, int heapsize, struct event *newevent)
#else /* not ANSI_DECLARATORS */
void eventheapinsert(heap, heapsize, newevent)
struct event **heap;
int heapsize;
struct event *newevent;
#endif /* not ANSI_DECLARATORS */

{
  REAL eventx, eventy;
  int eventnum;
  int parent;
  int notdone;

  eventx = newevent->xkey;
  eventy = newevent->ykey;
  eventnum = heapsize;
  notdone = eventnum > 0;
  while (notdone) {
    parent = (eventnum - 1) >> 1;
    if ((heap[parent]->ykey < eventy) ||
        ((heap[parent]->ykey == eventy)
         && (heap[parent]->xkey <= eventx))) {
      notdone = 0;
    } else {
      heap[eventnum] = heap[parent];
      heap[eventnum]->heapposition = eventnum;

      eventnum = parent;
      notdone = eventnum > 0;
    }
  }
  heap[eventnum] = newevent;
  newevent->heapposition = eventnum;
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void eventheapify(struct event **heap, int heapsize, int eventnum)
#else /* not ANSI_DECLARATORS */
void eventheapify(heap, heapsize, eventnum)
struct event **heap;
int heapsize;
int eventnum;
#endif /* not ANSI_DECLARATORS */

{
  struct event *thisevent;
  REAL eventx, eventy;
  int leftchild, rightchild;
  int smallest;
  int notdone;

  thisevent = heap[eventnum];
  eventx = thisevent->xkey;
  eventy = thisevent->ykey;
  leftchild = 2 * eventnum + 1;
  notdone = leftchild < heapsize;
  while (notdone) {
    if ((heap[leftchild]->ykey < eventy) ||
        ((heap[leftchild]->ykey == eventy)
         && (heap[leftchild]->xkey < eventx))) {
      smallest = leftchild;
    } else {
      smallest = eventnum;
    }
    rightchild = leftchild + 1;
    if (rightchild < heapsize) {
      if ((heap[rightchild]->ykey < heap[smallest]->ykey) ||
          ((heap[rightchild]->ykey == heap[smallest]->ykey)
           && (heap[rightchild]->xkey < heap[smallest]->xkey))) {
        smallest = rightchild;
      }
    }
    if (smallest == eventnum) {
      notdone = 0;
    } else {
      heap[eventnum] = heap[smallest];
      heap[eventnum]->heapposition = eventnum;
      heap[smallest] = thisevent;
      thisevent->heapposition = smallest;

      eventnum = smallest;
      leftchild = 2 * eventnum + 1;
      notdone = leftchild < heapsize;
    }
  }
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void eventheapdelete(struct event **heap, int heapsize, int eventnum)
#else /* not ANSI_DECLARATORS */
void eventheapdelete(heap, heapsize, eventnum)
struct event **heap;
int heapsize;
int eventnum;
#endif /* not ANSI_DECLARATORS */

{
  struct event *moveevent;
  REAL eventx, eventy;
  int parent;
  int notdone;

  moveevent = heap[heapsize - 1];
  if (eventnum > 0) {
    eventx = moveevent->xkey;
    eventy = moveevent->ykey;
    do {
      parent = (eventnum - 1) >> 1;
      if ((heap[parent]->ykey < eventy) ||
          ((heap[parent]->ykey == eventy)
           && (heap[parent]->xkey <= eventx))) {
        notdone = 0;
      } else {
        heap[eventnum] = heap[parent];
        heap[eventnum]->heapposition = eventnum;

        eventnum = parent;
        notdone = eventnum > 0;
      }
    } while (notdone);
  }
  heap[eventnum] = moveevent;
  moveevent->heapposition = eventnum;
  eventheapify(heap, heapsize - 1, eventnum);
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void createeventheap(struct mesh *m, struct event ***eventheap,
                     struct event **events, struct event **freeevents)
#else /* not ANSI_DECLARATORS */
void createeventheap(m, eventheap, events, freeevents)
struct mesh *m;
struct event ***eventheap;
struct event **events;
struct event **freeevents;
#endif /* not ANSI_DECLARATORS */

{
  vertex thisvertex;
  int maxevents;
  int i;

  maxevents = (3 * m->invertices) / 2;
  *eventheap = (struct event **) trimalloc(maxevents *
                                           (int) sizeof(struct event *));
  *events = (struct event *) trimalloc(maxevents * (int) sizeof(struct event));
  traversalinit(&m->vertices);
  for (i = 0; i < m->invertices; i++) {
    thisvertex = vertextraverse(m);
    (*events)[i].eventptr = (VOID *) thisvertex;
    (*events)[i].xkey = thisvertex[0];
    (*events)[i].ykey = thisvertex[1];
    eventheapinsert(*eventheap, i, *events + i);
  }
  *freeevents = (struct event *) NULL;
  for (i = maxevents - 1; i >= m->invertices; i--) {
    (*events)[i].eventptr = (VOID *) *freeevents;
    *freeevents = *events + i;
  }
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
int rightofhyperbola(struct mesh *m, struct otri *fronttri, vertex newsite)
#else /* not ANSI_DECLARATORS */
int rightofhyperbola(m, fronttri, newsite)
struct mesh *m;
struct otri *fronttri;
vertex newsite;
#endif /* not ANSI_DECLARATORS */

{
  vertex leftvertex, rightvertex;
  REAL dxa, dya, dxb, dyb;

  m->hyperbolacount++;

  dest(*fronttri, leftvertex);
  apex(*fronttri, rightvertex);
  if ((leftvertex[1] < rightvertex[1]) ||
      ((leftvertex[1] == rightvertex[1]) &&
       (leftvertex[0] < rightvertex[0]))) {
    if (newsite[0] >= rightvertex[0]) {
      return 1;
    }
  } else {
    if (newsite[0] <= leftvertex[0]) {
      return 0;
    }
  }
  dxa = leftvertex[0] - newsite[0];
  dya = leftvertex[1] - newsite[1];
  dxb = rightvertex[0] - newsite[0];
  dyb = rightvertex[1] - newsite[1];
  return dya * (dxb * dxb + dyb * dyb) > dyb * (dxa * dxa + dya * dya);
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
REAL circletop(struct mesh *m, vertex pa, vertex pb, vertex pc, REAL ccwabc)
#else /* not ANSI_DECLARATORS */
REAL circletop(m, pa, pb, pc, ccwabc)
struct mesh *m;
vertex pa;
vertex pb;
vertex pc;
REAL ccwabc;
#endif /* not ANSI_DECLARATORS */

{
  REAL xac, yac, xbc, ybc, xab, yab;
  REAL aclen2, bclen2, ablen2;

  m->circletopcount++;

  xac = pa[0] - pc[0];
  yac = pa[1] - pc[1];
  xbc = pb[0] - pc[0];
  ybc = pb[1] - pc[1];
  xab = pa[0] - pb[0];
  yab = pa[1] - pb[1];
  aclen2 = xac * xac + yac * yac;
  bclen2 = xbc * xbc + ybc * ybc;
  ablen2 = xab * xab + yab * yab;
  return pc[1] + (xac * bclen2 - xbc * aclen2 + sqrt(aclen2 * bclen2 * ablen2))
               / (2.0 * ccwabc);
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
void check4deadevent(struct otri *checktri, struct event **freeevents,
                     struct event **eventheap, int *heapsize)
#else /* not ANSI_DECLARATORS */
void check4deadevent(checktri, freeevents, eventheap, heapsize)
struct otri *checktri;
struct event **freeevents;
struct event **eventheap;
int *heapsize;
#endif /* not ANSI_DECLARATORS */

{
  struct event *deadevent;
  vertex eventvertex;
  int eventnum;

  org(*checktri, eventvertex);
  if (eventvertex != (vertex) NULL) {
    deadevent = (struct event *) eventvertex;
    eventnum = deadevent->heapposition;
    deadevent->eventptr = (VOID *) *freeevents;
    *freeevents = deadevent;
    eventheapdelete(eventheap, *heapsize, eventnum);
    (*heapsize)--;
    setorg(*checktri, NULL);
  }
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
struct splaynode *splay(struct mesh *m, struct splaynode *splaytree,
                        vertex searchpoint, struct otri *searchtri)
#else /* not ANSI_DECLARATORS */
struct splaynode *splay(m, splaytree, searchpoint, searchtri)
struct mesh *m;
struct splaynode *splaytree;
vertex searchpoint;
struct otri *searchtri;
#endif /* not ANSI_DECLARATORS */

{
  struct splaynode *child, *grandchild;
  struct splaynode *lefttree, *righttree;
  struct splaynode *leftright;
  vertex checkvertex;
  int rightofroot, rightofchild;

  if (splaytree == (struct splaynode *) NULL) {
    return (struct splaynode *) NULL;
  }
  dest(splaytree->keyedge, checkvertex);
  if (checkvertex == splaytree->keydest) {
    rightofroot = rightofhyperbola(m, &splaytree->keyedge, searchpoint);
    if (rightofroot) {
      otricopy(splaytree->keyedge, *searchtri);
      child = splaytree->rchild;
    } else {
      child = splaytree->lchild;
    }
    if (child == (struct splaynode *) NULL) {
      return splaytree;
    }
    dest(child->keyedge, checkvertex);
    if (checkvertex != child->keydest) {
      child = splay(m, child, searchpoint, searchtri);
      if (child == (struct splaynode *) NULL) {
        if (rightofroot) {
          splaytree->rchild = (struct splaynode *) NULL;
        } else {
          splaytree->lchild = (struct splaynode *) NULL;
        }
        return splaytree;
      }
    }
    rightofchild = rightofhyperbola(m, &child->keyedge, searchpoint);
    if (rightofchild) {
      otricopy(child->keyedge, *searchtri);
      grandchild = splay(m, child->rchild, searchpoint, searchtri);
      child->rchild = grandchild;
    } else {
      grandchild = splay(m, child->lchild, searchpoint, searchtri);
      child->lchild = grandchild;
    }
    if (grandchild == (struct splaynode *) NULL) {
      if (rightofroot) {
        splaytree->rchild = child->lchild;
        child->lchild = splaytree;
      } else {
        splaytree->lchild = child->rchild;
        child->rchild = splaytree;
      }
      return child;
    }
    if (rightofchild) {
      if (rightofroot) {
        splaytree->rchild = child->lchild;
        child->lchild = splaytree;
      } else {
        splaytree->lchild = grandchild->rchild;
        grandchild->rchild = splaytree;
      }
      child->rchild = grandchild->lchild;
      grandchild->lchild = child;
    } else {
      if (rightofroot) {
        splaytree->rchild = grandchild->lchild;
        grandchild->lchild = splaytree;
      } else {
        splaytree->lchild = child->rchild;
        child->rchild = splaytree;
      }
      child->lchild = grandchild->rchild;
      grandchild->rchild = child;
    }
    return grandchild;
  } else {
    lefttree = splay(m, splaytree->lchild, searchpoint, searchtri);
    righttree = splay(m, splaytree->rchild, searchpoint, searchtri);

    pooldealloc(&m->splaynodes, (VOID *) splaytree);
    if (lefttree == (struct splaynode *) NULL) {
      return righttree;
    } else if (righttree == (struct splaynode *) NULL) {
      return lefttree;
    } else if (lefttree->rchild == (struct splaynode *) NULL) {
      lefttree->rchild = righttree->lchild;
      righttree->lchild = lefttree;
      return righttree;
    } else if (righttree->lchild == (struct splaynode *) NULL) {
      righttree->lchild = lefttree->rchild;
      lefttree->rchild = righttree;
      return lefttree;
    } else {
/*      printf("Holy Toledo!!!\n"); */
      leftright = lefttree->rchild;
      while (leftright->rchild != (struct splaynode *) NULL) {
        leftright = leftright->rchild;
      }
      leftright->rchild = righttree;
      return lefttree;
    }
  }
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
struct splaynode *splayinsert(struct mesh *m, struct splaynode *splayroot,
                              struct otri *newkey, vertex searchpoint)
#else /* not ANSI_DECLARATORS */
struct splaynode *splayinsert(m, splayroot, newkey, searchpoint)
struct mesh *m;
struct splaynode *splayroot;
struct otri *newkey;
vertex searchpoint;
#endif /* not ANSI_DECLARATORS */

{
  struct splaynode *newsplaynode;

  newsplaynode = (struct splaynode *) poolalloc(&m->splaynodes);
  otricopy(*newkey, newsplaynode->keyedge);
  dest(*newkey, newsplaynode->keydest);
  if (splayroot == (struct splaynode *) NULL) {
    newsplaynode->lchild = (struct splaynode *) NULL;
    newsplaynode->rchild = (struct splaynode *) NULL;
  } else if (rightofhyperbola(m, &splayroot->keyedge, searchpoint)) {
    newsplaynode->lchild = splayroot;
    newsplaynode->rchild = splayroot->rchild;
    splayroot->rchild = (struct splaynode *) NULL;
  } else {
    newsplaynode->lchild = splayroot->lchild;
    newsplaynode->rchild = splayroot;
    splayroot->lchild = (struct splaynode *) NULL;
  }
  return newsplaynode;
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
struct splaynode *circletopinsert(struct mesh *m, struct behavior *b,
                                  struct splaynode *splayroot,
                                  struct otri *newkey,
                                  vertex pa, vertex pb, vertex pc, REAL topy)
#else /* not ANSI_DECLARATORS */
struct splaynode *circletopinsert(m, b, splayroot, newkey, pa, pb, pc, topy)
struct mesh *m;
struct behavior *b;
struct splaynode *splayroot;
struct otri *newkey;
vertex pa;
vertex pb;
vertex pc;
REAL topy;
#endif /* not ANSI_DECLARATORS */

{
  REAL ccwabc;
  REAL xac, yac, xbc, ybc;
  REAL aclen2, bclen2;
  REAL searchpoint[2];
  struct otri dummytri;

  ccwabc = counterclockwise(m, b, pa, pb, pc);
  xac = pa[0] - pc[0];
  yac = pa[1] - pc[1];
  xbc = pb[0] - pc[0];
  ybc = pb[1] - pc[1];
  aclen2 = xac * xac + yac * yac;
  bclen2 = xbc * xbc + ybc * ybc;
  searchpoint[0] = pc[0] - (yac * bclen2 - ybc * aclen2) / (2.0 * ccwabc);
  searchpoint[1] = topy;
  return splayinsert(m, splay(m, splayroot, (vertex) searchpoint, &dummytri),
                     newkey, (vertex) searchpoint);
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
struct splaynode *frontlocate(struct mesh *m, struct splaynode *splayroot,
                              struct otri *bottommost, vertex searchvertex,
                              struct otri *searchtri, int *farright)
#else /* not ANSI_DECLARATORS */
struct splaynode *frontlocate(m, splayroot, bottommost, searchvertex,
                              searchtri, farright)
struct mesh *m;
struct splaynode *splayroot;
struct otri *bottommost;
vertex searchvertex;
struct otri *searchtri;
int *farright;
#endif /* not ANSI_DECLARATORS */

{
  int farrightflag;
  triangle ptr;                       /* Temporary variable used by onext(). */

  otricopy(*bottommost, *searchtri);
  splayroot = splay(m, splayroot, searchvertex, searchtri);

  farrightflag = 0;
  while (!farrightflag && rightofhyperbola(m, searchtri, searchvertex)) {
    onextself(*searchtri);
    farrightflag = otriequal(*searchtri, *bottommost);
  }
  *farright = farrightflag;
  return splayroot;
}

#endif /* not REDUCED */

#ifndef REDUCED

#ifdef ANSI_DECLARATORS
long sweeplinedelaunay(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
long sweeplinedelaunay(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct event **eventheap;
  struct event *events;
  struct event *freeevents;
  struct event *nextevent;
  struct event *newevent;
  struct splaynode *splayroot;
  struct otri bottommost;
  struct otri searchtri;
  struct otri fliptri;
  struct otri lefttri, righttri, farlefttri, farrighttri;
  struct otri inserttri;
  vertex firstvertex, secondvertex;
  vertex nextvertex, lastvertex;
  vertex connectvertex;
  vertex leftvertex, midvertex, rightvertex;
  REAL lefttest, righttest;
  int heapsize;
  int check4events, farrightflag;
  triangle ptr;   /* Temporary variable used by sym(), onext(), and oprev(). */

  poolinit(&m->splaynodes, sizeof(struct splaynode), SPLAYNODEPERBLOCK,
           SPLAYNODEPERBLOCK, 0);
  splayroot = (struct splaynode *) NULL;

  if (b->verbose) {
    printf("  Placing vertices in event heap.\n");
  }
  createeventheap(m, &eventheap, &events, &freeevents);
  heapsize = m->invertices;

  if (b->verbose) {
    printf("  Forming triangulation.\n");
  }
  maketriangle(m, b, &lefttri);
  maketriangle(m, b, &righttri);
  bond(lefttri, righttri);
  lnextself(lefttri);
  lprevself(righttri);
  bond(lefttri, righttri);
  lnextself(lefttri);
  lprevself(righttri);
  bond(lefttri, righttri);
  firstvertex = (vertex) eventheap[0]->eventptr;
  eventheap[0]->eventptr = (VOID *) freeevents;
  freeevents = eventheap[0];
  eventheapdelete(eventheap, heapsize, 0);
  heapsize--;
  do {
    if (heapsize == 0) {
      printf("Error:  Input vertices are all identical.\n");
      triexit(1);
    }
    secondvertex = (vertex) eventheap[0]->eventptr;
    eventheap[0]->eventptr = (VOID *) freeevents;
    freeevents = eventheap[0];
    eventheapdelete(eventheap, heapsize, 0);
    heapsize--;
    if ((firstvertex[0] == secondvertex[0]) &&
        (firstvertex[1] == secondvertex[1])) {
      if (!b->quiet) {
        printf(
"Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.\n",
               secondvertex[0], secondvertex[1]);
      }
      setvertextype(secondvertex, UNDEADVERTEX);
      m->undeads++;
    }
  } while ((firstvertex[0] == secondvertex[0]) &&
           (firstvertex[1] == secondvertex[1]));
  setorg(lefttri, firstvertex);
  setdest(lefttri, secondvertex);
  setorg(righttri, secondvertex);
  setdest(righttri, firstvertex);
  lprev(lefttri, bottommost);
  lastvertex = secondvertex;
  while (heapsize > 0) {
    nextevent = eventheap[0];
    eventheapdelete(eventheap, heapsize, 0);
    heapsize--;
    check4events = 1;
    if (nextevent->xkey < m->xmin) {
      decode(nextevent->eventptr, fliptri);
      oprev(fliptri, farlefttri);
      check4deadevent(&farlefttri, &freeevents, eventheap, &heapsize);
      onext(fliptri, farrighttri);
      check4deadevent(&farrighttri, &freeevents, eventheap, &heapsize);

      if (otriequal(farlefttri, bottommost)) {
        lprev(fliptri, bottommost);
      }
      flip(m, b, &fliptri);
      setapex(fliptri, NULL);
      lprev(fliptri, lefttri);
      lnext(fliptri, righttri);
      sym(lefttri, farlefttri);

      if (randomnation(SAMPLERATE) == 0) {
        symself(fliptri);
        dest(fliptri, leftvertex);
        apex(fliptri, midvertex);
        org(fliptri, rightvertex);
        splayroot = circletopinsert(m, b, splayroot, &lefttri, leftvertex,
                                    midvertex, rightvertex, nextevent->ykey);
      }
    } else {
      nextvertex = (vertex) nextevent->eventptr;
      if ((nextvertex[0] == lastvertex[0]) &&
          (nextvertex[1] == lastvertex[1])) {
        if (!b->quiet) {
          printf(
"Warning:  A duplicate vertex at (%.12g, %.12g) appeared and was ignored.\n",
                 nextvertex[0], nextvertex[1]);
        }
        setvertextype(nextvertex, UNDEADVERTEX);
        m->undeads++;
        check4events = 0;
      } else {
        lastvertex = nextvertex;

        splayroot = frontlocate(m, splayroot, &bottommost, nextvertex,
                                &searchtri, &farrightflag);
/*
        otricopy(bottommost, searchtri);
        farrightflag = 0;
        while (!farrightflag && rightofhyperbola(m, &searchtri, nextvertex)) {
          onextself(searchtri);
          farrightflag = otriequal(searchtri, bottommost);
        }
*/

        check4deadevent(&searchtri, &freeevents, eventheap, &heapsize);

        otricopy(searchtri, farrighttri);
        sym(searchtri, farlefttri);
        maketriangle(m, b, &lefttri);
        maketriangle(m, b, &righttri);
        dest(farrighttri, connectvertex);
        setorg(lefttri, connectvertex);
        setdest(lefttri, nextvertex);
        setorg(righttri, nextvertex);
        setdest(righttri, connectvertex);
        bond(lefttri, righttri);
        lnextself(lefttri);
        lprevself(righttri);
        bond(lefttri, righttri);
        lnextself(lefttri);
        lprevself(righttri);
        bond(lefttri, farlefttri);
        bond(righttri, farrighttri);
        if (!farrightflag && otriequal(farrighttri, bottommost)) {
          otricopy(lefttri, bottommost);
        }

        if (randomnation(SAMPLERATE) == 0) {
          splayroot = splayinsert(m, splayroot, &lefttri, nextvertex);
        } else if (randomnation(SAMPLERATE) == 0) {
          lnext(righttri, inserttri);
          splayroot = splayinsert(m, splayroot, &inserttri, nextvertex);
        }
      }
    }
    nextevent->eventptr = (VOID *) freeevents;
    freeevents = nextevent;

    if (check4events) {
      apex(farlefttri, leftvertex);
      dest(lefttri, midvertex);
      apex(lefttri, rightvertex);
      lefttest = counterclockwise(m, b, leftvertex, midvertex, rightvertex);
      if (lefttest > 0.0) {
        newevent = freeevents;
        freeevents = (struct event *) freeevents->eventptr;
        newevent->xkey = m->xminextreme;
        newevent->ykey = circletop(m, leftvertex, midvertex, rightvertex,
                                   lefttest);
        newevent->eventptr = (VOID *) encode(lefttri);
        eventheapinsert(eventheap, heapsize, newevent);
        heapsize++;
        setorg(lefttri, newevent);
      }
      apex(righttri, leftvertex);
      org(righttri, midvertex);
      apex(farrighttri, rightvertex);
      righttest = counterclockwise(m, b, leftvertex, midvertex, rightvertex);
      if (righttest > 0.0) {
        newevent = freeevents;
        freeevents = (struct event *) freeevents->eventptr;
        newevent->xkey = m->xminextreme;
        newevent->ykey = circletop(m, leftvertex, midvertex, rightvertex,
                                   righttest);
        newevent->eventptr = (VOID *) encode(farrighttri);
        eventheapinsert(eventheap, heapsize, newevent);
        heapsize++;
        setorg(farrighttri, newevent);
      }
    }
  }

  pooldeinit(&m->splaynodes);
  lprevself(bottommost);
  return removeghosts(m, b, &bottommost);
}

#endif /* not REDUCED */

/**                                                                         **/
/**                                                                         **/
/********* Sweepline Delaunay triangulation ends here                *********/

/********* General mesh construction routines begin here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  delaunay()   Form a Delaunay triangulation.                              */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
long delaunay(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
long delaunay(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  long hulledges;

  m->eextras = 0;
  initializetrisubpools(m, b);

#ifdef REDUCED
  if (!b->quiet) {
    printf(
      "Constructing Delaunay triangulation by divide-and-conquer method.\n");
  }
  hulledges = divconqdelaunay(m, b);
#else /* not REDUCED */
  if (!b->quiet) {
    printf("Constructing Delaunay triangulation ");
    if (b->incremental) {
      printf("by incremental method.\n");
    } else if (b->sweepline) {
      printf("by sweepline method.\n");
    } else {
      printf("by divide-and-conquer method.\n");
    }
  }
  if (b->incremental) {
    hulledges = incrementaldelaunay(m, b);
  } else if (b->sweepline) {
    hulledges = sweeplinedelaunay(m, b);
  } else {
    hulledges = divconqdelaunay(m, b);
  }
#endif /* not REDUCED */

  if (m->triangles.items == 0) {
    /* The input vertices were all collinear, so there are no triangles. */
    return 0l;
  } else {
    return hulledges;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  reconstruct()   Reconstruct a triangulation from its .ele (and possibly  */
/*                  .poly) file.  Used when the -r switch is used.           */
/*                                                                           */
/*  Reads an .ele file and reconstructs the original mesh.  If the -p switch */
/*  is used, this procedure will also read a .poly file and reconstruct the  */
/*  subsegments of the original mesh.  If the -a switch is used, this        */
/*  procedure will also read an .area file and set a maximum area constraint */
/*  on each triangle.                                                        */
/*                                                                           */
/*  Vertices that are not corners of triangles, such as nodes on edges of    */
/*  subparametric elements, are discarded.                                   */
/*                                                                           */
/*  This routine finds the adjacencies between triangles (and subsegments)   */
/*  by forming one stack of triangles for each vertex.  Each triangle is on  */
/*  three different stacks simultaneously.  Each triangle's subsegment       */
/*  pointers are used to link the items in each stack.  This memory-saving   */
/*  feature makes the code harder to read.  The most important thing to keep */
/*  in mind is that each triangle is removed from a stack precisely when     */
/*  the corresponding pointer is adjusted to refer to a subsegment rather    */
/*  than the next triangle of the stack.                                     */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
int reconstruct(struct mesh *m, struct behavior *b, int *trianglelist,
                REAL *triangleattriblist, REAL *trianglearealist,
                int elements, int corners, int attribs,
                int *segmentlist,int *segmentmarkerlist, int numberofsegments)
#else /* not ANSI_DECLARATORS */
int reconstruct(m, b, trianglelist, triangleattriblist, trianglearealist,
                elements, corners, attribs, segmentlist, segmentmarkerlist,
                numberofsegments)
struct mesh *m;
struct behavior *b;
int *trianglelist;
REAL *triangleattriblist;
REAL *trianglearealist;
int elements;
int corners;
int attribs;
int *segmentlist;
int *segmentmarkerlist;
int numberofsegments;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
long reconstruct(struct mesh *m, struct behavior *b, char *elefilename,
                 char *areafilename, char *polyfilename, FILE *polyfile)
#else /* not ANSI_DECLARATORS */
long reconstruct(m, b, elefilename, areafilename, polyfilename, polyfile)
struct mesh *m;
struct behavior *b;
char *elefilename;
char *areafilename;
char *polyfilename;
FILE *polyfile;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  int vertexindex;
  int attribindex;
#else /* not TRILIBRARY */
  FILE *elefile;
  FILE *areafile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int areaelements;
#endif /* not TRILIBRARY */
  struct otri triangleloop;
  struct otri triangleleft;
  struct otri checktri;
  struct otri checkleft;
  struct otri checkneighbor;
  struct osub subsegloop;
  triangle *vertexarray;
  triangle *prevlink;
  triangle nexttri;
  vertex tdest, tapex;
  vertex checkdest, checkapex;
  vertex shorg;
  vertex killvertex;
  vertex segmentorg, segmentdest;
  REAL area;
  int corner[3];
  int end[2];
  int killvertexindex;
  int incorners;
  int segmentmarkers;
  int boundmarker;
  int aroundvertex;
  long hullsize;
  int notfound;
  long elementnumber, segmentnumber;
  int i, j;
  triangle ptr;                         /* Temporary variable used by sym(). */

#ifdef TRILIBRARY
  m->inelements = elements;
  incorners = corners;
  if (incorners < 3) {
    printf("Error:  Triangles must have at least 3 vertices.\n");
    triexit(1);
  }
  m->eextras = attribs;
#else /* not TRILIBRARY */
  /* Read the triangles from an .ele file. */
  if (!b->quiet) {
    printf("Opening %s.\n", elefilename);
  }
  elefile = fopen(elefilename, "r");
  if (elefile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", elefilename);
    triexit(1);
  }
  /* Read number of triangles, number of vertices per triangle, and */
  /*   number of triangle attributes from .ele file.                */
  stringptr = readline(inputline, elefile, elefilename);
  m->inelements = (int) strtol(stringptr, &stringptr, 0);
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    incorners = 3;
  } else {
    incorners = (int) strtol(stringptr, &stringptr, 0);
    if (incorners < 3) {
      printf("Error:  Triangles in %s must have at least 3 vertices.\n",
             elefilename);
      triexit(1);
    }
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    m->eextras = 0;
  } else {
    m->eextras = (int) strtol(stringptr, &stringptr, 0);
  }
#endif /* not TRILIBRARY */

  initializetrisubpools(m, b);

  /* Create the triangles. */
  for (elementnumber = 1; elementnumber <= m->inelements; elementnumber++) {
    maketriangle(m, b, &triangleloop);
    /* Mark the triangle as living. */
    triangleloop.tri[3] = (triangle) triangleloop.tri;
  }

  segmentmarkers = 0;
  if (b->poly) {
#ifdef TRILIBRARY
    m->insegments = numberofsegments;
    segmentmarkers = segmentmarkerlist != (int *) NULL;
#else /* not TRILIBRARY */
    /* Read number of segments and number of segment */
    /*   boundary markers from .poly file.           */
    stringptr = readline(inputline, polyfile, b->inpolyfilename);
    m->insegments = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr != '\0') {
      segmentmarkers = (int) strtol(stringptr, &stringptr, 0);
    }
#endif /* not TRILIBRARY */

    /* Create the subsegments. */
    for (segmentnumber = 1; segmentnumber <= m->insegments; segmentnumber++) {
      makesubseg(m, &subsegloop);
      /* Mark the subsegment as living. */
      subsegloop.ss[2] = (subseg) subsegloop.ss;
    }
  }

#ifdef TRILIBRARY
  vertexindex = 0;
  attribindex = 0;
#else /* not TRILIBRARY */
  if (b->vararea) {
    /* Open an .area file, check for consistency with the .ele file. */
    if (!b->quiet) {
      printf("Opening %s.\n", areafilename);
    }
    areafile = fopen(areafilename, "r");
    if (areafile == (FILE *) NULL) {
      printf("  Error:  Cannot access file %s.\n", areafilename);
      triexit(1);
    }
    stringptr = readline(inputline, areafile, areafilename);
    areaelements = (int) strtol(stringptr, &stringptr, 0);
    if (areaelements != m->inelements) {
      printf("Error:  %s and %s disagree on number of triangles.\n",
             elefilename, areafilename);
      triexit(1);
    }
  }
#endif /* not TRILIBRARY */

  if (!b->quiet) {
    printf("Reconstructing mesh.\n");
  }
  /* Allocate a temporary array that maps each vertex to some adjacent */
  /*   triangle.  I took care to allocate all the permanent memory for */
  /*   triangles and subsegments first.                                */
  vertexarray = (triangle *) trimalloc(m->vertices.items *
                                       (int) sizeof(triangle));
  /* Each vertex is initially unrepresented. */
  for (i = 0; i < m->vertices.items; i++) {
    vertexarray[i] = (triangle) m->dummytri;
  }

  if (b->verbose) {
    printf("  Assembling triangles.\n");
  }
  /* Read the triangles from the .ele file, and link */
  /*   together those that share an edge.            */
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  elementnumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
#ifdef TRILIBRARY
    /* Copy the triangle's three corners. */
    for (j = 0; j < 3; j++) {
      corner[j] = trianglelist[vertexindex++];
      if ((corner[j] < b->firstnumber) ||
          (corner[j] >= b->firstnumber + m->invertices)) {
        printf("Error:  Triangle %ld has an invalid vertex index.\n",
               elementnumber);
        triexit(1);
      }
    }
#else /* not TRILIBRARY */
    /* Read triangle number and the triangle's three corners. */
    stringptr = readline(inputline, elefile, elefilename);
    for (j = 0; j < 3; j++) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Triangle %ld is missing vertex %d in %s.\n",
               elementnumber, j + 1, elefilename);
        triexit(1);
      } else {
        corner[j] = (int) strtol(stringptr, &stringptr, 0);
        if ((corner[j] < b->firstnumber) ||
            (corner[j] >= b->firstnumber + m->invertices)) {
          printf("Error:  Triangle %ld has an invalid vertex index.\n",
                 elementnumber);
          triexit(1);
        }
      }
    }
#endif /* not TRILIBRARY */

    /* Find out about (and throw away) extra nodes. */
    for (j = 3; j < incorners; j++) {
#ifdef TRILIBRARY
      killvertexindex = trianglelist[vertexindex++];
#else /* not TRILIBRARY */
      stringptr = findfield(stringptr);
      if (*stringptr != '\0') {
        killvertexindex = (int) strtol(stringptr, &stringptr, 0);
#endif /* not TRILIBRARY */
        if ((killvertexindex >= b->firstnumber) &&
            (killvertexindex < b->firstnumber + m->invertices)) {
          /* Delete the non-corner vertex if it's not already deleted. */
          killvertex = getvertex(m, b, killvertexindex);
          if (vertextype(killvertex) != DEADVERTEX) {
            vertexdealloc(m, killvertex);
          }
        }
#ifndef TRILIBRARY
      }
#endif /* not TRILIBRARY */
    }

    /* Read the triangle's attributes. */
    for (j = 0; j < m->eextras; j++) {
#ifdef TRILIBRARY
      setelemattribute(triangleloop, j, triangleattriblist[attribindex++]);
#else /* not TRILIBRARY */
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        setelemattribute(triangleloop, j, 0);
      } else {
        setelemattribute(triangleloop, j,
                         (REAL) strtod(stringptr, &stringptr));
      }
#endif /* not TRILIBRARY */
    }

    if (b->vararea) {
#ifdef TRILIBRARY
      area = trianglearealist[elementnumber - b->firstnumber];
#else /* not TRILIBRARY */
      /* Read an area constraint from the .area file. */
      stringptr = readline(inputline, areafile, areafilename);
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        area = -1.0;                      /* No constraint on this triangle. */
      } else {
        area = (REAL) strtod(stringptr, &stringptr);
      }
#endif /* not TRILIBRARY */
      setareabound(triangleloop, area);
    }

    /* Set the triangle's vertices. */
    triangleloop.orient = 0;
    setorg(triangleloop, getvertex(m, b, corner[0]));
    setdest(triangleloop, getvertex(m, b, corner[1]));
    setapex(triangleloop, getvertex(m, b, corner[2]));
    /* Try linking the triangle to others that share these vertices. */
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      /* Take the number for the origin of triangleloop. */
      aroundvertex = corner[triangleloop.orient];
      /* Look for other triangles having this vertex. */
      nexttri = vertexarray[aroundvertex - b->firstnumber];
      /* Link the current triangle to the next one in the stack. */
      triangleloop.tri[6 + triangleloop.orient] = nexttri;
      /* Push the current triangle onto the stack. */
      vertexarray[aroundvertex - b->firstnumber] = encode(triangleloop);
      decode(nexttri, checktri);
      if (checktri.tri != m->dummytri) {
        dest(triangleloop, tdest);
        apex(triangleloop, tapex);
        /* Look for other triangles that share an edge. */
        do {
          dest(checktri, checkdest);
          apex(checktri, checkapex);
          if (tapex == checkdest) {
            /* The two triangles share an edge; bond them together. */
            lprev(triangleloop, triangleleft);
            bond(triangleleft, checktri);
          }
          if (tdest == checkapex) {
            /* The two triangles share an edge; bond them together. */
            lprev(checktri, checkleft);
            bond(triangleloop, checkleft);
          }
          /* Find the next triangle in the stack. */
          nexttri = checktri.tri[6 + checktri.orient];
          decode(nexttri, checktri);
        } while (checktri.tri != m->dummytri);
      }
    }
    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }

#ifdef TRILIBRARY
  vertexindex = 0;
#else /* not TRILIBRARY */
  fclose(elefile);
  if (b->vararea) {
    fclose(areafile);
  }
#endif /* not TRILIBRARY */

  hullsize = 0;                      /* Prepare to count the boundary edges. */
  if (b->poly) {
    if (b->verbose) {
      printf("  Marking segments in triangulation.\n");
    }
    /* Read the segments from the .poly file, and link them */
    /*   to their neighboring triangles.                    */
    boundmarker = 0;
    traversalinit(&m->subsegs);
    subsegloop.ss = subsegtraverse(m);
    segmentnumber = b->firstnumber;
    while (subsegloop.ss != (subseg *) NULL) {
#ifdef TRILIBRARY
      end[0] = segmentlist[vertexindex++];
      end[1] = segmentlist[vertexindex++];
      if (segmentmarkers) {
        boundmarker = segmentmarkerlist[segmentnumber - b->firstnumber];
      }
#else /* not TRILIBRARY */
      /* Read the endpoints of each segment, and possibly a boundary marker. */
      stringptr = readline(inputline, polyfile, b->inpolyfilename);
      /* Skip the first (segment number) field. */
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Segment %ld has no endpoints in %s.\n", segmentnumber,
               polyfilename);
        triexit(1);
      } else {
        end[0] = (int) strtol(stringptr, &stringptr, 0);
      }
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Segment %ld is missing its second endpoint in %s.\n",
               segmentnumber, polyfilename);
        triexit(1);
      } else {
        end[1] = (int) strtol(stringptr, &stringptr, 0);
      }
      if (segmentmarkers) {
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          boundmarker = 0;
        } else {
          boundmarker = (int) strtol(stringptr, &stringptr, 0);
        }
      }
#endif /* not TRILIBRARY */
      for (j = 0; j < 2; j++) {
        if ((end[j] < b->firstnumber) ||
            (end[j] >= b->firstnumber + m->invertices)) {
          printf("Error:  Segment %ld has an invalid vertex index.\n", 
                 segmentnumber);
          triexit(1);
        }
      }

      /* set the subsegment's vertices. */
      subsegloop.ssorient = 0;
      segmentorg = getvertex(m, b, end[0]);
      segmentdest = getvertex(m, b, end[1]);
      setsorg(subsegloop, segmentorg);
      setsdest(subsegloop, segmentdest);
      setsegorg(subsegloop, segmentorg);
      setsegdest(subsegloop, segmentdest);
      setmark(subsegloop, boundmarker);
      /* Try linking the subsegment to triangles that share these vertices. */
      for (subsegloop.ssorient = 0; subsegloop.ssorient < 2;
           subsegloop.ssorient++) {
        /* Take the number for the destination of subsegloop. */
        aroundvertex = end[1 - subsegloop.ssorient];
        /* Look for triangles having this vertex. */
        prevlink = &vertexarray[aroundvertex - b->firstnumber];
        nexttri = vertexarray[aroundvertex - b->firstnumber];
        decode(nexttri, checktri);
        sorg(subsegloop, shorg);
        notfound = 1;
        /* Look for triangles having this edge.  Note that I'm only       */
        /*   comparing each triangle's destination with the subsegment;   */
        /*   each triangle's apex is handled through a different vertex.  */
        /*   Because each triangle appears on three vertices' lists, each */
        /*   occurrence of a triangle on a list can (and does) represent  */
        /*   an edge.  In this way, most edges are represented twice, and */
        /*   every triangle-subsegment bond is represented once.          */
        while (notfound && (checktri.tri != m->dummytri)) {
          dest(checktri, checkdest);
          if (shorg == checkdest) {
            /* We have a match.  Remove this triangle from the list. */
            *prevlink = checktri.tri[6 + checktri.orient];
            /* Bond the subsegment to the triangle. */
            tsbond(checktri, subsegloop);
            /* Check if this is a boundary edge. */
            sym(checktri, checkneighbor);
            if (checkneighbor.tri == m->dummytri) {
              /* The next line doesn't insert a subsegment (because there's */
              /*   already one there), but it sets the boundary markers of  */
              /*   the existing subsegment and its vertices.                */
              insertsubseg(m, b, &checktri, 1);
              hullsize++;
            }
            notfound = 0;
          }
          /* Find the next triangle in the stack. */
          prevlink = &checktri.tri[6 + checktri.orient];
          nexttri = checktri.tri[6 + checktri.orient];
          decode(nexttri, checktri);
        }
      }
      subsegloop.ss = subsegtraverse(m);
      segmentnumber++;
    }
  }

  /* Mark the remaining edges as not being attached to any subsegment. */
  /* Also, count the (yet uncounted) boundary edges.                   */
  for (i = 0; i < m->vertices.items; i++) {
    /* Search the stack of triangles adjacent to a vertex. */
    nexttri = vertexarray[i];
    decode(nexttri, checktri);
    while (checktri.tri != m->dummytri) {
      /* Find the next triangle in the stack before this */
      /*   information gets overwritten.                 */
      nexttri = checktri.tri[6 + checktri.orient];
      /* No adjacent subsegment.  (This overwrites the stack info.) */
      tsdissolve(checktri);
      sym(checktri, checkneighbor);
      if (checkneighbor.tri == m->dummytri) {
        insertsubseg(m, b, &checktri, 1);
        hullsize++;
      }
      decode(nexttri, checktri);
    }
  }

  trifree((VOID *) vertexarray);
  return hullsize;
}

#endif /* not CDT_ONLY */

/**                                                                         **/
/**                                                                         **/
/********* General mesh construction routines end here               *********/

/********* Segment insertion begins here                             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  finddirection()   Find the first triangle on the path from one point     */
/*                    to another.                                            */
/*                                                                           */
/*  Finds the triangle that intersects a line segment drawn from the         */
/*  origin of `searchtri' to the point `searchpoint', and returns the result */
/*  in `searchtri'.  The origin of `searchtri' does not change, even though  */
/*  the triangle returned may differ from the one passed in.  This routine   */
/*  is used to find the direction to move in to get from one point to        */
/*  another.                                                                 */
/*                                                                           */
/*  The return value notes whether the destination or apex of the found      */
/*  triangle is collinear with the two points in question.                   */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
enum finddirectionresult finddirection(struct mesh *m, struct behavior *b,
                                       struct otri *searchtri,
                                       vertex searchpoint)
#else /* not ANSI_DECLARATORS */
enum finddirectionresult finddirection(m, b, searchtri, searchpoint)
struct mesh *m;
struct behavior *b;
struct otri *searchtri;
vertex searchpoint;
#endif /* not ANSI_DECLARATORS */

{
  struct otri checktri;
  vertex startvertex;
  vertex leftvertex, rightvertex;
  REAL leftccw, rightccw;
  int leftflag, rightflag;
  triangle ptr;           /* Temporary variable used by onext() and oprev(). */

  org(*searchtri, startvertex);
  dest(*searchtri, rightvertex);
  apex(*searchtri, leftvertex);
  /* Is `searchpoint' to the left? */
  leftccw = counterclockwise(m, b, searchpoint, startvertex, leftvertex);
  leftflag = leftccw > 0.0;
  /* Is `searchpoint' to the right? */
  rightccw = counterclockwise(m, b, startvertex, searchpoint, rightvertex);
  rightflag = rightccw > 0.0;
  if (leftflag && rightflag) {
    /* `searchtri' faces directly away from `searchpoint'.  We could go left */
    /*   or right.  Ask whether it's a triangle or a boundary on the left.   */
    onext(*searchtri, checktri);
    if (checktri.tri == m->dummytri) {
      leftflag = 0;
    } else {
      rightflag = 0;
    }
  }
  while (leftflag) {
    /* Turn left until satisfied. */
    onextself(*searchtri);
    if (searchtri->tri == m->dummytri) {
      printf("Internal error in finddirection():  Unable to find a\n");
      printf("  triangle leading from (%.12g, %.12g) to", startvertex[0],
             startvertex[1]);
      printf("  (%.12g, %.12g).\n", searchpoint[0], searchpoint[1]);
      internalerror();
    }
    apex(*searchtri, leftvertex);
    rightccw = leftccw;
    leftccw = counterclockwise(m, b, searchpoint, startvertex, leftvertex);
    leftflag = leftccw > 0.0;
  }
  while (rightflag) {
    /* Turn right until satisfied. */
    oprevself(*searchtri);
    if (searchtri->tri == m->dummytri) {
      printf("Internal error in finddirection():  Unable to find a\n");
      printf("  triangle leading from (%.12g, %.12g) to", startvertex[0],
             startvertex[1]);
      printf("  (%.12g, %.12g).\n", searchpoint[0], searchpoint[1]);
      internalerror();
    }
    dest(*searchtri, rightvertex);
    leftccw = rightccw;
    rightccw = counterclockwise(m, b, startvertex, searchpoint, rightvertex);
    rightflag = rightccw > 0.0;
  }
  if (leftccw == 0.0) {
    return LEFTCOLLINEAR;
  } else if (rightccw == 0.0) {
    return RIGHTCOLLINEAR;
  } else {
    return WITHIN;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  segmentintersection()   Find the intersection of an existing segment     */
/*                          and a segment that is being inserted.  Insert    */
/*                          a vertex at the intersection, splitting an       */
/*                          existing subsegment.                             */
/*                                                                           */
/*  The segment being inserted connects the apex of splittri to endpoint2.   */
/*  splitsubseg is the subsegment being split, and MUST adjoin splittri.     */
/*  Hence, endpoints of the subsegment being split are the origin and        */
/*  destination of splittri.                                                 */
/*                                                                           */
/*  On completion, splittri is a handle having the newly inserted            */
/*  intersection point as its origin, and endpoint1 as its destination.      */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void segmentintersection(struct mesh *m, struct behavior *b,
                         struct otri *splittri, struct osub *splitsubseg,
                         vertex endpoint2)
#else /* not ANSI_DECLARATORS */
void segmentintersection(m, b, splittri, splitsubseg, endpoint2)
struct mesh *m;
struct behavior *b;
struct otri *splittri;
struct osub *splitsubseg;
vertex endpoint2;
#endif /* not ANSI_DECLARATORS */

{
  struct osub opposubseg;
  vertex endpoint1;
  vertex torg, tdest;
  vertex leftvertex, rightvertex;
  vertex newvertex;
  enum insertvertexresult success;
  enum finddirectionresult collinear;
  REAL ex, ey;
  REAL tx, ty;
  REAL etx, ety;
  REAL split, denom;
  int i;
  triangle ptr;                       /* Temporary variable used by onext(). */
  subseg sptr;                        /* Temporary variable used by snext(). */

  /* Find the other three segment endpoints. */
  apex(*splittri, endpoint1);
  org(*splittri, torg);
  dest(*splittri, tdest);
  /* Segment intersection formulae; see the Antonio reference. */
  tx = tdest[0] - torg[0];
  ty = tdest[1] - torg[1];
  ex = endpoint2[0] - endpoint1[0];
  ey = endpoint2[1] - endpoint1[1];
  etx = torg[0] - endpoint2[0];
  ety = torg[1] - endpoint2[1];
  denom = ty * ex - tx * ey;
  if (denom == 0.0) {
    printf("Internal error in segmentintersection():");
    printf("  Attempt to find intersection of parallel segments.\n");
    internalerror();
  }
  split = (ey * etx - ex * ety) / denom;
  /* Create the new vertex. */
  newvertex = (vertex) poolalloc(&m->vertices);
  /* Interpolate its coordinate and attributes. */
  for (i = 0; i < 2 + m->nextras; i++) {
    newvertex[i] = torg[i] + split * (tdest[i] - torg[i]);
  }
  setvertexmark(newvertex, mark(*splitsubseg));
  setvertextype(newvertex, INPUTVERTEX);
  if (b->verbose > 1) {
    printf(
  "  Splitting subsegment (%.12g, %.12g) (%.12g, %.12g) at (%.12g, %.12g).\n",
           torg[0], torg[1], tdest[0], tdest[1], newvertex[0], newvertex[1]);
  }
  /* Insert the intersection vertex.  This should always succeed. */
  success = insertvertex(m, b, newvertex, splittri, splitsubseg, 0, 0);
  if (success != SUCCESSFULVERTEX) {
    printf("Internal error in segmentintersection():\n");
    printf("  Failure to split a segment.\n");
    internalerror();
  }
  /* Record a triangle whose origin is the new vertex. */
  setvertex2tri(newvertex, encode(*splittri));
  if (m->steinerleft > 0) {
    m->steinerleft--;
  }

  /* Divide the segment into two, and correct the segment endpoints. */
  ssymself(*splitsubseg);
  spivot(*splitsubseg, opposubseg);
  sdissolve(*splitsubseg);
  sdissolve(opposubseg);
  do {
    setsegorg(*splitsubseg, newvertex);
    snextself(*splitsubseg);
  } while (splitsubseg->ss != m->dummysub);
  do {
    setsegorg(opposubseg, newvertex);
    snextself(opposubseg);
  } while (opposubseg.ss != m->dummysub);

  /* Inserting the vertex may have caused edge flips.  We wish to rediscover */
  /*   the edge connecting endpoint1 to the new intersection vertex.         */
  collinear = finddirection(m, b, splittri, endpoint1);
  dest(*splittri, rightvertex);
  apex(*splittri, leftvertex);
  if ((leftvertex[0] == endpoint1[0]) && (leftvertex[1] == endpoint1[1])) {
    onextself(*splittri);
  } else if ((rightvertex[0] != endpoint1[0]) ||
             (rightvertex[1] != endpoint1[1])) {
    printf("Internal error in segmentintersection():\n");
    printf("  Topological inconsistency after splitting a segment.\n");
    internalerror();
  }
  /* `splittri' should have destination endpoint1. */
}

/*****************************************************************************/
/*                                                                           */
/*  scoutsegment()   Scout the first triangle on the path from one endpoint  */
/*                   to another, and check for completion (reaching the      */
/*                   second endpoint), a collinear vertex, or the            */
/*                   intersection of two segments.                           */
/*                                                                           */
/*  Returns one if the entire segment is successfully inserted, and zero if  */
/*  the job must be finished by conformingedge() or constrainededge().       */
/*                                                                           */
/*  If the first triangle on the path has the second endpoint as its         */
/*  destination or apex, a subsegment is inserted and the job is done.       */
/*                                                                           */
/*  If the first triangle on the path has a destination or apex that lies on */
/*  the segment, a subsegment is inserted connecting the first endpoint to   */
/*  the collinear vertex, and the search is continued from the collinear     */
/*  vertex.                                                                  */
/*                                                                           */
/*  If the first triangle on the path has a subsegment opposite its origin,  */
/*  then there is a segment that intersects the segment being inserted.      */
/*  Their intersection vertex is inserted, splitting the subsegment.         */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
int scoutsegment(struct mesh *m, struct behavior *b, struct otri *searchtri,
                 vertex endpoint2, int newmark)
#else /* not ANSI_DECLARATORS */
int scoutsegment(m, b, searchtri, endpoint2, newmark)
struct mesh *m;
struct behavior *b;
struct otri *searchtri;
vertex endpoint2;
int newmark;
#endif /* not ANSI_DECLARATORS */

{
  struct otri crosstri;
  struct osub crosssubseg;
  vertex leftvertex, rightvertex;
  enum finddirectionresult collinear;
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  collinear = finddirection(m, b, searchtri, endpoint2);
  dest(*searchtri, rightvertex);
  apex(*searchtri, leftvertex);
  if (((leftvertex[0] == endpoint2[0]) && (leftvertex[1] == endpoint2[1])) ||
      ((rightvertex[0] == endpoint2[0]) && (rightvertex[1] == endpoint2[1]))) {
    /* The segment is already an edge in the mesh. */
    if ((leftvertex[0] == endpoint2[0]) && (leftvertex[1] == endpoint2[1])) {
      lprevself(*searchtri);
    }
    /* Insert a subsegment, if there isn't already one there. */
    insertsubseg(m, b, searchtri, newmark);
    return 1;
  } else if (collinear == LEFTCOLLINEAR) {
    /* We've collided with a vertex between the segment's endpoints. */
    /* Make the collinear vertex be the triangle's origin. */
    lprevself(*searchtri);
    insertsubseg(m, b, searchtri, newmark);
    /* Insert the remainder of the segment. */
    return scoutsegment(m, b, searchtri, endpoint2, newmark);
  } else if (collinear == RIGHTCOLLINEAR) {
    /* We've collided with a vertex between the segment's endpoints. */
    insertsubseg(m, b, searchtri, newmark);
    /* Make the collinear vertex be the triangle's origin. */
    lnextself(*searchtri);
    /* Insert the remainder of the segment. */
    return scoutsegment(m, b, searchtri, endpoint2, newmark);
  } else {
    lnext(*searchtri, crosstri);
    tspivot(crosstri, crosssubseg);
    /* Check for a crossing segment. */
    if (crosssubseg.ss == m->dummysub) {
      return 0;
    } else {
      /* Insert a vertex at the intersection. */
      segmentintersection(m, b, &crosstri, &crosssubseg, endpoint2);
      otricopy(crosstri, *searchtri);
      insertsubseg(m, b, searchtri, newmark);
      /* Insert the remainder of the segment. */
      return scoutsegment(m, b, searchtri, endpoint2, newmark);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  conformingedge()   Force a segment into a conforming Delaunay            */
/*                     triangulation by inserting a vertex at its midpoint,  */
/*                     and recursively forcing in the two half-segments if   */
/*                     necessary.                                            */
/*                                                                           */
/*  Generates a sequence of subsegments connecting `endpoint1' to            */
/*  `endpoint2'.  `newmark' is the boundary marker of the segment, assigned  */
/*  to each new splitting vertex and subsegment.                             */
/*                                                                           */
/*  Note that conformingedge() does not always maintain the conforming       */
/*  Delaunay property.  Once inserted, segments are locked into place;       */
/*  vertices inserted later (to force other segments in) may render these    */
/*  fixed segments non-Delaunay.  The conforming Delaunay property will be   */
/*  restored by enforcequality() by splitting encroached subsegments.        */
/*                                                                           */
/*****************************************************************************/

#ifndef REDUCED
#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void conformingedge(struct mesh *m, struct behavior *b,
                    vertex endpoint1, vertex endpoint2, int newmark)
#else /* not ANSI_DECLARATORS */
void conformingedge(m, b, endpoint1, endpoint2, newmark)
struct mesh *m;
struct behavior *b;
vertex endpoint1;
vertex endpoint2;
int newmark;
#endif /* not ANSI_DECLARATORS */

{
  struct otri searchtri1, searchtri2;
  struct osub brokensubseg;
  vertex newvertex;
  vertex midvertex1, midvertex2;
  enum insertvertexresult success;
  int i;
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (b->verbose > 2) {
    printf("Forcing segment into triangulation by recursive splitting:\n");
    printf("  (%.12g, %.12g) (%.12g, %.12g)\n", endpoint1[0], endpoint1[1],
           endpoint2[0], endpoint2[1]);
  }
  /* Create a new vertex to insert in the middle of the segment. */
  newvertex = (vertex) poolalloc(&m->vertices);
  /* Interpolate coordinates and attributes. */
  for (i = 0; i < 2 + m->nextras; i++) {
    newvertex[i] = 0.5 * (endpoint1[i] + endpoint2[i]);
  }
  setvertexmark(newvertex, newmark);
  setvertextype(newvertex, SEGMENTVERTEX);
  /* No known triangle to search from. */
  searchtri1.tri = m->dummytri;
  /* Attempt to insert the new vertex. */
  success = insertvertex(m, b, newvertex, &searchtri1, (struct osub *) NULL,
                         0, 0);
  if (success == DUPLICATEVERTEX) {
    if (b->verbose > 2) {
      printf("  Segment intersects existing vertex (%.12g, %.12g).\n",
             newvertex[0], newvertex[1]);
    }
    /* Use the vertex that's already there. */
    vertexdealloc(m, newvertex);
    org(searchtri1, newvertex);
  } else {
    if (success == VIOLATINGVERTEX) {
      if (b->verbose > 2) {
        printf("  Two segments intersect at (%.12g, %.12g).\n",
               newvertex[0], newvertex[1]);
      }
      /* By fluke, we've landed right on another segment.  Split it. */
      tspivot(searchtri1, brokensubseg);
      success = insertvertex(m, b, newvertex, &searchtri1, &brokensubseg,
                             0, 0);
      if (success != SUCCESSFULVERTEX) {
        printf("Internal error in conformingedge():\n");
        printf("  Failure to split a segment.\n");
        internalerror();
      }
    }
    /* The vertex has been inserted successfully. */
    if (m->steinerleft > 0) {
      m->steinerleft--;
    }
  }
  otricopy(searchtri1, searchtri2);
  /* `searchtri1' and `searchtri2' are fastened at their origins to         */
  /*   `newvertex', and will be directed toward `endpoint1' and `endpoint2' */
  /*   respectively.  First, we must get `searchtri2' out of the way so it  */
  /*   won't be invalidated during the insertion of the first half of the   */
  /*   segment.                                                             */
  finddirection(m, b, &searchtri2, endpoint2);
  if (!scoutsegment(m, b, &searchtri1, endpoint1, newmark)) {
    /* The origin of searchtri1 may have changed if a collision with an */
    /*   intervening vertex on the segment occurred.                    */
    org(searchtri1, midvertex1);
    conformingedge(m, b, midvertex1, endpoint1, newmark);
  }
  if (!scoutsegment(m, b, &searchtri2, endpoint2, newmark)) {
    /* The origin of searchtri2 may have changed if a collision with an */
    /*   intervening vertex on the segment occurred.                    */
    org(searchtri2, midvertex2);
    conformingedge(m, b, midvertex2, endpoint2, newmark);
  }
}

#endif /* not CDT_ONLY */
#endif /* not REDUCED */

/*****************************************************************************/
/*                                                                           */
/*  delaunayfixup()   Enforce the Delaunay condition at an edge, fanning out */
/*                    recursively from an existing vertex.  Pay special      */
/*                    attention to stacking inverted triangles.              */
/*                                                                           */
/*  This is a support routine for inserting segments into a constrained      */
/*  Delaunay triangulation.                                                  */
/*                                                                           */
/*  The origin of fixuptri is treated as if it has just been inserted, and   */
/*  the local Delaunay condition needs to be enforced.  It is only enforced  */
/*  in one sector, however, that being the angular range defined by          */
/*  fixuptri.                                                                */
/*                                                                           */
/*  This routine also needs to make decisions regarding the "stacking" of    */
/*  triangles.  (Read the description of constrainededge() below before      */
/*  reading on here, so you understand the algorithm.)  If the position of   */
/*  the new vertex (the origin of fixuptri) indicates that the vertex before */
/*  it on the polygon is a reflex vertex, then "stack" the triangle by       */
/*  doing nothing.  (fixuptri is an inverted triangle, which is how stacked  */
/*  triangles are identified.)                                               */
/*                                                                           */
/*  Otherwise, check whether the vertex before that was a reflex vertex.     */
/*  If so, perform an edge flip, thereby eliminating an inverted triangle    */
/*  (popping it off the stack).  The edge flip may result in the creation    */
/*  of a new inverted triangle, depending on whether or not the new vertex   */
/*  is visible to the vertex three edges behind on the polygon.              */
/*                                                                           */
/*  If neither of the two vertices behind the new vertex are reflex          */
/*  vertices, fixuptri and fartri, the triangle opposite it, are not         */
/*  inverted; hence, ensure that the edge between them is locally Delaunay.  */
/*                                                                           */
/*  `leftside' indicates whether or not fixuptri is to the left of the       */
/*  segment being inserted.  (Imagine that the segment is pointing up from   */
/*  endpoint1 to endpoint2.)                                                 */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void delaunayfixup(struct mesh *m, struct behavior *b,
                   struct otri *fixuptri, int leftside)
#else /* not ANSI_DECLARATORS */
void delaunayfixup(m, b, fixuptri, leftside)
struct mesh *m;
struct behavior *b;
struct otri *fixuptri;
int leftside;
#endif /* not ANSI_DECLARATORS */

{
  struct otri neartri;
  struct otri fartri;
  struct osub faredge;
  vertex nearvertex, leftvertex, rightvertex, farvertex;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  lnext(*fixuptri, neartri);
  sym(neartri, fartri);
  /* Check if the edge opposite the origin of fixuptri can be flipped. */
  if (fartri.tri == m->dummytri) {
    return;
  }
  tspivot(neartri, faredge);
  if (faredge.ss != m->dummysub) {
    return;
  }
  /* Find all the relevant vertices. */
  apex(neartri, nearvertex);
  org(neartri, leftvertex);
  dest(neartri, rightvertex);
  apex(fartri, farvertex);
  /* Check whether the previous polygon vertex is a reflex vertex. */
  if (leftside) {
    if (counterclockwise(m, b, nearvertex, leftvertex, farvertex) <= 0.0) {
      /* leftvertex is a reflex vertex too.  Nothing can */
      /*   be done until a convex section is found.      */
      return;
    }
  } else {
    if (counterclockwise(m, b, farvertex, rightvertex, nearvertex) <= 0.0) {
      /* rightvertex is a reflex vertex too.  Nothing can */
      /*   be done until a convex section is found.       */
      return;
    }
  }
  if (counterclockwise(m, b, rightvertex, leftvertex, farvertex) > 0.0) {
    /* fartri is not an inverted triangle, and farvertex is not a reflex */
    /*   vertex.  As there are no reflex vertices, fixuptri isn't an     */
    /*   inverted triangle, either.  Hence, test the edge between the    */
    /*   triangles to ensure it is locally Delaunay.                     */
    if (incircle(m, b, leftvertex, farvertex, rightvertex, nearvertex) <=
        0.0) {
      return;
    }
    /* Not locally Delaunay; go on to an edge flip. */
  }        /* else fartri is inverted; remove it from the stack by flipping. */
  flip(m, b, &neartri);
  lprevself(*fixuptri);    /* Restore the origin of fixuptri after the flip. */
  /* Recursively process the two triangles that result from the flip. */
  delaunayfixup(m, b, fixuptri, leftside);
  delaunayfixup(m, b, &fartri, leftside);
}

/*****************************************************************************/
/*                                                                           */
/*  constrainededge()   Force a segment into a constrained Delaunay          */
/*                      triangulation by deleting the triangles it           */
/*                      intersects, and triangulating the polygons that      */
/*                      form on each side of it.                             */
/*                                                                           */
/*  Generates a single subsegment connecting `endpoint1' to `endpoint2'.     */
/*  The triangle `starttri' has `endpoint1' as its origin.  `newmark' is the */
/*  boundary marker of the segment.                                          */
/*                                                                           */
/*  To insert a segment, every triangle whose interior intersects the        */
/*  segment is deleted.  The union of these deleted triangles is a polygon   */
/*  (which is not necessarily monotone, but is close enough), which is       */
/*  divided into two polygons by the new segment.  This routine's task is    */
/*  to generate the Delaunay triangulation of these two polygons.            */
/*                                                                           */
/*  You might think of this routine's behavior as a two-step process.  The   */
/*  first step is to walk from endpoint1 to endpoint2, flipping each edge    */
/*  encountered.  This step creates a fan of edges connected to endpoint1,   */
/*  including the desired edge to endpoint2.  The second step enforces the   */
/*  Delaunay condition on each side of the segment in an incremental manner: */
/*  proceeding along the polygon from endpoint1 to endpoint2 (this is done   */
/*  independently on each side of the segment), each vertex is "enforced"    */
/*  as if it had just been inserted, but affecting only the previous         */
/*  vertices.  The result is the same as if the vertices had been inserted   */
/*  in the order they appear on the polygon, so the result is Delaunay.      */
/*                                                                           */
/*  In truth, constrainededge() interleaves these two steps.  The procedure  */
/*  walks from endpoint1 to endpoint2, and each time an edge is encountered  */
/*  and flipped, the newly exposed vertex (at the far end of the flipped     */
/*  edge) is "enforced" upon the previously flipped edges, usually affecting */
/*  only one side of the polygon (depending upon which side of the segment   */
/*  the vertex falls on).                                                    */
/*                                                                           */
/*  The algorithm is complicated by the need to handle polygons that are not */
/*  convex.  Although the polygon is not necessarily monotone, it can be     */
/*  triangulated in a manner similar to the stack-based algorithms for       */
/*  monotone polygons.  For each reflex vertex (local concavity) of the      */
/*  polygon, there will be an inverted triangle formed by one of the edge    */
/*  flips.  (An inverted triangle is one with negative area - that is, its   */
/*  vertices are arranged in clockwise order - and is best thought of as a   */
/*  wrinkle in the fabric of the mesh.)  Each inverted triangle can be       */
/*  thought of as a reflex vertex pushed on the stack, waiting to be fixed   */
/*  later.                                                                   */
/*                                                                           */
/*  A reflex vertex is popped from the stack when a vertex is inserted that  */
/*  is visible to the reflex vertex.  (However, if the vertex behind the     */
/*  reflex vertex is not visible to the reflex vertex, a new inverted        */
/*  triangle will take its place on the stack.)  These details are handled   */
/*  by the delaunayfixup() routine above.                                    */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void constrainededge(struct mesh *m, struct behavior *b,
                     struct otri *starttri, vertex endpoint2, int newmark)
#else /* not ANSI_DECLARATORS */
void constrainededge(m, b, starttri, endpoint2, newmark)
struct mesh *m;
struct behavior *b;
struct otri *starttri;
vertex endpoint2;
int newmark;
#endif /* not ANSI_DECLARATORS */

{
  struct otri fixuptri, fixuptri2;
  struct osub crosssubseg;
  vertex endpoint1;
  vertex farvertex;
  REAL area;
  int collision;
  int done;
  triangle ptr;             /* Temporary variable used by sym() and oprev(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  org(*starttri, endpoint1);
  lnext(*starttri, fixuptri);
  flip(m, b, &fixuptri);
  /* `collision' indicates whether we have found a vertex directly */
  /*   between endpoint1 and endpoint2.                            */
  collision = 0;
  done = 0;
  do {
    org(fixuptri, farvertex);
    /* `farvertex' is the extreme point of the polygon we are "digging" */
    /*   to get from endpoint1 to endpoint2.                           */
    if ((farvertex[0] == endpoint2[0]) && (farvertex[1] == endpoint2[1])) {
      oprev(fixuptri, fixuptri2);
      /* Enforce the Delaunay condition around endpoint2. */
      delaunayfixup(m, b, &fixuptri, 0);
      delaunayfixup(m, b, &fixuptri2, 1);
      done = 1;
    } else {
      /* Check whether farvertex is to the left or right of the segment */
      /*   being inserted, to decide which edge of fixuptri to dig      */
      /*   through next.                                                */
      area = counterclockwise(m, b, endpoint1, endpoint2, farvertex);
      if (area == 0.0) {
        /* We've collided with a vertex between endpoint1 and endpoint2. */
        collision = 1;
        oprev(fixuptri, fixuptri2);
        /* Enforce the Delaunay condition around farvertex. */
        delaunayfixup(m, b, &fixuptri, 0);
        delaunayfixup(m, b, &fixuptri2, 1);
        done = 1;
      } else {
        if (area > 0.0) {        /* farvertex is to the left of the segment. */
          oprev(fixuptri, fixuptri2);
          /* Enforce the Delaunay condition around farvertex, on the */
          /*   left side of the segment only.                        */
          delaunayfixup(m, b, &fixuptri2, 1);
          /* Flip the edge that crosses the segment.  After the edge is */
          /*   flipped, one of its endpoints is the fan vertex, and the */
          /*   destination of fixuptri is the fan vertex.               */
          lprevself(fixuptri);
        } else {                /* farvertex is to the right of the segment. */
          delaunayfixup(m, b, &fixuptri, 0);
          /* Flip the edge that crosses the segment.  After the edge is */
          /*   flipped, one of its endpoints is the fan vertex, and the */
          /*   destination of fixuptri is the fan vertex.               */
          oprevself(fixuptri);
        }
        /* Check for two intersecting segments. */
        tspivot(fixuptri, crosssubseg);
        if (crosssubseg.ss == m->dummysub) {
          flip(m, b, &fixuptri);    /* May create inverted triangle at left. */
        } else {
          /* We've collided with a segment between endpoint1 and endpoint2. */
          collision = 1;
          /* Insert a vertex at the intersection. */
          segmentintersection(m, b, &fixuptri, &crosssubseg, endpoint2);
          done = 1;
        }
      }
    }
  } while (!done);
  /* Insert a subsegment to make the segment permanent. */
  insertsubseg(m, b, &fixuptri, newmark);
  /* If there was a collision with an interceding vertex, install another */
  /*   segment connecting that vertex with endpoint2.                     */
  if (collision) {
    /* Insert the remainder of the segment. */
    if (!scoutsegment(m, b, &fixuptri, endpoint2, newmark)) {
      constrainededge(m, b, &fixuptri, endpoint2, newmark);
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*  insertsegment()   Insert a PSLG segment into a triangulation.            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void insertsegment(struct mesh *m, struct behavior *b,
                   vertex endpoint1, vertex endpoint2, int newmark)
#else /* not ANSI_DECLARATORS */
void insertsegment(m, b, endpoint1, endpoint2, newmark)
struct mesh *m;
struct behavior *b;
vertex endpoint1;
vertex endpoint2;
int newmark;
#endif /* not ANSI_DECLARATORS */

{
  struct otri searchtri1, searchtri2;
  triangle encodedtri;
  vertex checkvertex;
  triangle ptr;                         /* Temporary variable used by sym(). */

  if (b->verbose > 1) {
    printf("  Connecting (%.12g, %.12g) to (%.12g, %.12g).\n",
           endpoint1[0], endpoint1[1], endpoint2[0], endpoint2[1]);
  }

  /* Find a triangle whose origin is the segment's first endpoint. */
  checkvertex = (vertex) NULL;
  encodedtri = vertex2tri(endpoint1);
  if (encodedtri != (triangle) NULL) {
    decode(encodedtri, searchtri1);
    org(searchtri1, checkvertex);
  }
  if (checkvertex != endpoint1) {
    /* Find a boundary triangle to search from. */
    searchtri1.tri = m->dummytri;
    searchtri1.orient = 0;
    symself(searchtri1);
    /* Search for the segment's first endpoint by point location. */
    if (locate(m, b, endpoint1, &searchtri1) != ONVERTEX) {
      printf(
        "Internal error in insertsegment():  Unable to locate PSLG vertex\n");
      printf("  (%.12g, %.12g) in triangulation.\n",
             endpoint1[0], endpoint1[1]);
      internalerror();
    }
  }
  /* Remember this triangle to improve subsequent point location. */
  otricopy(searchtri1, m->recenttri);
  /* Scout the beginnings of a path from the first endpoint */
  /*   toward the second.                                   */
  if (scoutsegment(m, b, &searchtri1, endpoint2, newmark)) {
    /* The segment was easily inserted. */
    return;
  }
  /* The first endpoint may have changed if a collision with an intervening */
  /*   vertex on the segment occurred.                                      */
  org(searchtri1, endpoint1);

  /* Find a triangle whose origin is the segment's second endpoint. */
  checkvertex = (vertex) NULL;
  encodedtri = vertex2tri(endpoint2);
  if (encodedtri != (triangle) NULL) {
    decode(encodedtri, searchtri2);
    org(searchtri2, checkvertex);
  }
  if (checkvertex != endpoint2) {
    /* Find a boundary triangle to search from. */
    searchtri2.tri = m->dummytri;
    searchtri2.orient = 0;
    symself(searchtri2);
    /* Search for the segment's second endpoint by point location. */
    if (locate(m, b, endpoint2, &searchtri2) != ONVERTEX) {
      printf(
        "Internal error in insertsegment():  Unable to locate PSLG vertex\n");
      printf("  (%.12g, %.12g) in triangulation.\n",
             endpoint2[0], endpoint2[1]);
      internalerror();
    }
  }
  /* Remember this triangle to improve subsequent point location. */
  otricopy(searchtri2, m->recenttri);
  /* Scout the beginnings of a path from the second endpoint */
  /*   toward the first.                                     */
  if (scoutsegment(m, b, &searchtri2, endpoint1, newmark)) {
    /* The segment was easily inserted. */
    return;
  }
  /* The second endpoint may have changed if a collision with an intervening */
  /*   vertex on the segment occurred.                                       */
  org(searchtri2, endpoint2);

#ifndef REDUCED
#ifndef CDT_ONLY
  if (b->splitseg) {
    /* Insert vertices to force the segment into the triangulation. */
    conformingedge(m, b, endpoint1, endpoint2, newmark);
  } else {
#endif /* not CDT_ONLY */
#endif /* not REDUCED */
    /* Insert the segment directly into the triangulation. */
    constrainededge(m, b, &searchtri1, endpoint2, newmark);
#ifndef REDUCED
#ifndef CDT_ONLY
  }
#endif /* not CDT_ONLY */
#endif /* not REDUCED */
}

/*****************************************************************************/
/*                                                                           */
/*  markhull()   Cover the convex hull of a triangulation with subsegments.  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void markhull(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void markhull(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri hulltri;
  struct otri nexttri;
  struct otri starttri;
  triangle ptr;             /* Temporary variable used by sym() and oprev(). */

  /* Find a triangle handle on the hull. */
  hulltri.tri = m->dummytri;
  hulltri.orient = 0;
  symself(hulltri);
  /* Remember where we started so we know when to stop. */
  otricopy(hulltri, starttri);
  /* Go once counterclockwise around the convex hull. */
  do {
    /* Create a subsegment if there isn't already one here. */
    insertsubseg(m, b, &hulltri, 1);
    /* To find the next hull edge, go clockwise around the next vertex. */
    lnextself(hulltri);
    oprev(hulltri, nexttri);
    while (nexttri.tri != m->dummytri) {
      otricopy(nexttri, hulltri);
      oprev(hulltri, nexttri);
    }
  } while (!otriequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  formskeleton()   Create the segments of a triangulation, including PSLG  */
/*                   segments and edges on the convex hull.                  */
/*                                                                           */
/*  The PSLG segments are read from a .poly file.  The return value is the   */
/*  number of segments in the file.                                          */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void formskeleton(struct mesh *m, struct behavior *b, int *segmentlist,
                  int *segmentmarkerlist, int numberofsegments)
#else /* not ANSI_DECLARATORS */
void formskeleton(m, b, segmentlist, segmentmarkerlist, numberofsegments)
struct mesh *m;
struct behavior *b;
int *segmentlist;
int *segmentmarkerlist;
int numberofsegments;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void formskeleton(struct mesh *m, struct behavior *b,
                  FILE *polyfile, char *polyfilename)
#else /* not ANSI_DECLARATORS */
void formskeleton(m, b, polyfile, polyfilename)
struct mesh *m;
struct behavior *b;
FILE *polyfile;
char *polyfilename;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  char polyfilename[6];
  int index;
#else /* not TRILIBRARY */
  char inputline[INPUTLINESIZE];
  char *stringptr;
#endif /* not TRILIBRARY */
  vertex endpoint1, endpoint2;
  int segmentmarkers;
  int end1, end2;
  int boundmarker;
  int i;

  if (b->poly) {
    if (!b->quiet) {
      printf("Recovering segments in Delaunay triangulation.\n");
    }
#ifdef TRILIBRARY
    strcpy(polyfilename, "input");
    m->insegments = numberofsegments;
    segmentmarkers = segmentmarkerlist != (int *) NULL;
    index = 0;
#else /* not TRILIBRARY */
    /* Read the segments from a .poly file. */
    /* Read number of segments and number of boundary markers. */
    stringptr = readline(inputline, polyfile, polyfilename);
    m->insegments = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      segmentmarkers = 0;
    } else {
      segmentmarkers = (int) strtol(stringptr, &stringptr, 0);
    }
#endif /* not TRILIBRARY */
    /* If the input vertices are collinear, there is no triangulation, */
    /*   so don't try to insert segments.                              */
    if (m->triangles.items == 0) {
      return;
    }

    /* If segments are to be inserted, compute a mapping */
    /*   from vertices to triangles.                     */
    if (m->insegments > 0) {
      makevertexmap(m, b);
      if (b->verbose) {
        printf("  Recovering PSLG segments.\n");
      }
    }

    boundmarker = 0;
    /* Read and insert the segments. */
    for (i = 0; i < m->insegments; i++) {
#ifdef TRILIBRARY
      end1 = segmentlist[index++];
      end2 = segmentlist[index++];
      if (segmentmarkers) {
        boundmarker = segmentmarkerlist[i];
      }
#else /* not TRILIBRARY */
      stringptr = readline(inputline, polyfile, b->inpolyfilename);
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Segment %d has no endpoints in %s.\n",
               b->firstnumber + i, polyfilename);
        triexit(1);
      } else {
        end1 = (int) strtol(stringptr, &stringptr, 0);
      }
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Segment %d is missing its second endpoint in %s.\n",
               b->firstnumber + i, polyfilename);
        triexit(1);
      } else {
        end2 = (int) strtol(stringptr, &stringptr, 0);
      }
      if (segmentmarkers) {
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          boundmarker = 0;
        } else {
          boundmarker = (int) strtol(stringptr, &stringptr, 0);
        }
      }
#endif /* not TRILIBRARY */
      if ((end1 < b->firstnumber) ||
          (end1 >= b->firstnumber + m->invertices)) {
        if (!b->quiet) {
          printf("Warning:  Invalid first endpoint of segment %d in %s.\n",
                 b->firstnumber + i, polyfilename);
        }
      } else if ((end2 < b->firstnumber) ||
                 (end2 >= b->firstnumber + m->invertices)) {
        if (!b->quiet) {
          printf("Warning:  Invalid second endpoint of segment %d in %s.\n",
                 b->firstnumber + i, polyfilename);
        }
      } else {
        /* Find the vertices numbered `end1' and `end2'. */
        endpoint1 = getvertex(m, b, end1);
        endpoint2 = getvertex(m, b, end2);
        if ((endpoint1[0] == endpoint2[0]) && (endpoint1[1] == endpoint2[1])) {
          if (!b->quiet) {
            printf("Warning:  Endpoints of segment %d are coincident in %s.\n",
                   b->firstnumber + i, polyfilename);
          }
        } else {
          insertsegment(m, b, endpoint1, endpoint2, boundmarker);
        }
      }
    }
  } else {
    m->insegments = 0;
  }
  if (b->convex || !b->poly) {
    /* Enclose the convex hull with subsegments. */
    if (b->verbose) {
      printf("  Enclosing convex hull with segments.\n");
    }
    markhull(m, b);
  }
}

/**                                                                         **/
/**                                                                         **/
/********* Segment insertion ends here                               *********/

/********* Carving out holes and concavities begins here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  infecthull()   Virally infect all of the triangles of the convex hull    */
/*                 that are not protected by subsegments.  Where there are   */
/*                 subsegments, set boundary markers as appropriate.         */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void infecthull(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void infecthull(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri hulltri;
  struct otri nexttri;
  struct otri starttri;
  struct osub hullsubseg;
  triangle **deadtriangle;
  vertex horg, hdest;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (b->verbose) {
    printf("  Marking concavities (external triangles) for elimination.\n");
  }
  /* Find a triangle handle on the hull. */
  hulltri.tri = m->dummytri;
  hulltri.orient = 0;
  symself(hulltri);
  /* Remember where we started so we know when to stop. */
  otricopy(hulltri, starttri);
  /* Go once counterclockwise around the convex hull. */
  do {
    /* Ignore triangles that are already infected. */
    if (!infected(hulltri)) {
      /* Is the triangle protected by a subsegment? */
      tspivot(hulltri, hullsubseg);
      if (hullsubseg.ss == m->dummysub) {
        /* The triangle is not protected; infect it. */
        if (!infected(hulltri)) {
          infect(hulltri);
          deadtriangle = (triangle **) poolalloc(&m->viri);
          *deadtriangle = hulltri.tri;
        }
      } else {
        /* The triangle is protected; set boundary markers if appropriate. */
        if (mark(hullsubseg) == 0) {
          setmark(hullsubseg, 1);
          org(hulltri, horg);
          dest(hulltri, hdest);
          if (vertexmark(horg) == 0) {
            setvertexmark(horg, 1);
          }
          if (vertexmark(hdest) == 0) {
            setvertexmark(hdest, 1);
          }
        }
      }
    }
    /* To find the next hull edge, go clockwise around the next vertex. */
    lnextself(hulltri);
    oprev(hulltri, nexttri);
    while (nexttri.tri != m->dummytri) {
      otricopy(nexttri, hulltri);
      oprev(hulltri, nexttri);
    }
  } while (!otriequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  plague()   Spread the virus from all infected triangles to any neighbors */
/*             not protected by subsegments.  Delete all infected triangles. */
/*                                                                           */
/*  This is the procedure that actually creates holes and concavities.       */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase identifies all   */
/*  the triangles that will die, and marks them as infected.  They are       */
/*  marked to ensure that each triangle is added to the virus pool only      */
/*  once, so the procedure will terminate.                                   */
/*                                                                           */
/*  The second phase actually eliminates the infected triangles.  It also    */
/*  eliminates orphaned vertices.                                            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void plague(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void plague(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri testtri;
  struct otri neighbor;
  triangle **virusloop;
  triangle **deadtriangle;
  struct osub neighborsubseg;
  vertex testvertex;
  vertex norg, ndest;
  vertex deadorg, deaddest, deadapex;
  int killorg;
  triangle ptr;             /* Temporary variable used by sym() and onext(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (b->verbose) {
    printf("  Marking neighbors of marked triangles.\n");
  }
  /* Loop through all the infected triangles, spreading the virus to */
  /*   their neighbors, then to their neighbors' neighbors.          */
  traversalinit(&m->viri);
  virusloop = (triangle **) traverse(&m->viri);
  while (virusloop != (triangle **) NULL) {
    testtri.tri = *virusloop;
    /* A triangle is marked as infected by messing with one of its pointers */
    /*   to subsegments, setting it to an illegal value.  Hence, we have to */
    /*   temporarily uninfect this triangle so that we can examine its      */
    /*   adjacent subsegments.                                              */
    uninfect(testtri);
    if (b->verbose > 2) {
      /* Assign the triangle an orientation for convenience in */
      /*   checking its vertices.                              */
      testtri.orient = 0;
      org(testtri, deadorg);
      dest(testtri, deaddest);
      apex(testtri, deadapex);
      printf("    Checking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             deadorg[0], deadorg[1], deaddest[0], deaddest[1],
             deadapex[0], deadapex[1]);
    }
    /* Check each of the triangle's three neighbors. */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      /* Find the neighbor. */
      sym(testtri, neighbor);
      /* Check for a subsegment between the triangle and its neighbor. */
      tspivot(testtri, neighborsubseg);
      /* Check if the neighbor is nonexistent or already infected. */
      if ((neighbor.tri == m->dummytri) || infected(neighbor)) {
        if (neighborsubseg.ss != m->dummysub) {
          /* There is a subsegment separating the triangle from its      */
          /*   neighbor, but both triangles are dying, so the subsegment */
          /*   dies too.                                                 */
          subsegdealloc(m, neighborsubseg.ss);
          if (neighbor.tri != m->dummytri) {
            /* Make sure the subsegment doesn't get deallocated again */
            /*   later when the infected neighbor is visited.         */
            uninfect(neighbor);
            tsdissolve(neighbor);
            infect(neighbor);
          }
        }
      } else {                   /* The neighbor exists and is not infected. */
        if (neighborsubseg.ss == m->dummysub) {
          /* There is no subsegment protecting the neighbor, so */
          /*   the neighbor becomes infected.                   */
          if (b->verbose > 2) {
            org(neighbor, deadorg);
            dest(neighbor, deaddest);
            apex(neighbor, deadapex);
            printf(
              "    Marking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
                   deadorg[0], deadorg[1], deaddest[0], deaddest[1],
                   deadapex[0], deadapex[1]);
          }
          infect(neighbor);
          /* Ensure that the neighbor's neighbors will be infected. */
          deadtriangle = (triangle **) poolalloc(&m->viri);
          *deadtriangle = neighbor.tri;
        } else {               /* The neighbor is protected by a subsegment. */
          /* Remove this triangle from the subsegment. */
          stdissolve(neighborsubseg);
          /* The subsegment becomes a boundary.  Set markers accordingly. */
          if (mark(neighborsubseg) == 0) {
            setmark(neighborsubseg, 1);
          }
          org(neighbor, norg);
          dest(neighbor, ndest);
          if (vertexmark(norg) == 0) {
            setvertexmark(norg, 1);
          }
          if (vertexmark(ndest) == 0) {
            setvertexmark(ndest, 1);
          }
        }
      }
    }
    /* Remark the triangle as infected, so it doesn't get added to the */
    /*   virus pool again.                                             */
    infect(testtri);
    virusloop = (triangle **) traverse(&m->viri);
  }

  if (b->verbose) {
    printf("  Deleting marked triangles.\n");
  }

  traversalinit(&m->viri);
  virusloop = (triangle **) traverse(&m->viri);
  while (virusloop != (triangle **) NULL) {
    testtri.tri = *virusloop;

    /* Check each of the three corners of the triangle for elimination. */
    /*   This is done by walking around each vertex, checking if it is  */
    /*   still connected to at least one live triangle.                 */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      org(testtri, testvertex);
      /* Check if the vertex has already been tested. */
      if (testvertex != (vertex) NULL) {
        killorg = 1;
        /* Mark the corner of the triangle as having been tested. */
        setorg(testtri, NULL);
        /* Walk counterclockwise about the vertex. */
        onext(testtri, neighbor);
        /* Stop upon reaching a boundary or the starting triangle. */
        while ((neighbor.tri != m->dummytri) &&
               (!otriequal(neighbor, testtri))) {
          if (infected(neighbor)) {
            /* Mark the corner of this triangle as having been tested. */
            setorg(neighbor, NULL);
          } else {
            /* A live triangle.  The vertex survives. */
            killorg = 0;
          }
          /* Walk counterclockwise about the vertex. */
          onextself(neighbor);
        }
        /* If we reached a boundary, we must walk clockwise as well. */
        if (neighbor.tri == m->dummytri) {
          /* Walk clockwise about the vertex. */
          oprev(testtri, neighbor);
          /* Stop upon reaching a boundary. */
          while (neighbor.tri != m->dummytri) {
            if (infected(neighbor)) {
            /* Mark the corner of this triangle as having been tested. */
              setorg(neighbor, NULL);
            } else {
              /* A live triangle.  The vertex survives. */
              killorg = 0;
            }
            /* Walk clockwise about the vertex. */
            oprevself(neighbor);
          }
        }
        if (killorg) {
          if (b->verbose > 1) {
            printf("    Deleting vertex (%.12g, %.12g)\n",
                   testvertex[0], testvertex[1]);
          }
          setvertextype(testvertex, UNDEADVERTEX);
          m->undeads++;
        }
      }
    }

    /* Record changes in the number of boundary edges, and disconnect */
    /*   dead triangles from their neighbors.                         */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      sym(testtri, neighbor);
      if (neighbor.tri == m->dummytri) {
        /* There is no neighboring triangle on this edge, so this edge    */
        /*   is a boundary edge.  This triangle is being deleted, so this */
        /*   boundary edge is deleted.                                    */
        m->hullsize--;
      } else {
        /* Disconnect the triangle from its neighbor. */
        dissolve(neighbor);
        /* There is a neighboring triangle on this edge, so this edge */
        /*   becomes a boundary edge when this triangle is deleted.   */
        m->hullsize++;
      }
    }
    /* Return the dead triangle to the pool of triangles. */
    triangledealloc(m, testtri.tri);
    virusloop = (triangle **) traverse(&m->viri);
  }
  /* Empty the virus pool. */
  poolrestart(&m->viri);
}

/*****************************************************************************/
/*                                                                           */
/*  regionplague()   Spread regional attributes and/or area constraints      */
/*                   (from a .poly file) throughout the mesh.                */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase spreads an       */
/*  attribute and/or an area constraint through a (segment-bounded) region.  */
/*  The triangles are marked to ensure that each triangle is added to the    */
/*  virus pool only once, so the procedure will terminate.                   */
/*                                                                           */
/*  The second phase uninfects all infected triangles, returning them to     */
/*  normal.                                                                  */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void regionplague(struct mesh *m, struct behavior *b,
                  REAL attribute, REAL area)
#else /* not ANSI_DECLARATORS */
void regionplague(m, b, attribute, area)
struct mesh *m;
struct behavior *b;
REAL attribute;
REAL area;
#endif /* not ANSI_DECLARATORS */

{
  struct otri testtri;
  struct otri neighbor;
  triangle **virusloop;
  triangle **regiontri;
  struct osub neighborsubseg;
  vertex regionorg, regiondest, regionapex;
  triangle ptr;             /* Temporary variable used by sym() and onext(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (b->verbose > 1) {
    printf("  Marking neighbors of marked triangles.\n");
  }
  /* Loop through all the infected triangles, spreading the attribute      */
  /*   and/or area constraint to their neighbors, then to their neighbors' */
  /*   neighbors.                                                          */
  traversalinit(&m->viri);
  virusloop = (triangle **) traverse(&m->viri);
  while (virusloop != (triangle **) NULL) {
    testtri.tri = *virusloop;
    /* A triangle is marked as infected by messing with one of its pointers */
    /*   to subsegments, setting it to an illegal value.  Hence, we have to */
    /*   temporarily uninfect this triangle so that we can examine its      */
    /*   adjacent subsegments.                                              */
    uninfect(testtri);
    if (b->regionattrib) {
      /* Set an attribute. */
      setelemattribute(testtri, m->eextras, attribute);
    }
    if (b->vararea) {
      /* Set an area constraint. */
      setareabound(testtri, area);
    }
    if (b->verbose > 2) {
      /* Assign the triangle an orientation for convenience in */
      /*   checking its vertices.                              */
      testtri.orient = 0;
      org(testtri, regionorg);
      dest(testtri, regiondest);
      apex(testtri, regionapex);
      printf("    Checking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
             regionorg[0], regionorg[1], regiondest[0], regiondest[1],
             regionapex[0], regionapex[1]);
    }
    /* Check each of the triangle's three neighbors. */
    for (testtri.orient = 0; testtri.orient < 3; testtri.orient++) {
      /* Find the neighbor. */
      sym(testtri, neighbor);
      /* Check for a subsegment between the triangle and its neighbor. */
      tspivot(testtri, neighborsubseg);
      /* Make sure the neighbor exists, is not already infected, and */
      /*   isn't protected by a subsegment.                          */
      if ((neighbor.tri != m->dummytri) && !infected(neighbor)
          && (neighborsubseg.ss == m->dummysub)) {
        if (b->verbose > 2) {
          org(neighbor, regionorg);
          dest(neighbor, regiondest);
          apex(neighbor, regionapex);
          printf("    Marking (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
                 regionorg[0], regionorg[1], regiondest[0], regiondest[1],
                 regionapex[0], regionapex[1]);
        }
        /* Infect the neighbor. */
        infect(neighbor);
        /* Ensure that the neighbor's neighbors will be infected. */
        regiontri = (triangle **) poolalloc(&m->viri);
        *regiontri = neighbor.tri;
      }
    }
    /* Remark the triangle as infected, so it doesn't get added to the */
    /*   virus pool again.                                             */
    infect(testtri);
    virusloop = (triangle **) traverse(&m->viri);
  }

  /* Uninfect all triangles. */
  if (b->verbose > 1) {
    printf("  Unmarking marked triangles.\n");
  }
  traversalinit(&m->viri);
  virusloop = (triangle **) traverse(&m->viri);
  while (virusloop != (triangle **) NULL) {
    testtri.tri = *virusloop;
    uninfect(testtri);
    virusloop = (triangle **) traverse(&m->viri);
  }
  /* Empty the virus pool. */
  poolrestart(&m->viri);
}

/*****************************************************************************/
/*                                                                           */
/*  carveholes()   Find the holes and infect them.  Find the area            */
/*                 constraints and infect them.  Infect the convex hull.     */
/*                 Spread the infection and kill triangles.  Spread the      */
/*                 area constraints.                                         */
/*                                                                           */
/*  This routine mainly calls other routines to carry out all these          */
/*  functions.                                                               */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void carveholes(struct mesh *m, struct behavior *b, REAL *holelist, int holes,
                REAL *regionlist, int regions)
#else /* not ANSI_DECLARATORS */
void carveholes(m, b, holelist, holes, regionlist, regions)
struct mesh *m;
struct behavior *b;
REAL *holelist;
int holes;
REAL *regionlist;
int regions;
#endif /* not ANSI_DECLARATORS */

{
  struct otri searchtri;
  struct otri triangleloop;
  struct otri *regiontris;
  triangle **holetri;
  triangle **regiontri;
  vertex searchorg, searchdest;
  enum locateresult intersect;
  int i;
  triangle ptr;                         /* Temporary variable used by sym(). */

  if (!(b->quiet || (b->noholes && b->convex))) {
    printf("Removing unwanted triangles.\n");
    if (b->verbose && (holes > 0)) {
      printf("  Marking holes for elimination.\n");
    }
  }

  if (regions > 0) {
    /* Allocate storage for the triangles in which region points fall. */
    regiontris = (struct otri *) trimalloc(regions *
                                           (int) sizeof(struct otri));
  } else {
    regiontris = (struct otri *) NULL;
  }

  if (((holes > 0) && !b->noholes) || !b->convex || (regions > 0)) {
    /* Initialize a pool of viri to be used for holes, concavities, */
    /*   regional attributes, and/or regional area constraints.     */
    poolinit(&m->viri, sizeof(triangle *), VIRUSPERBLOCK, VIRUSPERBLOCK, 0);
  }

  if (!b->convex) {
    /* Mark as infected any unprotected triangles on the boundary. */
    /*   This is one way by which concavities are created.         */
    infecthull(m, b);
  }

  if ((holes > 0) && !b->noholes) {
    /* Infect each triangle in which a hole lies. */
    for (i = 0; i < 2 * holes; i += 2) {
      /* Ignore holes that aren't within the bounds of the mesh. */
      if ((holelist[i] >= m->xmin) && (holelist[i] <= m->xmax)
          && (holelist[i + 1] >= m->ymin) && (holelist[i + 1] <= m->ymax)) {
        /* Start searching from some triangle on the outer boundary. */
        searchtri.tri = m->dummytri;
        searchtri.orient = 0;
        symself(searchtri);
        /* Ensure that the hole is to the left of this boundary edge; */
        /*   otherwise, locate() will falsely report that the hole    */
        /*   falls within the starting triangle.                      */
        org(searchtri, searchorg);
        dest(searchtri, searchdest);
        if (counterclockwise(m, b, searchorg, searchdest, &holelist[i]) >
            0.0) {
          /* Find a triangle that contains the hole. */
          intersect = locate(m, b, &holelist[i], &searchtri);
          if ((intersect != OUTSIDE) && (!infected(searchtri))) {
            /* Infect the triangle.  This is done by marking the triangle  */
            /*   as infected and including the triangle in the virus pool. */
            infect(searchtri);
            holetri = (triangle **) poolalloc(&m->viri);
            *holetri = searchtri.tri;
          }
        }
      }
    }
  }

  /* Now, we have to find all the regions BEFORE we carve the holes, because */
  /*   locate() won't work when the triangulation is no longer convex.       */
  /*   (Incidentally, this is the reason why regional attributes and area    */
  /*   constraints can't be used when refining a preexisting mesh, which     */
  /*   might not be convex; they can only be used with a freshly             */
  /*   triangulated PSLG.)                                                   */
  if (regions > 0) {
    /* Find the starting triangle for each region. */
    for (i = 0; i < regions; i++) {
      regiontris[i].tri = m->dummytri;
      /* Ignore region points that aren't within the bounds of the mesh. */
      if ((regionlist[4 * i] >= m->xmin) && (regionlist[4 * i] <= m->xmax) &&
          (regionlist[4 * i + 1] >= m->ymin) &&
          (regionlist[4 * i + 1] <= m->ymax)) {
        /* Start searching from some triangle on the outer boundary. */
        searchtri.tri = m->dummytri;
        searchtri.orient = 0;
        symself(searchtri);
        /* Ensure that the region point is to the left of this boundary */
        /*   edge; otherwise, locate() will falsely report that the     */
        /*   region point falls within the starting triangle.           */
        org(searchtri, searchorg);
        dest(searchtri, searchdest);
        if (counterclockwise(m, b, searchorg, searchdest, &regionlist[4 * i]) >
            0.0) {
          /* Find a triangle that contains the region point. */
          intersect = locate(m, b, &regionlist[4 * i], &searchtri);
          if ((intersect != OUTSIDE) && (!infected(searchtri))) {
            /* Record the triangle for processing after the */
            /*   holes have been carved.                    */
            otricopy(searchtri, regiontris[i]);
          }
        }
      }
    }
  }

  if (m->viri.items > 0) {
    /* Carve the holes and concavities. */
    plague(m, b);
  }
  /* The virus pool should be empty now. */

  if (regions > 0) {
    if (!b->quiet) {
      if (b->regionattrib) {
        if (b->vararea) {
          printf("Spreading regional attributes and area constraints.\n");
        } else {
          printf("Spreading regional attributes.\n");
        }
      } else { 
        printf("Spreading regional area constraints.\n");
      }
    }
    if (b->regionattrib && !b->refine) {
      /* Assign every triangle a regional attribute of zero. */
      traversalinit(&m->triangles);
      triangleloop.orient = 0;
      triangleloop.tri = triangletraverse(m);
      while (triangleloop.tri != (triangle *) NULL) {
        setelemattribute(triangleloop, m->eextras, 0.0);
        triangleloop.tri = triangletraverse(m);
      }
    }
    for (i = 0; i < regions; i++) {
      if (regiontris[i].tri != m->dummytri) {
        /* Make sure the triangle under consideration still exists. */
        /*   It may have been eaten by the virus.                   */
        if (!deadtri(regiontris[i].tri)) {
          /* Put one triangle in the virus pool. */
          infect(regiontris[i]);
          regiontri = (triangle **) poolalloc(&m->viri);
          *regiontri = regiontris[i].tri;
          /* Apply one region's attribute and/or area constraint. */
          regionplague(m, b, regionlist[4 * i + 2], regionlist[4 * i + 3]);
          /* The virus pool should be empty now. */
        }
      }
    }
    if (b->regionattrib && !b->refine) {
      /* Note the fact that each triangle has an additional attribute. */
      m->eextras++;
    }
  }

  /* Free up memory. */
  if (((holes > 0) && !b->noholes) || !b->convex || (regions > 0)) {
    pooldeinit(&m->viri);
  }
  if (regions > 0) {
    trifree((VOID *) regiontris);
  }
}

/**                                                                         **/
/**                                                                         **/
/********* Carving out holes and concavities ends here               *********/

/********* Mesh quality maintenance begins here                      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  tallyencs()   Traverse the entire list of subsegments, and check each    */
/*                to see if it is encroached.  If so, add it to the list.    */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void tallyencs(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void tallyencs(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct osub subsegloop;
  int dummy;

  traversalinit(&m->subsegs);
  subsegloop.ssorient = 0;
  subsegloop.ss = subsegtraverse(m);
  while (subsegloop.ss != (subseg *) NULL) {
    /* If the segment is encroached, add it to the list. */
    dummy = checkseg4encroach(m, b, &subsegloop);
    subsegloop.ss = subsegtraverse(m);
  }
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  precisionerror()  Print an error message for precision problems.         */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

void precisionerror()
{
  printf("Try increasing the area criterion and/or reducing the minimum\n");
  printf("  allowable angle so that tiny triangles are not created.\n");
#ifdef SINGLE
  printf("Alternatively, try recompiling me with double precision\n");
  printf("  arithmetic (by removing \"#define SINGLE\" from the\n");
  printf("  source file or \"-DSINGLE\" from the makefile).\n");
#endif /* SINGLE */
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  splitencsegs()   Split all the encroached subsegments.                   */
/*                                                                           */
/*  Each encroached subsegment is repaired by splitting it - inserting a     */
/*  vertex at or near its midpoint.  Newly inserted vertices may encroach    */
/*  upon other subsegments; these are also repaired.                         */
/*                                                                           */
/*  `triflaws' is a flag that specifies whether one should take note of new  */
/*  bad triangles that result from inserting vertices to repair encroached   */
/*  subsegments.                                                             */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void splitencsegs(struct mesh *m, struct behavior *b, int triflaws)
#else /* not ANSI_DECLARATORS */
void splitencsegs(m, b, triflaws)
struct mesh *m;
struct behavior *b;
int triflaws;
#endif /* not ANSI_DECLARATORS */

{
  struct otri enctri;
  struct otri testtri;
  struct osub testsh;
  struct osub currentenc;
  struct badsubseg *encloop;
  vertex eorg, edest, eapex;
  vertex newvertex;
  enum insertvertexresult success;
  REAL segmentlength, nearestpoweroftwo;
  REAL split;
  REAL multiplier, divisor;
  int acuteorg, acuteorg2, acutedest, acutedest2;
  int dummy;
  int i;
  triangle ptr;                     /* Temporary variable used by stpivot(). */
  subseg sptr;                        /* Temporary variable used by snext(). */

  /* Note that steinerleft == -1 if an unlimited number */
  /*   of Steiner points is allowed.                    */
  while ((m->badsubsegs.items > 0) && (m->steinerleft != 0)) {
    traversalinit(&m->badsubsegs);
    encloop = badsubsegtraverse(m);
    while ((encloop != (struct badsubseg *) NULL) && (m->steinerleft != 0)) {
      sdecode(encloop->encsubseg, currentenc);
      sorg(currentenc, eorg);
      sdest(currentenc, edest);
      /* Make sure that this segment is still the same segment it was   */
      /*   when it was determined to be encroached.  If the segment was */
      /*   enqueued multiple times (because several newly inserted      */
      /*   vertices encroached it), it may have already been split.     */
      if (!deadsubseg(currentenc.ss) &&
          (eorg == encloop->subsegorg) && (edest == encloop->subsegdest)) {
        /* To decide where to split a segment, we need to know if the   */
        /*   segment shares an endpoint with an adjacent segment.       */
        /*   The concern is that, if we simply split every encroached   */
        /*   segment in its center, two adjacent segments with a small  */
        /*   angle between them might lead to an infinite loop; each    */
        /*   vertex added to split one segment will encroach upon the   */
        /*   other segment, which must then be split with a vertex that */
        /*   will encroach upon the first segment, and so on forever.   */
        /* To avoid this, imagine a set of concentric circles, whose    */
        /*   radii are powers of two, about each segment endpoint.      */
        /*   These concentric circles determine where the segment is    */
        /*   split.  (If both endpoints are shared with adjacent        */
        /*   segments, split the segment in the middle, and apply the   */
        /*   concentric circles for later splittings.)                  */

        /* Is the origin shared with another segment? */
        stpivot(currentenc, enctri);
        lnext(enctri, testtri);
        tspivot(testtri, testsh);
        acuteorg = testsh.ss != m->dummysub;
        /* Is the destination shared with another segment? */
        lnextself(testtri);
        tspivot(testtri, testsh);
        acutedest = testsh.ss != m->dummysub;

        /* If we're using Chew's algorithm (rather than Ruppert's) */
        /*   to define encroachment, delete free vertices from the */
        /*   subsegment's diametral circle.                        */
        if (!b->conformdel && !acuteorg && !acutedest) {
          apex(enctri, eapex);
          while ((vertextype(eapex) == FREEVERTEX) &&
                 ((eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                  (eorg[1] - eapex[1]) * (edest[1] - eapex[1]) < 0.0)) {
            deletevertex(m, b, &testtri);
            stpivot(currentenc, enctri);
            apex(enctri, eapex);
            lprev(enctri, testtri);
          }
        }

        /* Now, check the other side of the segment, if there's a triangle */
        /*   there.                                                        */
        sym(enctri, testtri);
        if (testtri.tri != m->dummytri) {
          /* Is the destination shared with another segment? */
          lnextself(testtri);
          tspivot(testtri, testsh);
          acutedest2 = testsh.ss != m->dummysub;
          acutedest = acutedest || acutedest2;
          /* Is the origin shared with another segment? */
          lnextself(testtri);
          tspivot(testtri, testsh);
          acuteorg2 = testsh.ss != m->dummysub;
          acuteorg = acuteorg || acuteorg2;

          /* Delete free vertices from the subsegment's diametral circle. */
          if (!b->conformdel && !acuteorg2 && !acutedest2) {
            org(testtri, eapex);
            while ((vertextype(eapex) == FREEVERTEX) &&
                   ((eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                    (eorg[1] - eapex[1]) * (edest[1] - eapex[1]) < 0.0)) {
              deletevertex(m, b, &testtri);
              sym(enctri, testtri);
              apex(testtri, eapex);
              lprevself(testtri);
            }
          }
        }

        /* Use the concentric circles if exactly one endpoint is shared */
        /*   with another adjacent segment.                             */
        if (acuteorg || acutedest) {
          segmentlength = sqrt((edest[0] - eorg[0]) * (edest[0] - eorg[0]) +
                               (edest[1] - eorg[1]) * (edest[1] - eorg[1]));
          /* Find the power of two that most evenly splits the segment.  */
          /*   The worst case is a 2:1 ratio between subsegment lengths. */
          nearestpoweroftwo = 1.0;
          while (segmentlength > 3.0 * nearestpoweroftwo) {
            nearestpoweroftwo *= 2.0;
          }
          while (segmentlength < 1.5 * nearestpoweroftwo) {
            nearestpoweroftwo *= 0.5;
          }
          /* Where do we split the segment? */
          split = nearestpoweroftwo / segmentlength;
          if (acutedest) {
            split = 1.0 - split;
          }
        } else {
          /* If we're not worried about adjacent segments, split */
          /*   this segment in the middle.                       */
          split = 0.5;
        }

        /* Create the new vertex. */
        newvertex = (vertex) poolalloc(&m->vertices);
        /* Interpolate its coordinate and attributes. */
        for (i = 0; i < 2 + m->nextras; i++) {
          newvertex[i] = eorg[i] + split * (edest[i] - eorg[i]);
        }

        if (!b->noexact) {
          /* Roundoff in the above calculation may yield a `newvertex'   */
          /*   that is not precisely collinear with `eorg' and `edest'.  */
          /*   Improve collinearity by one step of iterative refinement. */
          multiplier = counterclockwise(m, b, eorg, edest, newvertex);
          divisor = ((eorg[0] - edest[0]) * (eorg[0] - edest[0]) +
                     (eorg[1] - edest[1]) * (eorg[1] - edest[1]));
          if ((multiplier != 0.0) && (divisor != 0.0)) {
            multiplier = multiplier / divisor;
            /* Watch out for NANs. */
            if (multiplier == multiplier) {
              newvertex[0] += multiplier * (edest[1] - eorg[1]);
              newvertex[1] += multiplier * (eorg[0] - edest[0]);
            }
          }
        }

        setvertexmark(newvertex, mark(currentenc));
        setvertextype(newvertex, SEGMENTVERTEX);
        if (b->verbose > 1) {
          printf(
  "  Splitting subsegment (%.12g, %.12g) (%.12g, %.12g) at (%.12g, %.12g).\n",
                 eorg[0], eorg[1], edest[0], edest[1],
                 newvertex[0], newvertex[1]);
        }
        /* Check whether the new vertex lies on an endpoint. */
        if (((newvertex[0] == eorg[0]) && (newvertex[1] == eorg[1])) ||
            ((newvertex[0] == edest[0]) && (newvertex[1] == edest[1]))) {
          printf("Error:  Ran out of precision at (%.12g, %.12g).\n",
                 newvertex[0], newvertex[1]);
          printf("I attempted to split a segment to a smaller size than\n");
          printf("  can be accommodated by the finite precision of\n");
          printf("  floating point arithmetic.\n");
          precisionerror();
          triexit(1);
        }
        /* Insert the splitting vertex.  This should always succeed. */
        success = insertvertex(m, b, newvertex, &enctri, &currentenc,
                               1, triflaws);
        if ((success != SUCCESSFULVERTEX) && (success != ENCROACHINGVERTEX)) {
          printf("Internal error in splitencsegs():\n");
          printf("  Failure to split a segment.\n");
          internalerror();
        }
        if (m->steinerleft > 0) {
          m->steinerleft--;
        }
        /* Check the two new subsegments to see if they're encroached. */
        dummy = checkseg4encroach(m, b, &currentenc);
        snextself(currentenc);
        dummy = checkseg4encroach(m, b, &currentenc);
      }

      badsubsegdealloc(m, encloop);
      encloop = badsubsegtraverse(m);
    }
  }
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  tallyfaces()   Test every triangle in the mesh for quality measures.     */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void tallyfaces(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void tallyfaces(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop;

  if (b->verbose) {
    printf("  Making a list of bad triangles.\n");
  }
  traversalinit(&m->triangles);
  triangleloop.orient = 0;
  triangleloop.tri = triangletraverse(m);
  while (triangleloop.tri != (triangle *) NULL) {
    /* If the triangle is bad, enqueue it. */
    testtriangle(m, b, &triangleloop);
    triangleloop.tri = triangletraverse(m);
  }
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  splittriangle()   Inserts a vertex at the circumcenter of a triangle.    */
/*                    Deletes the newly inserted vertex if it encroaches     */
/*                    upon a segment.                                        */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void splittriangle(struct mesh *m, struct behavior *b,
                   struct badtriang *badtri)
#else /* not ANSI_DECLARATORS */
void splittriangle(m, b, badtri)
struct mesh *m;
struct behavior *b;
struct badtriang *badtri;
#endif /* not ANSI_DECLARATORS */

{
  struct otri badotri;
  vertex borg, bdest, bapex;
  vertex newvertex;
  REAL xi, eta;
  enum insertvertexresult success;
  int errorflag;
  int i;

  decode(badtri->poortri, badotri);
  org(badotri, borg);
  dest(badotri, bdest);
  apex(badotri, bapex);
  /* Make sure that this triangle is still the same triangle it was      */
  /*   when it was tested and determined to be of bad quality.           */
  /*   Subsequent transformations may have made it a different triangle. */
  if (!deadtri(badotri.tri) && (borg == badtri->triangorg) &&
      (bdest == badtri->triangdest) && (bapex == badtri->triangapex)) {
    if (b->verbose > 1) {
      printf("  Splitting this triangle at its circumcenter:\n");
      printf("    (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n", borg[0],
             borg[1], bdest[0], bdest[1], bapex[0], bapex[1]);
    }

    errorflag = 0;
    /* Create a new vertex at the triangle's circumcenter. */
    newvertex = (vertex) poolalloc(&m->vertices);
    findcircumcenter(m, b, borg, bdest, bapex, newvertex, &xi, &eta, 1);

    /* Check whether the new vertex lies on a triangle vertex. */
    if (((newvertex[0] == borg[0]) && (newvertex[1] == borg[1])) ||
        ((newvertex[0] == bdest[0]) && (newvertex[1] == bdest[1])) ||
        ((newvertex[0] == bapex[0]) && (newvertex[1] == bapex[1]))) {
      if (!b->quiet) {
        printf(
             "Warning:  New vertex (%.12g, %.12g) falls on existing vertex.\n",
               newvertex[0], newvertex[1]);
        errorflag = 1;
      }
      vertexdealloc(m, newvertex);
    } else {
      for (i = 2; i < 2 + m->nextras; i++) {
        /* Interpolate the vertex attributes at the circumcenter. */
        newvertex[i] = borg[i] + xi * (bdest[i] - borg[i])
                              + eta * (bapex[i] - borg[i]);
      }
      /* The new vertex must be in the interior, and therefore is a */
      /*   free vertex with a marker of zero.                       */
      setvertexmark(newvertex, 0);
      setvertextype(newvertex, FREEVERTEX);

      /* Ensure that the handle `badotri' does not represent the longest  */
      /*   edge of the triangle.  This ensures that the circumcenter must */
      /*   fall to the left of this edge, so point location will work.    */
      /*   (If the angle org-apex-dest exceeds 90 degrees, then the       */
      /*   circumcenter lies outside the org-dest edge, and eta is        */
      /*   negative.  Roundoff error might prevent eta from being         */
      /*   negative when it should be, so I test eta against xi.)         */
      if (eta < xi) {
        lprevself(badotri);
      }

      /* Insert the circumcenter, searching from the edge of the triangle, */
      /*   and maintain the Delaunay property of the triangulation.        */
      success = insertvertex(m, b, newvertex, &badotri, (struct osub *) NULL,
                             1, 1);
      if (success == SUCCESSFULVERTEX) {
        if (m->steinerleft > 0) {
          m->steinerleft--;
        }
      } else if (success == ENCROACHINGVERTEX) {
        /* If the newly inserted vertex encroaches upon a subsegment, */
        /*   delete the new vertex.                                   */
        undovertex(m, b);
        if (b->verbose > 1) {
          printf("  Rejecting (%.12g, %.12g).\n", newvertex[0], newvertex[1]);
        }
        vertexdealloc(m, newvertex);
      } else if (success == VIOLATINGVERTEX) {
        /* Failed to insert the new vertex, but some subsegment was */
        /*   marked as being encroached.                            */
        vertexdealloc(m, newvertex);
      } else {                                 /* success == DUPLICATEVERTEX */
        /* Couldn't insert the new vertex because a vertex is already there. */
        if (!b->quiet) {
          printf(
            "Warning:  New vertex (%.12g, %.12g) falls on existing vertex.\n",
                 newvertex[0], newvertex[1]);
          errorflag = 1;
        }
        vertexdealloc(m, newvertex);
      }
    }
    if (errorflag) {
      if (b->verbose) {
        printf("  The new vertex is at the circumcenter of triangle\n");
        printf("    (%.12g, %.12g) (%.12g, %.12g) (%.12g, %.12g)\n",
               borg[0], borg[1], bdest[0], bdest[1], bapex[0], bapex[1]);
      }
      printf("This probably means that I am trying to refine triangles\n");
      printf("  to a smaller size than can be accommodated by the finite\n");
      printf("  precision of floating point arithmetic.  (You can be\n");
      printf("  sure of this if I fail to terminate.)\n");
      precisionerror();
    }
  }
}

#endif /* not CDT_ONLY */

/*****************************************************************************/
/*                                                                           */
/*  enforcequality()   Remove all the encroached subsegments and bad         */
/*                     triangles from the triangulation.                     */
/*                                                                           */
/*****************************************************************************/

#ifndef CDT_ONLY

#ifdef ANSI_DECLARATORS
void enforcequality(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void enforcequality(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct badtriang *badtri;
  int i;

  if (!b->quiet) {
    printf("Adding Steiner points to enforce quality.\n");
  }
  /* Initialize the pool of encroached subsegments. */
  poolinit(&m->badsubsegs, sizeof(struct badsubseg), BADSUBSEGPERBLOCK,
           BADSUBSEGPERBLOCK, 0);
  if (b->verbose) {
    printf("  Looking for encroached subsegments.\n");
  }
  /* Test all segments to see if they're encroached. */
  tallyencs(m, b);
  if (b->verbose && (m->badsubsegs.items > 0)) {
    printf("  Splitting encroached subsegments.\n");
  }
  /* Fix encroached subsegments without noting bad triangles. */
  splitencsegs(m, b, 0);
  /* At this point, if we haven't run out of Steiner points, the */
  /*   triangulation should be (conforming) Delaunay.            */

  /* Next, we worry about enforcing triangle quality. */
  if ((b->minangle > 0.0) || b->vararea || b->fixedarea || b->usertest) {
    /* Initialize the pool of bad triangles. */
    poolinit(&m->badtriangles, sizeof(struct badtriang), BADTRIPERBLOCK,
             BADTRIPERBLOCK, 0);
    /* Initialize the queues of bad triangles. */
    for (i = 0; i < 4096; i++) {
      m->queuefront[i] = (struct badtriang *) NULL;
    }
    m->firstnonemptyq = -1;
    /* Test all triangles to see if they're bad. */
    tallyfaces(m, b);
    /* Initialize the pool of recently flipped triangles. */
    poolinit(&m->flipstackers, sizeof(struct flipstacker), FLIPSTACKERPERBLOCK,
             FLIPSTACKERPERBLOCK, 0);
    m->checkquality = 1;
    if (b->verbose) {
      printf("  Splitting bad triangles.\n");
    }
    while ((m->badtriangles.items > 0) && (m->steinerleft != 0)) {
      /* Fix one bad triangle by inserting a vertex at its circumcenter. */
      badtri = dequeuebadtriang(m);
      splittriangle(m, b, badtri);
      if (m->badsubsegs.items > 0) {
        /* Put bad triangle back in queue for another try later. */
        enqueuebadtriang(m, b, badtri);
        /* Fix any encroached subsegments that resulted. */
        /*   Record any new bad triangles that result.   */
        splitencsegs(m, b, 1);
      } else {
        /* Return the bad triangle to the pool. */
        pooldealloc(&m->badtriangles, (VOID *) badtri);
      }
    }
  }
  /* At this point, if the "-D" switch was selected and we haven't run out  */
  /*   of Steiner points, the triangulation should be (conforming) Delaunay */
  /*   and have no low-quality triangles.                                   */

  /* Might we have run out of Steiner points too soon? */
  if (!b->quiet && b->conformdel && (m->badsubsegs.items > 0) &&
      (m->steinerleft == 0)) {
    printf("\nWarning:  I ran out of Steiner points, but the mesh has\n");
    if (m->badsubsegs.items == 1) {
      printf("  one encroached subsegment, and therefore might not be truly\n"
             );
    } else {
      printf("  %ld encroached subsegments, and therefore might not be truly\n"
             , m->badsubsegs.items);
    }
    printf("  Delaunay.  If the Delaunay property is important to you,\n");
    printf("  try increasing the number of Steiner points (controlled by\n");
    printf("  the -S switch) slightly and try again.\n\n");
  }
}

#endif /* not CDT_ONLY */

/**                                                                         **/
/**                                                                         **/
/********* Mesh quality maintenance ends here                        *********/

/*****************************************************************************/
/*                                                                           */
/*  highorder()   Create extra nodes for quadratic subparametric elements.   */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void highorder(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void highorder(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop, trisym;
  struct osub checkmark;
  vertex newvertex;
  vertex torg, tdest;
  int i;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

  if (!b->quiet) {
    printf("Adding vertices for second-order triangles.\n");
  }
  /* The following line ensures that dead items in the pool of nodes    */
  /*   cannot be allocated for the extra nodes associated with high     */
  /*   order elements.  This ensures that the primary nodes (at the     */
  /*   corners of elements) will occur earlier in the output files, and */
  /*   have lower indices, than the extra nodes.                        */
  m->vertices.deaditemstack = (VOID *) NULL;

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each triangle.  If there isn't another triangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent triangle, operate on the edge only if the current triangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (triangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
        org(triangleloop, torg);
        dest(triangleloop, tdest);
        /* Create a new node in the middle of the edge.  Interpolate */
        /*   its attributes.                                         */
        newvertex = (vertex) poolalloc(&m->vertices);
        for (i = 0; i < 2 + m->nextras; i++) {
          newvertex[i] = 0.5 * (torg[i] + tdest[i]);
        }
        /* Set the new node's marker to zero or one, depending on */
        /*   whether it lies on a boundary.                       */
        setvertexmark(newvertex, trisym.tri == m->dummytri);
        setvertextype(newvertex,
                      trisym.tri == m->dummytri ? FREEVERTEX : SEGMENTVERTEX);
        if (b->usesegments) {
          tspivot(triangleloop, checkmark);
          /* If this edge is a segment, transfer the marker to the new node. */
          if (checkmark.ss != m->dummysub) {
            setvertexmark(newvertex, mark(checkmark));
            setvertextype(newvertex, SEGMENTVERTEX);
          }
        }
        if (b->verbose > 1) {
          printf("  Creating (%.12g, %.12g).\n", newvertex[0], newvertex[1]);
        }
        /* Record the new node in the (one or two) adjacent elements. */
        triangleloop.tri[m->highorderindex + triangleloop.orient] =
                (triangle) newvertex;
        if (trisym.tri != m->dummytri) {
          trisym.tri[m->highorderindex + trisym.orient] = (triangle) newvertex;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }
}

/********* File I/O routines begin here                              *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  readline()   Read a nonempty line from a file.                           */
/*                                                                           */
/*  A line is considered "nonempty" if it contains something that looks like */
/*  a number.  Comments (prefaced by `#') are ignored.                       */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
char *readline(char *string, FILE *infile, char *infilename)
#else /* not ANSI_DECLARATORS */
char *readline(string, infile, infilename)
char *string;
FILE *infile;
char *infilename;
#endif /* not ANSI_DECLARATORS */

{
  char *result;

  /* Search for something that looks like a number. */
  do {
    result = fgets(string, INPUTLINESIZE, infile);
    if (result == (char *) NULL) {
      printf("  Error:  Unexpected end of file in %s.\n", infilename);
      triexit(1);
    }
    /* Skip anything that doesn't look like a number, a comment, */
    /*   or the end of a line.                                   */
    while ((*result != '\0') && (*result != '#')
           && (*result != '.') && (*result != '+') && (*result != '-')
           && ((*result < '0') || (*result > '9'))) {
      result++;
    }
  /* If it's a comment or end of line, read another line and try again. */
  } while ((*result == '#') || (*result == '\0'));
  return result;
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  findfield()   Find the next field of a string.                           */
/*                                                                           */
/*  Jumps past the current field by searching for whitespace, then jumps     */
/*  past the whitespace to find the next field.                              */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
char *findfield(char *string)
#else /* not ANSI_DECLARATORS */
char *findfield(string)
char *string;
#endif /* not ANSI_DECLARATORS */

{
  char *result;

  result = string;
  /* Skip the current field.  Stop upon reaching whitespace. */
  while ((*result != '\0') && (*result != '#')
         && (*result != ' ') && (*result != '\t')) {
    result++;
  }
  /* Now skip the whitespace and anything else that doesn't look like a */
  /*   number, a comment, or the end of a line.                         */
  while ((*result != '\0') && (*result != '#')
         && (*result != '.') && (*result != '+') && (*result != '-')
         && ((*result < '0') || (*result > '9'))) {
    result++;
  }
  /* Check for a comment (prefixed with `#'). */
  if (*result == '#') {
    *result = '\0';
  }
  return result;
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  readnodes()   Read the vertices from a file, which may be a .node or     */
/*                .poly file.                                                */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
void readnodes(struct mesh *m, struct behavior *b, char *nodefilename,
               char *polyfilename, FILE **polyfile)
#else /* not ANSI_DECLARATORS */
void readnodes(m, b, nodefilename, polyfilename, polyfile)
struct mesh *m;
struct behavior *b;
char *nodefilename;
char *polyfilename;
FILE **polyfile;
#endif /* not ANSI_DECLARATORS */

{
  FILE *infile;
  vertex vertexloop;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  char *infilename;
  REAL x, y;
  int firstnode;
  int nodemarkers;
  int currentmarker;
  int i, j;

  if (b->poly) {
    /* Read the vertices from a .poly file. */
    if (!b->quiet) {
      printf("Opening %s.\n", polyfilename);
    }
    *polyfile = fopen(polyfilename, "r");
    if (*polyfile == (FILE *) NULL) {
      printf("  Error:  Cannot access file %s.\n", polyfilename);
      triexit(1);
    }
    /* Read number of vertices, number of dimensions, number of vertex */
    /*   attributes, and number of boundary markers.                   */
    stringptr = readline(inputline, *polyfile, polyfilename);
    m->invertices = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      m->mesh_dim = 2;
    } else {
      m->mesh_dim = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      m->nextras = 0;
    } else {
      m->nextras = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      nodemarkers = 0;
    } else {
      nodemarkers = (int) strtol(stringptr, &stringptr, 0);
    }
    if (m->invertices > 0) {
      infile = *polyfile;
      infilename = polyfilename;
      m->readnodefile = 0;
    } else {
      /* If the .poly file claims there are zero vertices, that means that */
      /*   the vertices should be read from a separate .node file.         */
      m->readnodefile = 1;
      infilename = nodefilename;
    }
  } else {
    m->readnodefile = 1;
    infilename = nodefilename;
    *polyfile = (FILE *) NULL;
  }

  if (m->readnodefile) {
    /* Read the vertices from a .node file. */
    if (!b->quiet) {
      printf("Opening %s.\n", nodefilename);
    }
    infile = fopen(nodefilename, "r");
    if (infile == (FILE *) NULL) {
      printf("  Error:  Cannot access file %s.\n", nodefilename);
      triexit(1);
    }
    /* Read number of vertices, number of dimensions, number of vertex */
    /*   attributes, and number of boundary markers.                   */
    stringptr = readline(inputline, infile, nodefilename);
    m->invertices = (int) strtol(stringptr, &stringptr, 0);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      m->mesh_dim = 2;
    } else {
      m->mesh_dim = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      m->nextras = 0;
    } else {
      m->nextras = (int) strtol(stringptr, &stringptr, 0);
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      nodemarkers = 0;
    } else {
      nodemarkers = (int) strtol(stringptr, &stringptr, 0);
    }
  }

  if (m->invertices < 3) {
    printf("Error:  Input must have at least three input vertices.\n");
    triexit(1);
  }
  if (m->mesh_dim != 2) {
    printf("Error:  Triangle only works with two-dimensional meshes.\n");
    triexit(1);
  }
  if (m->nextras == 0) {
    b->weighted = 0;
  }

  initializevertexpool(m, b);

  /* Read the vertices. */
  for (i = 0; i < m->invertices; i++) {
    vertexloop = (vertex) poolalloc(&m->vertices);
    stringptr = readline(inputline, infile, infilename);
    if (i == 0) {
      firstnode = (int) strtol(stringptr, &stringptr, 0);
      if ((firstnode == 0) || (firstnode == 1)) {
        b->firstnumber = firstnode;
      }
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Vertex %d has no x coordinate.\n", b->firstnumber + i);
      triexit(1);
    }
    x = (REAL) strtod(stringptr, &stringptr);
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Vertex %d has no y coordinate.\n", b->firstnumber + i);
      triexit(1);
    }
    y = (REAL) strtod(stringptr, &stringptr);
    vertexloop[0] = x;
    vertexloop[1] = y;
    /* Read the vertex attributes. */
    for (j = 2; j < 2 + m->nextras; j++) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        vertexloop[j] = 0.0;
      } else {
        vertexloop[j] = (REAL) strtod(stringptr, &stringptr);
      }
    }
    if (nodemarkers) {
      /* Read a vertex marker. */
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        setvertexmark(vertexloop, 0);
      } else {
        currentmarker = (int) strtol(stringptr, &stringptr, 0);
        setvertexmark(vertexloop, currentmarker);
      }
    } else {
      /* If no markers are specified in the file, they default to zero. */
      setvertexmark(vertexloop, 0);
    }
    setvertextype(vertexloop, INPUTVERTEX);
    /* Determine the smallest and largest x and y coordinates. */
    if (i == 0) {
      m->xmin = m->xmax = x;
      m->ymin = m->ymax = y;
    } else {
      m->xmin = (x < m->xmin) ? x : m->xmin;
      m->xmax = (x > m->xmax) ? x : m->xmax;
      m->ymin = (y < m->ymin) ? y : m->ymin;
      m->ymax = (y > m->ymax) ? y : m->ymax;
    }
  }
  if (m->readnodefile) {
    fclose(infile);
  }

  /* Nonexistent x value used as a flag to mark circle events in sweepline */
  /*   Delaunay algorithm.                                                 */
  m->xminextreme = 10 * m->xmin - 9 * m->xmax;
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  transfernodes()   Read the vertices from memory.                         */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void transfernodes(struct mesh *m, struct behavior *b, REAL *pointlist,
                   REAL *pointattriblist, int *pointmarkerlist,
                   int numberofpoints, int numberofpointattribs)
#else /* not ANSI_DECLARATORS */
void transfernodes(m, b, pointlist, pointattriblist, pointmarkerlist,
                   numberofpoints, numberofpointattribs)
struct mesh *m;
struct behavior *b;
REAL *pointlist;
REAL *pointattriblist;
int *pointmarkerlist;
int numberofpoints;
int numberofpointattribs;
#endif /* not ANSI_DECLARATORS */

{
  vertex vertexloop;
  REAL x, y;
  int i, j;
  int coordindex;
  int attribindex;

  m->invertices = numberofpoints;
  m->mesh_dim = 2;
  m->nextras = numberofpointattribs;
  m->readnodefile = 0;
  if (m->invertices < 3) {
    printf("Error:  Input must have at least three input vertices.\n");
    triexit(1);
  }
  if (m->nextras == 0) {
    b->weighted = 0;
  }

  initializevertexpool(m, b);

  /* Read the vertices. */
  coordindex = 0;
  attribindex = 0;
  for (i = 0; i < m->invertices; i++) {
    vertexloop = (vertex) poolalloc(&m->vertices);
    /* Read the vertex coordinates. */
    x = vertexloop[0] = pointlist[coordindex++];
    y = vertexloop[1] = pointlist[coordindex++];
    /* Read the vertex attributes. */
    for (j = 0; j < numberofpointattribs; j++) {
      vertexloop[2 + j] = pointattriblist[attribindex++];
    }
    if (pointmarkerlist != (int *) NULL) {
      /* Read a vertex marker. */
      setvertexmark(vertexloop, pointmarkerlist[i]);
    } else {
      /* If no markers are specified, they default to zero. */
      setvertexmark(vertexloop, 0);
    }
    setvertextype(vertexloop, INPUTVERTEX);
    /* Determine the smallest and largest x and y coordinates. */
    if (i == 0) {
      m->xmin = m->xmax = x;
      m->ymin = m->ymax = y;
    } else {
      m->xmin = (x < m->xmin) ? x : m->xmin;
      m->xmax = (x > m->xmax) ? x : m->xmax;
      m->ymin = (y < m->ymin) ? y : m->ymin;
      m->ymax = (y > m->ymax) ? y : m->ymax;
    }
  }

  /* Nonexistent x value used as a flag to mark circle events in sweepline */
  /*   Delaunay algorithm.                                                 */
  m->xminextreme = 10 * m->xmin - 9 * m->xmax;
}

#endif /* TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  readholes()   Read the holes, and possibly regional attributes and area  */
/*                constraints, from a .poly file.                            */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
void readholes(struct mesh *m, struct behavior *b,
               FILE *polyfile, char *polyfilename, REAL **hlist, int *holes,
               REAL **rlist, int *regions)
#else /* not ANSI_DECLARATORS */
void readholes(m, b, polyfile, polyfilename, hlist, holes, rlist, regions)
struct mesh *m;
struct behavior *b;
FILE *polyfile;
char *polyfilename;
REAL **hlist;
int *holes;
REAL **rlist;
int *regions;
#endif /* not ANSI_DECLARATORS */

{
  REAL *holelist;
  REAL *regionlist;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int index;
  int i;

  /* Read the holes. */
  stringptr = readline(inputline, polyfile, polyfilename);
  *holes = (int) strtol(stringptr, &stringptr, 0);
  if (*holes > 0) {
    holelist = (REAL *) trimalloc(2 * *holes * (int) sizeof(REAL));
    *hlist = holelist;
    for (i = 0; i < 2 * *holes; i += 2) {
      stringptr = readline(inputline, polyfile, polyfilename);
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Hole %d has no x coordinate.\n",
               b->firstnumber + (i >> 1));
        triexit(1);
      } else {
        holelist[i] = (REAL) strtod(stringptr, &stringptr);
      }
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Hole %d has no y coordinate.\n",
               b->firstnumber + (i >> 1));
        triexit(1);
      } else {
        holelist[i + 1] = (REAL) strtod(stringptr, &stringptr);
      }
    }
  } else {
    *hlist = (REAL *) NULL;
  }

#ifndef CDT_ONLY
  if ((b->regionattrib || b->vararea) && !b->refine) {
    /* Read the area constraints. */
    stringptr = readline(inputline, polyfile, polyfilename);
    *regions = (int) strtol(stringptr, &stringptr, 0);
    if (*regions > 0) {
      regionlist = (REAL *) trimalloc(4 * *regions * (int) sizeof(REAL));
      *rlist = regionlist;
      index = 0;
      for (i = 0; i < *regions; i++) {
        stringptr = readline(inputline, polyfile, polyfilename);
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no x coordinate.\n",
                 b->firstnumber + i);
          triexit(1);
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Region %d has no y coordinate.\n",
                 b->firstnumber + i);
          triexit(1);
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          printf(
            "Error:  Region %d has no region attribute or area constraint.\n",
                 b->firstnumber + i);
          triexit(1);
        } else {
          regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
        }
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          regionlist[index] = regionlist[index - 1];
        } else {
          regionlist[index] = (REAL) strtod(stringptr, &stringptr);
        }
        index++;
      }
    }
  } else {
    /* Set `*regions' to zero to avoid an accidental free() later. */
    *regions = 0;
    *rlist = (REAL *) NULL;
  }
#endif /* not CDT_ONLY */

  fclose(polyfile);
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  finishfile()   Write the command line to the output file so the user     */
/*                 can remember how the file was generated.  Close the file. */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
void finishfile(FILE *outfile, int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void finishfile(outfile, argc, argv)
FILE *outfile;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

{
  int i;

  fprintf(outfile, "# Generated by");
  for (i = 0; i < argc; i++) {
    fprintf(outfile, " ");
    fputs(argv[i], outfile);
  }
  fprintf(outfile, "\n");
  fclose(outfile);
}

#endif /* not TRILIBRARY */

/*****************************************************************************/
/*                                                                           */
/*  writenodes()   Number the vertices and write them to a .node file.       */
/*                                                                           */
/*  To save memory, the vertex numbers are written over the boundary markers */
/*  after the vertices are written to a file.                                */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writenodes(struct mesh *m, struct behavior *b, REAL **pointlist,
                REAL **pointattriblist, int **pointmarkerlist)
#else /* not ANSI_DECLARATORS */
void writenodes(m, b, pointlist, pointattriblist, pointmarkerlist)
struct mesh *m;
struct behavior *b;
REAL **pointlist;
REAL **pointattriblist;
int **pointmarkerlist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writenodes(struct mesh *m, struct behavior *b, char *nodefilename,
                int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writenodes(m, b, nodefilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *nodefilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  REAL *plist;
  REAL *palist;
  int *pmlist;
  int coordindex;
  int attribindex;
#else /* not TRILIBRARY */
  FILE *outfile;
#endif /* not TRILIBRARY */
  vertex vertexloop;
  long outvertices;
  int vertexnumber;
  int i;

  if (b->jettison) {
    outvertices = m->vertices.items - m->undeads;
  } else {
    outvertices = m->vertices.items;
  }

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing vertices.\n");
  }
  /* Allocate memory for output vertices if necessary. */
  if (*pointlist == (REAL *) NULL) {
    *pointlist = (REAL *) trimalloc((int) (outvertices * 2 * sizeof(REAL)));
  }
  /* Allocate memory for output vertex attributes if necessary. */
  if ((m->nextras > 0) && (*pointattriblist == (REAL *) NULL)) {
    *pointattriblist = (REAL *) trimalloc((int) (outvertices * m->nextras *
                                                 sizeof(REAL)));
  }
  /* Allocate memory for output vertex markers if necessary. */
  if (!b->nobound && (*pointmarkerlist == (int *) NULL)) {
    *pointmarkerlist = (int *) trimalloc((int) (outvertices * sizeof(int)));
  }
  plist = *pointlist;
  palist = *pointattriblist;
  pmlist = *pointmarkerlist;
  coordindex = 0;
  attribindex = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", nodefilename);
  }
  outfile = fopen(nodefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", nodefilename);
    triexit(1);
  }
  /* Number of vertices, number of dimensions, number of vertex attributes, */
  /*   and number of boundary markers (zero or one).                        */
  fprintf(outfile, "%ld  %d  %d  %d\n", outvertices, m->mesh_dim,
          m->nextras, 1 - b->nobound);
#endif /* not TRILIBRARY */

  traversalinit(&m->vertices);
  vertexnumber = b->firstnumber;
  vertexloop = vertextraverse(m);
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
#ifdef TRILIBRARY
      /* X and y coordinates. */
      plist[coordindex++] = vertexloop[0];
      plist[coordindex++] = vertexloop[1];
      /* Vertex attributes. */
      for (i = 0; i < m->nextras; i++) {
        palist[attribindex++] = vertexloop[2 + i];
      }
      if (!b->nobound) {
        /* Copy the boundary marker. */
        pmlist[vertexnumber - b->firstnumber] = vertexmark(vertexloop);
      }
#else /* not TRILIBRARY */
      /* Vertex number, x and y coordinates. */
      fprintf(outfile, "%4d    %.17g  %.17g", vertexnumber, vertexloop[0],
              vertexloop[1]);
      for (i = 0; i < m->nextras; i++) {
        /* Write an attribute. */
        fprintf(outfile, "  %.17g", vertexloop[i + 2]);
      }
      if (b->nobound) {
        fprintf(outfile, "\n");
      } else {
        /* Write the boundary marker. */
        fprintf(outfile, "    %d\n", vertexmark(vertexloop));
      }
#endif /* not TRILIBRARY */

      setvertexmark(vertexloop, vertexnumber);
      vertexnumber++;
    }
    vertexloop = vertextraverse(m);
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

/*****************************************************************************/
/*                                                                           */
/*  numbernodes()   Number the vertices.                                     */
/*                                                                           */
/*  Each vertex is assigned a marker equal to its number.                    */
/*                                                                           */
/*  Used when writenodes() is not called because no .node file is written.   */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void numbernodes(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void numbernodes(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  vertex vertexloop;
  int vertexnumber;

  traversalinit(&m->vertices);
  vertexnumber = b->firstnumber;
  vertexloop = vertextraverse(m);
  while (vertexloop != (vertex) NULL) {
    setvertexmark(vertexloop, vertexnumber);
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      vertexnumber++;
    }
    vertexloop = vertextraverse(m);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  writeelements()   Write the triangles to an .ele file.                   */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writeelements(struct mesh *m, struct behavior *b,
                   int **trianglelist, REAL **triangleattriblist)
#else /* not ANSI_DECLARATORS */
void writeelements(m, b, trianglelist, triangleattriblist)
struct mesh *m;
struct behavior *b;
int **trianglelist;
REAL **triangleattriblist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writeelements(struct mesh *m, struct behavior *b, char *elefilename,
                   int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writeelements(m, b, elefilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *elefilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  int *tlist;
  REAL *talist;
  int vertexindex;
  int attribindex;
#else /* not TRILIBRARY */
  FILE *outfile;
#endif /* not TRILIBRARY */
  struct otri triangleloop;
  vertex p1, p2, p3;
  vertex mid1, mid2, mid3;
  long elementnumber;
  int i;

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing triangles.\n");
  }
  /* Allocate memory for output triangles if necessary. */
  if (*trianglelist == (int *) NULL) {
    *trianglelist = (int *) trimalloc((int) (m->triangles.items *
                                             ((b->order + 1) * (b->order + 2) /
                                              2) * sizeof(int)));
  }
  /* Allocate memory for output triangle attributes if necessary. */
  if ((m->eextras > 0) && (*triangleattriblist == (REAL *) NULL)) {
    *triangleattriblist = (REAL *) trimalloc((int) (m->triangles.items *
                                                    m->eextras *
                                                    sizeof(REAL)));
  }
  tlist = *trianglelist;
  talist = *triangleattriblist;
  vertexindex = 0;
  attribindex = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", elefilename);
  }
  outfile = fopen(elefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", elefilename);
    triexit(1);
  }
  /* Number of triangles, vertices per triangle, attributes per triangle. */
  fprintf(outfile, "%ld  %d  %d\n", m->triangles.items,
          (b->order + 1) * (b->order + 2) / 2, m->eextras);
#endif /* not TRILIBRARY */

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  elementnumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, p1);
    dest(triangleloop, p2);
    apex(triangleloop, p3);
    if (b->order == 1) {
#ifdef TRILIBRARY
      tlist[vertexindex++] = vertexmark(p1);
      tlist[vertexindex++] = vertexmark(p2);
      tlist[vertexindex++] = vertexmark(p3);
#else /* not TRILIBRARY */
      /* Triangle number, indices for three vertices. */
      fprintf(outfile, "%4ld    %4d  %4d  %4d", elementnumber,
              vertexmark(p1), vertexmark(p2), vertexmark(p3));
#endif /* not TRILIBRARY */
    } else {
      mid1 = (vertex) triangleloop.tri[m->highorderindex + 1];
      mid2 = (vertex) triangleloop.tri[m->highorderindex + 2];
      mid3 = (vertex) triangleloop.tri[m->highorderindex];
#ifdef TRILIBRARY
      tlist[vertexindex++] = vertexmark(p1);
      tlist[vertexindex++] = vertexmark(p2);
      tlist[vertexindex++] = vertexmark(p3);
      tlist[vertexindex++] = vertexmark(mid1);
      tlist[vertexindex++] = vertexmark(mid2);
      tlist[vertexindex++] = vertexmark(mid3);
#else /* not TRILIBRARY */
      /* Triangle number, indices for six vertices. */
      fprintf(outfile, "%4ld    %4d  %4d  %4d  %4d  %4d  %4d", elementnumber,
              vertexmark(p1), vertexmark(p2), vertexmark(p3), vertexmark(mid1),
              vertexmark(mid2), vertexmark(mid3));
#endif /* not TRILIBRARY */
    }

#ifdef TRILIBRARY
    for (i = 0; i < m->eextras; i++) {
      talist[attribindex++] = elemattribute(triangleloop, i);
    }
#else /* not TRILIBRARY */
    for (i = 0; i < m->eextras; i++) {
      fprintf(outfile, "  %.17g", elemattribute(triangleloop, i));
    }
    fprintf(outfile, "\n");
#endif /* not TRILIBRARY */

    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

/*****************************************************************************/
/*                                                                           */
/*  writepoly()   Write the segments and holes to a .poly file.              */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writepoly(struct mesh *m, struct behavior *b,
               int **segmentlist, int **segmentmarkerlist)
#else /* not ANSI_DECLARATORS */
void writepoly(m, b, segmentlist, segmentmarkerlist)
struct mesh *m;
struct behavior *b;
int **segmentlist;
int **segmentmarkerlist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writepoly(struct mesh *m, struct behavior *b, char *polyfilename,
               REAL *holelist, int holes, REAL *regionlist, int regions,
               int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writepoly(m, b, polyfilename, holelist, holes, regionlist, regions,
               argc, argv)
struct mesh *m;
struct behavior *b;
char *polyfilename;
REAL *holelist;
int holes;
REAL *regionlist;
int regions;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  int *slist;
  int *smlist;
  int index;
#else /* not TRILIBRARY */
  FILE *outfile;
  long holenumber, regionnumber;
#endif /* not TRILIBRARY */
  struct osub subsegloop;
  vertex endpoint1, endpoint2;
  long subsegnumber;

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing segments.\n");
  }
  /* Allocate memory for output segments if necessary. */
  if (*segmentlist == (int *) NULL) {
    *segmentlist = (int *) trimalloc((int) (m->subsegs.items * 2 *
                                            sizeof(int)));
  }
  /* Allocate memory for output segment markers if necessary. */
  if (!b->nobound && (*segmentmarkerlist == (int *) NULL)) {
    *segmentmarkerlist = (int *) trimalloc((int) (m->subsegs.items *
                                                  sizeof(int)));
  }
  slist = *segmentlist;
  smlist = *segmentmarkerlist;
  index = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", polyfilename);
  }
  outfile = fopen(polyfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", polyfilename);
    triexit(1);
  }
  /* The zero indicates that the vertices are in a separate .node file. */
  /*   Followed by number of dimensions, number of vertex attributes,   */
  /*   and number of boundary markers (zero or one).                    */
  fprintf(outfile, "%d  %d  %d  %d\n", 0, m->mesh_dim, m->nextras,
          1 - b->nobound);
  /* Number of segments, number of boundary markers (zero or one). */
  fprintf(outfile, "%ld  %d\n", m->subsegs.items, 1 - b->nobound);
#endif /* not TRILIBRARY */

  traversalinit(&m->subsegs);
  subsegloop.ss = subsegtraverse(m);
  subsegloop.ssorient = 0;
  subsegnumber = b->firstnumber;
  while (subsegloop.ss != (subseg *) NULL) {
    sorg(subsegloop, endpoint1);
    sdest(subsegloop, endpoint2);
#ifdef TRILIBRARY
    /* Copy indices of the segment's two endpoints. */
    slist[index++] = vertexmark(endpoint1);
    slist[index++] = vertexmark(endpoint2);
    if (!b->nobound) {
      /* Copy the boundary marker. */
      smlist[subsegnumber - b->firstnumber] = mark(subsegloop);
    }
#else /* not TRILIBRARY */
    /* Segment number, indices of its two endpoints, and possibly a marker. */
    if (b->nobound) {
      fprintf(outfile, "%4ld    %4d  %4d\n", subsegnumber,
              vertexmark(endpoint1), vertexmark(endpoint2));
    } else {
      fprintf(outfile, "%4ld    %4d  %4d    %4d\n", subsegnumber,
              vertexmark(endpoint1), vertexmark(endpoint2), mark(subsegloop));
    }
#endif /* not TRILIBRARY */

    subsegloop.ss = subsegtraverse(m);
    subsegnumber++;
  }

#ifndef TRILIBRARY
#ifndef CDT_ONLY
  fprintf(outfile, "%d\n", holes);
  if (holes > 0) {
    for (holenumber = 0; holenumber < holes; holenumber++) {
      /* Hole number, x and y coordinates. */
      fprintf(outfile, "%4ld   %.17g  %.17g\n", b->firstnumber + holenumber,
              holelist[2 * holenumber], holelist[2 * holenumber + 1]);
    }
  }
  if (regions > 0) {
    fprintf(outfile, "%d\n", regions);
    for (regionnumber = 0; regionnumber < regions; regionnumber++) {
      /* Region number, x and y coordinates, attribute, maximum area. */
      fprintf(outfile, "%4ld   %.17g  %.17g  %.17g  %.17g\n",
              b->firstnumber + regionnumber,
              regionlist[4 * regionnumber], regionlist[4 * regionnumber + 1],
              regionlist[4 * regionnumber + 2],
              regionlist[4 * regionnumber + 3]);
    }
  }
#endif /* not CDT_ONLY */

  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

/*****************************************************************************/
/*                                                                           */
/*  writeedges()   Write the edges to an .edge file.                         */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writeedges(struct mesh *m, struct behavior *b,
                int **edgelist, int **edgemarkerlist)
#else /* not ANSI_DECLARATORS */
void writeedges(m, b, edgelist, edgemarkerlist)
struct mesh *m;
struct behavior *b;
int **edgelist;
int **edgemarkerlist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writeedges(struct mesh *m, struct behavior *b, char *edgefilename,
                int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writeedges(m, b, edgefilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *edgefilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  int *elist;
  int *emlist;
  int index;
#else /* not TRILIBRARY */
  FILE *outfile;
#endif /* not TRILIBRARY */
  struct otri triangleloop, trisym;
  struct osub checkmark;
  vertex p1, p2;
  long edgenumber;
  triangle ptr;                         /* Temporary variable used by sym(). */
  subseg sptr;                      /* Temporary variable used by tspivot(). */

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing edges.\n");
  }
  /* Allocate memory for edges if necessary. */
  if (*edgelist == (int *) NULL) {
    *edgelist = (int *) trimalloc((int) (m->edges * 2 * sizeof(int)));
  }
  /* Allocate memory for edge markers if necessary. */
  if (!b->nobound && (*edgemarkerlist == (int *) NULL)) {
    *edgemarkerlist = (int *) trimalloc((int) (m->edges * sizeof(int)));
  }
  elist = *edgelist;
  emlist = *edgemarkerlist;
  index = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", edgefilename);
  }
  outfile = fopen(edgefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", edgefilename);
    triexit(1);
  }
  /* Number of edges, number of boundary markers (zero or one). */
  fprintf(outfile, "%ld  %d\n", m->edges, 1 - b->nobound);
#endif /* not TRILIBRARY */

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  edgenumber = b->firstnumber;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each triangle.  If there isn't another triangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent triangle, operate on the edge only if the current triangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (triangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
        org(triangleloop, p1);
        dest(triangleloop, p2);
#ifdef TRILIBRARY
        elist[index++] = vertexmark(p1);
        elist[index++] = vertexmark(p2);
#endif /* TRILIBRARY */
        if (b->nobound) {
#ifndef TRILIBRARY
          /* Edge number, indices of two endpoints. */
          fprintf(outfile, "%4ld   %d  %d\n", edgenumber,
                  vertexmark(p1), vertexmark(p2));
#endif /* not TRILIBRARY */
        } else {
          /* Edge number, indices of two endpoints, and a boundary marker. */
          /*   If there's no subsegment, the boundary marker is zero.      */
          if (b->usesegments) {
            tspivot(triangleloop, checkmark);
            if (checkmark.ss == m->dummysub) {
#ifdef TRILIBRARY
              emlist[edgenumber - b->firstnumber] = 0;
#else /* not TRILIBRARY */
              fprintf(outfile, "%4ld   %d  %d  %d\n", edgenumber,
                      vertexmark(p1), vertexmark(p2), 0);
#endif /* not TRILIBRARY */
            } else {
#ifdef TRILIBRARY
              emlist[edgenumber - b->firstnumber] = mark(checkmark);
#else /* not TRILIBRARY */
              fprintf(outfile, "%4ld   %d  %d  %d\n", edgenumber,
                      vertexmark(p1), vertexmark(p2), mark(checkmark));
#endif /* not TRILIBRARY */
            }
          } else {
#ifdef TRILIBRARY
            emlist[edgenumber - b->firstnumber] = trisym.tri == m->dummytri;
#else /* not TRILIBRARY */
            fprintf(outfile, "%4ld   %d  %d  %d\n", edgenumber,
                    vertexmark(p1), vertexmark(p2), trisym.tri == m->dummytri);
#endif /* not TRILIBRARY */
          }
        }
        edgenumber++;
      }
    }
    triangleloop.tri = triangletraverse(m);
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

/*****************************************************************************/
/*                                                                           */
/*  writevoronoi()   Write the Voronoi diagram to a .v.node and .v.edge      */
/*                   file.                                                   */
/*                                                                           */
/*  The Voronoi diagram is the geometric dual of the Delaunay triangulation. */
/*  Hence, the Voronoi vertices are listed by traversing the Delaunay        */
/*  triangles, and the Voronoi edges are listed by traversing the Delaunay   */
/*  edges.                                                                   */
/*                                                                           */
/*  WARNING:  In order to assign numbers to the Voronoi vertices, this       */
/*  procedure messes up the subsegments or the extra nodes of every          */
/*  element.  Hence, you should call this procedure last.                    */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writevoronoi(struct mesh *m, struct behavior *b, REAL **vpointlist,
                  REAL **vpointattriblist, int **vpointmarkerlist,
                  int **vedgelist, int **vedgemarkerlist, REAL **vnormlist)
#else /* not ANSI_DECLARATORS */
void writevoronoi(m, b, vpointlist, vpointattriblist, vpointmarkerlist,
                  vedgelist, vedgemarkerlist, vnormlist)
struct mesh *m;
struct behavior *b;
REAL **vpointlist;
REAL **vpointattriblist;
int **vpointmarkerlist;
int **vedgelist;
int **vedgemarkerlist;
REAL **vnormlist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writevoronoi(struct mesh *m, struct behavior *b, char *vnodefilename,
                  char *vedgefilename, int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writevoronoi(m, b, vnodefilename, vedgefilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *vnodefilename;
char *vedgefilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  REAL *plist;
  REAL *palist;
  int *elist;
  REAL *normlist;
  int coordindex;
  int attribindex;
#else /* not TRILIBRARY */
  FILE *outfile;
#endif /* not TRILIBRARY */
  struct otri triangleloop, trisym;
  vertex torg, tdest, tapex;
  REAL circumcenter[2];
  REAL xi, eta;
  long vnodenumber, vedgenumber;
  int p1, p2;
  int i;
  triangle ptr;                         /* Temporary variable used by sym(). */

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing Voronoi vertices.\n");
  }
  /* Allocate memory for Voronoi vertices if necessary. */
  if (*vpointlist == (REAL *) NULL) {
    *vpointlist = (REAL *) trimalloc((int) (m->triangles.items * 2 *
                                            sizeof(REAL)));
  }
  /* Allocate memory for Voronoi vertex attributes if necessary. */
  if (*vpointattriblist == (REAL *) NULL) {
    *vpointattriblist = (REAL *) trimalloc((int) (m->triangles.items *
                                                  m->nextras * sizeof(REAL)));
  }
  *vpointmarkerlist = (int *) NULL;
  plist = *vpointlist;
  palist = *vpointattriblist;
  coordindex = 0;
  attribindex = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", vnodefilename);
  }
  outfile = fopen(vnodefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", vnodefilename);
    triexit(1);
  }
  /* Number of triangles, two dimensions, number of vertex attributes, */
  /*   no markers.                                                     */
  fprintf(outfile, "%ld  %d  %d  %d\n", m->triangles.items, 2, m->nextras, 0);
#endif /* not TRILIBRARY */

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  vnodenumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, torg);
    dest(triangleloop, tdest);
    apex(triangleloop, tapex);
    findcircumcenter(m, b, torg, tdest, tapex, circumcenter, &xi, &eta, 0);
#ifdef TRILIBRARY
    /* X and y coordinates. */
    plist[coordindex++] = circumcenter[0];
    plist[coordindex++] = circumcenter[1];
    for (i = 2; i < 2 + m->nextras; i++) {
      /* Interpolate the vertex attributes at the circumcenter. */
      palist[attribindex++] = torg[i] + xi * (tdest[i] - torg[i])
                                     + eta * (tapex[i] - torg[i]);
    }
#else /* not TRILIBRARY */
    /* Voronoi vertex number, x and y coordinates. */
    fprintf(outfile, "%4ld    %.17g  %.17g", vnodenumber, circumcenter[0],
            circumcenter[1]);
    for (i = 2; i < 2 + m->nextras; i++) {
      /* Interpolate the vertex attributes at the circumcenter. */
      fprintf(outfile, "  %.17g", torg[i] + xi * (tdest[i] - torg[i])
                                         + eta * (tapex[i] - torg[i]));
    }
    fprintf(outfile, "\n");
#endif /* not TRILIBRARY */

    * (int *) (triangleloop.tri + 6) = (int) vnodenumber;
    triangleloop.tri = triangletraverse(m);
    vnodenumber++;
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing Voronoi edges.\n");
  }
  /* Allocate memory for output Voronoi edges if necessary. */
  if (*vedgelist == (int *) NULL) {
    *vedgelist = (int *) trimalloc((int) (m->edges * 2 * sizeof(int)));
  }
  *vedgemarkerlist = (int *) NULL;
  /* Allocate memory for output Voronoi norms if necessary. */
  if (*vnormlist == (REAL *) NULL) {
    *vnormlist = (REAL *) trimalloc((int) (m->edges * 2 * sizeof(REAL)));
  }
  elist = *vedgelist;
  normlist = *vnormlist;
  coordindex = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", vedgefilename);
  }
  outfile = fopen(vedgefilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", vedgefilename);
    triexit(1);
  }
  /* Number of edges, zero boundary markers. */
  fprintf(outfile, "%ld  %d\n", m->edges, 0);
#endif /* not TRILIBRARY */

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  vedgenumber = b->firstnumber;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each triangle.  If there isn't another triangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent triangle, operate on the edge only if the current triangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (triangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
         triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
        /* Find the number of this triangle (and Voronoi vertex). */
        p1 = * (int *) (triangleloop.tri + 6);
        if (trisym.tri == m->dummytri) {
          org(triangleloop, torg);
          dest(triangleloop, tdest);
#ifdef TRILIBRARY
          /* Copy an infinite ray.  Index of one endpoint, and -1. */
          elist[coordindex] = p1;
          normlist[coordindex++] = tdest[1] - torg[1];
          elist[coordindex] = -1;
          normlist[coordindex++] = torg[0] - tdest[0];
#else /* not TRILIBRARY */
          /* Write an infinite ray.  Edge number, index of one endpoint, -1, */
          /*   and x and y coordinates of a vector representing the          */
          /*   direction of the ray.                                         */
          fprintf(outfile, "%4ld   %d  %d   %.17g  %.17g\n", vedgenumber,
                  p1, -1, tdest[1] - torg[1], torg[0] - tdest[0]);
#endif /* not TRILIBRARY */
        } else {
          /* Find the number of the adjacent triangle (and Voronoi vertex). */
          p2 = * (int *) (trisym.tri + 6);
          /* Finite edge.  Write indices of two endpoints. */
#ifdef TRILIBRARY
          elist[coordindex] = p1;
          normlist[coordindex++] = 0.0;
          elist[coordindex] = p2;
          normlist[coordindex++] = 0.0;
#else /* not TRILIBRARY */
          fprintf(outfile, "%4ld   %d  %d\n", vedgenumber, p1, p2);
#endif /* not TRILIBRARY */
        }
        vedgenumber++;
      }
    }
    triangleloop.tri = triangletraverse(m);
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writeneighbors(struct mesh *m, struct behavior *b, int **neighborlist)
#else /* not ANSI_DECLARATORS */
void writeneighbors(m, b, neighborlist)
struct mesh *m;
struct behavior *b;
int **neighborlist;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
void writeneighbors(struct mesh *m, struct behavior *b, char *neighborfilename,
                    int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writeneighbors(m, b, neighborfilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *neighborfilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
#ifdef TRILIBRARY
  int *nlist;
  int index;
#else /* not TRILIBRARY */
  FILE *outfile;
#endif /* not TRILIBRARY */
  struct otri triangleloop, trisym;
  long elementnumber;
  int neighbor1, neighbor2, neighbor3;
  triangle ptr;                         /* Temporary variable used by sym(). */

#ifdef TRILIBRARY
  if (!b->quiet) {
    printf("Writing neighbors.\n");
  }
  /* Allocate memory for neighbors if necessary. */
  if (*neighborlist == (int *) NULL) {
    *neighborlist = (int *) trimalloc((int) (m->triangles.items * 3 *
                                             sizeof(int)));
  }
  nlist = *neighborlist;
  index = 0;
#else /* not TRILIBRARY */
  if (!b->quiet) {
    printf("Writing %s.\n", neighborfilename);
  }
  outfile = fopen(neighborfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", neighborfilename);
    triexit(1);
  }
  /* Number of triangles, three neighbors per triangle. */
  fprintf(outfile, "%ld  %d\n", m->triangles.items, 3);
#endif /* not TRILIBRARY */

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  elementnumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
    * (int *) (triangleloop.tri + 6) = (int) elementnumber;
    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }
  * (int *) (m->dummytri + 6) = -1;

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  elementnumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
    triangleloop.orient = 1;
    sym(triangleloop, trisym);
    neighbor1 = * (int *) (trisym.tri + 6);
    triangleloop.orient = 2;
    sym(triangleloop, trisym);
    neighbor2 = * (int *) (trisym.tri + 6);
    triangleloop.orient = 0;
    sym(triangleloop, trisym);
    neighbor3 = * (int *) (trisym.tri + 6);
#ifdef TRILIBRARY
    nlist[index++] = neighbor1;
    nlist[index++] = neighbor2;
    nlist[index++] = neighbor3;
#else /* not TRILIBRARY */
    /* Triangle number, neighboring triangle numbers. */
    fprintf(outfile, "%4ld    %d  %d  %d\n", elementnumber,
            neighbor1, neighbor2, neighbor3);
#endif /* not TRILIBRARY */

    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }

#ifndef TRILIBRARY
  finishfile(outfile, argc, argv);
#endif /* not TRILIBRARY */
}

/*****************************************************************************/
/*                                                                           */
/*  writeoff()   Write the triangulation to an .off file.                    */
/*                                                                           */
/*  OFF stands for the Object File Format, a format used by the Geometry     */
/*  Center's Geomview package.                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef TRILIBRARY

#ifdef ANSI_DECLARATORS
void writeoff(struct mesh *m, struct behavior *b, char *offfilename,
              int argc, char **argv)
#else /* not ANSI_DECLARATORS */
void writeoff(m, b, offfilename, argc, argv)
struct mesh *m;
struct behavior *b;
char *offfilename;
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

{
  FILE *outfile;
  struct otri triangleloop;
  vertex vertexloop;
  vertex p1, p2, p3;
  long outvertices;

  if (!b->quiet) {
    printf("Writing %s.\n", offfilename);
  }

  if (b->jettison) {
    outvertices = m->vertices.items - m->undeads;
  } else {
    outvertices = m->vertices.items;
  }

  outfile = fopen(offfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("  Error:  Cannot create file %s.\n", offfilename);
    triexit(1);
  }
  /* Number of vertices, triangles, and edges. */
  fprintf(outfile, "OFF\n%ld  %ld  %ld\n", outvertices, m->triangles.items,
          m->edges);

  /* Write the vertices. */
  traversalinit(&m->vertices);
  vertexloop = vertextraverse(m);
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      /* The "0.0" is here because the OFF format uses 3D coordinates. */
      fprintf(outfile, " %.17g  %.17g  %.17g\n", vertexloop[0], vertexloop[1],
              0.0);
    }
    vertexloop = vertextraverse(m);
  }

  /* Write the triangles. */
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, p1);
    dest(triangleloop, p2);
    apex(triangleloop, p3);
    /* The "3" means a three-vertex polygon. */
    fprintf(outfile, " 3   %4d  %4d  %4d\n", vertexmark(p1) - b->firstnumber,
            vertexmark(p2) - b->firstnumber, vertexmark(p3) - b->firstnumber);
    triangleloop.tri = triangletraverse(m);
  }
  finishfile(outfile, argc, argv);
}

#endif /* not TRILIBRARY */

/**                                                                         **/
/**                                                                         **/
/********* File I/O routines end here                                *********/

/*****************************************************************************/
/*                                                                           */
/*  quality_statistics()   Print statistics about the quality of the mesh.   */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void quality_statistics(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void quality_statistics(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  struct otri triangleloop;
  vertex p[3];
  REAL cossquaretable[8];
  REAL ratiotable[16];
  REAL dx[3], dy[3];
  REAL edgelength[3];
  REAL dotproduct;
  REAL cossquare;
  REAL triarea;
  REAL shortest, longest;
  REAL trilongest2;
  REAL smallestarea, biggestarea;
  REAL triminaltitude2;
  REAL minaltitude;
  REAL triaspect2;
  REAL worstaspect;
  REAL smallestangle, biggestangle;
  REAL radconst, degconst;
  int angletable[18];
  int aspecttable[16];
  int aspectindex;
  int tendegree;
  int acutebiggest;
  int i, ii, j, k;

  printf("Mesh quality statistics:\n\n");
  radconst = PI / 18.0;
  degconst = 180.0 / PI;
  for (i = 0; i < 8; i++) {
    cossquaretable[i] = cos(radconst * (REAL) (i + 1));
    cossquaretable[i] = cossquaretable[i] * cossquaretable[i];
  }
  for (i = 0; i < 18; i++) {
    angletable[i] = 0;
  }

  ratiotable[0]  =      1.5;      ratiotable[1]  =     2.0;
  ratiotable[2]  =      2.5;      ratiotable[3]  =     3.0;
  ratiotable[4]  =      4.0;      ratiotable[5]  =     6.0;
  ratiotable[6]  =     10.0;      ratiotable[7]  =    15.0;
  ratiotable[8]  =     25.0;      ratiotable[9]  =    50.0;
  ratiotable[10] =    100.0;      ratiotable[11] =   300.0;
  ratiotable[12] =   1000.0;      ratiotable[13] = 10000.0;
  ratiotable[14] = 100000.0;      ratiotable[15] =     0.0;
  for (i = 0; i < 16; i++) {
    aspecttable[i] = 0;
  }

  worstaspect = 0.0;
  minaltitude = m->xmax - m->xmin + m->ymax - m->ymin;
  minaltitude = minaltitude * minaltitude;
  shortest = minaltitude;
  longest = 0.0;
  smallestarea = minaltitude;
  biggestarea = 0.0;
  worstaspect = 0.0;
  smallestangle = 0.0;
  biggestangle = 2.0;
  acutebiggest = 1;

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, p[0]);
    dest(triangleloop, p[1]);
    apex(triangleloop, p[2]);
    trilongest2 = 0.0;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dx[i] = p[j][0] - p[k][0];
      dy[i] = p[j][1] - p[k][1];
      edgelength[i] = dx[i] * dx[i] + dy[i] * dy[i];
      if (edgelength[i] > trilongest2) {
        trilongest2 = edgelength[i];
      }
      if (edgelength[i] > longest) {
        longest = edgelength[i];
      }
      if (edgelength[i] < shortest) {
        shortest = edgelength[i];
      }
    }

    triarea = counterclockwise(m, b, p[0], p[1], p[2]);
    if (triarea < smallestarea) {
      smallestarea = triarea;
    }
    if (triarea > biggestarea) {
      biggestarea = triarea;
    }
    triminaltitude2 = triarea * triarea / trilongest2;
    if (triminaltitude2 < minaltitude) {
      minaltitude = triminaltitude2;
    }
    triaspect2 = trilongest2 / triminaltitude2;
    if (triaspect2 > worstaspect) {
      worstaspect = triaspect2;
    }
    aspectindex = 0;
    while ((triaspect2 > ratiotable[aspectindex] * ratiotable[aspectindex])
           && (aspectindex < 15)) {
      aspectindex++;
    }
    aspecttable[aspectindex]++;

    for (i = 0; i < 3; i++) {
      j = plus1mod3[i];
      k = minus1mod3[i];
      dotproduct = dx[j] * dx[k] + dy[j] * dy[k];
      cossquare = dotproduct * dotproduct / (edgelength[j] * edgelength[k]);
      tendegree = 8;
      for (ii = 7; ii >= 0; ii--) {
        if (cossquare > cossquaretable[ii]) {
          tendegree = ii;
        }
      }
      if (dotproduct <= 0.0) {
        angletable[tendegree]++;
        if (cossquare > smallestangle) {
          smallestangle = cossquare;
        }
        if (acutebiggest && (cossquare < biggestangle)) {
          biggestangle = cossquare;
        }
      } else {
        angletable[17 - tendegree]++;
        if (acutebiggest || (cossquare > biggestangle)) {
          biggestangle = cossquare;
          acutebiggest = 0;
        }
      }
    }
    triangleloop.tri = triangletraverse(m);
  }

  shortest = sqrt(shortest);
  longest = sqrt(longest);
  minaltitude = sqrt(minaltitude);
  worstaspect = sqrt(worstaspect);
  smallestarea *= 0.5;
  biggestarea *= 0.5;
  if (smallestangle >= 1.0) {
    smallestangle = 0.0;
  } else {
    smallestangle = degconst * acos(sqrt(smallestangle));
  }
  if (biggestangle >= 1.0) {
    biggestangle = 180.0;
  } else {
    if (acutebiggest) {
      biggestangle = degconst * acos(sqrt(biggestangle));
    } else {
      biggestangle = 180.0 - degconst * acos(sqrt(biggestangle));
    }
  }

  printf("  Smallest area: %16.5g   |  Largest area: %16.5g\n",
         smallestarea, biggestarea);
  printf("  Shortest edge: %16.5g   |  Longest edge: %16.5g\n",
         shortest, longest);
  printf("  Shortest altitude: %12.5g   |  Largest aspect ratio: %8.5g\n\n",
         minaltitude, worstaspect);

  printf("  Triangle aspect ratio histogram:\n");
  printf("  1.1547 - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
         ratiotable[0], aspecttable[0], ratiotable[7], ratiotable[8],
         aspecttable[8]);
  for (i = 1; i < 7; i++) {
    printf("  %6.6g - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
           ratiotable[i - 1], ratiotable[i], aspecttable[i],
           ratiotable[i + 7], ratiotable[i + 8], aspecttable[i + 8]);
  }
  printf("  %6.6g - %-6.6g    :  %8d    | %6.6g -            :  %8d\n",
         ratiotable[6], ratiotable[7], aspecttable[7], ratiotable[14],
         aspecttable[15]);
  printf("  (Aspect ratio is longest edge divided by shortest altitude)\n\n");

  printf("  Smallest angle: %15.5g   |  Largest angle: %15.5g\n\n",
         smallestangle, biggestangle);

  printf("  Angle histogram:\n");
  for (i = 0; i < 9; i++) {
    printf("    %3d - %3d degrees:  %8d    |    %3d - %3d degrees:  %8d\n",
           i * 10, i * 10 + 10, angletable[i],
           i * 10 + 90, i * 10 + 100, angletable[i + 9]);
  }
  printf("\n");
}

/*****************************************************************************/
/*                                                                           */
/*  statistics()   Print all sorts of cool facts.                            */
/*                                                                           */
/*****************************************************************************/

#ifdef ANSI_DECLARATORS
void statistics(struct mesh *m, struct behavior *b)
#else /* not ANSI_DECLARATORS */
void statistics(m, b)
struct mesh *m;
struct behavior *b;
#endif /* not ANSI_DECLARATORS */

{
  printf("\nStatistics:\n\n");
  printf("  Input vertices: %d\n", m->invertices);
  if (b->refine) {
    printf("  Input triangles: %d\n", m->inelements);
  }
  if (b->poly) {
    printf("  Input segments: %d\n", m->insegments);
    if (!b->refine) {
      printf("  Input holes: %d\n", m->holes);
    }
  }

  printf("\n  Mesh vertices: %ld\n", m->vertices.items - m->undeads);
  printf("  Mesh triangles: %ld\n", m->triangles.items);
  printf("  Mesh edges: %ld\n", m->edges);
  printf("  Mesh exterior boundary edges: %ld\n", m->hullsize);
  if (b->poly || b->refine) {
    printf("  Mesh interior boundary edges: %ld\n",
           m->subsegs.items - m->hullsize);
    printf("  Mesh subsegments (constrained edges): %ld\n",
           m->subsegs.items);
  }
  printf("\n");

  if (b->verbose) {
    quality_statistics(m, b);
    printf("Memory allocation statistics:\n\n");
    printf("  Maximum number of vertices: %ld\n", m->vertices.maxitems);
    printf("  Maximum number of triangles: %ld\n", m->triangles.maxitems);
    if (m->subsegs.maxitems > 0) {
      printf("  Maximum number of subsegments: %ld\n", m->subsegs.maxitems);
    }
    if (m->viri.maxitems > 0) {
      printf("  Maximum number of viri: %ld\n", m->viri.maxitems);
    }
    if (m->badsubsegs.maxitems > 0) {
      printf("  Maximum number of encroached subsegments: %ld\n",
             m->badsubsegs.maxitems);
    }
    if (m->badtriangles.maxitems > 0) {
      printf("  Maximum number of bad triangles: %ld\n",
             m->badtriangles.maxitems);
    }
    if (m->flipstackers.maxitems > 0) {
      printf("  Maximum number of stacked triangle flips: %ld\n",
             m->flipstackers.maxitems);
    }
    if (m->splaynodes.maxitems > 0) {
      printf("  Maximum number of splay tree nodes: %ld\n",
             m->splaynodes.maxitems);
    }
    printf("  Approximate heap memory use (bytes): %ld\n\n",
           m->vertices.maxitems * m->vertices.itembytes +
           m->triangles.maxitems * m->triangles.itembytes +
           m->subsegs.maxitems * m->subsegs.itembytes +
           m->viri.maxitems * m->viri.itembytes +
           m->badsubsegs.maxitems * m->badsubsegs.itembytes +
           m->badtriangles.maxitems * m->badtriangles.itembytes +
           m->flipstackers.maxitems * m->flipstackers.itembytes +
           m->splaynodes.maxitems * m->splaynodes.itembytes);

    printf("Algorithmic statistics:\n\n");
    if (!b->weighted) {
      printf("  Number of incircle tests: %ld\n", m->incirclecount);
    } else {
      printf("  Number of 3D orientation tests: %ld\n", m->orient3dcount);
    }
    printf("  Number of 2D orientation tests: %ld\n", m->counterclockcount);
    if (m->hyperbolacount > 0) {
      printf("  Number of right-of-hyperbola tests: %ld\n",
             m->hyperbolacount);
    }
    if (m->circletopcount > 0) {
      printf("  Number of circle top computations: %ld\n",
             m->circletopcount);
    }
    if (m->circumcentercount > 0) {
      printf("  Number of triangle circumcenter computations: %ld\n",
             m->circumcentercount);
    }
    printf("\n");
  }
}

/*****************************************************************************/
/*                                                                           */
/*  main() or triangulate()   Gosh, do everything.                           */
/*                                                                           */
/*  The sequence is roughly as follows.  Many of these steps can be skipped, */
/*  depending on the command line switches.                                  */
/*                                                                           */
/*  - Initialize constants and parse the command line.                       */
/*  - Read the vertices from a file and either                               */
/*    - triangulate them (no -r), or                                         */
/*    - read an old mesh from files and reconstruct it (-r).                 */
/*  - Insert the PSLG segments (-p), and possibly segments on the convex     */
/*      hull (-c).                                                           */
/*  - Read the holes (-p), regional attributes (-pA), and regional area      */
/*      constraints (-pa).  Carve the holes and concavities, and spread the  */
/*      regional attributes and area constraints.                            */
/*  - Enforce the constraints on minimum angle (-q) and maximum area (-a).   */
/*      Also enforce the conforming Delaunay property (-q and -a).           */
/*  - Compute the number of edges in the resulting mesh.                     */
/*  - Promote the mesh's linear triangles to higher order elements (-o).     */
/*  - Write the output files and print the statistics.                       */
/*  - Check the consistency and Delaunay property of the mesh (-C).          */
/*                                                                           */
/*****************************************************************************/

#ifdef TRILIBRARY

#ifdef ANSI_DECLARATORS
void triangulate(char *triswitches, struct triangulateio *in,
                 struct triangulateio *out, struct triangulateio *vorout)
#else /* not ANSI_DECLARATORS */
void triangulate(triswitches, in, out, vorout)
char *triswitches;
struct triangulateio *in;
struct triangulateio *out;
struct triangulateio *vorout;
#endif /* not ANSI_DECLARATORS */

#else /* not TRILIBRARY */

#ifdef ANSI_DECLARATORS
int main(int argc, char **argv)
#else /* not ANSI_DECLARATORS */
int main(argc, argv)
int argc;
char **argv;
#endif /* not ANSI_DECLARATORS */

#endif /* not TRILIBRARY */

{
  struct mesh m;
  struct behavior b;
  REAL *holearray;                                        /* Array of holes. */
  REAL *regionarray;   /* Array of regional attributes and area constraints. */
#ifndef TRILIBRARY
  FILE *polyfile;
#endif /* not TRILIBRARY */
#ifndef NO_TIMER
  /* Variables for timing the performance of Triangle.  The types are */
  /*   defined in sys/time.h.                                         */
  struct timeval tv0, tv1, tv2, tv3, tv4, tv5, tv6;
  struct timezone tz;
#endif /* not NO_TIMER */

#ifndef NO_TIMER
  gettimeofday(&tv0, &tz);
#endif /* not NO_TIMER */

  triangleinit(&m);
#ifdef TRILIBRARY
  parsecommandline(1, &triswitches, &b);
#else /* not TRILIBRARY */
  parsecommandline(argc, argv, &b);
#endif /* not TRILIBRARY */
  m.steinerleft = b.steiner;

#ifdef TRILIBRARY
  transfernodes(&m, &b, in->pointlist, in->pointattributelist,
                in->pointmarkerlist, in->numberofpoints,
                in->numberofpointattributes);
#else /* not TRILIBRARY */
  readnodes(&m, &b, b.innodefilename, b.inpolyfilename, &polyfile);
#endif /* not TRILIBRARY */

#ifndef NO_TIMER
  if (!b.quiet) {
    gettimeofday(&tv1, &tz);
  }
#endif /* not NO_TIMER */

#ifdef CDT_ONLY
  m.hullsize = delaunay(&m, &b);                /* Triangulate the vertices. */
#else /* not CDT_ONLY */
  if (b.refine) {
    /* Read and reconstruct a mesh. */
#ifdef TRILIBRARY
    m.hullsize = reconstruct(&m, &b, in->trianglelist,
                             in->triangleattributelist, in->trianglearealist,
                             in->numberoftriangles, in->numberofcorners,
                             in->numberoftriangleattributes,
                             in->segmentlist, in->segmentmarkerlist,
                             in->numberofsegments);
#else /* not TRILIBRARY */
    m.hullsize = reconstruct(&m, &b, b.inelefilename, b.areafilename,
                             b.inpolyfilename, polyfile);
#endif /* not TRILIBRARY */
  } else {
    m.hullsize = delaunay(&m, &b);              /* Triangulate the vertices. */
  }
#endif /* not CDT_ONLY */

#ifndef NO_TIMER
  if (!b.quiet) {
    gettimeofday(&tv2, &tz);
    if (b.refine) {
      printf("Mesh reconstruction");
    } else {
      printf("Delaunay");
    }
    printf(" milliseconds:  %ld\n", 1000l * (tv2.tv_sec - tv1.tv_sec) +
           (tv2.tv_usec - tv1.tv_usec) / 1000l);
  }
#endif /* not NO_TIMER */

  /* Ensure that no vertex can be mistaken for a triangular bounding */
  /*   box vertex in insertvertex().                                 */
  m.infvertex1 = (vertex) NULL;
  m.infvertex2 = (vertex) NULL;
  m.infvertex3 = (vertex) NULL;

  if (b.usesegments) {
    m.checksegments = 1;                /* Segments will be introduced next. */
    if (!b.refine) {
      /* Insert PSLG segments and/or convex hull segments. */
#ifdef TRILIBRARY
      formskeleton(&m, &b, in->segmentlist,
                   in->segmentmarkerlist, in->numberofsegments);
#else /* not TRILIBRARY */
      formskeleton(&m, &b, polyfile, b.inpolyfilename);
#endif /* not TRILIBRARY */
    }
  }

#ifndef NO_TIMER
  if (!b.quiet) {
    gettimeofday(&tv3, &tz);
    if (b.usesegments && !b.refine) {
      printf("Segment milliseconds:  %ld\n",
             1000l * (tv3.tv_sec - tv2.tv_sec) +
             (tv3.tv_usec - tv2.tv_usec) / 1000l);
    }
  }
#endif /* not NO_TIMER */

  if (b.poly && (m.triangles.items > 0)) {
#ifdef TRILIBRARY
    holearray = in->holelist;
    m.holes = in->numberofholes;
    regionarray = in->regionlist;
    m.regions = in->numberofregions;
#else /* not TRILIBRARY */
    readholes(&m, &b, polyfile, b.inpolyfilename, &holearray, &m.holes,
              &regionarray, &m.regions);
#endif /* not TRILIBRARY */
    if (!b.refine) {
      /* Carve out holes and concavities. */
      carveholes(&m, &b, holearray, m.holes, regionarray, m.regions);
    }
  } else {
    /* Without a PSLG, there can be no holes or regional attributes   */
    /*   or area constraints.  The following are set to zero to avoid */
    /*   an accidental free() later.                                  */
    m.holes = 0;
    m.regions = 0;
  }

#ifndef NO_TIMER
  if (!b.quiet) {
    gettimeofday(&tv4, &tz);
    if (b.poly && !b.refine) {
      printf("Hole milliseconds:  %ld\n", 1000l * (tv4.tv_sec - tv3.tv_sec) +
             (tv4.tv_usec - tv3.tv_usec) / 1000l);
    }
  }
#endif /* not NO_TIMER */

#ifndef CDT_ONLY
  if (b.quality && (m.triangles.items > 0)) {
    enforcequality(&m, &b);           /* Enforce angle and area constraints. */
  }
#endif /* not CDT_ONLY */

#ifndef NO_TIMER
  if (!b.quiet) {
    gettimeofday(&tv5, &tz);
#ifndef CDT_ONLY
    if (b.quality) {
      printf("Quality milliseconds:  %ld\n",
             1000l * (tv5.tv_sec - tv4.tv_sec) +
             (tv5.tv_usec - tv4.tv_usec) / 1000l);
    }
#endif /* not CDT_ONLY */
  }
#endif /* not NO_TIMER */

  /* Calculate the number of edges. */
  m.edges = (3l * m.triangles.items + m.hullsize) / 2l;

  if (b.order > 1) {
    highorder(&m, &b);       /* Promote elements to higher polynomial order. */
  }
  if (!b.quiet) {
    printf("\n");
  }

#ifdef TRILIBRARY
  if (b.jettison) {
    out->numberofpoints = m.vertices.items - m.undeads;
  } else {
    out->numberofpoints = m.vertices.items;
  }
  out->numberofpointattributes = m.nextras;
  out->numberoftriangles = m.triangles.items;
  out->numberofcorners = (b.order + 1) * (b.order + 2) / 2;
  out->numberoftriangleattributes = m.eextras;
  out->numberofedges = m.edges;
  if (b.usesegments) {
    out->numberofsegments = m.subsegs.items;
  } else {
    out->numberofsegments = m.hullsize;
  }
  if (vorout != (struct triangulateio *) NULL) {
    vorout->numberofpoints = m.triangles.items;
    vorout->numberofpointattributes = m.nextras;
    vorout->numberofedges = m.edges;
  }
#endif /* TRILIBRARY */
  /* If not using iteration numbers, don't write a .node file if one was */
  /*   read, because the original one would be overwritten!              */
  if (b.nonodewritten || (b.noiterationnum && m.readnodefile)) {
    if (!b.quiet) {
#ifdef TRILIBRARY
      printf("NOT writing vertices.\n");
#else /* not TRILIBRARY */
      printf("NOT writing a .node file.\n");
#endif /* not TRILIBRARY */
    }
    numbernodes(&m, &b);         /* We must remember to number the vertices. */
  } else {
    /* writenodes() numbers the vertices too. */
#ifdef TRILIBRARY
    writenodes(&m, &b, &out->pointlist, &out->pointattributelist,
               &out->pointmarkerlist);
#else /* not TRILIBRARY */
    writenodes(&m, &b, b.outnodefilename, argc, argv);
#endif /* TRILIBRARY */
  }
  if (b.noelewritten) {
    if (!b.quiet) {
#ifdef TRILIBRARY
      printf("NOT writing triangles.\n");
#else /* not TRILIBRARY */
      printf("NOT writing an .ele file.\n");
#endif /* not TRILIBRARY */
    }
  } else {
#ifdef TRILIBRARY
    writeelements(&m, &b, &out->trianglelist, &out->triangleattributelist);
#else /* not TRILIBRARY */
    writeelements(&m, &b, b.outelefilename, argc, argv);
#endif /* not TRILIBRARY */
  }
  /* The -c switch (convex switch) causes a PSLG to be written */
  /*   even if none was read.                                  */
  if (b.poly || b.convex) {
    /* If not using iteration numbers, don't overwrite the .poly file. */
    if (b.nopolywritten || b.noiterationnum) {
      if (!b.quiet) {
#ifdef TRILIBRARY
        printf("NOT writing segments.\n");
#else /* not TRILIBRARY */
        printf("NOT writing a .poly file.\n");
#endif /* not TRILIBRARY */
      }
    } else {
#ifdef TRILIBRARY
      writepoly(&m, &b, &out->segmentlist, &out->segmentmarkerlist);
      out->numberofholes = m.holes;
      out->numberofregions = m.regions;
      if (b.poly) {
        out->holelist = in->holelist;
        out->regionlist = in->regionlist;
      } else {
        out->holelist = (REAL *) NULL;
        out->regionlist = (REAL *) NULL;
      }
#else /* not TRILIBRARY */
      writepoly(&m, &b, b.outpolyfilename, holearray, m.holes, regionarray,
                m.regions, argc, argv);
#endif /* not TRILIBRARY */
    }
  }
#ifndef TRILIBRARY
#ifndef CDT_ONLY
  if (m.regions > 0) {
    trifree((VOID *) regionarray);
  }
#endif /* not CDT_ONLY */
  if (m.holes > 0) {
    trifree((VOID *) holearray);
  }
  if (b.geomview) {
    writeoff(&m, &b, b.offfilename, argc, argv);
  }
#endif /* not TRILIBRARY */
  if (b.edgesout) {
#ifdef TRILIBRARY
    writeedges(&m, &b, &out->edgelist, &out->edgemarkerlist);
#else /* not TRILIBRARY */
    writeedges(&m, &b, b.edgefilename, argc, argv);
#endif /* not TRILIBRARY */
  }
  if (b.voronoi) {
#ifdef TRILIBRARY
    writevoronoi(&m, &b, &vorout->pointlist, &vorout->pointattributelist,
                 &vorout->pointmarkerlist, &vorout->edgelist,
                 &vorout->edgemarkerlist, &vorout->normlist);
#else /* not TRILIBRARY */
    writevoronoi(&m, &b, b.vnodefilename, b.vedgefilename, argc, argv);
#endif /* not TRILIBRARY */
  }
  if (b.neighbors) {
#ifdef TRILIBRARY
    writeneighbors(&m, &b, &out->neighborlist);
#else /* not TRILIBRARY */
    writeneighbors(&m, &b, b.neighborfilename, argc, argv);
#endif /* not TRILIBRARY */
  }

  if (!b.quiet) {
#ifndef NO_TIMER
    gettimeofday(&tv6, &tz);
    printf("\nOutput milliseconds:  %ld\n",
           1000l * (tv6.tv_sec - tv5.tv_sec) +
           (tv6.tv_usec - tv5.tv_usec) / 1000l);
    printf("Total running milliseconds:  %ld\n",
           1000l * (tv6.tv_sec - tv0.tv_sec) +
           (tv6.tv_usec - tv0.tv_usec) / 1000l);
#endif /* not NO_TIMER */

    statistics(&m, &b);
  }

#ifndef REDUCED
  if (b.docheck) {
    checkmesh(&m, &b);
    checkdelaunay(&m, &b);
  }
#endif /* not REDUCED */

  triangledeinit(&m, &b);
#ifndef TRILIBRARY
  return 0;
#endif /* not TRILIBRARY */
}
