/*****************************************************************************/
/*                                                                           */
/*    ,d88^^o 888                                   o    o                   */
/*    8888    888o^88,  o88^^o Y88b    o    /      d8b  d8b      o88^^8o     */
/*    "Y88b   888  888 d888   b Y88b  d8b  /      d888bdY88b    d888  88b    */
/*     "Y88b, 888  888 8888   8  Y888/Y88b/      / Y88Y Y888b   8888oo888    */
/*    o  8888 888  888 q888   p   Y8/  Y8/      /   YY   Y888b  q888         */
/*    "oo88P" 888  888  "88oo"     Y    Y      /          Y888b  "88oooo"    */
/*                                                                           */
/*  A Display Program for Meshes and More.                                   */
/*  (showme.c)                                                               */
/*                                                                           */
/*  Version 1.6                                                              */
/*  July 28, 2005                                                            */
/*                                                                           */
/*  Copyright 1996, 1998, 2005                                               */
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
/*      http://www.cs.cmu.edu/~quake/showme.html                             */
/*                                                                           */
/*  Show Me was created as part of the Archimedes project in the School of   */
/*    Computer Science at Carnegie Mellon University.  Archimedes is a       */
/*    system for compiling parallel finite element solvers.  For further     */
/*    information, see Anja Feldmann, Omar Ghattas, John R. Gilbert, Gary L. */
/*    Miller, David R. O'Hallaron, Eric J. Schwabe, Jonathan R. Shewchuk,    */
/*    and Shang-Hua Teng.  "Automated Parallel Solution of Unstructured PDE  */
/*    Problems."  To appear in Communications of the ACM, we hope.           */
/*                                                                           */
/*  If you make any improvements to this code, please please please let me   */
/*    know, so that I may obtain the improvements.  Even if you don't change */
/*    the code, I'd still love to hear what it's being used for.             */
/*                                                                           */
/*  Disclaimer:  Neither I nor Carnegie Mellon warrant this code in any way  */
/*    whatsoever.  Use at your own risk.                                     */
/*                                                                           */
/*****************************************************************************/

/* For single precision (which will save some memory and reduce paging),     */
/*   write "#define SINGLE" below.                                           */
/*                                                                           */
/* For double precision (which will allow you to display triangulations of   */
/*   a finer resolution), leave SINGLE undefined.                            */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif

/* Maximum number of characters in a file name (including the null).         */

#define FILENAMESIZE 2048

/* Maximum number of characters in a line read from a file (including the    */
/*   null).                                                                  */

#define INPUTLINESIZE 1024

#define STARTWIDTH 414
#define STARTHEIGHT 414
#define MINWIDTH 50
#define MINHEIGHT 50
#define BUTTONHEIGHT 21
#define BUTTONROWS 3
#define PANELHEIGHT (BUTTONHEIGHT * BUTTONROWS)
#define MAXCOLORS 64

#define IMAGE_TYPES 7
#define NOTHING -1
#define NODE 0
#define POLY 1
#define ELE 2
#define EDGE 3
#define PART 4
#define ADJ 5
#define VORO 6

#define STARTEXPLOSION 0.5

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

/* A necessary forward declaration.                                          */

int load_image();

Display *display;
int screen;
Window rootwindow;
Window mainwindow;
Window quitwin;
Window leftwin;
Window rightwin;
Window upwin;
Window downwin;
Window resetwin;
Window pswin;
Window epswin;
Window expwin;
Window exppluswin;
Window expminuswin;
Window widthpluswin;
Window widthminuswin;
Window versionpluswin;
Window versionminuswin;
Window fillwin;
Window nodewin[2];
Window polywin[2];
Window elewin[2];
Window edgewin[2];
Window partwin[2];
Window adjwin[2];
Window voronoiwin[2];

int windowdepth;
XEvent event;
Colormap rootmap;
XFontStruct *font;
int width, height;
int black, white;
int showme_foreground;
GC fontgc;
GC blackfontgc;
GC linegc;
GC trianglegc;
int colors[MAXCOLORS];
XColor rgb[MAXCOLORS];
int color;

int start_image, current_image;
int start_inc, current_inc;
int loweriteration;
int line_width;
int loaded[2][IMAGE_TYPES];
REAL xlo[2][IMAGE_TYPES], ylo[2][IMAGE_TYPES];
REAL xhi[2][IMAGE_TYPES], yhi[2][IMAGE_TYPES];
REAL xscale, yscale;
REAL xoffset, yoffset;
int zoom;

int nodes[2], node_dim[2];
REAL *nodeptr[2];
int polynodes[2], poly_dim[2], polyedges[2], polyholes[2];
REAL *polynodeptr[2], *polyholeptr[2];
int *polyedgeptr[2];
int elems[2], ele_corners[2];
int *eleptr[2];
int edges[2];
int *edgeptr[2];
REAL *normptr[2];
int subdomains[2];
int *partpart[2];
REAL *partcenter[2], *partshift[2];
int adjsubdomains[2];
int *adjptr[2];
int vnodes[2], vnode_dim[2];
REAL *vnodeptr[2];
int vedges[2];
int *vedgeptr[2];
REAL *vnormptr[2];
int firstnumber[2];

int quiet, fillelem, bw_ps, explode;
REAL explosion;

char filename[FILENAMESIZE];
char nodefilename[2][FILENAMESIZE];
char polyfilename[2][FILENAMESIZE];
char elefilename[2][FILENAMESIZE];
char edgefilename[2][FILENAMESIZE];
char partfilename[2][FILENAMESIZE];
char adjfilename[2][FILENAMESIZE];
char vnodefilename[2][FILENAMESIZE];
char vedgefilename[2][FILENAMESIZE];

char *colorname[] = {"aquamarine", "red", "green yellow", "magenta",
                     "yellow", "green", "orange", "blue",
                     "white", "sandy brown", "cyan", "moccasin",
                     "cadet blue", "coral", "cornflower blue", "sky blue",
                     "firebrick", "forest green", "gold", "goldenrod",
                     "gray", "hot pink", "chartreuse", "pale violet red",
                     "indian red", "khaki", "lavender", "light blue",
                     "light gray", "light steel blue", "lime green", "azure",
                     "maroon", "medium aquamarine", "dodger blue", "honeydew",
                     "medium orchid", "medium sea green", "moccasin",
                     "medium slate blue", "medium spring green",
                     "medium turquoise", "medium violet red",
                     "orange red", "chocolate", "light goldenrod",
                     "orchid", "pale green", "pink", "plum",
                     "purple", "salmon", "sea green",
                     "sienna", "slate blue", "spring green",
                     "steel blue", "tan", "thistle", "turquoise",
                     "violet", "violet red", "wheat",
                     "yellow green"};

void syntax()
{
  printf("showme [-bfw_Qh] input_file\n");
  printf("    -b  Black and white PostScript (default is color).\n");
  printf("    -f  Fill triangles of partitioned mesh with color.\n");
  printf("    -w  Set line width to some specified number.\n");
  printf("    -Q  Quiet:  No terminal output except errors.\n");
  printf("    -h  Help:  Detailed instructions for Show Me.\n");
  exit(0);
}

void info()
{
  printf("Show Me\n");
  printf("A Display Program for Meshes and More.\n");
  printf("Version 1.6\n\n");
  printf(
"Copyright 1996 Jonathan Richard Shewchuk  (bugs/comments to jrs@cs.cmu.edu)\n"
);
  printf("School of Computer Science / Carnegie Mellon University\n");
  printf("5000 Forbes Avenue / Pittsburgh, Pennsylvania  15213-3891\n");
  printf(
"Created as part of the Archimedes project (tools for parallel FEM).\n");
  printf(
"Supported in part by NSF Grant CMS-9318163 and an NSERC 1967 Scholarship.\n");
  printf("There is no warranty whatsoever.  Use at your own risk.\n");
#ifdef SINGLE
  printf("This executable is compiled for single precision arithmetic.\n\n\n");
#else
  printf("This executable is compiled for double precision arithmetic.\n\n\n");
#endif
  printf(
"Show Me graphically displays the contents of geometric files, especially\n");
  printf(
"those generated by Triangle, my two-dimensional quality mesh generator and\n"
);
  printf(
"Delaunay triangulator.  Show Me can also write images in PostScript form.\n");
  printf(
"Show Me is also useful for checking the consistency of the files you create\n"
);
  printf(
"as input to Triangle; Show Me does these checks more thoroughly than\n");
  printf("Triangle does.  The command syntax is:\n\n");
  printf("showme [-bfw_Qh] input_file\n\n");
  printf(
"The underscore indicates that a number should follow the -w switch.\n");
  printf(
"input_file may be one of several types of file.  It must have extension\n");
  printf(
".node, .poly, .ele, .edge, .part, or .adj.  If no extension is provided,\n");
  printf(
"Show Me will assume the extension .ele.  A .node file represents a set of\n");
  printf(
"points; a .poly file represents a Planar Straight Line Graph; an .ele file\n"
);
  printf(
"(coupled with a .node file) represents the elements of a mesh or the\n");
  printf(
"triangles of a triangulation; an .edge file (coupled with a .node file)\n");
  printf(
"represents a set of edges; a .part file specifies a partition of a mesh;\n");
  printf(
"and a .adj file represents the adjacency graph defined by a partition.\n");
  printf("\n");
  printf("Command Line Switches:\n");
  printf("\n");
  printf(
"    -b  Makes all PostScript output black and white.  If this switch is not\n"
);
  printf(
"        selected, color PostScript is used for partitioned meshes and\n");
  printf("        adjacency graphs (.part and .adj files).\n");
  printf(
"    -f  On color displays and in color PostScript, displays partitioned\n");
  printf(
"        meshes by filling triangles with color, rather than by coloring the\n"
);
  printf(
"        edges.  This switch will result in a clearer picture if all\n");
  printf(
"        triangles are reasonably large, and a less clear picture if small\n");
  printf(
"        triangles are present.  (There is also a button to toggle this\n");
  printf("        behavior.)\n");
  printf(
"    -w  Followed by an integer, specifies the line width used in all\n");
  printf(
"        images.  (There are also buttons to change the line width.)\n");
  printf(
"    -Q  Quiet:  Suppresses all explanation of what Show Me is doing, unless\n"
);
  printf("        an error occurs.\n");
  printf("    -h  Help:  Displays these instructions.\n");
  printf("\n");
  printf("Controls:\n");
  printf("\n");
  printf(
"  To zoom in on an image, point at the location where you want a closer\n");
  printf(
"  look, and click the left mouse button.  To zoom out, click the right\n");
  printf(
"  mouse button.  In either case, the point you click on will be centered in\n"
);
  printf(
"  the window.  If you want to know the coordinates of a point, click the\n");
  printf(
"  middle mouse button; the coordinates will be printed on the terminal you\n"
);
  printf("  invoked Show Me from.\n\n");
  printf(
"  If you resize the window, the image will grow or shrink to match.\n");
  printf("\n");
  printf(
"  There is a panel of control buttons at the bottom of the Show Me window:\n"
);
  printf("\n");
  printf("  Quit:  Shuts down Show Me.\n");
  printf("  <, >, ^, v:  Moves the image in the indicated direction.\n");
  printf(
"  Reset: Unzooms and centers the image in the window.  When you switch from\n"
);
  printf(
"    one image to another, the viewing region does not change, so you may\n");
  printf(
"    need to reset the new image to make it fully visible.  This often is\n");
  printf(
"    the case when switching between Delaunay triangulations and their\n");
  printf(
"    corresponding Voronoi diagrams, as Voronoi vertices can be far from the\n"
);
  printf("    initial point set.\n");
  printf(
"  Width+, -:  Increases or decreases the width of all lines and points.\n");
  printf(
"  Exp, +, -:  These buttons appear only when you are viewing a partitioned\n"
);
  printf(
"    mesh (.part file).  `Exp' toggles between an exploded and non-exploded\n"
);
  printf(
"    image of the mesh.  The non-exploded image will not show the partition\n"
);
  printf(
"    on a black and white monitor.  `+' and `-' allow you to adjust the\n");
  printf(
"    spacing between pieces of the mesh to better distinguish them.\n");
  printf(
"  Fill:  This button appears only when you are viewing a partitioned mesh\n");
  printf(
"    (.part file).  It toggles between color-filled triangles and colored\n");
  printf(
"    edges (as the -f switch does).  Filled triangles look better when all\n");
  printf(
"    triangles are reasonably large; colored edges look better when there\n");
  printf("    are very small triangles present.\n");
  printf(
"  PS:  Creates a PostScript file containing the image you are viewing.  If\n"
);
  printf(
"    the -b switch is selected, all PostScript output will be black and\n");
  printf(
"    white; otherwise, .part.ps and .adj.ps files will be color, independent\n"
);
  printf(
"    of whether you are using a color monitor.  Normally the output will\n");
  printf(
"    preserve the properties of the image you see on the screen, including\n");
  printf(
"    zoom and line width; however, if black and white output is selected (-b\n"
);
  printf(
"    switch), partitioned meshes will always be drawn exploded.  The output\n"
);
  printf(
"    file name depends on the image being viewed.  If you want several\n");
  printf(
"    different snapshots (zooming in on different parts) of the same object,\n"
);
  printf(
"    you'll have to rename each file after Show Me creates it so that it\n");
  printf("    isn't overwritten by the next snapshot.\n");
  printf(
"  EPS:  Creates an encapsulated PostScript file, suitable for inclusion in\n"
);
  printf(
"    documents.  Otherwise, this button is just like the PS button.  (The\n");
  printf(
"    only difference is that .eps files lack a `showpage' command at the\n");
  printf("    end.)\n\n");
  printf(
"  There are two nearly-identical rows of buttons that load different images\n"
);
  printf("  from disk.  Each row contains the following buttons:\n\n");
  printf("  node:  Loads a .node file.\n");
  printf(
"  poly:  Loads a .poly file (and possibly an associated .node file).\n");
  printf("  ele:  Loads an .ele file (and associated .node file).\n");
  printf("  edge:  Loads an .edge file (and associated .node file).\n");
  printf(
"  part:  Loads a .part file (and associated .node and .ele files).\n");
  printf(
"  adj:  Loads an .adj file (and associated .node, .ele, and .part files).\n");
  printf("  voro:  Loads a .v.node and .v.edge file for a Voronoi diagram.\n");
  printf("\n");
  printf(
"  Each row represents a different iteration number of the geometry files.\n");
  printf(
"  For a full explanation of iteration numbers, read the instructions for\n");
  printf(
"  Triangle.  Briefly, iteration numbers are used to allow a user to easily\n"
);
  printf(
"  represent a sequence of related triangulations.  Iteration numbers are\n");
  printf(
"  used in the names of geometry files; for instance, mymesh.3.ele is a\n");
  printf(
"  triangle file with iteration number three, and mymesh.ele has an implicit\n"
);
  printf("  iteration number of zero.\n\n");
  printf(
"  The control buttons at the right end of each row display the two\n");
  printf(
"  iterations currently under view.  These buttons can be clicked to\n");
  printf(
"  increase or decrease the iteration numbers, and thus conveniently view\n");
  printf("  a sequence of meshes.\n\n");
  printf(
"  Show Me keeps each file in memory after loading it, but you can force\n");
  printf(
"  Show Me to reread a set of files (for one iteration number) by reclicking\n"
);
  printf(
"  the button that corresponds to the current image.  This is convenient if\n"
);
  printf("  you have changed a geometry file.\n\n");
  printf("File Formats:\n\n");
  printf(
"  All files may contain comments prefixed by the character '#'.  Points,\n");
  printf(
"  segments, holes, triangles, edges, and subdomains must be numbered\n");
  printf(
"  consecutively, starting from either 1 or 0.  Whichever you choose, all\n");
  printf(
"  input files must be consistent (for any single iteration number); if the\n"
);
  printf(
"  nodes are numbered from 1, so must be all other objects.  Show Me\n");
  printf(
"  automatically detects your choice while reading a .node (or .poly) file.\n"
);
  printf("  Examples of these file formats are given below.\n\n");
  printf("  .node files:\n");
  printf(
"    First line:  <# of points> <dimension (must be 2)> <# of attributes>\n");
  printf(
"                                           <# of boundary markers (0 or 1)>\n"
);
  printf(
"    Remaining lines:  <point #> <x> <y> [attributes] [boundary marker]\n");
  printf("\n");
  printf(
"    The attributes, which are typically floating-point values of physical\n");
  printf(
"    quantities (such as mass or conductivity) associated with the nodes of\n"
);
  printf(
"    a finite element mesh, are ignored by Show Me.  Show Me also ignores\n");
  printf(
"    boundary markers.  See the instructions for Triangle to find out what\n");
  printf("    attributes and boundary markers are.\n\n");
  printf("  .poly files:\n");
  printf(
"    First line:  <# of points> <dimension (must be 2)> <# of attributes>\n");
  printf(
"                                           <# of boundary markers (0 or 1)>\n"
);
  printf(
"    Following lines:  <point #> <x> <y> [attributes] [boundary marker]\n");
  printf("    One line:  <# of segments> <# of boundary markers (0 or 1)>\n");
  printf(
"    Following lines:  <segment #> <endpoint> <endpoint> [boundary marker]\n");
  printf("    One line:  <# of holes>\n");
  printf("    Following lines:  <hole #> <x> <y>\n");
  printf("    [Optional additional lines that are ignored]\n\n");
  printf(
"    A .poly file represents a Planar Straight Line Graph (PSLG), an idea\n");
  printf(
"    familiar to computational geometers.  By definition, a PSLG is just a\n");
  printf(
"    list of points and edges.  A .poly file also contains some additional\n");
  printf("    information.\n\n");
  printf(
"    The first section lists all the points, and is identical to the format\n"
);
  printf(
"    of .node files.  <# of points> may be set to zero to indicate that the\n"
);
  printf(
"    points are listed in a separate .node file; .poly files produced by\n");
  printf(
"    Triangle always have this format.  When Show Me reads such a file, it\n");
  printf("    also reads the corresponding .node file.\n\n");
  printf(
"    The second section lists the segments.  Segments are edges whose\n");
  printf(
"    presence in a triangulation produced from the PSLG is enforced.  Each\n");
  printf(
"    segment is specified by listing the indices of its two endpoints.  This\n"
);
  printf(
"    means that its endpoints must be included in the point list.  Each\n");
  printf(
"    segment, like each point, may have a boundary marker, which is ignored\n"
);
  printf("    by Show Me.\n\n");
  printf(
"    The third section lists holes and concavities that are desired in any\n");
  printf(
"    triangulation generated from the PSLG.  Holes are specified by\n");
  printf("    identifying a point inside each hole.\n\n");
  printf("  .ele files:\n");
  printf(
"    First line:  <# of triangles> <points per triangle> <# of attributes>\n");
  printf(
"    Remaining lines:  <triangle #> <point> <point> <point> ... [attributes]\n"
);
  printf("\n");
  printf(
"    Points are indices into the corresponding .node file.  Show Me ignores\n"
);
  printf(
"    all but the first three points of each triangle; these should be the\n");
  printf(
"    corners listed in counterclockwise order around the triangle.  The\n");
  printf("    attributes are ignored by Show Me.\n\n");
  printf("  .edge files:\n");
  printf("    First line:  <# of edges> <# of boundary markers (0 or 1)>\n");
  printf(
"    Following lines:  <edge #> <endpoint> <endpoint> [boundary marker]\n");
  printf("\n");
  printf(
"    Endpoints are indices into the corresponding .node file.  The boundary\n"
);
  printf("    markers are ignored by Show Me.\n\n");
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
  printf("  .part files:\n");
  printf("    First line:  <# of triangles> <# of subdomains>\n");
  printf("    Remaining lines:  <triangle #> <subdomain #>\n\n");
  printf(
"    The set of triangles is partitioned by a .part file; each triangle is\n");
  printf("    mapped to a subdomain.\n\n");
  printf("  .adj files:\n");
  printf("    First line:  <# of subdomains>\n");
  printf("    Remaining lines:  <adjacency matrix entry>\n\n");
  printf(
"    An .adj file represents adjacencies between subdomains (presumably\n");
  printf("    computed by a partitioner).  The first line is followed by\n");
  printf(
"    (subdomains X subdomains) lines, each containing one entry of the\n");
  printf(
"    adjacency matrix.  A nonzero entry indicates that two subdomains are\n");
  printf("    adjacent (share a point).\n\n");
  printf("Example:\n\n");
  printf(
"  Here is a sample file `box.poly' describing a square with a square hole:\n"
);
  printf("\n");
  printf(
"    # A box with eight points in 2D, no attributes, no boundary marker.\n");
  printf("    8 2 0 0\n");
  printf("    # Outer box has these vertices:\n");
  printf("     1   0 0\n");
  printf("     2   0 3\n");
  printf("     3   3 0\n");
  printf("     4   3 3\n");
  printf("    # Inner square has these vertices:\n");
  printf("     5   1 1\n");
  printf("     6   1 2\n");
  printf("     7   2 1\n");
  printf("     8   2 2\n");
  printf("    # Five segments without boundary markers.\n");
  printf("    5 0\n");
  printf("     1   1 2          # Left side of outer box.\n");
  printf("     2   5 7          # Segments 2 through 5 enclose the hole.\n");
  printf("     3   7 8\n");
  printf("     4   8 6\n");
  printf("     5   6 5\n");
  printf("    # One hole in the middle of the inner square.\n");
  printf("    1\n");
  printf("     1   1.5 1.5\n\n");
  printf(
"  After this PSLG is triangulated by Triangle, the resulting triangulation\n"
);
  printf(
"  consists of a .node and .ele file.  Here is the former, `box.1.node',\n");
  printf("  which duplicates the points of the PSLG:\n\n");
  printf("    8  2  0  0\n");
  printf("       1    0  0\n");
  printf("       2    0  3\n");
  printf("       3    3  0\n");
  printf("       4    3  3\n");
  printf("       5    1  1\n");
  printf("       6    1  2\n");
  printf("       7    2  1\n");
  printf("       8    2  2\n");
  printf("    # Generated by triangle -pcBev box\n");
  printf("\n");
  printf("  Here is the triangulation file, `box.1.ele'.\n");
  printf("\n");
  printf("    8  3  0\n");
  printf("       1       1     5     6\n");
  printf("       2       5     1     3\n");
  printf("       3       2     6     8\n");
  printf("       4       6     2     1\n");
  printf("       5       7     3     4\n");
  printf("       6       3     7     5\n");
  printf("       7       8     4     2\n");
  printf("       8       4     8     7\n");
  printf("    # Generated by triangle -pcBev box\n\n");
  printf("  Here is the edge file for the triangulation, `box.1.edge'.\n\n");
  printf("    16  0\n");
  printf("       1   1  5\n");
  printf("       2   5  6\n");
  printf("       3   6  1\n");
  printf("       4   1  3\n");
  printf("       5   3  5\n");
  printf("       6   2  6\n");
  printf("       7   6  8\n");
  printf("       8   8  2\n");
  printf("       9   2  1\n");
  printf("      10   7  3\n");
  printf("      11   3  4\n");
  printf("      12   4  7\n");
  printf("      13   7  5\n");
  printf("      14   8  4\n");
  printf("      15   4  2\n");
  printf("      16   8  7\n");
  printf("    # Generated by triangle -pcBev box\n");
  printf("\n");
  printf(
"  Here's a file `box.1.part' that partitions the mesh into four subdomains.\n"
);
  printf("\n");
  printf("    8  4\n");
  printf("       1    3\n");
  printf("       2    3\n");
  printf("       3    4\n");
  printf("       4    4\n");
  printf("       5    1\n");
  printf("       6    1\n");
  printf("       7    2\n");
  printf("       8    2\n");
  printf("    # Generated by slice -s4 box.1\n\n");
  printf(
"  Here's a file `box.1.adj' that represents the resulting adjacencies.\n");
  printf("\n");
  printf("    4\n");
  printf("      9\n");
  printf("      2\n");
  printf("      2\n");
  printf("      0\n");
  printf("      2\n");
  printf("      9\n");
  printf("      0\n");
  printf("      2\n");
  printf("      2\n");
  printf("      0\n");
  printf("      9\n");
  printf("      2\n");
  printf("      0\n");
  printf("      2\n");
  printf("      2\n");
  printf("      9\n");
  printf("\n");
  printf("Display Speed:\n");
  printf("\n");
  printf(
"  It is worthwhile to note that .edge files typically plot and print twice\n"
);
  printf(
"  as quickly as .ele files, because .ele files cause each internal edge to\n"
);
  printf(
"  be drawn twice.  For the same reason, PostScript files created from edge\n"
);
  printf("  sets are smaller than those created from triangulations.\n\n");
  printf("Show Me on the Web:\n\n");
  printf(
"  To see an illustrated, updated version of these instructions, check out\n");
  printf("\n");
  printf("    http://www.cs.cmu.edu/~quake/showme.html\n");
  printf("\n");
  printf("A Brief Plea:\n");
  printf("\n");
  printf(
"  If you use Show Me (or Triangle), and especially if you use it to\n");
  printf(
"  accomplish real work, I would like very much to hear from you.  A short\n");
  printf(
"  letter or email (to jrs@cs.cmu.edu) describing how you use Show Me (and\n");
  printf(
"  its sister programs) will mean a lot to me.  The more people I know\n");
  printf(
"  are using my programs, the more easily I can justify spending time on\n");
  printf(
"  improvements, which in turn will benefit you.  Also, I can put you\n");
  printf(
"  on a list to receive email whenever new versions are available.\n");
  printf("\n");
  printf(
"  If you use a PostScript file generated by Show Me in a publication,\n");
  printf("  please include an acknowledgment as well.\n\n");
  exit(0);
}

void set_filenames(filename, lowermeshnumber)
char *filename;
int lowermeshnumber;
{
  char numberstring[100];
  int i;

  for (i = 0; i < 2; i++) {
    strcpy(nodefilename[i], filename);
    strcpy(polyfilename[i], filename);
    strcpy(elefilename[i], filename);
    strcpy(edgefilename[i], filename);
    strcpy(partfilename[i], filename);
    strcpy(adjfilename[i], filename);
    strcpy(vnodefilename[i], filename);
    strcpy(vedgefilename[i], filename);

    if (lowermeshnumber + i > 0) {
      sprintf(numberstring, ".%d", lowermeshnumber + i);
      strcat(nodefilename[i], numberstring);
      strcat(polyfilename[i], numberstring);
      strcat(elefilename[i], numberstring);
      strcat(edgefilename[i], numberstring);
      strcat(partfilename[i], numberstring);
      strcat(adjfilename[i], numberstring);
      strcat(vnodefilename[i], numberstring);
      strcat(vedgefilename[i], numberstring);
    }

    strcat(nodefilename[i], ".node");
    strcat(polyfilename[i], ".poly");
    strcat(elefilename[i], ".ele");
    strcat(edgefilename[i], ".edge");
    strcat(partfilename[i], ".part");
    strcat(adjfilename[i], ".adj");
    strcat(vnodefilename[i], ".v.node");
    strcat(vedgefilename[i], ".v.edge");
  }
}

void parsecommandline(argc, argv)
int argc;
char **argv;
{
  int increment;
  int meshnumber;
  int i, j;

  quiet = 0;
  fillelem = 0;
  line_width = 1;
  bw_ps = 0;
  start_image = ELE;
  filename[0] = '\0';
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      for (j = 1; argv[i][j] != '\0'; j++) {
        if (argv[i][j] == 'f') {
          fillelem = 1;
	}
        if (argv[i][j] == 'w') {
          if ((argv[i][j + 1] >= '1') && (argv[i][j + 1] <= '9')) {
            line_width = 0;
            while ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) {
              j++;
              line_width = line_width * 10 + (int) (argv[i][j] - '0');
            }
            if (line_width > 100) {
              printf("Error:  Line width cannot exceed 100.\n");
              line_width = 1;
	    }
	  }
	}
        if (argv[i][j] == 'b') {
          bw_ps = 1;
	}
        if (argv[i][j] == 'Q') {
          quiet = 1;
	}
        if ((argv[i][j] == 'h') || (argv[i][j] == 'H') ||
            (argv[i][j] == '?')) {
          info();
	}
      }
    } else {
      strcpy(filename, argv[i]);
    }
  }
  if (filename[0] == '\0') {
    syntax();
  }
  if (!strcmp(&filename[strlen(filename) - 5], ".node")) {
    filename[strlen(filename) - 5] = '\0';
    start_image = NODE;
  }
  if (!strcmp(&filename[strlen(filename) - 5], ".poly")) {
    filename[strlen(filename) - 5] = '\0';
    start_image = POLY;
  }
  if (!strcmp(&filename[strlen(filename) - 4], ".ele")) {
    filename[strlen(filename) - 4] = '\0';
    start_image = ELE;
  }
  if (!strcmp(&filename[strlen(filename) - 5], ".edge")) {
    filename[strlen(filename) - 5] = '\0';
    start_image = EDGE;
  }
  if (!strcmp(&filename[strlen(filename) - 5], ".part")) {
    filename[strlen(filename) - 5] = '\0';
    start_image = PART;
  }
  if (!strcmp(&filename[strlen(filename) - 4], ".adj")) {
    filename[strlen(filename) - 4] = '\0';
    start_image = ADJ;
  }

  increment = 0;
  j = 1;
  while (filename[j] != '\0') {
    if ((filename[j] == '.') && (filename[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((filename[j] >= '0') && (filename[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (filename[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (filename[j] != '\0');
  }
  if (increment > 0) {
    filename[increment - 1] = '\0';
  }

  if (meshnumber == 0) {
    start_inc = 0;
    loweriteration = 0;
  } else {
    start_inc = 1;
    loweriteration = meshnumber - 1;
  }
  set_filenames(filename, loweriteration);
}

void free_inc(inc)
int inc;
{
  if (loaded[inc][NODE]) {
    free(nodeptr[inc]);
  }
  if (loaded[inc][POLY]) {
    if (polynodes[inc] > 0) {
      free(polynodeptr[inc]);
    }
    free(polyedgeptr[inc]);
    free(polyholeptr[inc]);
  }
  if (loaded[inc][ELE]) {
    free(eleptr[inc]);
  }
  if (loaded[inc][PART]) {
    free(partpart[inc]);
    free(partcenter[inc]);
    free(partshift[inc]);
  }
  if (loaded[inc][EDGE]) {
    free(edgeptr[inc]);
    free(normptr[inc]);
  }
  if (loaded[inc][ADJ]) {
    free(adjptr[inc]);
  }
  if (loaded[inc][VORO]) {
    free(vnodeptr[inc]);
    free(vedgeptr[inc]);
    free(vnormptr[inc]);
  }
}

void move_inc(inc)
int inc;
{
  int i;

  free_inc(1 - inc);
  for (i = 0; i < IMAGE_TYPES; i++) {
    loaded[1 - inc][i] = loaded[inc][i];
    loaded[inc][i] = 0;
    xlo[1 - inc][i] = xlo[inc][i];
    ylo[1 - inc][i] = ylo[inc][i];
    xhi[1 - inc][i] = xhi[inc][i];
    yhi[1 - inc][i] = yhi[inc][i];
  }
  nodes[1 - inc] = nodes[inc];
  node_dim[1 - inc] = node_dim[inc];
  nodeptr[1 - inc] = nodeptr[inc];
  polynodes[1 - inc] = polynodes[inc];
  poly_dim[1 - inc] = poly_dim[inc];
  polyedges[1 - inc] = polyedges[inc];
  polyholes[1 - inc] = polyholes[inc];
  polynodeptr[1 - inc] = polynodeptr[inc];
  polyedgeptr[1 - inc] = polyedgeptr[inc];
  polyholeptr[1 - inc] = polyholeptr[inc];
  elems[1 - inc] = elems[inc];
  ele_corners[1 - inc] = ele_corners[inc];
  eleptr[1 - inc] = eleptr[inc];
  edges[1 - inc] = edges[inc];
  edgeptr[1 - inc] = edgeptr[inc];
  normptr[1 - inc] = normptr[inc];
  subdomains[1 - inc] = subdomains[inc];
  partpart[1 - inc] = partpart[inc];
  partcenter[1 - inc] = partcenter[inc];
  partshift[1 - inc] = partshift[inc];
  adjsubdomains[1 - inc] = adjsubdomains[inc];
  adjptr[1 - inc] = adjptr[inc];
  vnodes[1 - inc] = vnodes[inc];
  vnode_dim[1 - inc] = vnode_dim[inc];
  vnodeptr[1 - inc] = vnodeptr[inc];
  vedges[1 - inc] = vedges[inc];
  vedgeptr[1 - inc] = vedgeptr[inc];
  vnormptr[1 - inc] = vnormptr[inc];
  firstnumber[1 - inc] = firstnumber[inc];
  firstnumber[inc] = -1;
}

void unload_inc(inc)
int inc;
{
  int i;

  current_image = NOTHING;
  for (i = 0; i < IMAGE_TYPES; i++) {
    loaded[inc][i] = 0;
    firstnumber[inc] = -1;
  }
}

void showme_init()
{
  current_image = NOTHING;
  current_inc = 0;
  explosion = STARTEXPLOSION;
  unload_inc(0);
  unload_inc(1);
}

char *readline(string, infile, infilename)
char *string;
FILE *infile;
char *infilename;
{
  char *result;

  do {
    result = fgets(string, INPUTLINESIZE, infile);
    if (result == (char *) NULL) {
      printf("  Error:  Unexpected end of file in %s.\n",
             infilename);
      exit(1);
    }
    while ((*result != '\0') && (*result != '#')
           && (*result != '.') && (*result != '+') && (*result != '-')
           && ((*result < '0') || (*result > '9'))) {
      result++;
    }
  } while ((*result == '#') || (*result == '\0'));
  return result;
}

char *findfield(string)
char *string;
{
  char *result;

  result = string;
  while ((*result != '\0') && (*result != '#')
         && (*result != ' ') && (*result != '\t')) {
    result++;
  }
  while ((*result != '\0') && (*result != '#')
         && (*result != '.') && (*result != '+') && (*result != '-')
         && ((*result < '0') || (*result > '9'))) {
    result++;
  }
  if (*result == '#') {
    *result = '\0';
  }
  return result;
}

int load_node(fname, firstnumber, nodes, dim, ptr, xmin, ymin, xmax, ymax)
char *fname;
int *firstnumber;
int *nodes;
int *dim;
REAL **ptr;
REAL *xmin;
REAL *ymin;
REAL *xmax;
REAL *ymax;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int extras;
  int nodemarks;
  int index;
  int nodenumber;
  int i, j;
  int smallerr;
  REAL x, y;

  *xmin = *ymin = 0.0;
  *xmax = *ymax = 1.0;
  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", fname);
    return 1;
  }
  stringptr = readline(inputline, infile, fname);
  *nodes = (int) strtol (stringptr, &stringptr, 0);
  if (*nodes < 3) {
    printf("  Error:  %s contains %d points.\n", fname, *nodes);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    *dim = 2;
  } else {
    *dim = (int) strtol (stringptr, &stringptr, 0);
  }
  if (*dim < 1) {
    printf("  Error:  %s has dimensionality %d.\n", fname, *dim);
    return 1;
  }
  if (*dim != 2) {
    printf("  I only understand two-dimensional meshes.\n");
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    extras = 0;
  } else {
    extras = (int) strtol (stringptr, &stringptr, 0);
  }
  if (extras < 0) {
    printf("  Error:  %s has negative value for number of attributes.\n",
           fname);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    nodemarks = 0;
  } else {
    nodemarks = (int) strtol (stringptr, &stringptr, 0);
  }
  if (nodemarks < 0) {
    printf("  Warning:  %s has negative value for number of point markers.\n",
           fname);
  }
  if (nodemarks > 1) {
    printf(
   "  Warning:  %s has value greater than one for number of point markers.\n",
           fname);
  }
  *ptr = (REAL *) malloc((*nodes + 1) * *dim * sizeof(REAL));
  if (*ptr == (REAL *) NULL) {
    printf("  Out of memory.\n");
    return 1;
  }
  index = *dim;
  smallerr = 1;
  for (i = 0; i < *nodes; i++) {
    stringptr = readline(inputline, infile, fname);
    nodenumber = (int) strtol (stringptr, &stringptr, 0);
    if ((i == 0) && (*firstnumber == -1)) {
      if (nodenumber == 0) {
        *firstnumber = 0;
      } else {
        *firstnumber = 1;
      }
    }
    if ((nodenumber != *firstnumber + i) && (smallerr)) {
      printf("  Warning:  Points in %s are not numbered correctly\n", fname);
      printf("            (starting with point %d).\n", *firstnumber + i);
      smallerr = 0;
    }
    for (j = 0; j < *dim; j++) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Point %d is missing a coordinate in %s.\n",
               *firstnumber + i, fname);
        free(*ptr);
        return 1;
      }
      (*ptr)[index++] = (REAL) strtod(stringptr, &stringptr);
    }
  }
  fclose(infile);
  index = *dim;
  *xmin = *xmax = (*ptr)[index];
  *ymin = *ymax = (*ptr)[index + 1];
  for (i = 2; i <= *nodes; i++) {
    index += *dim;
    x = (*ptr)[index];
    y = (*ptr)[index + 1];
    if (x < *xmin) {
      *xmin = x;
    }
    if (y < *ymin) {
      *ymin = y;
    }
    if (x > *xmax) {
      *xmax = x;
    }
    if (y > *ymax) {
      *ymax = y;
    }
  }
  if (*xmin == *xmax) {
    *xmin -= 0.5;
    *xmax += 0.5;
  }
  if (*ymin == *ymax) {
    *ymin -= 0.5;
    *ymax += 0.5;
  }
  return 0;
}

int load_poly(inc, fname, firstnumber, pnodes, dim, edges, holes, nodeptr,
              edgeptr, holeptr, xmin, ymin, xmax, ymax)
int inc;
char *fname;
int *firstnumber;
int *pnodes;
int *dim;
int *edges;
int *holes;
REAL **nodeptr;
int **edgeptr;
REAL **holeptr;
REAL *xmin;
REAL *ymin;
REAL *xmax;
REAL *ymax;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int extras;
  int nodemarks;
  int segmentmarks;
  int index;
  int nodenumber, edgenumber, holenumber;
  int maxnode;
  int i, j;
  int smallerr;
  REAL x, y;

  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", fname);
    return 1;
  }
  stringptr = readline(inputline, infile, fname);
  *pnodes = (int) strtol (stringptr, &stringptr, 0);
  if (*pnodes == 0) {
    if (!loaded[inc][NODE]) {
      if (load_image(inc, NODE)) {
        return 1;
      }
    }
    maxnode = nodes[inc];
    *xmin = xlo[inc][NODE];
    *ymin = ylo[inc][NODE];
    *xmax = xhi[inc][NODE];
    *ymax = yhi[inc][NODE];
  } else {
    if (*pnodes < 1) {
      printf("  Error:  %s contains %d points.\n", fname, *pnodes);
      return 1;
    }
    maxnode = *pnodes;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    *dim = 2;
  } else {
    *dim = (int) strtol (stringptr, &stringptr, 0);
  }
  if (*dim < 1) {
    printf("  Error:  %s has dimensionality %d.\n", fname, *dim);
    return 1;
  }
  if (*dim != 2) {
    printf("  I only understand two-dimensional .poly files.\n");
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    extras = 0;
  } else {
    extras = (int) strtol (stringptr, &stringptr, 0);
  }
  if (extras < 0) {
    printf("  Error:  %s has negative value for number of attributes.\n",
           fname);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    nodemarks = 0;
  } else {
    nodemarks = (int) strtol (stringptr, &stringptr, 0);
  }
  if (nodemarks < 0) {
    printf("  Warning:  %s has negative value for number of point markers.\n",
           fname);
  }
  if (nodemarks > 1) {
    printf(
   "  Warning:  %s has value greater than one for number of point markers.\n",
           fname);
  }
  if (*pnodes > 0) {
    *nodeptr = (REAL *) malloc((*pnodes + 1) * *dim * sizeof(REAL));
    if (*nodeptr == (REAL *) NULL) {
      printf("  Out of memory.\n");
      return 1;
    }
    index = *dim;
    smallerr = 1;
    for (i = 0; i < *pnodes; i++) {
      stringptr = readline(inputline, infile, fname);
      nodenumber = (int) strtol (stringptr, &stringptr, 0);
      if ((i == 0) && (*firstnumber == -1)) {
        if (nodenumber == 0) {
          *firstnumber = 0;
        } else {
          *firstnumber = 1;
        }
      }
      if ((nodenumber != *firstnumber + i) && (smallerr)) {
        printf("  Warning:  Points in %s are not numbered correctly.\n",
               fname);
        printf("            (starting with point %d).\n", *firstnumber + i);
        smallerr = 0;
      }
      for (j = 0; j < *dim; j++) {
        stringptr = findfield(stringptr);
        if (*stringptr == '\0') {
          printf("Error:  Point %d is missing a coordinate in %s.\n",
                 *firstnumber + i, fname);
          free(*nodeptr);
          return 1;
        }
        (*nodeptr)[index++] = (REAL) strtod(stringptr, &stringptr);
      }
    }
  }
  stringptr = readline(inputline, infile, fname);
  *edges = (int) strtol (stringptr, &stringptr, 0);
  if (*edges < 0) {
    printf("  Error:  %s contains %d segments.\n", fname, *edges);
    free(*nodeptr);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    segmentmarks = 0;
  } else {
    segmentmarks = (int) strtol (stringptr, &stringptr, 0);
  }
  if (segmentmarks < 0) {
    printf("  Error:  %s has negative value for number of segment markers.\n",
           fname);
    free(*nodeptr);
    return 1;
  }
  if (segmentmarks > 1) {
    printf(
   "  Error:  %s has value greater than one for number of segment markers.\n",
           fname);
    free(*nodeptr);
    return 1;
  }
  *edgeptr = (int *) malloc(((*edges + 1) << 1) * sizeof(int));
  if (*edgeptr == (int *) NULL) {
    printf("  Out of memory.\n");
    free(*nodeptr);
    return 1;
  }
  index = 2;
  smallerr = 1;
  for (i = *firstnumber; i < *firstnumber + *edges; i++) {
    stringptr = readline(inputline, infile, fname);
    edgenumber = (int) strtol (stringptr, &stringptr, 0);
    if ((edgenumber != i) && (smallerr)) {
      printf("  Warning:  Segments in %s are not numbered correctly.\n",
             fname);
      printf("            (starting with segment %d).\n", i);
      smallerr = 0;
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Segment %d is missing its endpoints in %s.\n", i, fname);
      free(*nodeptr);
      free(*edgeptr);
      return 1;
    }
    (*edgeptr)[index] = (int) strtol (stringptr, &stringptr, 0) + 1 -
                        *firstnumber;
    if (((*edgeptr)[index] < 1) || ((*edgeptr)[index] > maxnode)) {
      printf("Error:  Segment %d has invalid endpoint in %s.\n", i, fname);
      return 1;
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Segment %d is missing an endpoint in %s.\n", i, fname);
      free(*nodeptr);
      free(*edgeptr);
      return 1;
    }
    (*edgeptr)[index + 1] = (int) strtol (stringptr, &stringptr, 0) + 1 -
                            *firstnumber;
    if (((*edgeptr)[index + 1] < 1) || ((*edgeptr)[index + 1] > maxnode)) {
      printf("Error:  Segment %d has invalid endpoint in %s.\n", i, fname);
      return 1;
    }
    index += 2;
  }
  stringptr = readline(inputline, infile, fname);
  *holes = (int) strtol (stringptr, &stringptr, 0);
  if (*holes < 0) {
    printf("  Error:  %s contains %d holes.\n", fname, *holes);
    free(*nodeptr);
    free(*edgeptr);
    return 1;
  }
  *holeptr = (REAL *) malloc((*holes + 1) * *dim * sizeof(REAL));
  if (*holeptr == (REAL *) NULL) {
    printf("  Out of memory.\n");
    free(*nodeptr);
    free(*edgeptr);
    return 1;
  }
  index = *dim;
  smallerr = 1;
  for (i = *firstnumber; i < *firstnumber + *holes; i++) {
    stringptr = readline(inputline, infile, fname);
    holenumber = (int) strtol (stringptr, &stringptr, 0);
    if ((holenumber != i) && (smallerr)) {
      printf("  Warning:  Holes in %s are not numbered correctly.\n", fname);
      printf("            (starting with hole %d).\n", i);
      smallerr = 0;
    }
    for (j = 0; j < *dim; j++) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Hole %d is missing a coordinate in %s.\n", i,
               fname);
        free(*nodeptr);
        free(*edgeptr);
        free(*holeptr);
        return 1;
      }
      (*holeptr)[index++] = (REAL) strtod(stringptr, &stringptr);
    }
  }
  fclose(infile);
  if (*pnodes > 0) {
    index = *dim;
    *xmin = *xmax = (*nodeptr)[index];
    *ymin = *ymax = (*nodeptr)[index + 1];
    for (i = 2; i <= *pnodes; i++) {
      index += *dim;
      x = (*nodeptr)[index];
      y = (*nodeptr)[index + 1];
      if (x < *xmin) {
        *xmin = x;
      }
      if (y < *ymin) {
        *ymin = y;
      }
      if (x > *xmax) {
        *xmax = x;
      }
      if (y > *ymax) {
        *ymax = y;
      }
    }
  }
  index = *dim;
  for (i = 1; i <= *holes; i++) {
    x = (*holeptr)[index];
    y = (*holeptr)[index + 1];
    if (x < *xmin) {
      *xmin = x;
    }
    if (y < *ymin) {
      *ymin = y;
    }
    if (x > *xmax) {
      *xmax = x;
    }
    if (y > *ymax) {
      *ymax = y;
    }
    index += *dim;
  }
  return 0;
}

int load_ele(fname, firstnumber, nodes, elems, corners, ptr)
char *fname;
int firstnumber;
int nodes;
int *elems;
int *corners;
int **ptr;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int extras;
  int index;
  int elemnumber;
  int i, j;
  int smallerr;

  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", fname);
    return 1;
  }
  stringptr = readline(inputline, infile, fname);
  *elems = (int) strtol (stringptr, &stringptr, 0);
  if (*elems < 1) {
    printf("  Error:  %s contains %d triangles.\n", fname, *elems);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    *corners = 3;
  } else {
    *corners = (int) strtol (stringptr, &stringptr, 0);
  }
  if (*corners < 3) {
    printf("  Error:  Triangles in %s have only %d corners.\n", fname,
           *corners);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    extras = 0;
  } else {
    extras = (int) strtol (stringptr, &stringptr, 0);
  }
  if (extras < 0) {
    printf("  Error:  %s has negative value for extra fields.\n", fname);
    return 1;
  }
  *ptr = (int *) malloc((*elems + 1) * 3 * sizeof(int));
  if (*ptr == (int *) NULL) {
    printf("  Out of memory.\n");
    return 1;
  }
  index = 3;
  smallerr = 1;
  for (i = firstnumber; i < firstnumber + *elems; i++) {
    stringptr = readline(inputline, infile, fname);
    elemnumber = (int) strtol (stringptr, &stringptr, 0);
    if ((elemnumber != i) && (smallerr)) {
      printf("  Warning:  Triangles in %s are not numbered correctly.\n",
             fname);
      printf("            (starting with triangle %d).\n", i);
      smallerr = 0;
    }
    for (j = 0; j < 3; j++) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Triangle %d is missing a corner in %s.\n", i, fname);
        free(*ptr);
        return 1;
      }
      (*ptr)[index] = (int) strtol (stringptr, &stringptr, 0) + 1 -
                      firstnumber;
      if (((*ptr)[index] < 1) || ((*ptr)[index] > nodes)) {
        printf("Error:  Triangle %d has invalid corner in %s.\n", i, fname);
        return 1;
      }
      index++;
    }
  }
  fclose(infile);
  return 0;
}

int load_edge(fname, firstnumber, nodes, edges, edgeptr, normptr)
char *fname;
int firstnumber;
int nodes;
int *edges;
int **edgeptr;
REAL **normptr;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int index;
  int edgenumber;
  int edgemarks;
  int i;
  int smallerr;

  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
    if (infile == (FILE *) NULL) {
      printf("  Error:  Cannot access file %s.\n", fname);
      return 1;
    }
  stringptr = readline(inputline, infile, fname);
  *edges = (int) strtol (stringptr, &stringptr, 0);
  if (*edges < 1) {
    printf("  Error:  %s contains %d edges.\n", fname, *edges);
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    edgemarks = 0;
  } else {
    edgemarks = (int) strtol (stringptr, &stringptr, 0);
  }
  if (edgemarks < 0) {
    printf("  Error:  %s has negative value for number of edge markers.\n",
           fname);
    return 1;
  }
  if (edgemarks > 1) {
    printf(
   "  Error:  %s has value greater than one for number of edge markers.\n",
           fname);
    return 1;
  }
  *edgeptr = (int *) malloc(((*edges + 1) << 1) * sizeof(int));
  if (*edgeptr == (int *) NULL) {
    printf("  Out of memory.\n");
    return 1;
  }
  *normptr = (REAL *) malloc(((*edges + 1) << 1) * sizeof(REAL));
  if (*normptr == (REAL *) NULL) {
    printf("  Out of memory.\n");
    free(*edgeptr);
    return 1;
  }
  index = 2;
  smallerr = 1;
  for (i = firstnumber; i < firstnumber + *edges; i++) {
    stringptr = readline(inputline, infile, fname);
    edgenumber = (int) strtol (stringptr, &stringptr, 0);
    if ((edgenumber != i) && (smallerr)) {
      printf("  Warning:  Edges in %s are not numbered correctly.\n", fname);
      printf("            (starting with edge %d).\n", i);
      smallerr = 0;
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Edge %d is missing its endpoints in %s.\n", i, fname);
      free(*edgeptr);
      free(*normptr);
      return 1;
    }
    (*edgeptr)[index] = (int) strtol (stringptr, &stringptr, 0) + 1 -
                        firstnumber;
    if (((*edgeptr)[index] < 1) || ((*edgeptr)[index] > nodes)) {
      printf("Error:  Edge %d has invalid endpoint in %s.\n", i, fname);
      return 1;
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Edge %d is missing an endpoint in %s.\n", i, fname);
      free(*edgeptr);
      free(*normptr);
      return 1;
    }
    (*edgeptr)[index + 1] = (int) strtol (stringptr, &stringptr, 0);
    if ((*edgeptr)[index + 1] == -1) {
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Edge %d is missing its direction in %s.\n", i, fname);
        free(*edgeptr);
        free(*normptr);
        return 1;
      }
      (*normptr)[index] = (REAL) strtod(stringptr, &stringptr);
      stringptr = findfield(stringptr);
      if (*stringptr == '\0') {
        printf("Error:  Edge %d is missing a direction coordinate in %s.\n",
               i, fname);
        free(*edgeptr);
        free(*normptr);
        return 1;
      }
      (*normptr)[index + 1] = (REAL) strtod(stringptr, &stringptr);
    } else {
      (*edgeptr)[index + 1] += 1 - firstnumber;
      if (((*edgeptr)[index + 1] < 1) || ((*edgeptr)[index + 1] > nodes)) {
        printf("Error:  Edge %d has invalid endpoint in %s.\n", i, fname);
        return 1;
      }
    }
    index += 2;
  }
  fclose(infile);
  return 0;
}

int load_part(fname, dim, firstnumber, elems, nodeptr, eleptr, parts,
              partition, partcenter, partshift)
char *fname;
int dim;
int firstnumber;
int elems;
REAL *nodeptr;
int *eleptr;
int *parts;
int **partition;
REAL **partcenter;
REAL **partshift;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int partelems;
  int index;
  int elemnumber;
  int i, j;
  int smallerr;
  int *partsize;

  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", fname);
    return 1;
  }
  stringptr = readline(inputline, infile, fname);
  partelems = (int) strtol (stringptr, &stringptr, 0);
  if (partelems != elems) {
    printf(
      "  Error:  .ele and .part files do not agree on number of triangles.\n");
    return 1;
  }
  stringptr = findfield(stringptr);
  if (*stringptr == '\0') {
    *parts = 1;
  } else {
    *parts = (int) strtol (stringptr, &stringptr, 0);
  }
  if (*parts < 1) {
    printf("  Error:  %s specifies %d subdomains.\n", fname, *parts);
    return 1;
  }
  *partition = (int *) malloc((elems + 1) * sizeof(int));
  if (*partition == (int *) NULL) {
    printf("  Out of memory.\n");
    return 1;
  }
  smallerr = 1;
  for (i = firstnumber; i < firstnumber + partelems; i++) {
    stringptr = readline(inputline, infile, fname);
    elemnumber = (int) strtol (stringptr, &stringptr, 0);
    if ((elemnumber != i) && (smallerr)) {
      printf("  Warning:  Triangles in %s are not numbered correctly.\n",
             fname);
      printf("            (starting with triangle %d).\n", i);
      smallerr = 0;
    }
    stringptr = findfield(stringptr);
    if (*stringptr == '\0') {
      printf("Error:  Triangle %d has no subdomain in %s.\n", i, fname);
      free(*partition);
      return 1;
    }
    j = i + 1 - firstnumber;
    (*partition)[j] = (int) strtol (stringptr, &stringptr, 0) - firstnumber;
    if (((*partition)[j] >= *parts) || ((*partition)[j] < 0)) {
      printf("  Error:  Triangle %d of %s has an invalid subdomain.\n",
             i, fname);
      free(*partition);
      return 1;
    }
  }
  fclose(infile);
  *partcenter = (REAL *) malloc(((*parts + 1) << 1) * sizeof(REAL));
  if (*partcenter == (REAL *) NULL) {
    printf("Error:  Out of memory.\n");
    free(*partition);
    return 1;
  }
  *partshift = (REAL *) malloc((*parts << 1) * sizeof(REAL));
  if (*partshift == (REAL *) NULL) {
    printf("Error:  Out of memory.\n");
    free(*partition);
    free(*partcenter);
    return 1;
  }
  partsize = (int *) malloc((*parts + 1) * sizeof(int));
  if (partsize == (int *) NULL) {
    printf("Error:  Out of memory.\n");
    free(*partition);
    free(*partcenter);
    free(*partshift);
    return 1;
  }
  for (i = 0; i <= *parts; i++) {
    partsize[i] = 0;
    (*partcenter)[i << 1] = 0.0;
    (*partcenter)[(i << 1) + 1] = 0.0;
  }
  index = 3;
  for (i = 1; i <= elems; i++) {
    partsize[(*partition)[i]] += 3;
    for (j = 0; j < 3; j++) {
      (*partcenter)[(*partition)[i] << 1] +=
                nodeptr[eleptr[index] * dim];
      (*partcenter)[((*partition)[i] << 1) + 1] +=
                nodeptr[eleptr[index++] * dim + 1];
    }
  }
  for (i = 0; i < *parts; i++) {
    (*partcenter)[i << 1] /= (REAL) partsize[i];
    (*partcenter)[(i << 1) + 1] /= (REAL) partsize[i];
    (*partcenter)[*parts << 1] += (*partcenter)[i << 1];
    (*partcenter)[(*parts << 1) + 1] += (*partcenter)[(i << 1) + 1];
  }
  (*partcenter)[*parts << 1] /= (REAL) *parts;
  (*partcenter)[(*parts << 1) + 1] /= (REAL) *parts;
  free(partsize);
  return 0;
}

int load_adj(fname, subdomains, ptr)
char *fname;
int *subdomains;
int **ptr;
{
  FILE *infile;
  char inputline[INPUTLINESIZE];
  char *stringptr;
  int i, j;

  if (!quiet) {
    printf("Opening %s.\n", fname);
  }
  infile = fopen(fname, "r");
  if (infile == (FILE *) NULL) {
    printf("  Error:  Cannot access file %s.\n", fname);
    return 1;
  }
  stringptr = readline(inputline, infile, fname);
  *subdomains = (int) strtol (stringptr, &stringptr, 0);
  if (*subdomains < 1) {
    printf("  Error:  %s contains %d subdomains.\n", fname, *subdomains);
    return 1;
  }
  *ptr = (int *) malloc(*subdomains * *subdomains * sizeof(int));
  if (*ptr == (int *) NULL) {
    printf("  Out of memory.\n");
    return 1;
  }
  for (i = 0; i < *subdomains; i++) {
    for (j = 0; j < *subdomains; j++) {
      stringptr = readline(inputline, infile, fname);
      (*ptr)[i * *subdomains + j] = (int) strtol (stringptr, &stringptr, 0);
    }
  }
  return 0;
}

void findpartshift(parts, explosion, partcenter, partshift)
int parts;
REAL explosion;
REAL *partcenter;
REAL *partshift;
{
  int i;

  for (i = 0; i < parts; i++) {
    partshift[i << 1] = explosion *
          (partcenter[i << 1] - partcenter[parts << 1]);
    partshift[(i << 1) + 1] = explosion *
          (partcenter[(i << 1) + 1] - partcenter[(parts << 1) + 1]);
  }
}

int load_image(inc, image)
int inc;
int image;
{
  int error;

  switch (image) {
    case NODE:
      error = load_node(nodefilename[inc], &firstnumber[inc], &nodes[inc],
                        &node_dim[inc], &nodeptr[inc], &xlo[inc][NODE],
                        &ylo[inc][NODE], &xhi[inc][NODE], &yhi[inc][NODE]);
      break;
    case POLY:
      error = load_poly(inc, polyfilename[inc], &firstnumber[inc],
                        &polynodes[inc], &poly_dim[inc], &polyedges[inc],
                        &polyholes[inc], &polynodeptr[inc], &polyedgeptr[inc],
                        &polyholeptr[inc],
                        &xlo[inc][POLY], &ylo[inc][POLY],
                        &xhi[inc][POLY], &yhi[inc][POLY]);
      break;
    case ELE:
      error = load_ele(elefilename[inc], firstnumber[inc], nodes[inc],
                       &elems[inc], &ele_corners[inc], &eleptr[inc]);
      xlo[inc][ELE] = xlo[inc][NODE];
      ylo[inc][ELE] = ylo[inc][NODE];
      xhi[inc][ELE] = xhi[inc][NODE];
      yhi[inc][ELE] = yhi[inc][NODE];
      break;
    case EDGE:
      error = load_edge(edgefilename[inc], firstnumber[inc], nodes[inc],
                        &edges[inc], &edgeptr[inc], &normptr[inc]);
      xlo[inc][EDGE] = xlo[inc][NODE];
      ylo[inc][EDGE] = ylo[inc][NODE];
      xhi[inc][EDGE] = xhi[inc][NODE];
      yhi[inc][EDGE] = yhi[inc][NODE];
      break;
    case PART:
      error = load_part(partfilename[inc], node_dim[inc], firstnumber[inc],
                        elems[inc], nodeptr[inc], eleptr[inc],
                        &subdomains[inc], &partpart[inc], &partcenter[inc],
                        &partshift[inc]);
      if (!error) {
        findpartshift(subdomains[inc], explosion, partcenter[inc],
                      partshift[inc]);
      }
      xlo[inc][PART] = xlo[inc][NODE];
      ylo[inc][PART] = ylo[inc][NODE];
      xhi[inc][PART] = xhi[inc][NODE];
      yhi[inc][PART] = yhi[inc][NODE];
      break;
    case ADJ:
      error = load_adj(adjfilename[inc], &adjsubdomains[inc], &adjptr[inc]);
      xlo[inc][ADJ] = xlo[inc][NODE];
      ylo[inc][ADJ] = ylo[inc][NODE];
      xhi[inc][ADJ] = xhi[inc][NODE];
      yhi[inc][ADJ] = yhi[inc][NODE];
      break;
    case VORO:
      error = load_node(vnodefilename[inc], &firstnumber[inc], &vnodes[inc],
                        &vnode_dim[inc], &vnodeptr[inc], &xlo[inc][VORO],
                        &ylo[inc][VORO], &xhi[inc][VORO], &yhi[inc][VORO]);
      if (!error) {
        error = load_edge(vedgefilename[inc], firstnumber[inc], vnodes[inc],
                          &vedges[inc], &vedgeptr[inc], &vnormptr[inc]);
      }
      break;
    default:
      error = 1;
  }
  if (!error) {
    loaded[inc][image] = 1;
  }
  return error;
}

void choose_image(inc, image)
int inc;
int image;
{
  if (!loaded[inc][image]) {
    if ((image == ELE) || (image == EDGE) || (image == PART)
        || (image == ADJ)) {
      if (!loaded[inc][NODE]) {
        if (load_image(inc, NODE)) {
          return;
        }
      }
    }
    if ((image == PART) || (image == ADJ)) {
      if (!loaded[inc][ELE]) {
        if (load_image(inc, ELE)) {
          return;
        }
      }
    }
    if (image == ADJ) {
      if (!loaded[inc][PART]) {
        if (load_image(inc, PART)) {
          return;
        }
      }
    }
    if (load_image(inc, image)) {
      return;
    }
  }
  current_inc = inc;
  current_image = image;
}

Window make_button(name, x, y, width)
char *name;
int x;
int y;
int width;
{
  XSetWindowAttributes attr;
  XSizeHints hints;
  Window button;

  attr.background_pixel = black;
  attr.border_pixel = white;
  attr.backing_store = NotUseful;
  attr.event_mask = ExposureMask | ButtonReleaseMask | ButtonPressMask;
  attr.bit_gravity = SouthWestGravity;
  attr.win_gravity = SouthWestGravity;
  attr.save_under = False;
  button = XCreateWindow(display, mainwindow, x, y, width, BUTTONHEIGHT - 4,
                         2, 0, InputOutput, CopyFromParent,
                         CWBackPixel | CWBorderPixel | CWEventMask |
                         CWBitGravity | CWWinGravity | CWBackingStore |
                         CWSaveUnder, &attr);
  hints.width = width;
  hints.height = BUTTONHEIGHT - 4;
  hints.min_width = 0;
  hints.min_height = BUTTONHEIGHT - 4;
  hints.max_width = width;
  hints.max_height = BUTTONHEIGHT - 4;
  hints.width_inc = 1;
  hints.height_inc = 1;
  hints.flags = PMinSize | PMaxSize | PSize | PResizeInc;
  XSetStandardProperties(display, button, name, "showme", None, (char **) NULL,
                         0, &hints);
  return button;
}

void make_buttons(y)
int y;
{
  int i;

  for (i = 1; i >= 0; i--) {
    nodewin[i] = make_button("node", 0, y + (1 - i) * BUTTONHEIGHT, 42);
    XMapWindow(display, nodewin[i]);
    polywin[i] = make_button("poly", 44, y + (1 - i) * BUTTONHEIGHT, 42);
    XMapWindow(display, polywin[i]);
    elewin[i] = make_button("ele", 88, y + (1 - i) * BUTTONHEIGHT, 33);
    XMapWindow(display, elewin[i]);
    edgewin[i] = make_button("edge", 123, y + (1 - i) * BUTTONHEIGHT, 42);
    XMapWindow(display, edgewin[i]);
    partwin[i] = make_button("part", 167, y + (1 - i) * BUTTONHEIGHT, 42);
    XMapWindow(display, partwin[i]);
    adjwin[i] = make_button("adj", 211, y + (1 - i) * BUTTONHEIGHT, 33);
    XMapWindow(display, adjwin[i]);
    voronoiwin[i] = make_button("voro", 246, y + (1 - i) * BUTTONHEIGHT, 42);
    XMapWindow(display, voronoiwin[i]);
  }
  versionpluswin = make_button("    +", 290, y, 52);
  XMapWindow(display, versionpluswin);
  versionminuswin = make_button("    -", 290, y + BUTTONHEIGHT, 52);
  XMapWindow(display, versionminuswin);

  quitwin = make_button("Quit", 0, y + 2 * BUTTONHEIGHT, 42);
  XMapWindow(display, quitwin);
  leftwin = make_button("<", 44, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, leftwin);
  rightwin = make_button(">", 60, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, rightwin);
  upwin = make_button("^", 76, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, upwin);
  downwin = make_button("v", 92, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, downwin);
  resetwin = make_button("Reset", 108, y + 2 * BUTTONHEIGHT, 52);
  XMapWindow(display, resetwin);
  widthpluswin = make_button("Width+", 162, y + 2 * BUTTONHEIGHT, 61);
  XMapWindow(display, widthpluswin);
  widthminuswin = make_button("-", 225, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, widthminuswin);
  expwin = make_button("Exp", 241, y + 2 * BUTTONHEIGHT, 33);
  XMapWindow(display, expwin);
  exppluswin = make_button("+", 276, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, exppluswin);
  expminuswin = make_button("-", 292, y + 2 * BUTTONHEIGHT, 14);
  XMapWindow(display, expminuswin);
  fillwin = make_button("Fill", 308, y + 2 * BUTTONHEIGHT, 41);
  XMapWindow(display, fillwin);
  pswin = make_button("PS", 351, y + 2 * BUTTONHEIGHT, 24);
  XMapWindow(display, pswin);
  epswin = make_button("EPS", 377, y + 2 * BUTTONHEIGHT, 33);
  XMapWindow(display, epswin);
}

void fill_button(button)
Window button;
{
  int x, y;
  unsigned int w, h, d, b;
  Window rootw;

  XGetGeometry(display, button, &rootw, &x, &y, &w, &h, &d, &b);
  XFillRectangle(display, button, fontgc, 0, 0, w, h);
}

void draw_buttons()
{
  char numberstring[32];
  char buttonstring[6];
  int i;

  for (i = 1; i >= 0; i--) {
    if ((current_image == NODE) && (current_inc == i)) {
      fill_button(nodewin[i]);
      XDrawString(display, nodewin[i], blackfontgc, 2, 13, "node", 4);
    } else {
      XClearWindow(display, nodewin[i]);
      XDrawString(display, nodewin[i], fontgc, 2, 13, "node", 4);
    }
    if ((current_image == POLY) && (current_inc == i)) {
      fill_button(polywin[i]);
      XDrawString(display, polywin[i], blackfontgc, 2, 13, "poly", 4);
    } else {
      XClearWindow(display, polywin[i]);
      XDrawString(display, polywin[i], fontgc, 2, 13, "poly", 4);
    }
    if ((current_image == ELE) && (current_inc == i)) {
      fill_button(elewin[i]);
      XDrawString(display, elewin[i], blackfontgc, 2, 13, "ele", 3);
    } else {
      XClearWindow(display, elewin[i]);
      XDrawString(display, elewin[i], fontgc, 2, 13, "ele", 3);
    }
    if ((current_image == EDGE) && (current_inc == i)) {
      fill_button(edgewin[i]);
      XDrawString(display, edgewin[i], blackfontgc, 2, 13, "edge", 4);
    } else {
      XClearWindow(display, edgewin[i]);
      XDrawString(display, edgewin[i], fontgc, 2, 13, "edge", 4);
    }
    if ((current_image == PART) && (current_inc == i)) {
      fill_button(partwin[i]);
      XDrawString(display, partwin[i], blackfontgc, 2, 13, "part", 4);
    } else {
      XClearWindow(display, partwin[i]);
      XDrawString(display, partwin[i], fontgc, 2, 13, "part", 4);
    }
    if ((current_image == ADJ) && (current_inc == i)) {
      fill_button(adjwin[i]);
      XDrawString(display, adjwin[i], blackfontgc, 2, 13, "adj", 3);
    } else {
      XClearWindow(display, adjwin[i]);
      XDrawString(display, adjwin[i], fontgc, 2, 13, "adj", 3);
    }
    if ((current_image == VORO) && (current_inc == i)) {
      fill_button(voronoiwin[i]);
      XDrawString(display, voronoiwin[i], blackfontgc, 2, 13, "voro", 4);
    } else {
      XClearWindow(display, voronoiwin[i]);
      XDrawString(display, voronoiwin[i], fontgc, 2, 13, "voro", 4);
    }
  }

  XClearWindow(display, versionpluswin);
  sprintf(numberstring, "%d", loweriteration + 1);
  sprintf(buttonstring, "%-4.4s+", numberstring);
  XDrawString(display, versionpluswin, fontgc, 2, 13, buttonstring, 5);
  XClearWindow(display, versionminuswin);
  sprintf(numberstring, "%d", loweriteration);
  if (loweriteration == 0) {
    sprintf(buttonstring, "%-4.4s", numberstring);
  } else {
    sprintf(buttonstring, "%-4.4s-", numberstring);
  }
  XDrawString(display, versionminuswin, fontgc, 2, 13, buttonstring, 5);

  XClearWindow(display, quitwin);
  XDrawString(display, quitwin, fontgc, 2, 13, "Quit", 4);
  XClearWindow(display, leftwin);
  XDrawString(display, leftwin, fontgc, 2, 13, "<", 1);
  XClearWindow(display, rightwin);
  XDrawString(display, rightwin, fontgc, 2, 13, ">", 1);
  XClearWindow(display, upwin);
  XDrawString(display, upwin, fontgc, 2, 13, "^", 1);
  XClearWindow(display, downwin);
  XDrawString(display, downwin, fontgc, 2, 13, "v", 1);
  XClearWindow(display, resetwin);
  XDrawString(display, resetwin, fontgc, 2, 13, "Reset", 6);
  XClearWindow(display, widthpluswin);
  if (line_width < 100) {
    XDrawString(display, widthpluswin, fontgc, 2, 13, "Width+", 6);
  } else {
    XDrawString(display, widthpluswin, fontgc, 2, 13, "Width ", 6);
  }
  XClearWindow(display, widthminuswin);
  if (line_width > 1) {
    XDrawString(display, widthminuswin, fontgc, 2, 13, "-", 1);
  }
  XClearWindow(display, expwin);
  XClearWindow(display, exppluswin);
  XClearWindow(display, expminuswin);
  XClearWindow(display, fillwin);
  if (current_image == PART) {
    if (explode) {
      fill_button(expwin);
      XDrawString(display, expwin, blackfontgc, 2, 13, "Exp", 3);
    } else {
      XDrawString(display, expwin, fontgc, 2, 13, "Exp", 3);
    }
    XDrawString(display, exppluswin, fontgc, 2, 13, "+", 1);
    XDrawString(display, expminuswin, fontgc, 2, 13, "-", 1);
    if (fillelem) {
      fill_button(fillwin);
      XDrawString(display, fillwin, blackfontgc, 2, 13, "Fill", 4);
    } else {
      XDrawString(display, fillwin, fontgc, 2, 13, "Fill", 4);
    }
  }
  XClearWindow(display, pswin);
  XDrawString(display, pswin, fontgc, 2, 13, "PS", 2);
  XClearWindow(display, epswin);
  XDrawString(display, epswin, fontgc, 2, 13, "EPS", 3);
}

void showme_window(argc, argv)
int argc;
char **argv;
{
  XSetWindowAttributes attr;
  XSizeHints hints;
  XGCValues fontvalues, linevalues;
  XColor alloc_color, exact_color;
  int i;

  display = XOpenDisplay((char *) NULL);
  if (!display) {
    printf("Error:  Cannot open display.\n");
    exit(1);
  }
  screen = DefaultScreen(display);
  rootwindow = DefaultRootWindow(display);
  black = BlackPixel(display, screen);
  white = WhitePixel(display, screen);
  windowdepth = DefaultDepth(display, screen);
  rootmap = DefaultColormap(display, screen);
  width = STARTWIDTH;
  height = STARTHEIGHT;
  attr.background_pixel = black;
  attr.border_pixel = white;
  attr.backing_store = NotUseful;
  attr.event_mask = ExposureMask | ButtonReleaseMask | ButtonPressMask |
                    StructureNotifyMask;
  attr.bit_gravity = NorthWestGravity;
  attr.win_gravity = NorthWestGravity;
  attr.save_under = False;
  mainwindow = XCreateWindow(display, rootwindow, 0, 0, width,
                             height + PANELHEIGHT, 3, 0,
                             InputOutput, CopyFromParent,
                             CWBackPixel | CWBorderPixel | CWEventMask |
                             CWBitGravity | CWWinGravity | CWBackingStore |
                             CWSaveUnder, &attr);
  hints.width = width;
  hints.height = height + PANELHEIGHT;
  hints.min_width = MINWIDTH;
  hints.min_height = MINHEIGHT + PANELHEIGHT;
  hints.width_inc = 1;
  hints.height_inc = 1;
  hints.flags = PMinSize | PSize | PResizeInc;
  XSetStandardProperties(display, mainwindow, "Show Me", "showme", None,
                         argv, argc, &hints);
  XChangeProperty(display, mainwindow, XA_WM_CLASS, XA_STRING, 8,
                  PropModeReplace, "showme\0Archimedes", 18);
  XClearWindow(display, mainwindow);
  XMapWindow(display, mainwindow);
  if ((windowdepth > 1) &&
      XAllocNamedColor(display, rootmap, "yellow", &alloc_color,
                       &exact_color)) {
    color = 1;
    explode = bw_ps;
    fontvalues.foreground = alloc_color.pixel;
    linevalues.foreground = alloc_color.pixel;
    showme_foreground = alloc_color.pixel;
    for (i = 0; i < 64; i++) {
      if (XAllocNamedColor(display, rootmap, colorname[i], &alloc_color,
                           &rgb[i])) {
        colors[i] = alloc_color.pixel;
      } else {
        colors[i] = white;
        rgb[i].red = alloc_color.red;
        rgb[i].green = alloc_color.green;
        rgb[i].blue = alloc_color.blue;
        if (!quiet) {
          printf("Warning:  I could not allocate %s.\n", colorname[i]);
        }
      }
    }
  } else {
    color = 0;
    fillelem = 0;
    explode = 1;
    fontvalues.foreground = white;
    linevalues.foreground = white;
    showme_foreground = white;
  }
  font = XLoadQueryFont(display, "9x15");
  fontvalues.background = black;
  fontvalues.font = font->fid;
  fontvalues.fill_style = FillSolid;
  fontvalues.line_width = 2;
  fontgc = XCreateGC(display, rootwindow, GCForeground | GCBackground |
                      GCFont | GCLineWidth | GCFillStyle, &fontvalues);
  fontvalues.foreground = black;
  blackfontgc = XCreateGC(display, rootwindow, GCForeground | GCBackground |
                         GCFont | GCLineWidth | GCFillStyle, &fontvalues);
  linevalues.background = black;
  linevalues.line_width = line_width;
  linevalues.cap_style = CapRound;
  linevalues.join_style = JoinRound;
  linevalues.fill_style = FillSolid;
  linegc = XCreateGC(display, rootwindow, GCForeground | GCBackground |
                     GCLineWidth | GCCapStyle | GCJoinStyle | GCFillStyle,
                     &linevalues);
  trianglegc = XCreateGC(display, rootwindow, GCForeground | GCBackground |
                         GCLineWidth | GCCapStyle | GCJoinStyle | GCFillStyle,
                         &linevalues);
  make_buttons(height);
  XFlush(display);
}

void draw_node(nodes, dim, ptr, xscale, yscale, xoffset, yoffset)
int nodes;
int dim;
REAL *ptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i;
  int index;

  index = dim;
  for (i = 1; i <= nodes; i++) {
    XFillRectangle(display, mainwindow, linegc,
                   (int) (ptr[index] * xscale + xoffset) - (line_width >> 1),
                   (int) (ptr[index + 1] * yscale + yoffset) -
                         (line_width >> 1), line_width, line_width);
    index += dim;
  }
}

void draw_poly(nodes, dim, edges, holes, nodeptr, edgeptr, holeptr,
               xscale, yscale, xoffset, yoffset)
int nodes;
int dim;
int edges;
int holes;
REAL *nodeptr;
int *edgeptr;
REAL *holeptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i;
  int index;
  REAL *point1, *point2;
  int x1, y1, x2, y2;

  index = dim;
  for (i = 1; i <= nodes; i++) {
    XFillRectangle(display, mainwindow, linegc,
                   (int) (nodeptr[index] * xscale + xoffset) -
                         (line_width >> 1),
                   (int) (nodeptr[index + 1] * yscale + yoffset) -
                         (line_width >> 1), line_width, line_width);
    index += dim;
  }
  index = 2;
  for (i = 1; i <= edges; i++) {
    point1 = &nodeptr[edgeptr[index++] * dim];
    point2 = &nodeptr[edgeptr[index++] * dim];
    XDrawLine(display, mainwindow, linegc,
              (int) (point1[0] * xscale + xoffset),
              (int) (point1[1] * yscale + yoffset),
              (int) (point2[0] * xscale + xoffset),
              (int) (point2[1] * yscale + yoffset));
  }
  index = dim;
  if (color) {
    XSetForeground(display, linegc, colors[0]);
  }
  for (i = 1; i <= holes; i++) {
    x1 = (int) (holeptr[index] * xscale + xoffset) - 3;
    y1 = (int) (holeptr[index + 1] * yscale + yoffset) - 3;
    x2 = x1 + 6;
    y2 = y1 + 6;
    XDrawLine(display, mainwindow, linegc, x1, y1, x2, y2);
    XDrawLine(display, mainwindow, linegc, x1, y2, x2, y1);
    index += dim;
  }
  XSetForeground(display, linegc, showme_foreground);
}

void draw_ele(inc, elems, corners, ptr, partition, shift,
              xscale, yscale, xoffset, yoffset)
int inc;
int elems;
int corners;
int *ptr;
int *partition;
REAL *shift;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i, j;
  int index;
  REAL shiftx, shifty;
  REAL *prevpoint, *nowpoint;
  XPoint *vertices;

  if (color && fillelem && (partition != (int *) NULL)) {
    vertices = (XPoint *) malloc(3 * sizeof(XPoint));
    if (vertices == (XPoint *) NULL) {
      printf("Error:  Out of memory.\n");
      exit(1);
    }
  }
  index = 3;
  for (i = 1; i <= elems; i++) {
    if ((partition != (int *) NULL) && explode) {
      shiftx = shift[partition[i] << 1];
      shifty = shift[(partition[i] << 1) + 1];
    }
    if (color && (partition != (int *) NULL)) {
      if (fillelem) {
        XSetForeground(display, trianglegc, colors[partition[i] & 63]);
      } else {
        XSetForeground(display, linegc, colors[partition[i] & 63]);
      }
    }
    if (color && fillelem && (partition != (int *) NULL)) {
      if ((partition != (int *) NULL) && explode) {
        for (j = 0; j < 3; j++) {
          nowpoint = &nodeptr[inc][ptr[index + j] * node_dim[inc]];
          vertices[j].x = (nowpoint[0] + shiftx) * xscale + xoffset;
          vertices[j].y = (nowpoint[1] + shifty) * yscale + yoffset;
        }
      } else {
        for (j = 0; j < 3; j++) {
          nowpoint = &nodeptr[inc][ptr[index + j] * node_dim[inc]];
          vertices[j].x = nowpoint[0] * xscale + xoffset;
          vertices[j].y = nowpoint[1] * yscale + yoffset;
        }
      }
      XFillPolygon(display, mainwindow, trianglegc, vertices, 3,
                   Convex, CoordModeOrigin);
    }
    prevpoint = &nodeptr[inc][ptr[index + 2] * node_dim[inc]];
    if ((partition != (int *) NULL) && explode) {
      for (j = 0; j < 3; j++) {
        nowpoint = &nodeptr[inc][ptr[index++] * node_dim[inc]];
        XDrawLine(display, mainwindow, linegc,
                  (int) ((prevpoint[0] + shiftx) * xscale + xoffset),
                  (int) ((prevpoint[1] + shifty) * yscale + yoffset),
                  (int) ((nowpoint[0] + shiftx) * xscale + xoffset),
                  (int) ((nowpoint[1] + shifty) * yscale + yoffset));
        prevpoint = nowpoint;
      }
    } else {
      for (j = 0; j < 3; j++) {
        nowpoint = &nodeptr[inc][ptr[index++] * node_dim[inc]];
        XDrawLine(display, mainwindow, linegc,
                  (int) (prevpoint[0] * xscale + xoffset),
                  (int) (prevpoint[1] * yscale + yoffset),
                  (int) (nowpoint[0] * xscale + xoffset),
                  (int) (nowpoint[1] * yscale + yoffset));
        prevpoint = nowpoint;
      }
    }
  }
  if (color && fillelem && (partition != (int *) NULL)) {
    free(vertices);
  }
  XSetForeground(display, linegc, showme_foreground);
}

void draw_edge(nodes, dim, edges, nodeptr, edgeptr, normptr,
               xscale, yscale, xoffset, yoffset)
int nodes;
int dim;
int edges;
REAL *nodeptr;
int *edgeptr;
REAL *normptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i;
  int index;
  REAL *point1, *point2;
  REAL normx, normy;
  REAL normmult, normmultx, normmulty;
  REAL windowxmin, windowymin, windowxmax, windowymax;

  index = 2;
  for (i = 1; i <= edges; i++) {
    point1 = &nodeptr[edgeptr[index++] * dim];
    if (edgeptr[index] == -1) {
      normx = normptr[index - 1];
      normy = normptr[index++];
      normmultx = 0.0;
      if (normx > 0) {
        windowxmax = (width - 1 - xoffset) / xscale;
        normmultx = (windowxmax - point1[0]) / normx;
      } else if (normx < 0) {
        windowxmin = -xoffset / xscale;
        normmultx = (windowxmin - point1[0]) / normx;
      }
      normmulty = 0.0;
      if (normy > 0) {
        windowymax = -yoffset / yscale;
        normmulty = (windowymax - point1[1]) / normy;
      } else if (normy < 0) {
        windowymin = (height - 1 - yoffset) / yscale;
        normmulty = (windowymin - point1[1]) / normy;
      }
      if (normmultx == 0.0) {
        normmult = normmulty;
      } else if (normmulty == 0.0) {
        normmult = normmultx;
      } else if (normmultx < normmulty) {
        normmult = normmultx;
      } else {
        normmult = normmulty;
      }
      if (normmult > 0.0) {
        XDrawLine(display, mainwindow, linegc,
                  (int) (point1[0] * xscale + xoffset),
                  (int) (point1[1] * yscale + yoffset),
                  (int) ((point1[0] + normmult * normx) * xscale + xoffset),
                  (int) ((point1[1] + normmult * normy) * yscale + yoffset));
      }
    } else {
      point2 = &nodeptr[edgeptr[index++] * dim];
      XDrawLine(display, mainwindow, linegc,
                (int) (point1[0] * xscale + xoffset),
                (int) (point1[1] * yscale + yoffset),
                (int) (point2[0] * xscale + xoffset),
                (int) (point2[1] * yscale + yoffset));
    }
  }
}

void draw_adj(dim, subdomains, ptr, center, xscale, yscale,
              xoffset, yoffset)
int dim;
int subdomains;
int *ptr;
REAL *center;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i, j;
  REAL *point1, *point2;

  for (i = 0; i < subdomains; i++) {
    for (j = i + 1; j < subdomains; j++) {
      if (ptr[i * subdomains + j]) {
        point1 = &center[i * dim];
        point2 = &center[j * dim];
        XDrawLine(display, mainwindow, linegc,
                  (int) (point1[0] * xscale + xoffset),
                  (int) (point1[1] * yscale + yoffset),
                  (int) (point2[0] * xscale + xoffset),
                  (int) (point2[1] * yscale + yoffset));
      }
    }
  }
  for (i = 0; i < subdomains; i++) {
    point1 = &center[i * dim];
    if (color) {
      XSetForeground(display, linegc, colors[i & 63]);
    }
    XFillArc(display, mainwindow, linegc,
             (int) (point1[0] * xscale + xoffset) - 5 - (line_width >> 1),
             (int) (point1[1] * yscale + yoffset) - 5 - (line_width >> 1),
             line_width + 10, line_width + 10, 0, 23040);
  }
  XSetForeground(display, linegc, showme_foreground);
}

void draw(inc, image, xmin, ymin, xmax, ymax)
int inc;
int image;
REAL xmin;
REAL ymin;
REAL xmax;
REAL ymax;
{
  draw_buttons();
  XClearWindow(display, mainwindow);
  if (image == NOTHING) {
    return;
  }
  if (!loaded[inc][image]) {
    return;
  }
  if ((image == PART) && explode) {
    xmin += (xmin - partcenter[inc][subdomains[inc] << 1]) * explosion;
    xmax += (xmax - partcenter[inc][subdomains[inc] << 1]) * explosion;
    ymin += (ymin - partcenter[inc][(subdomains[inc] << 1) + 1]) * explosion;
    ymax += (ymax - partcenter[inc][(subdomains[inc] << 1) + 1]) * explosion;
  }
  xscale = (REAL) (width - line_width - 4) / (xmax - xmin);
  yscale = (REAL) (height - line_width - 4) / (ymax - ymin);
  if (xscale > yscale) {
    xscale = yscale;
  } else {
    yscale = xscale;
  }
  xoffset = 0.5 * ((REAL) width - xscale * (xmax - xmin)) -
            xscale * xmin;
  yoffset = (REAL) height - 0.5 * ((REAL) height - yscale * (ymax - ymin)) +
            yscale * ymin;
  yscale = - yscale;
  switch(image) {
    case NODE:
      draw_node(nodes[inc], node_dim[inc], nodeptr[inc],
                xscale, yscale, xoffset, yoffset);
      break;
    case POLY:
      if (polynodes[inc] > 0) {
        draw_poly(polynodes[inc], poly_dim[inc], polyedges[inc],
                  polyholes[inc], polynodeptr[inc], polyedgeptr[inc],
                  polyholeptr[inc],
                  xscale, yscale, xoffset, yoffset);
      } else {
        draw_poly(nodes[inc], node_dim[inc], polyedges[inc],
                  polyholes[inc], nodeptr[inc], polyedgeptr[inc],
                  polyholeptr[inc],
                  xscale, yscale, xoffset, yoffset);
      }
      break;
    case ELE:
      draw_ele(inc, elems[inc], ele_corners[inc], eleptr[inc],
               (int *) NULL, (REAL *) NULL,
               xscale, yscale, xoffset, yoffset);
      break;
    case EDGE:
      draw_edge(nodes[inc], node_dim[inc], edges[inc],
                nodeptr[inc], edgeptr[inc], normptr[inc],
                xscale, yscale, xoffset, yoffset);
      break;
    case PART:
      draw_ele(inc, elems[inc], ele_corners[inc], eleptr[inc],
               partpart[inc], partshift[inc],
               xscale, yscale, xoffset, yoffset);
      break;
    case ADJ:
      draw_adj(node_dim[inc], adjsubdomains[inc], adjptr[inc], partcenter[inc],
               xscale, yscale, xoffset, yoffset);
      break;
    case VORO:
      if (loaded[inc][NODE]) {
        draw_node(nodes[inc], node_dim[inc], nodeptr[inc],
                  xscale, yscale, xoffset, yoffset);
      }
      draw_edge(vnodes[inc], vnode_dim[inc], vedges[inc],
                vnodeptr[inc], vedgeptr[inc], vnormptr[inc],
                xscale, yscale, xoffset, yoffset);
      break;
    default:
      break;
  }
}

void addps(instring, outstring, eps)
char *instring;
char *outstring;
int eps;
{
  strcpy(outstring, instring);
  if (eps) {
    strcat(outstring, ".eps");
  } else {
    strcat(outstring, ".ps");
  }
}

int print_head(fname, file, llcornerx, llcornery, eps)
char *fname;
FILE **file;
int llcornerx;
int llcornery;
int eps;
{
  if (!quiet) {
    printf("Writing %s\n", fname);
  }
  *file = fopen(fname, "w");
  if (*file == (FILE *) NULL) {
    printf("  Error:  Could not open %s\n", fname);
    return 1;
  }
  if (eps) {
    fprintf(*file, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  } else {
    fprintf(*file, "%%!PS-Adobe-2.0\n");
  }
  fprintf(*file, "%%%%BoundingBox: %d %d %d %d\n", llcornerx, llcornery,
          612 - llcornerx, 792 - llcornery);
  fprintf(*file, "%%%%Creator: Show Me\n");
  fprintf(*file, "%%%%EndComments\n\n");
  fprintf(*file, "/m {moveto} bind def\n");
  fprintf(*file, "/l {lineto} bind def\n");
  fprintf(*file, "/s {setrgbcolor} bind def\n");
  fprintf(*file, "/g {gsave fill grestore} bind def\n");
  fprintf(*file, "/k {stroke} bind def\n\n");
  fprintf(*file, "1 setlinecap\n");
  fprintf(*file, "1 setlinejoin\n");
  fprintf(*file, "%d setlinewidth\n", line_width);
  fprintf(*file, "%d %d m\n", llcornerx, llcornery);
  fprintf(*file, "%d %d l\n", 612 - llcornerx, llcornery);
  fprintf(*file, "%d %d l\n", 612 - llcornerx, 792 - llcornery);
  fprintf(*file, "%d %d l\n", llcornerx, 792 - llcornery);
  fprintf(*file, "closepath\nclip\nnewpath\n");
  return 0;
}

void print_node(nodefile, nodes, dim, ptr, xscale, yscale,
                xoffset, yoffset)
FILE *nodefile;
int nodes;
int dim;
REAL *ptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i;
  int index;

  index = dim;
  for (i = 1; i <= nodes; i++) {
    fprintf(nodefile, "%d %d %d 0 360 arc\nfill\n",
            (int) (ptr[index] * xscale + xoffset),
            (int) (ptr[index + 1] * yscale + yoffset),
            1 + (line_width >> 1));
    index += dim;
  }
}

void print_poly(polyfile, nodes, dim, edges, holes, nodeptr, edgeptr, holeptr,
                xscale, yscale, xoffset, yoffset)
FILE *polyfile;
int nodes;
int dim;
int edges;
int holes;
REAL *nodeptr;
int *edgeptr;
REAL *holeptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
{
  int i;
  int index;
  REAL *point1, *point2;

  index = dim;
  for (i = 1; i <= nodes; i++) {
    fprintf(polyfile, "%d %d %d 0 360 arc\nfill\n",
            (int) (nodeptr[index] * xscale + xoffset),
            (int) (nodeptr[index + 1] * yscale + yoffset),
            1 + (line_width >> 1));
    index += dim;
  }
  index = 2;
  for (i = 1; i <= edges; i++) {
    point1 = &nodeptr[edgeptr[index++] * dim];
    point2 = &nodeptr[edgeptr[index++] * dim];
    fprintf(polyfile, "%d %d m\n",
            (int) (point1[0] * xscale + xoffset),
            (int) (point1[1] * yscale + yoffset));
    fprintf(polyfile, "%d %d l\nk\n",
            (int) (point2[0] * xscale + xoffset),
            (int) (point2[1] * yscale + yoffset));
  }
}

void print_ele(elefile, nodes, dim, elems, corners, nodeptr, eleptr,
               partition, shift,
               xscale, yscale, xoffset, yoffset, llcornerx, llcornery)
FILE *elefile;
int nodes;
int dim;
int elems;
int corners;
REAL *nodeptr;
int *eleptr;
int *partition;
REAL *shift;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
int llcornerx;
int llcornery;
{
  int i, j;
  int index, colorindex;
  REAL shiftx, shifty;
  REAL *nowpoint;

  index = 3;
  if ((partition != (int *) NULL) && !bw_ps) {
    fprintf(elefile, "0 0 0 s\n");
    fprintf(elefile, "%d %d m\n", llcornerx, llcornery);
    fprintf(elefile, "%d %d l\n", 612 - llcornerx, llcornery);
    fprintf(elefile, "%d %d l\n", 612 - llcornerx, 792 - llcornery);
    fprintf(elefile, "%d %d l\n", llcornerx, 792 - llcornery);
    fprintf(elefile, "fill\n");
  }
  for (i = 1; i <= elems; i++) {
    if ((partition != (int *) NULL) && !bw_ps) {
      colorindex = partition[i] & 63;
      fprintf(elefile, "%6.3f %6.3f %6.3f s\n",
              (REAL) rgb[colorindex].red / 65535.0,
              (REAL) rgb[colorindex].green / 65535.0,
              (REAL) rgb[colorindex].blue / 65535.0);
    }
    nowpoint = &nodeptr[eleptr[index + 2] * dim];
    if ((partition != (int *) NULL) && (explode || bw_ps)) {
      shiftx = shift[partition[i] << 1];
      shifty = shift[(partition[i] << 1) + 1];
      fprintf(elefile, "%d %d m\n",
              (int) ((nowpoint[0] + shiftx) * xscale + xoffset),
              (int) ((nowpoint[1] + shifty) * yscale + yoffset));
      for (j = 0; j < 3; j++) {
        nowpoint = &nodeptr[eleptr[index++] * dim];
        fprintf(elefile, "%d %d l\n",
                (int) ((nowpoint[0] + shiftx) * xscale + xoffset),
                (int) ((nowpoint[1] + shifty) * yscale + yoffset));
      }
    } else {
      fprintf(elefile, "%d %d m\n",
              (int) (nowpoint[0] * xscale + xoffset),
              (int) (nowpoint[1] * yscale + yoffset));
      for (j = 0; j < 3; j++) {
        nowpoint = &nodeptr[eleptr[index++] * dim];
        fprintf(elefile, "%d %d l\n",
                (int) (nowpoint[0] * xscale + xoffset),
                (int) (nowpoint[1] * yscale + yoffset));
      }
    }
    if (fillelem && (partition != (int *) NULL) && !bw_ps) {
      fprintf(elefile, "g\n1 1 0 s\n");
    }
    fprintf(elefile, "k\n");
  }
}

void print_edge(edgefile, nodes, dim, edges, nodeptr, edgeptr, normptr,
                xscale, yscale, xoffset, yoffset, llcornerx, llcornery)
FILE *edgefile;
int nodes;
int dim;
int edges;
REAL *nodeptr;
int *edgeptr;
REAL *normptr;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
int llcornerx;
int llcornery;
{
  int i;
  int index;
  REAL *point1, *point2;
  REAL normx, normy;
  REAL normmult, normmultx, normmulty;
  REAL windowxmin, windowymin, windowxmax, windowymax;

  index = 2;
  for (i = 1; i <= edges; i++) {
    point1 = &nodeptr[edgeptr[index++] * dim];
    if (edgeptr[index] == -1) {
      normx = normptr[index - 1];
      normy = normptr[index++];
      normmultx = 0.0;
      if (normx > 0) {
        windowxmax = ((REAL) (612 - llcornerx) - xoffset) / xscale;
        normmultx = (windowxmax - point1[0]) / normx;
      } else if (normx < 0) {
        windowxmin = ((REAL) llcornerx - xoffset) / xscale;
        normmultx = (windowxmin - point1[0]) / normx;
      }
      normmulty = 0.0;
      if (normy > 0) {
        windowymax = ((REAL) (792 - llcornery) - yoffset) / yscale;
        normmulty = (windowymax - point1[1]) / normy;
      } else if (normy < 0) {
        windowymin = ((REAL) llcornery - yoffset) / yscale;
        normmulty = (windowymin - point1[1]) / normy;
      }
      if (normmultx == 0.0) {
        normmult = normmulty;
      } else if (normmulty == 0.0) {
        normmult = normmultx;
      } else if (normmultx < normmulty) {
        normmult = normmultx;
      } else {
        normmult = normmulty;
      }
      if (normmult > 0.0) {
        fprintf(edgefile, "%d %d m\n",
                (int) (point1[0] * xscale + xoffset),
                (int) (point1[1] * yscale + yoffset));
        fprintf(edgefile, "%d %d l\nk\n",
                (int) ((point1[0] + normmult * normx) * xscale + xoffset),
                (int) ((point1[1] + normmult * normy) * yscale + yoffset));
      }
    } else {
      point2 = &nodeptr[edgeptr[index++] * dim];
      fprintf(edgefile, "%d %d m\n",
              (int) (point1[0] * xscale + xoffset),
              (int) (point1[1] * yscale + yoffset));
      fprintf(edgefile, "%d %d l\nk\n",
              (int) (point2[0] * xscale + xoffset),
              (int) (point2[1] * yscale + yoffset));
    }
  }
}

void print_adj(adjfile, dim, subdomains, ptr, center, xscale, yscale,
               xoffset, yoffset, llcornerx, llcornery)
FILE *adjfile;
int dim;
int subdomains;
int *ptr;
REAL *center;
REAL xscale;
REAL yscale;
REAL xoffset;
REAL yoffset;
int llcornerx;
int llcornery;
{
  int i, j;
  REAL *point1, *point2;
  int colorindex;

  if (!bw_ps) {
    fprintf(adjfile, "0 0 0 s\n");
    fprintf(adjfile, "%d %d m\n", llcornerx, llcornery);
    fprintf(adjfile, "%d %d l\n", 612 - llcornerx, llcornery);
    fprintf(adjfile, "%d %d l\n", 612 - llcornerx, 792 - llcornery);
    fprintf(adjfile, "%d %d l\n", llcornerx, 792 - llcornery);
    fprintf(adjfile, "fill\n");
    fprintf(adjfile, "1 1 0 s\n");
  }
  for (i = 0; i < subdomains; i++) {
    for (j = i + 1; j < subdomains; j++) {
      if (ptr[i * subdomains + j]) {
        point1 = &center[i * dim];
        point2 = &center[j * dim];
        fprintf(adjfile, "%d %d m\n",
                (int) (point1[0] * xscale + xoffset),
                (int) (point1[1] * yscale + yoffset));
        fprintf(adjfile, "%d %d l\nk\n",
                (int) (point2[0] * xscale + xoffset),
                (int) (point2[1] * yscale + yoffset));
      }
    }
  }
  for (i = 0; i < subdomains; i++) {
    point1 = &center[i * dim];
    if (!bw_ps) {
      colorindex = i & 63;
      fprintf(adjfile, "%6.3f %6.3f %6.3f s\n",
              (REAL) rgb[colorindex].red / 65535.0,
              (REAL) rgb[colorindex].green / 65535.0,
              (REAL) rgb[colorindex].blue / 65535.0);
      fprintf(adjfile, "%d %d %d 0 360 arc\nfill\n",
              (int) (point1[0] * xscale + xoffset),
              (int) (point1[1] * yscale + yoffset),
              5 + (line_width >> 1));
    } else {
      fprintf(adjfile, "%d %d %d 0 360 arc\nfill\n",
              (int) (point1[0] * xscale + xoffset),
              (int) (point1[1] * yscale + yoffset),
              3 + (line_width >> 1));
    }
  }
}

void print(inc, image, xmin, ymin, xmax, ymax, eps)
int inc;
int image;
REAL xmin;
REAL ymin;
REAL xmax;
REAL ymax;
int eps;
{
  REAL xxscale, yyscale, xxoffset, yyoffset;
  char psfilename[FILENAMESIZE];
  int llcornerx, llcornery;
  FILE *psfile;

  if (image == NOTHING) {
    return;
  }
  if (!loaded[inc][image]) {
    return;
  }
  if ((image == PART) && (explode || bw_ps)) {
    xmin += (xmin - partcenter[inc][subdomains[inc] << 1]) * explosion;
    xmax += (xmax - partcenter[inc][subdomains[inc] << 1]) * explosion;
    ymin += (ymin - partcenter[inc][(subdomains[inc] << 1) + 1]) * explosion;
    ymax += (ymax - partcenter[inc][(subdomains[inc] << 1) + 1]) * explosion;
  }
  xxscale = (460.0 - (REAL) line_width) / (xmax - xmin);
  yyscale = (640.0 - (REAL) line_width) / (ymax - ymin);
  if (xxscale > yyscale) {
    xxscale = yyscale;
    llcornerx = (604 - (int) (yyscale * (xmax - xmin)) - line_width) >> 1;
    llcornery = 72;
  } else {
    yyscale = xxscale;
    llcornerx = 72;
    llcornery = (784 - (int) (xxscale * (ymax - ymin)) - line_width) >> 1;
  }
  xxoffset = 0.5 * (612.0 - xxscale * (xmax - xmin)) - xxscale * xmin +
             (line_width >> 1);
  yyoffset = 0.5 * (792.0 - yyscale * (ymax - ymin)) - yyscale * ymin +
             (line_width >> 1);
  switch(image) {
    case NODE:
      addps(nodefilename[inc], psfilename, eps);
      break;
    case POLY:
      addps(polyfilename[inc], psfilename, eps);
      break;
    case ELE:
      addps(elefilename[inc], psfilename, eps);
      break;
    case EDGE:
      addps(edgefilename[inc], psfilename, eps);
      break;
    case PART:
      addps(partfilename[inc], psfilename, eps);
      break;
    case ADJ:
      addps(adjfilename[inc], psfilename, eps);
      break;
    case VORO:
      addps(vedgefilename[inc], psfilename, eps);
      break;
    default:
      break;
  }
  if (print_head(psfilename, &psfile, llcornerx, llcornery, eps)) {
    return;
  }
  switch(image) {
    case NODE:
      print_node(psfile, nodes[inc], node_dim[inc], nodeptr[inc],
                 xxscale, yyscale, xxoffset, yyoffset);
      break;
    case POLY:
      if (polynodes[inc] > 0) {
        print_poly(psfile, polynodes[inc], poly_dim[inc], polyedges[inc],
                   polyholes[inc], polynodeptr[inc], polyedgeptr[inc],
                   polyholeptr[inc], xxscale, yyscale, xxoffset, yyoffset);
      } else {
        print_poly(psfile, nodes[inc], node_dim[inc], polyedges[inc],
                   polyholes[inc], nodeptr[inc], polyedgeptr[inc],
                   polyholeptr[inc], xxscale, yyscale, xxoffset, yyoffset);
      }
      break;
    case ELE:
      print_ele(psfile, nodes[inc], node_dim[inc], elems[inc],
                ele_corners[inc], nodeptr[inc], eleptr[inc],
                (int *) NULL, (REAL *) NULL,
                xxscale, yyscale, xxoffset, yyoffset, llcornerx, llcornery);
      break;
    case EDGE:
      print_edge(psfile, nodes[inc], node_dim[inc], edges[inc],
                 nodeptr[inc], edgeptr[inc], normptr[inc],
                 xxscale, yyscale, xxoffset, yyoffset, llcornerx, llcornery);
      break;
    case PART:
      print_ele(psfile, nodes[inc], node_dim[inc], elems[inc],
                ele_corners[inc], nodeptr[inc], eleptr[inc],
                partpart[inc], partshift[inc],
                xxscale, yyscale, xxoffset, yyoffset, llcornerx, llcornery);
      break;
    case ADJ:
      print_adj(psfile, node_dim[inc], adjsubdomains[inc], adjptr[inc],
                partcenter[inc],
                xxscale, yyscale, xxoffset, yyoffset, llcornerx, llcornery);
      break;
    case VORO:
      print_edge(psfile, vnodes[inc], vnode_dim[inc], vedges[inc],
                 vnodeptr[inc], vedgeptr[inc], vnormptr[inc],
                 xxscale, yyscale, xxoffset, yyoffset, llcornerx, llcornery);
      break;
    default:
      break;
  }
  if (!eps) {
    fprintf(psfile, "showpage\n");
  }
  fclose(psfile);
}

int main(argc, argv)
int argc;
char **argv;
{
  REAL xmin, ymin, xmax, ymax;
  REAL xptr, yptr, xspan, yspan;
  int past_image;
  int new_image;
  int new_inc;

  parsecommandline(argc, argv);
  showme_init();
  choose_image(start_inc, start_image);
  showme_window(argc, argv);

  if (current_image != NOTHING) {
    xmin = xlo[current_inc][current_image];
    ymin = ylo[current_inc][current_image];
    xmax = xhi[current_inc][current_image];
    ymax = yhi[current_inc][current_image];
    zoom = 0;
  }

  XMaskEvent(display, ExposureMask, &event);
  while (1) {
    switch (event.type) {
      case ButtonRelease:
        if (event.xany.window == quitwin) {
          XDestroyWindow(display, mainwindow);
          XCloseDisplay(display);
          return 0;
        } else if (event.xany.window == leftwin) {
          xspan = 0.25 * (xmax - xmin);
          xmin += xspan;
          xmax += xspan;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else if (event.xany.window == rightwin) {
          xspan = 0.25 * (xmax - xmin);
          xmin -= xspan;
          xmax -= xspan;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else if (event.xany.window == upwin) {
          yspan = 0.25 * (ymax - ymin);
          ymin -= yspan;
          ymax -= yspan;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else if (event.xany.window == downwin) {
          yspan = 0.25 * (ymax - ymin);
          ymin += yspan;
          ymax += yspan;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else if (event.xany.window == resetwin) {
          xmin = xlo[current_inc][current_image];
          ymin = ylo[current_inc][current_image];
          xmax = xhi[current_inc][current_image];
          ymax = yhi[current_inc][current_image];
          zoom = 0;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else if (event.xany.window == widthpluswin) {
          if (line_width < 100) {
            line_width++;
            XSetLineAttributes(display, linegc, line_width, LineSolid,
                               CapRound, JoinRound);
            XSetLineAttributes(display, trianglegc, line_width, LineSolid,
                               CapRound, JoinRound);
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == widthminuswin) {
          if (line_width > 1) {
            line_width--;
            XSetLineAttributes(display, linegc, line_width, LineSolid,
                               CapRound, JoinRound);
            XSetLineAttributes(display, trianglegc, line_width, LineSolid,
                               CapRound, JoinRound);
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == expwin) {
          if ((current_image == PART) && loaded[current_inc][PART]) {
            explode = !explode;
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == exppluswin) {
          if ((current_image == PART) && loaded[PART] && explode) {
            explosion += 0.125;
            findpartshift(subdomains[current_inc], explosion,
                          partcenter[current_inc], partshift[current_inc]);
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == expminuswin) {
          if ((current_image == PART) && loaded[PART] && explode &&
              (explosion >= 0.125)) {
            explosion -= 0.125;
            findpartshift(subdomains[current_inc], explosion,
                          partcenter[current_inc], partshift[current_inc]);
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == fillwin) {
          if ((current_image == PART) && loaded[PART]) {
            fillelem = !fillelem;
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          }
        } else if (event.xany.window == pswin) {
          fill_button(pswin);
          XFlush(display);
          print(current_inc, current_image, xmin, ymin, xmax, ymax, 0);
          XClearWindow(display, pswin);
          XDrawString(display, pswin, fontgc, 2, 13, "PS", 2);
        } else if (event.xany.window == epswin) {
          fill_button(epswin);
          XFlush(display);
          print(current_inc, current_image, xmin, ymin, xmax, ymax, 1);
          XClearWindow(display, epswin);
          XDrawString(display, epswin, fontgc, 2, 13, "EPS", 3);
        } else if (event.xany.window == versionpluswin) {
          move_inc(1);
          loweriteration++;
          set_filenames(filename, loweriteration);
          if (current_inc == 1) {
            current_inc = 0;
          } else {
            current_image = NOTHING;
            XClearWindow(display, mainwindow);
          }
          draw_buttons();
        } else if (event.xany.window == versionminuswin) {
          if (loweriteration > 0) {
            move_inc(0);
            loweriteration--;
            set_filenames(filename, loweriteration);
            if (current_inc == 0) {
              current_inc = 1;
            } else {
              current_image = NOTHING;
              XClearWindow(display, mainwindow);
            }
            draw_buttons();
          }
        } else if ((event.xany.window == nodewin[0]) ||
                   (event.xany.window == polywin[0]) ||
                   (event.xany.window == elewin[0]) ||
                   (event.xany.window == edgewin[0]) ||
                   (event.xany.window == partwin[0]) ||
                   (event.xany.window == adjwin[0]) ||
                   (event.xany.window == voronoiwin[0]) ||
                   (event.xany.window == nodewin[1]) ||
                   (event.xany.window == polywin[1]) ||
                   (event.xany.window == elewin[1]) ||
                   (event.xany.window == edgewin[1]) ||
                   (event.xany.window == partwin[1]) ||
                   (event.xany.window == adjwin[1]) ||
                   (event.xany.window == voronoiwin[1])) {
          if (event.xany.window == nodewin[0]) {
            new_inc = 0;
            new_image = NODE;
          }
          if (event.xany.window == polywin[0]) {
            new_inc = 0;
            new_image = POLY;
          }
          if (event.xany.window == elewin[0]) {
            new_inc = 0;
            new_image = ELE;
          }
          if (event.xany.window == edgewin[0]) {
            new_inc = 0;
            new_image = EDGE;
          }
          if (event.xany.window == partwin[0]) {
            new_inc = 0;
            new_image = PART;
          }
          if (event.xany.window == adjwin[0]) {
            new_inc = 0;
            new_image = ADJ;
          }
          if (event.xany.window == voronoiwin[0]) {
            new_inc = 0;
            new_image = VORO;
          }
          if (event.xany.window == nodewin[1]) {
            new_inc = 1;
            new_image = NODE;
          }
          if (event.xany.window == polywin[1]) {
            new_inc = 1;
            new_image = POLY;
          }
          if (event.xany.window == elewin[1]) {
            new_inc = 1;
            new_image = ELE;
          }
          if (event.xany.window == edgewin[1]) {
            new_inc = 1;
            new_image = EDGE;
          }
          if (event.xany.window == partwin[1]) {
            new_inc = 1;
            new_image = PART;
          }
          if (event.xany.window == adjwin[1]) {
            new_inc = 1;
            new_image = ADJ;
          }
          if (event.xany.window == voronoiwin[1]) {
            new_inc = 1;
            new_image = VORO;
          }
          past_image = current_image;
          if ((current_inc == new_inc) && (current_image == new_image)) {
            free_inc(new_inc);
            unload_inc(new_inc);
          }
          choose_image(new_inc, new_image);
          if ((past_image == NOTHING) && (current_image != NOTHING)) {
            xmin = xlo[current_inc][current_image];
            ymin = ylo[current_inc][current_image];
            xmax = xhi[current_inc][current_image];
            ymax = yhi[current_inc][current_image];
            zoom = 0;
          }
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        } else {
          xptr = ((REAL) event.xbutton.x - xoffset) / xscale;
          yptr = ((REAL) event.xbutton.y - yoffset) / yscale;
          if ((current_image == PART) && loaded[PART] && explode) {
            xptr = (xptr + partcenter[current_inc]
                                     [subdomains[current_inc] << 1]
                    * explosion) / (1.0 + explosion);
            yptr = (yptr + partcenter[current_inc]
                                     [(subdomains[current_inc] << 1) + 1]
                    * explosion) / (1.0 + explosion);
          }
          if ((event.xbutton.button == Button1)
              || (event.xbutton.button == Button3)) {
            if (event.xbutton.button == Button1) {
              xspan = 0.25 * (xmax - xmin);
              yspan = 0.25 * (ymax - ymin);
              zoom++;
            } else {
              xspan = xmax - xmin;
              yspan = ymax - ymin;
              zoom--;
            }
            xmin = xptr - xspan;
            ymin = yptr - yspan;
            xmax = xptr + xspan;
            ymax = yptr + yspan;
            draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          } else if (event.xbutton.button == Button2) {
            printf("x = %.4g, y = %.4g\n", xptr, yptr);
          }
        }
        break;
      case DestroyNotify:
        XDestroyWindow(display, mainwindow);
        XCloseDisplay(display);
        return 0;
      case ConfigureNotify:
        if ((width != event.xconfigure.width) ||
            (height != event.xconfigure.height - PANELHEIGHT)) {
          width = event.xconfigure.width;
          height = event.xconfigure.height - PANELHEIGHT;
          draw(current_inc, current_image, xmin, ymin, xmax, ymax);
          while (XCheckMaskEvent(display, ExposureMask, &event));
        }
        break;
      case Expose:
        draw(current_inc, current_image, xmin, ymin, xmax, ymax);
        while (XCheckMaskEvent(display, ExposureMask, &event));
        break;
      default:
        break;
    }
    XNextEvent(display, &event);
  }
}
