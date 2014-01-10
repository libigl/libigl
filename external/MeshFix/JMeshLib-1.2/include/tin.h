/****************************************************************************
* JMeshLib                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef _TIN_H
#define _TIN_H

#include "jmesh.h"

//! Triangulation

//! This class represents a manifold and oriented triangle mesh.
//! Vertices, Edges and Triangles are stored in the Lists V, E and T
//! respectively. Methods boundaries(), handles() and shells() may
//! be used to retrieve the respective topological entities.
//! Navigation of the mesh is based on the topological relationships
//! stored in each Vertex, Edge and Triangle.
//! Some methods would require a global update only to maintain
//! consistent values of the protected fields n_boundaries, n_handles
//! and n_shells. On the other hand, the same methods would work in
//! constant time if these values would not need to be updated.
//! To keep complexity as low as possible,
//! we make use of the 'dirty bits' d_boundaries, d_handles and 
//! d_shells to mark that the respective entities must be updated.
//! The access functions boundaries(), handles() and shells()
//! check the status of the respective dirty bit and do a global
//! update (i.e., eulerUpdate()) only if necessary.
//! The complexity of the methods is provided as a function of a
//! generic 'N', which is O(V.numels()) = O(E.numels()) = O(T.numels()).

class Triangulation
{
 protected:

 int n_boundaries;		//!< Number of boundary loops
 int n_handles;			//!< Number of handles
 int n_shells;			//!< Number of connected components

 bool d_boundaries;		//!< Dirty bit for n_boundaries
 bool d_handles;		//!< Dirty bit for n_handles
 bool d_shells;			//!< Dirty bit for n_shells

 public :
 
 List V;			//!< Vertex set
 List E;			//!< Edge set
 List T;			//!< Triangle set


 // Constructors
 
 //! Empty triangulation. Should be used only prior to a call to load().
 Triangulation();

 //! Pre-defined triangulation. Currently, only "tetrahedron" is recognized.
 Triangulation(const char *);

 //! Clones an existing Trianglation.
 Triangulation(const Triangulation *, const bool clone_info =false);

 //! Clones an existing connected component.

 //! Creates a new Triangulation out of a connected component of an existing
 //! Triangulation. 't' is a triangle of the connected component that must
 //! be copied. If 'keep_ref' is TRUE, each element of the existing mesh
 //! keeps a pointer to the corresponding new element in the 'info' field.
 Triangulation(const Triangle *t, const bool keep_ref =false);

 //! Destructor. Frees the memory allocated for all the mesh elements.
 //! Warning! This method uses the freeNodes() method of the class List,
 //! which is not guaranteed to work correctly on systems other than
 //! Linux (see the documentation of List for details).
 //! Assuming that the method works correctly, however, the calling
 //! function is responsible of freeing the memory that was possibly
 //! allocated for objects pointed to by the 'info' field of mesh
 //! elements. Clearly, this must be done before calling the destructor.
 ~Triangulation();


 //! Get the number of boundary loops of the triangle mesh. O(1) or O(N).
 int boundaries() {if (d_boundaries) eulerUpdate(); return n_boundaries;}

 //! Get the number of handles of the triangle mesh. O(1) or O(N).
 int handles() {if (d_handles) eulerUpdate(); return n_handles;}

 //! Get the number of connected components of the triangle mesh. O(1) or O(N).
 int shells() {if (d_shells) eulerUpdate(); return n_shells;}


 /////////////////////////////////////////////////////////////////////////////
 //
 // Input/Output methods (Implemented in "MESH_STRUCTURE/io.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Initialize the triangle mesh from the file 'filename'.

 //! The file format is automatically deduced from the magic number
 //! or the filename extension. If 'doupdate' is FALSE, the
 //! global update for the topological entities is prevented.
 //! Currently, the following file formats are supported:
 //! Open Inventor (IV), VRML 1.0 and 2.0 (WRL), Object File Format (OFF),
 //! IMATI Ver-Tri (VER, TRI).
 //! A non-zero value is returned in case of error. Specifically,
 //! IO_CANTOPEN means that the file couldn't be opened for reading.
 //! IO_FORMAT means that the file format was not recognized by the loader.
 //! IO_UNKNOWN represents all the other errors.
 //! The calling function is responsible of verifying that the mesh is
 //! empty before calling this method.

 int load(const char *filename, const bool update=1);

 int cutAndStitch();	//!< Convert to manifold
 bool CreateIndexedTriangle(ExtVertex **, int, int, int);
 int loadIV(const char *);		//!< Loads IV
 int loadVRML1(const char *);		//!< Loads VRML 1.0
 int loadOFF(const char *);		//!< Loads OFF
 int loadPLY(const char *);		//!< Loads PLY
 int loadVerTri(const char *);		//!< Loads VER-TRI
 int loadVRML2(const char *);		//!< Loads VRML 2.0
 int loadOBJ(const char *);		//!< Loads OBJ
 int loadSTL(const char *);		//!< Loads STL

 protected:
 void closeLoadingSession(FILE *, int, ExtVertex **, bool);
 void coordBackApproximation();

 public:

 //! Save the triangle mesh to file 'filename'.

 //! The file format is deduced from the filename extension
 //! (wrl = vrml 1.0), (iv = OpenInventor), (off = Object
 //! file format), (ply = PLY format), (tri = IMATI Ver-Tri).
 //! If 'back_approx' is set, vertex coordinates are approximated
 //! to reflect the limited precision of floating point
 //! representation in ASCII files. This should be used when
 //! coherence is necessary between in-memory and saved data.
 //! A non-zero return value is returned if errors occur.
 
 int save(const char *filename, bool back_approx=0);

 int saveIV(const char *);		//!< Saves IV
 int saveOFF(const char *);		//!< Saves OFF 1.0
 int saveOBJ(const char *);		//!< Saves OBJ
 int saveSTL(const char *);		//!< Saves STL
 int savePLY(const char *, bool ascii = 1); //!< Saves PLY 1.0 (ascii or binary)
 int saveVerTri(const char *);		//!< Saves Ver-Tri

 //! Saves the triangle mesh to a VRML 1.0 file.
 //! The value of 'mode' specifies whether to use additional
 //! information attached to mesh elements in order to assign
 //! them a proper color.
 //! IO_CSAVE_OVERALL assigns a unique color for the entire mesh (default).
 //! IO_CSAVE_PERFACE assigns a color to each triangle depending on the value
 //! of its 'info' field.
 //! IO_CSAVE_PERVERTEX	assigns a color to each vertex depending on the value
 //! of its 'info' field.
 //! IO_CSAVE_PERFACE_INDEXED assigns one of five base colors to each triangle
 //! depending on the value of its 'mask' field.
 //! IO_CSAVE_PERVERTEX_INDEXED assigns one of five base colors to each vertex
 //! depending on the value of its 'mask' field.
 int saveVRML1(const char *, const int mode=0);


 //! Append another triangle mesh to the existing one.

 //! This method works exactly as the 'load()' method, except for the fact
 //! that it does not assume that the mesh is empty.
 int append(const char *filename, const bool doupdate=1);


 /////////////////////////////////////////////////////////////////////////////
 //
 // Primitive Construction (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Creates an Edge connecting two existing mesh vertices.

 //! Returns the newly created edge. If an edge connecting the two vertices
 //! already exists in the mesh, then no new edge is created and the old one
 //! is returned.
 Edge     *CreateEdge(Vertex *v1, Vertex *v2);


 //! Creates an Edge connecting two existing mesh Extended vertices.

 //! Returns the newly created edge. If an edge connecting the two vertices
 //! already exists in the mesh, then no new edge is created and the old one
 //! is returned.
 //! If 'check' is FALSE, the check for previously existing edges is skipped.
 Edge     *CreateEdge(ExtVertex *v1, ExtVertex *v2, const bool check=1);


 //! Creates a properly oriented Triangle bounded by three existing mesh edges.

 //! Returns the newly created Triangle. If e1, e2 and e3
 //! are not suitable for creating a properly oriented and
 //! manifold triangle, the creation fails and NULL is returned.
 Triangle *CreateTriangle(Edge *e1, Edge *e2, Edge *e3);


 //! Creates an arbitrarily oriented Triangle bounded by three existing mesh edges.

 //! Returns the newly created Triangle. If either e1, e2 or e3
 //! has already two incident triangles, the creation fails and NULL is returned.
 //! This method assumes that e1, e2 and e3 are incident to exactly three vertices.
 Triangle *CreateUnorientedTriangle(Edge *, Edge *, Edge *);


 //! Creates a new Edge 'e' and an oriented Triangle bounded by 'e', 'e1' and 'e2'.

 //! The newly created triangle is returned, unless 'e1' and 'e2' do not share a
 //! vertex or they are not boundary edges. In this cases, NULL is returned.
 Triangle *EulerEdgeTriangle(Edge *e1, Edge *e2);

 /////////////////////////////////////////////////////////////////////////////
 //
 // Primitive Destruction (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////


 //! Unlinks a triangle from the mesh. O(1).

 //! Resulting isolated vertices and edges are unlinked too.
 //! If necessary, this method duplicates non-manifold vertices that can
 //! occur due to the removal of the triangle.
 //! The unlinked triangle, along with the other possible unlinked elements,
 //! must be removed from the List T through removeUnlinkedElements().
 void unlinkTriangle(Triangle *);


 //! Unlinks a triangle from the mesh. O(1).

 //! No check is performed on the resulting topology, which may be inconsistent.
 //! The unlinked triangle, along with the other possible unlinked elements,
 //! must be removed from the List T through removeUnlinkedElements().
 void unlinkTriangleNoManifold(Triangle *);


 //! Removes a triangle from the mesh. O(N).

 //! This is equivalent to an unlinkTriangle(t) followed by a
 //! removeUnlinkedElements().
 void removeTriangle(Triangle *t) {unlinkTriangle(t); removeUnlinkedElements();}

 //! Removes all the unlinked triangles from List T. Returns the number of removed triangles. O(N).
 int removeTriangles();

 //! Removes all the unlinked edges from List E. Returns the number of removed edges. O(N).
 int removeEdges();

 //! Removes all the unlinked vertices from List V. Returns the number of removed vertices. O(N).
 int removeVertices();

 //! Removes all the unlinked elements from the lists. Returns the number of removed elements. O(N).
 int removeUnlinkedElements() {return removeTriangles()+removeEdges()+removeVertices();}


 /////////////////////////////////////////////////////////////////////////////
 //
 // Methods acting on selections (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Deselects all the triangles. O(N).
 void deselectTriangles();

 //! Removes all the selected triangles. O(N).
 void removeSelectedTriangles();

 //! Selects all the triangles having at least one boundary vertex. O(N).
 //! Returns the number of selected triangles.
 int selectBoundaryTriangles();

 //! Enlarges the current selection of one triangle in width. O(N).

 //! Each triangle sharing at least one vertex with a currently selected
 //! triangle becomes selected.
 //! Returns the number of newly selected triangles.
 int growSelection();

 //! Shrinks the current selection of one triangle in width. O(N).

 //! Each triangle sharing at least one vertex with a currently unselected
 //! triangle becomes unselected.
 void shrinkSelection();

 //! Inverts the selection status of all the triangles. O(N).

 //! If 't0' is not NULL, then only the connected component containing 't0'
 //! is inverted.
 void invertSelection(Triangle *t0 =NULL);

 //! If 't0' is selected, deselect everything but the selected triangles connected to 't0'
 void reselectSelection(Triangle *t0);

 //! Creates a new Triangulation out of an existing selection containing 't0'. O(output).

 //! If necessary, non-manifold vertices are properly duplicated.
 //! If 'keep_ref' is set to TRUE, then elements of the original mesh point
 //! (through their info field) to corresponding elements of the newly created copy.

 Triangulation *createSubMeshFromSelection(Triangle *t0 = NULL, bool keep_ref = 0);

 //! Marks all the triangles within distance L from 'p' as selected. O(output).

 //! A triangle is considered to be within distance L from 'p' only if all
 //! its three vertices are so.
 //! Point 'p' is assumed to belong to triangle 't0', which is required to
 //! limit the complexity to the size of the selected region.
 //! Returns the number of selected triangles.
 int  selectSphericalRegion(Triangle *t0, const double L, const Point *p);


 //! Marks all the triangles within distance L from 'p' as deselected. O(output).

 //! A triangle is considered to be within distance L from 'p' only if all
 //! its three vertices are so.
 //! Point 'p' is assumed to belong to triangle 't0', which is required to
 //! limit the complexity to the size of the selected region.
 //! Returns the number of deselected triangles.
 int  deselectSphericalRegion(Triangle *t0, const double L, const Point *p);


 //! Deselects all the triangles farther than L from 'p'. O(N).

 //! A triangle is considered to be farther than L from 'p' if at least
 //! one of its three vertices is so.
 //! Point 'p' is assumed to belong to triangle 't0'. Passing 't0' avoids
 //! the non robust and expensive computation of point-in-triangle.
 void reselectSphericalRegion(Triangle *t0, const double L, const Point *p);

 //! Re-triangulates the currently selected region using a Delaunay-like approach. O(SlogS).

 //! A common plane is computed as the average of the planes of the triangles selected;
 //! then, the vertices of the region are projected on the plane and edges are iteratively
 //! swapped up to convergence (which is guaranteed on planar and simple domains).
 //! Finally, the vertices are moved back to their original positions. This operation is
 //! particularly useful to improve the quality of nearly flat regions. The selection must
 //! be simple and its Gauss map must be entirely contained in a semi-sphere.
 //! Returns TRUE on success, FALSE otherwise.
 bool retriangulateSelectedRegion();


 //! TRUE iff the set of selected triangles in 'l' is simply connected. O(l->numels()).
 bool isSelectionSimple(List *l);

 //! Unmarks all the elements but leaves the selection status of triangles as is. O(N).
 void unmarkEverythingButSelections();


 //! Selects all the triangles of the connected component containing t0. O(N).

 //! If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
 //! Returns the number of selected triangles.
 int selectConnectedComponent(Triangle *t0, bool stop_on_sharp=0);


 //! Deselects all the triangles of the connected component containing t0. O(N).

 //! If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
 //! Returns the number of deselected triangles.
 int deselectConnectedComponent(Triangle *t0, bool stop_on_sharp=0);

 //! Append to the current mesh a copy of all the elements of 't'.

 //! The newly created elements form a new selection.
 void append(Triangulation *t);

 /////////////////////////////////////////////////////////////////////////////
 //
 // Region manipulation (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Make a list of triangles within distance L from 'p'. O(output).

 //! Starting from 't0', which is assumed to contain 'p', add a triangle at a
 //! time to the list as long as all the vertices stay within distance L from 'p'.
 List     *getRegion(Triangle *t0, const double L, const Point *p);

 //! Removes triangles within distance L from 'p'. O(N).

 //! Starting from 't0', which is assumed to contain 'p', remove a triangle at a
 //! time as long as all its vertices stay within distance L from 'p'.
 void      removeRegion(Triangle *t0, const double L, const Point *p);

 //! Get the vertex next to 'v' on the boundary of the region. O(1).
 Vertex   *nextVertexOnRegionBoundary(Vertex *v) const;

 //! Retrieve internal vertices of a region. O(l->numels()).

 //! This method returns a list containing an edge of the region's boundary
 //! as its first element, and all the internal vertices as the remaining elements.
 List     *getRegionInternalVertices(List *l);

 //! Transform the vertices of the shell containing 't0' using the matrix m. O(S).
 void      transformShell(Triangle *t0, const Matrix4x4& m);

 //! Remove all the triangles belonging to the shell containing 't0'. O(N).
 void 	   removeShell(Triangle *t0);


 /////////////////////////////////////////////////////////////////////////////
 //
 // Global Operations (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Tag sharp edges based on threshold angle. O(N).

 //! Tag as sharp all the edges in which the normals of the two incident
 //! triangles form an angle greater than 't'.
 void   sharpEdgeTagging(const double t);

 //! Unmark all the elements. O(N).
 void   unmarkEverything();

 //! Bounding box longest edge. 'b' and 't' are set as the longest diagonal end-points. O(N).
 double getBoundingBox(Point& b, Point& t) const;

 //! Bounding box longest diagonal. O(N).
 double bboxLongestDiagonal() {Point a, b; getBoundingBox(a, b); return a.distance(b);}

 //! Approximate bounding ball radius. O(N).
 double getBoundingBallRadius() const;

 //! Total area of the mesh. O(N).
 double area() const;

 //! Total volume of the mesh assuming that boundaries() = 0. O(N).
 double volume() const;

 //! Scale the mesh to make it fit within a cube [0,0,0]-[s,s,s]. O(N).
 void   normalize(const double s =1.0);

 //! Transform the mesh geometry using the transformation matrix m. O(N).
 void   transform(const Matrix4x4& m);

 //! Randomly move vertices along their normals. O(N).

 //! Displacement is bounded by 'p'% of the bounding ball radius.
 void   addNormalNoise(const double p);

 //! Iteratively swaps edges to minimize the Delaunay minimum angle. O(N).

 //! Edges tagged as sharp are constrained not to swap.
 //! On generically curved manifolds this process is not guaranteed to converge.
 //! This method returns TRUE if convergence is reached, FALSE otherwise.
 bool   iterativeEdgeSwaps();


 /////////////////////////////////////////////////////////////////////////////
 //
 // Surface topology manipulation (Implemented in "MESH_STRUCTURE/tin.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 //! Invert all the triangle and edge orientations. O(N).
 void      flipNormals();

 //! Invert the orientation of triangles and edges belonging to the shell containing 't0'. O(S).
 void      flipNormals(Triangle *t0);

 //! Return the triangle with the maximum 'z' coordinate in the shell containing 't0'. O(N).

 //! Useful for orienting meshes bounding solids.
 Triangle *topTriangle(Triangle *t0);

 //! Updates the values of n_boundaries, n_handles and n_shells. O(N).

 //! The relative dirty bits are set to zero.
 void      eulerUpdate();

 //! Duplicates edges and vertices to make the mesh homeomorphic to a disk. O(N).
 void      openToDisk();


 /////////////////////////////////////////////////////////////////////////////
 //
 // Handling wrong topology (Implemented in "MESH_STRUCTURE/checkAndRepair.C")
 //
 /////////////////////////////////////////////////////////////////////////////

 int       removeSmallestComponents();		// (in "MESH_STRUCTURE/checkAndRepair.C")
 int       forceNormalConsistence();		// (in "MESH_STRUCTURE/checkAndRepair.C")
 int       forceNormalConsistence(Triangle *);	// (in "MESH_STRUCTURE/checkAndRepair.C")
 int       duplicateNonManifoldVertices();	// (in "MESH_STRUCTURE/checkAndRepair.C")
 const char *checkConnectivity();			// (in "MESH_STRUCTURE/checkAndRepair.C")
 int       removeDuplicatedTriangles();		// (in "MESH_STRUCTURE/checkAndRepair.C")
 int       checkAndRepair();			// (in "MESH_STRUCTURE/checkAndRepair.C")


 // Degenerate geometry manipulation (Implemented in "MESH_STRUCTURE/checkAndRepair.C")

 int     mergeCoincidentEdges();
 int     removeDegenerateTriangles();
 int     removeOverlappingTriangles();
 int     selectTinyHandles(double);
 Vertex *checkGeometry();


 // Triangulation methods

 int     StarTriangulateHole(Edge *);
 int     TriangulateHole(Edge *, Point *);
 int     TriangulateHole(Edge *, List *);
 Vertex *watsonInsert(Point *, List *, int);
 int     retriangulateVT(Vertex *);
 Vertex *splitEdge(Edge *, Point *, bool =0);
 Vertex *splitTriangle(Triangle *, Point *, bool =0);

 // Debug and work-in-progress

 void printReport();
};

#define FOREACHTRIANGLE(Tt, n) for (n = T.head(), Tt = (n)?((Triangle *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Triangle *)n->data):NULL)
#define FOREACHEDGE(Tt, n) for (n = E.head(), Tt = (n)?((Edge *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Edge *)n->data):NULL)
#define FOREACHVERTEX(Tt, n) for (n = V.head(), Tt = (n)?((Vertex *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Vertex *)n->data):NULL)

#define MARK_VISIT(a)   ((a)->mask |= ((unsigned char)1))
#define IS_VISITED(a)   ((a)->mask &  ((unsigned char)1))
#define UNMARK_VISIT(a) ((a)->mask &= (~((unsigned char)1)))

#define MARK_VISIT2(a)   ((a)->mask |= ((unsigned char)2))
#define IS_VISITED2(a)   ((a)->mask &  ((unsigned char)2))
#define UNMARK_VISIT2(a) ((a)->mask &= (~((unsigned char)2)))

#define MARK_BIT(a,b)   ((a)->mask |= ((unsigned char)(1<<b)))
#define IS_BIT(a,b)     ((a)->mask &  ((unsigned char)(1<<b)))
#define UNMARK_BIT(a,b) ((a)->mask &= (~((unsigned char)(1<<b))))

#define TAG_SHARPEDGE(a)   (MARK_BIT((a),7))
#define IS_SHARPEDGE(a)    (IS_BIT((a),7))
#define UNTAG_SHARPEDGE(a) (UNMARK_BIT((a),7))


// Errors from loading

#define IO_CANTOPEN	10
#define IO_FORMAT	20
#define IO_UNKNOWN	30

#define IO_CSAVE_OVERALL		0
#define IO_CSAVE_PERFACE		1
#define IO_CSAVE_PERVERTEX		2
#define IO_CSAVE_PERFACE_INDEXED	3
#define IO_CSAVE_PERVERTEX_INDEXED	4

#endif //_TIN_H

