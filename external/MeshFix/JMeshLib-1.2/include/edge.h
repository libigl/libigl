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

#ifndef _EDGE_H
#define _EDGE_H

#include "j_mesh.h"
#include "vertex.h"

//! Edge of a Triangulation.

//! This  class  represents an edge of a triangulation. An Edge is the main
//! part of the Triangulation data structure. Each edge has an  orientation
//! (i.e. from v1 to v2) which forces the order in which incident triangles
//! (t1 and t2) are stored in the class. When looking the edge so  that  it
//! points  "upwards", if the normal of t1 points towards the observer then
//! t1 must be on the left of the  edge.  The  field  mask  is  useful  for
//! assigning up to 256 different states to the edge.


class Edge
{
 public :

 Vertex *v1,*v2;		//!< End-points
 class Triangle *t1,*t2;	//!< Incident triangles
 unsigned char mask;		//!< bit-mask for marking purposes
 void *info;			//!< Further information

 Edge(Vertex *s, Vertex *d);	//!< Constructor
 ~Edge();			//!< Destructor

 //! TRUE iff edge is properly linked to a Triangulation.
 bool isLinked()	const 	{return (v1 != NULL);}

 //! TRUE iff 'v' is an end-point of the edge.
 bool hasVertex(const Vertex *v) const {return (v1==v || v2==v);}

 //! TRUE iff 't' is incident to the edge.
 bool hasTriangle(const Triangle *t) const {return (t1==t || t2==t);}

 //! TRUE if both 'va' and 'vb' are vertices of the edge.
 bool hasVertices(const Vertex *va, const Vertex *vb) const {return ((v1==va && v2==vb) || (v2==va && v1==vb));}

 //! Euclidean length of the edge.
 double length() const 	{return v1->distance(v2);}

 //! Squared length of the edge.
 double squaredLength() const 	{return v1->squaredDistance(v2);}

 //! Convert to vector v2-v1.
 Point toVector() const 	{return (*v2)-(*v1);}

 //! Convert to normalized vector (v2-v1)/|v2-v1|.
 Point toUnitVector() const;

 //! Return the edge's mid-point.
 Point getMidPoint()	const {return ((*v1)+(*v2))/2.0;}

 //! Invert the edge's orientation.
 void invert() {p_swap((void **)(&v1), (void **)(&v2)); p_swap((void **)(&t1), (void **)(&t2));}

 //! Triangle on the left of the edge when looking from 'v'. NULL if 'v' is not a vertex of the edge.
 Triangle *leftTriangle(const Vertex *v) const {return ((v1 == v)?(t1):((v2 == v)?(t2):(NULL)));}

 //! Triangle on the right of the edge when looking from 'v'. NULL if 'v' is not a vertex of the edge.
 Triangle *rightTriangle(const Vertex *v) const {return ((v1 == v)?(t2):((v2 == v)?(t1):(NULL)));}

 //! Vertex opposite to 'v'. NULL if 'v' is not a vertex of the edge.
 Vertex *oppositeVertex(const Vertex *v) const {return ((v1 == v)?(v2):((v2 == v)?(v1):(NULL)));}

 //! Incident triangle opposite to 't'. NULL if 't' is not incident to the edge.
 Triangle *oppositeTriangle(const Triangle *t) const {return ((t1 == t)?(t2):((t2 == t)?(t1):(NULL)));}

 //! Replace vertex 'a' with vertex 'b' in the edge and return TRUE. If 'a' is not a vertex of the edge return FALSE.
 bool replaceVertex(const Vertex *a, Vertex *b) {if (v1==a) v1=b; else if (v2==a) v2=b; else return 0; return 1;}

 //! Replace incident triangle 'a' with 'b' and return TRUE. If 'a' is not incident to the edge return FALSE.
 bool replaceTriangle(const Triangle *a, Triangle *b) {if (t1==a) t1=b; else if (t2==a) t2=b; else return 0; return 1;}

 //! Vertex shared with edge 'b'. NULL if this and 'b' do not share any vertex.
 Vertex *commonVertex(const Edge *b) const {return ((v1 == b->v1 || v1 == b->v2)?(v1):((v2 == b->v1 || v2 == b->v2)?(v2):(NULL)));}

 //! TRUE iff edge is on the boundary (i.e., one of the two incident triangles is NULL).
 bool isOnBoundary() const {return (t1 == NULL || t2 == NULL);}

 //! TRUE iff edge is isolated (i.e., both the two incident triangles are NULL).
 bool isIsolated() const {return (t1 == NULL && t2 == NULL);}

 //! If the edge is on boundary return its only incident triangle. NULL otherwise.
 Triangle *getBoundaryTriangle() const {return (t2 == NULL)?(t1):((t1 == NULL)?(t2):(NULL));}

 //! Print the coordinates of the end-ponts to the file handler pointed to by 'f' (stdout by default).
 void printEdge(FILE *f =stdout) const {v1->printPoint(f); v2->printPoint(f);}

 //! Return the normal at the edge as the average of the normals of the two incident triangles.

 //! A null (0,0,0) vector is returned if the edge is on boundary.
 Point getNormal() const;

 //! Combinatorial edge-swap.

 //! Vertices of the edge are replaced with vertices of the two incident triangles which are opposite to this edge.
 //! Connectivity information is updated properly.
 //! If the edge is on boundary or if the edge after the swap already exists return FALSE and do not change anything.
 //! Return TRUE on success.
 //! If 'fast' is set, no topological check is performed.
 bool swap(const bool fast=0);

 //! Edge collapse.

 //! This method collapses the edge and  updates  the  connectivity  of  the
 //! neighboring  elements consistently. The edge will be transformed into a
 //! vertex with the coordinates of 'p'.
 //! This method returns TRUE on success, FALSE otherwise.
 //! Failure occurs when the collapse would produce an invalid connectivity graph.
 //! Caution! If the collapse succeeds the
 //! edge,  its  incident  triangles and the second vertex are unlinked, but
 //! they are still present in the lists of the Triangulation.
 //! The calling function is responsible of removing them from the lists using
 //! the method removeUnlinkedElements().
 bool collapse(const Point& p);

 //! Edge collapse.

 //! This method collapses the edge and  updates  the  connectivity  of  the
 //! neighboring  elements consistently. The edge will be transformed into a
 //! vertex placed at the edge's mid-point.
 //! This method returns TRUE on success, FALSE otherwise.
 //! Failure occurs when the collapse would produce an invalid connectivity graph.
 //! Caution! If the collapse succeeds the
 //! edge,  its  incident  triangles and the second vertex are unlinked, but
 //! they are still present in the lists of the Triangulation.
 //! The calling function is responsible of removing them from the lists using
 //! the method removeUnlinkedElements().
 bool collapse();

 //! Merge with another boundary edge.

 //! If both this and 'e' are boundary edges, the edge 'e' is identified with
 //! this one, and the connectivity of the neighboring elements is updated consistently.
 //! This method returns TRUE on success, FALSE otherwise.
 //! Failure occurs when the merge would produce an invalid connectivity graph (i.e., non orientable).
 //! Caution! If the merge succeeds the edge 'e' and its two  end-points
 //! are  unlinked,  but  they are still present in the lists of the
 //! Triangulation. It's responsibility of the calling  function  to  remove
 //! them from the lists using the method removeUnlinkedElements().
 bool merge(Edge *e);

 //! Stitching primitive.

 //! If there is a copy of this edge incident to one of the end-points,
 //! identify it with this edge, and update the connectivity properly.
 //! This method returns TRUE on success, FALSE otherwise.
 //! Caution! If the stitch succeeds, the duplicated edge
 //! is unlinked, but it is still present in the lists of the
 //! Triangulation. It's responsibility of the calling function to remove
 //! it from the lists using the method removeEdges().
 bool stitch();

 //! Dihedral angle at the edge.
 double dihedralAngle() const;

 //! Angle between the normals of the two incident triangles.

 //! Angle between the normals of the two incident triangles. If
 //! the edge is on boundary or one or both the incident triangles are
 //! degenerate, return -1.
 double curvature() const;

 //! Return the minimum among the six angles of the two incident triangles.2PI if on boundary.
 double delaunayMinAngle() const;
};

//! Edge comparison based on length to be used with jqsort() or abstractHeap.
int edgeCompare(const void *a, const void *b);

//! Lexycographic edge comparison to be used with jqsort() or abstractHeap.
int lexEdgeCompare(const void *a, const void *b);

#endif //_EDGE_H

