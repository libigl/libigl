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

#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include "jmesh.h"

//! Triangle of a Triangulation.

//! This  class represents a triangle of a triangulation. Each Triangle has
//! an orientation (clockwise or counter-clockwise) due  to  the  order  in
//! which  its  edges are stored in the class. When looking at the triangle
//! so that (e1, e2, e3) are sorted counter-clockwise, the  normal  at  the
//! triangle  points  towards  the  observer.  The field mask is useful for
//! assigning up to 256 different states to the edge.


class Triangle
{
 public :
 
 Edge *e1, *e2, *e3; 		//!< Edges of the triangle
 void *info;			//!< Further information
 unsigned char mask;		//!< bit-mask for marking purposes

 Triangle(Edge *, Edge *, Edge *);		//!< Constructor

 bool isLinked() const {return (e1!=NULL);}	//!< TRUE if properly linked

 //! Inverts the orientation of the triangle
 void invert() {p_swap((void **)(&e2), (void **)(&e3));}

 Vertex *v1() const {return e1->commonVertex(e2);}	//!< First vertex
 Vertex *v2() const {return e2->commonVertex(e3);}	//!< Second vertex
 Vertex *v3() const {return e3->commonVertex(e1);}	//!< Third vertex

 //! First adjacent triangle. NULL if boundary.
 Triangle *t1() const {return e1->oppositeTriangle(this);}

 //! Second adjacent triangle. NULL if boundary.
 Triangle *t2() const {return e2->oppositeTriangle(this);}

 //! Third adjacent triangle. NULL if boundary.
 Triangle *t3() const {return e3->oppositeTriangle(this);}

 //! TRUE iff 'e' is an edge of the triangle.
 bool hasEdge(const Edge *e) const {return (e==e1 || e==e2 || e==e3);}

 //! TRUE iff 'v' is a vertex of the triangle.
 bool hasVertex(const Vertex *v) const {return (e1->hasVertex(v) || e2->hasVertex(v) || e3->hasVertex(v));}

 //! Triangle's edge opposite to 'v'. NULL if 'v' is not a vertex of the triangle.
 Edge *oppositeEdge(const Vertex *v) const
  {return ((!e1->hasVertex(v))?(e1):((!e2->hasVertex(v))?(e2):((!e3->hasVertex(v))?(e3):(NULL))));}

 //! Adjacent triangle opposite to 'v'. NULL if 'v' is not a vertex of the triangle.
 Triangle *oppositeTriangle(const Vertex *v) const
  {return ((!e1->hasVertex(v))?(t1()):((!e2->hasVertex(v))?(t2()):((!e3->hasVertex(v))?(t3()):(NULL))));}

 //! Triangle's vertex opposite to 'e'. NULL if 'e' is not an edge of the triangle.
 Vertex *oppositeVertex(const Edge *e) const
  {return (e==e1)?(v2()):((e==e2)?(v3()):((e==e3)?(v1()):(NULL)));}

 //! Triangle adjacent to the next edge of 'e'. NULL if 'e' is not an edge of the triangle.
 Triangle *rightTriangle(const Edge *e) const
  {return (e==e1)?(t2()):((e==e2)?(t3()):((e==e3)?(t1()):(NULL)));}

 //! Triangle adjacent to the previous edge of 'e'. NULL if 'e' is not an edge of the triangle.
 Triangle *leftTriangle(const Edge *e) const
  {return (e==e1)?(t3()):((e==e2)?(t1()):((e==e3)?(t2()):(NULL)));}

 //! Edge next to 'e' in the ordering or the triangle. NULL if 'e' is not an edge of the triangle.
 Edge *nextEdge(const Edge *e) const {return ((e==e1)?(e2):((e==e2)?(e3):((e==e3)?(e1):(NULL))));}

 //! Edge preceeding 'e' in the ordering or the triangle. NULL if 'e' is not an edge of the triangle.
 Edge *prevEdge(const Edge *e) const {return ((e==e1)?(e3):((e==e2)?(e1):((e==e3)?(e2):(NULL))));}

 //! Vertex next to 'v' in the ordering or the triangle. NULL if 'v' is not a vertex of the triangle.
 Vertex *nextVertex(const Vertex *v) const
  {return (!e1->hasVertex(v))?(v3()):((!e2->hasVertex(v))?(v1()):((!e3->hasVertex(v))?(v2()):(NULL)));}

 //! Vertex preceeding 'v' in the ordering or the triangle. NULL if 'v' is not a vertex of the triangle.
 Vertex *prevVertex(const Vertex *v) const
  {return (!e1->hasVertex(v))?(v1()):((!e2->hasVertex(v))?(v2()):((!e3->hasVertex(v))?(v3()):(NULL)));}

 //! If this triangle shares an edge with 'b', then such an edge is returned. NULL otherwise.
 Edge *commonEdge(const Triangle *b) const
	{return ((e1 == b->e1 || e1 == b->e2 || e1 == b->e3)?(e1):\
		(((e2 == b->e1 || e2 == b->e2 || e2 == b->e3)?(e2):\
		(((e3 == b->e1 || e3 == b->e2 || e3 == b->e3)?(e3):(NULL))))));}

 //! Replace edge 'a' with edge 'b' in the triangle and return TRUE. If 'a' is not an edge of the triangle return FALSE.
 bool replaceEdge(const Edge *a, Edge *b)
  {if (e1==a) e1=b; else if (e2==a) e2=b; else if (e3==a) e3=b; else return 0; return 1;}

 //! Return a normal vector with direction (v1-v2) cross (v2-v3). If triangle is degenerate return a null vector.
 Point  getNormal() const;

 //! Return a vector orthogonal to the plane of the triangle. If triangle is degenerate return a null vector.
 Point  getVector() const;

 //! TRUE if the oriantation is consistent with the one of 't' OR if this and 't' do not share any edge. 
 bool checkAdjNor(const Triangle *t) const;

 //! Area of the triangle (Heron's formula). 
 double area() const;

 //! Perimeter of the triangle. 
 double perimeter() const;

 //! Return the triangle's barycenter.
 Point  getCenter() const;

 //! Return the center of the triangle's bounding sphere.
 Point  getCircleCenter() const;

 //! TRUE iff 'p' is inside the triangle's bounding sphere.
 bool inSphere(const Point *p) const;

 //! Angle at vertex 'v'. Return -1 if 'v' is not a vertex of the triangle.
 double getAngle(const Vertex *v) const;

 //! Angle between the normal vector of this and the one of 't'. Return -1 if one or both the triangles are degenerate.
 double getDAngle(const Triangle *t) const;

 //! Distance of 'p' from the plane of the triangle. Return -1 if triangle is degenerate.
 double distanceFromPoint(const Point *p) const;

 //! Squared distance of 'p' from the plane of the triangle. Return -1 if triangle is degenerate.
 double squaredDistanceFromPoint(const Point *p) const;

 //! Distance of 'p' from the closest point of the triangle. Return -1 if triangle is degenerate.
 //! If 'c' is not NULL, its coordinates are set to the ones of the closest point.
 double pointTriangleDistance(const Point *p, Point *c =NULL) const;

 //! Squared distance of 'p' from the closest point of the triangle. Return -1 if triangle is degenerate.
 double pointTriangleSquaredDistance(const Point *p) const;

 //! Projection of 'p' on the plane of the triangle. Return INFINITE_POINT if triangle is degenerate.
 Point project(const Point *p) const;

 //! TRUE iff the projection of 'p' on the plane of the triangle is inside the triangle.
 bool isInside(const Point *p) const;

 //! Triangle-edge intersection.

 //! Checks whether the triangle intersects the edge 'e'. If it does,
 //! the method returns TRUE and 'p' is initialized with the point of
 //! intersection. If there is no intersection or if the edge lies
 //! on the plane of the triangle, FALSE is returned. 't' is a threshold
 //! distance within which a point is considered to be on the face.
 bool intersectsEdge(const Edge *e, Point *p, const double =0.0) const;

 //! Cap-like degeneracy check.

 //! If the angle at a vertex is flat, return the edge opposite to such vertex.
 //! Return NULL otherwise.
 Edge *isCap() const;

 //! Needle-like degeneracy check.

 //! If an edge of the triangle has zero length, return that edge.
 //! If the angle at a vertex is zero, return the edge opposite to such vertex.
 //! Return NULL otherwise.
 Edge *isNeedle() const;

 //! Degeneracy check. Return TRUE iff isCap() or isNeedle().
 bool isDegenerate() const;

 //! Return TRUE iff one of the adjacent triangles overlaps with this one.
 bool overlaps() const;

 //! Print the coordinates of the three vertices to the file handler pointed to by 'f' (stdout by default).
 void printTriangle(FILE *f =stdout) const;
};

#endif // _TRIANGLE_H

