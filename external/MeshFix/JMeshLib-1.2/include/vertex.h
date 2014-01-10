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

#ifndef _VERTEX_H
#define _VERTEX_H

#include "j_mesh.h"
#include "list.h"
#include "point.h"

//! Vertex of a Triangulation

//! This class represents a vertex of a manifold and oriented triangulation.
//! The base-class Point
//! describes  geometrical  and  additional  attributes  of the vertex. The
//! field 'e0' is sufficient to retrieve  all  the  neighboring  elements  in
//! optimal  time,  while  the field 'mask' is useful for assigning up to 256
//! different states to the vertex.


class Vertex : public Point
{
 public :
 class Edge *e0;			//!< One of the incident edges
 unsigned char mask;			//!< bit-mask for marking purposes

 //! Creates a new vertex with coordinates (0,0,0).
 Vertex();

 //! Creates a new vertex with the same coordinates (x,y,z).
 Vertex(const coord& x, const coord& y, const coord& z);

 //! Creates a new vertex with the same coordinates as 'p'. The info field is not copied.
 Vertex(const Point *p);

 //! Creates a new vertex with the same coordinates as 'p'. The info field is not copied.
 Vertex(const Point& p);
 ~Vertex();				//!< Destructor

 bool isLinked() const {return (e0!=0);} //!< TRUE iff vertex is not isolated

 //! List of adjacent vertices.

 //! Returns the list of vertices which are linked to this through an  edge.
 //! The  list is counter-clockwise ordered. In the case of an internal ver-
 //! tex the list starts from the opposite vertex of e0. If the vertex is on
 //! the  boundary,  the  list starts from the opposite vertex of the clock-
 //! wise-most boundary edge.
 List *VV() const;

 //! List of incident edges.

 //! Returns the list of edges incident at this vertex. The list is counter-clockwise
 //! ordered.  In  the case of an internal vertex the list starts from 'e0'
 //! If the vertex is on the boundary, the list starts from its clockwise-most
 //! incident boundary edge.
 List *VE() const;
 
 //! List of incident triangles.

 //! Returns  the  list  of triangles incident at this. The list is counter-
 //! clockwise ordered. In the case of an internal vertex  the  list  starts
 //! from  the  triangle  on  the left of e0, when looking from this. If the
 //! vertex is on the boundary, the  list  starts  from  the  clockwise-most
 //! boundary triangle.
 List *VT() const;

 //! Returns the edge connecting this vertex to 'v'. NULL if such an edge does not exist.
 class Edge *getEdge(const Vertex *v) const;
 int valence() const;				//!< Returns the number of incident edges
 int isOnBoundary() const;			//!< TRUE iff vertex is on boundary
 
 //! Returns the edge following this vertex on the boundary.
 
 //! This edge is the counterclockwise-most incident edge.
 //! Returns NULL if this vertex is not on the boundary.
 Edge *nextBoundaryEdge() const;

 //! Returns the edge preceeding this vertex on the boundary.
 
 //! This edge is the clockwise-most incident edge.
 //! Returns NULL if this vertex is not on the boundary.
 Edge *prevBoundaryEdge() const;

 //! Returns the vertex following this one on the boundary.
 
 //! If the vertex is on the boundary, this is equivalent to nextBoundaryEdge()->oppositeVertex(v)
 //! otherwise returns NULL.
 Vertex *nextOnBoundary() const;

 //! Returns the vertex preceeding this one on the boundary.
 
 //! If the vertex is on the boundary, this is equivalent to prevBoundaryEdge()->oppositeVertex(v)
 //! otherwise returns NULL.
 Vertex *prevOnBoundary() const;

 //! Normal at the vertex computed as the sum of incident triangle normals weighted on their incidence angle.
 Point getNormal() const;

 //! Returns  the angle between the two incident boundary edges. If the vertex is not on the boundary, returns -1.
 double getBoundaryAngle() const;
 
 //! Discriminant Angle for triangulating 3D polygons.
 
 //! This method is useful when patching holes, and represents  a  heuristic
 //! for choosing which vertex of the hole's boundary must be patched first.
 //! Several cases are considered, including degenerate ones. Let 'v1' and 'v2'
 //! be  the two boundary vertices which are linked to this one through a boundary
 //! edge. If 'v1' and 'v2' coincide, the method returns a negative  number.
 //! If  'v1'  ,  'this'  ,  'v2'  form  a flat angle, the method returns 3PI (270
 //! degrees). If the angle formed by 'v1' ,  'this'  ,  'v2'  is  0,  the  method
 //! returns 0. If the vertex is not on the boundary, the method returns the
 //! limit number DBL_MAX. In all the other cases the method returns the sum
 //! of three angles A + D1 + D2, where A is the angle formed by v1 , this ,
 //! v2 , while D1 is the angle between the  normal  of  the  clockwise-most
 //! incident  boundary  triangle and the normal of the triangle v1 , this ,
 //! v2; D2 is the analogous for the counterclockwise-most incident boundary
 //! triangle.
 double getAngleForTriangulation() const;

 //! Discriminant Angle for triangulating flat (or nearly flat) polygons.

 //! This  method  returns the angle between the two incident boundary edges
 //! when projected onto the plane whose normal is 'n'. This  angle
 //! may be more than PI, because it represents the aperture of the non-tri-
 //! angulated region around the vertex when projected on the plane. If  the
 //! vertex  is  not  on  the  boundary,  the method returns the limit value
 //! DBL_MAX.
 double getAngleOnAveragePlane(Point *n) const;
 
 double totalAngle() const;		//!< Sum of incident angles. Returns -1 if on boundary.

 //! Excess angle. Returns DBL_MAX if on boundary.
 double gaussianCurvature() const
  {double t=totalAngle(); return (t>=0)?(t):(DBL_MAX);}

 double totalDihedralAngle() const;	//!< Sum of signed dihedral angles. Returns DBL_MAX if on boundary.
 double voronoiArea() const;		//!< A third of the total area of incident triangles.

 //! Zips the gap starting from here.
 int zip(const bool =1);

 Edge *inverseCollapse(Vertex *, Vertex *, Vertex *);
 Edge *inverseCollapse(Vertex *, Edge *, Edge *, Edge *, Edge *, Edge *, class Triangle *, class Triangle *);
 int  getTopology(const double&) const;
};

//! Scans the nodes 'n' of a list 'l' of vertices 'v'.
#define FOREACHVVVERTEX(l, v, n) for (n = l->head(), v = (n)?((Vertex *)n->data):NULL; n != NULL; n=n->next(), v = (n)?((Vertex *)n->data):NULL)

//! Scans the nodes 'n' of a list 'l' of edges 'e'.
#define FOREACHVEEDGE(l, e, n) for (n = l->head(), e = (n)?((Edge *)n->data):NULL; n != NULL; n=n->next(), e = (n)?((Edge *)n->data):NULL)

//! Scans the nodes 'n' of a list 'l' of triangles 't'.
#define FOREACHVTTRIANGLE(l, t, n) for (n = l->head(), t = (n)?((Triangle *)n->data):NULL; n != NULL; n=n->next(), t = (n)?((Triangle *)n->data):NULL)


//! Extended vertex for temporary use during connectivity creation.

//! This class is used to allow the reconstruction of the connectivity
//! in linear time (average case) and to handle badly oriented input files.
//! It provides a complete VE relation.

class ExtVertex
{
 public :
 Vertex *v;
 List VE;

 ExtVertex(Vertex *a) {v=a;}
};

#endif //_VERTEX_H

