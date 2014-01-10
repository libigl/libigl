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

#ifndef _POINT_H
#define _POINT_H

#include "j_mesh.h"

//! Geometric point definition

//! This class represents a point in the Euclidean 3D space. It can be used
//! to  represent  3D vectors originating at (0,0,0) and terminating at the
//! corresponding point. Several methods of  this  class  are  intended  to
//! manipulate  vectors  rather  than  points;  for  example, a call of the
//! method normalize is an actual normalization if the object is a  vector,
//! but  it  has  to  be intended as a projection on the unit sphere if the
//! object is intended to be a point. An object of type Point is a  triplet
//! (x,y,z)  of  coordinates  endowed with a pointer 'info' to possible additional
//! information. Each coordinate is a number of type 'coord' which, by
//! default,  is  a standard double. Operations on points include addition,
//! subtraction, cross and dot product, and many others. This class  implements
//! several useful operations using vector arithmethic. For example,
//! the simple piece of code "A = B*C;" assignes to A the value of the  dot
//! product of B and C.
//! Nearly zero or nearly flat angles are automatically snapped to
//! exactly zero and exactly flat angles if the difference is smaller
//! than the global variable _acos_tolerance. This is the very basic application
//! of our version of the epsilon geometry for robust computation.


class Point
{
 public :
 coord x,y,z;					//!< Coordinates
 void *info;					//!< Further information

 //! Creates a new point with coordinates (0,0,0).
 Point() {x = y = z = 0; info = NULL;}

 //! Creates a new point with the same coordinates as 's'. The info field is not copied.
 Point(const Point *s) {x = s->x; y = s->y; z = s->z; info = NULL;}

 //! Creates a new point with the same coordinates as 's'. The info field is not copied.
 Point(const Point& s) {x = s.x; y = s.y; z = s.z; info = NULL;}

 //! Creates a new point with coordinates (a,b,c).
 Point(const coord& a, const coord& b, const coord& c) {x = a; y = b; z = c; info = NULL;}

 //! Set the coordinates to (a,b,c).
 void	setValue(const coord& a, const coord& b, const coord& c) {x = a; y = b; z = c;}

 //! Set the coordinates as those of 'p'
 void	setValue(const Point& p) {x = p.x; y = p.y; z = p.z;}

 //! Set the coordinates as those of '*p'
 void	setValue(const Point *p) {x = p->x; y = p->y; z = p->z;}

 //! Returns the vector difference
 Point 	operator-(const Point& p) const {return Point(x-p.x, y-p.y, z-p.z);}

 //! Returns the vector sum
 Point 	operator+(const Point& p) const {return Point(x+p.x, y+p.y, z+p.z);}

 //! Sums another point
 void 	operator+=(const Point& p) {x+=p.x; y+=p.y; z+=p.z;}

 //! Subtracts another point
 void 	operator-=(const Point& p) {x-=p.x; y-=p.y; z-=p.z;}

 //! Returns the Cross Product
 Point 	operator&(const Point& p) const {return Point(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);}

 //! Returns the Dot Product
 double operator*(const Point& p) const {return (x*p.x+y*p.y+z*p.z);}

 //! Returns the product with a scalar
 Point  operator*(const double& d) const {return Point(x*d,y*d,z*d);}

 //! Multiplies by a scalar
 void 	operator*=(const double& m) {x*=m; y*=m; z*=m;}

 //! Divides by a scalar
 void 	operator/=(const double& m) {x/=m; y/=m; z/=m;}

 //! Returns the vector divided by the scalar
 Point 	operator/(const double& d) const {return Point(x/d,y/d,z/d);}

 //! TRUE iff coordinates are equal
 bool  	operator==(const Point& p) const {return (x==p.x && y==p.y && z==p.z);}

 //! FALSE iff coordinates are equal
 bool  	operator!=(const Point& p) const {return (x!=p.x || y!=p.y || z!=p.z);}

 //! Returns the inverse vector
 Point 	inverse() const {return Point(-x,-y,-z);}

 //! Inverts the vector
 void 	invert() {x=-x; y=-y; z=-z;}

 //! TRUE if vector is (0,0,0)
 bool  	isNull() const {return (x==0 && y==0 && z==0);}

 //! Distance from origin
 double length() const {return sqrt(x*x + y*y + z*z);}

 //! Squared distance from origin
 double squaredLength() const {return (x*x + y*y + z*z);}

 //! Divides the vector by its length. If isNull() the application exits with an error.
 void 	normalize();

 //! Rotates the vector around 'axis' by 'ang' radians ccw.
 void  	rotate(const Point& axis, const double& ang);

 //! Projects the vector on the plane with normal 'n' passing through the origin.
 void   project(const Point *n);

 //! TRUE iff 'a', this vector and 'b' are not collinear
 bool 	notAligned(const Point *a, const Point *b) const;

 //! Distance from 'b'
 double distance(const Point& b) const {return (((*(this))-(b)).length());}

 //! Distance from '*b'
 double distance(const Point *b) const {return (((*(this))-(*b)).length());}

 //! Squared distance from '*b'
 double squaredDistance(const Point *b) const {return (((*(this))-(*b)).squaredLength());}

 //! Distance from straight line through 'a' and 'b'
 double distanceFromLine(const Point *a, const Point *b) const;

 //! Distance from straight line through 'a' and 'b'. *cc is set to the closest line point.
 double distanceFromLine(const Point *a, const Point *b, Point *cc) const;

 double distanceFromEdge(const Point *a, const Point *b) const; //!< Distance from segment a-b

 //! Distance from segment a-b. *cc is set to the closest edge point.
 double distanceFromEdge(const Point *a, const Point *b, Point *cc) const;

 //! Distance between the straight lines through (this) - l1_p2 and l2_p1 - l2_p2.
 double distanceLineLine(const Point *l1_p2, const Point *l2_p1, const Point *l2_p2) const;

 //!< Angle between this vector and 'v' in radians.
 double getAngle(const Point& v) const;

 //! Angle defined by <a, *this, b> in radians.
 double getAngle(const Point& a, const Point& b) const {return (a-(*this)).getAngle(b-(*this));}

 //! Angle defined by <*a, *this, *b> in radians.
 double getAngle(const Point *a, const Point *b) const {return ((*a)-(*this)).getAngle((*b)-(*this));}

 //! Returns the solution of the linear system Ax = d, where A is a 3x3 matrix whose rows are row1, row2 and row3, d = this
 Point  linearSystem(const Point& row1, const Point& row2, const Point& row3);

 //! Side test.

 //! When looking from the direction pointed to by this vector, this method returns 1 if the points 'p1',
 //! 'p2' and 'p3' turn right, -1 if they turn left, 0 if they are aligned.
 //! Notice that in this latter case the three point do not need to be linearly dependent.
 int 	side3D(const Point *p1, const Point *p2, const Point *p3) const;

 //! Sets the point as the intersection of a segment and a plane.

 //! Initializes the coordinates with the intersection of the segment  p1-p2
 //! and the plane passing through 'source' with normal 'normal'. If the segment
 //! lies entirely on the plane, this method returns 2 and  the  coordinates
 //! are  initialized  with  those  of  'p1'. If there is no intersection, the
 //! method returns 0 and the coordinates are not  modified.  Otherwise  the
 //! method returns 1.
 int    intersectionWithPlane(const Point *p1, const Point *p2, const Point *source, const Point *normal);

 //! Sets the point as the intersection of a segment and a plane.

 //! Initializes the coordinates with the intersection of the segment  p1-p2
 //! and the plane of equation ax+by+cz+d = 0. If the segment
 //! lies entirely on the plane, this method returns 2 and  the  coordinates
 //! are  initialized  with  those  of  'p1'. If there is no intersection, the
 //! method returns 0 and the coordinates are not  modified.  Otherwise  the
 //! method returns 1.
 int    intersectionWithPlane(const Point *p1, const Point *p2, const double& a, const double& b, const double& c, const double& d);

 //! Line-line closest point computation.

 //! Computes the closest points of the line passing through this and this2,
 //! and the line passing through p1 and p2. The computed points are used to
 //! initialize the  coordinates  of  cpOnThis  and  cpOnOther.  The  method
 //! returns 0 if the lines are parallel, 1 otherwise.
 int    closestPoints(const Point *this2, const Point *p1, const Point *p2, Point *cpOnThis, Point *cpOnOther) const;

 //! Returns the projection of the point on the straight line though 'a' and 'b'.
 Point  projection(const Point *a, const Point *b) const;

 //! Prints the coordinates of the point to a file handler. stdout is the default.
 void 	printPoint(FILE *fp =stdout) const {fprintf(fp,"%f %f %f,\n",x,y,z);}		// Debug
};

//! Lexycographic comparison to be used with jqsort() or abstractHeap.
int xyzCompare(const void *p1, const void *p2);

//! Static point with DBL_MAX coordinates.
extern const Point INFINITE_POINT;

//! Checks whether a point is INFINITE_POINT.
#define IS_FINITE_POINT(p) ((p).x < DBL_MAX)

#endif // _POINT_H

