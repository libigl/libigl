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

#include "point.h"
#include <stdlib.h>
#include <limits.h>
#include <errno.h>


const Point INFINITE_POINT(DBL_MAX, DBL_MAX, DBL_MAX);


//////// Lexicographic Point comparison for qsort //////////

int xyzCompare(const void *a, const void *b)
{
 coord c;

 if ((c=(((Point *)a)->x - ((Point *)b)->x)) < 0) return -1;
 if (c > 0) return 1;
 if ((c=(((Point *)a)->y - ((Point *)b)->y)) < 0) return -1;
 if (c > 0) return 1;
 if ((c=(((Point *)a)->z - ((Point *)b)->z)) < 0) return -1;
 if (c > 0) return 1;

 return 0;
}


//////////////// Normalization /////////////////////////

void Point::normalize()
{
 double l = length();

 if (l == 0) JMesh::error("normalize : Trying to normalize a null vector !\n");

 x/=l;
 y/=l;
 z/=l;
}


//////////////////// Point rotation ////////////////////
/////////// 'ang' radians CCW around 'axis' ////////////

void Point::rotate(const Point& a, const double& ang)
{
 double l, q[4], m[3][3];
 if ((l = a.length())==0.0) return;
 l = sin(ang/2.0)/l;

 q[0] = a.x*l;
 q[1] = a.y*l;
 q[2] = a.z*l;
 q[3] = cos(ang/2.0);

 m[0][0] = 1 - 2.0*(q[1]*q[1] + q[2]*q[2]);
 m[0][1] =     2.0*(q[0]*q[1] + q[2]*q[3]);
 m[0][2] =     2.0*(q[2]*q[0] - q[1]*q[3]);

 m[1][0] =     2.0*(q[0]*q[1] - q[2]*q[3]);
 m[1][1] = 1 - 2.0*(q[2]*q[2] + q[0]*q[0]);
 m[1][2] =     2.0*(q[1]*q[2] + q[0]*q[3]);

 m[2][0] =     2.0*(q[2]*q[0] + q[1]*q[3]);
 m[2][1] =     2.0*(q[1]*q[2] - q[0]*q[3]);
 m[2][2] = 1 - 2.0*(q[1]*q[1] + q[0]*q[0]);

 q[0] = x; q[1] = y; q[2] = z;
 x = m[0][0]*q[0] + m[1][0]*q[1] + m[2][0]*q[2];
 y = m[0][1]*q[0] + m[1][1]*q[1] + m[2][1]*q[2];
 z = m[0][2]*q[0] + m[1][2]*q[1] + m[2][2]*q[2];
}


///// Project the point on the plane whose normal is 'nor' /////

void Point::project(const Point *nor)
{
 Point pr = (*this)-((*nor)*((*this)*(*nor)));
 x = pr.x; y = pr.y; z = pr.z;
}


////////////// Alignment check /////////////

bool Point::notAligned(const Point *A, const Point *B) const
{
 Point p1 = ((*this)-(*A)), p2 = ((*this)-(*B));
 if (p1.length()*p2.length() == 0.0) return 0;
 double a = p1.getAngle(p2);

 return (a != 0.0 && a != M_PI);
}


/////////// Distance from the line passing through A and B ////////

double Point::distanceFromLine(const Point *A, const Point *B) const
{
 Point BA = (*B)-(*A);
 double lba = BA.length();

 if (lba == 0.0) JMesh::error("distanceFromLine : Degenerate line passed !\n");

 return ((((*this)-(*A))&BA).length())/(lba);
}


/////////////////// Distance from a line ///////////////////////
//// 'cc' is initialized as the point of the line whose     ////
//// distance from 'this' is minimum.                       ////

double Point::distanceFromLine(const Point *A, const Point *B, Point *cc) const
{
 Point AB = (*A)-(*B);
 Point AP = (*A)-(*this);
 Point BP = (*B)-(*this);

 if (AP.isNull())
 {
  cc->x = A->x; cc->y = A->y; cc->z = A->z; 
  return 0.0;
 }
 else if (BP.isNull())
 {
  cc->x = B->x; cc->y = B->y; cc->z = B->z; 
  return 0.0;
 }

 double t = (AB*AB);
 if (t == 0.0) JMesh::error("distanceFromLine : Degenerate line passed !\n");
 else t = (AP*AB)/(-t);
 cc->x = t*AB.x + A->x;
 cc->y = t*AB.y + A->y;
 cc->z = t*AB.z + A->z;
 return distanceFromLine(A,B);
}


////////////// Projection on the line passing through A and B ///////////

Point Point::projection(const Point *A, const Point *B) const
{
 Point BA = (*B)-(*A);
 double l = BA*BA;
 if (l == 0.0) JMesh::error("projection : Degenerate line passed !\n");

 return ((*A)+(BA*((BA*((*this)-(*A)))/(l))));
}


////////////// Distance from a segment /////////////////

double Point::distanceFromEdge(const Point *A, const Point *B) const
{
 Point AP = (*A)-(*this); double apl = AP.length();
 Point BP = (*B)-(*this); double bpl = BP.length();

 if (apl == 0 || bpl == 0.0) return 0.0;

 Point AB = (*A)-(*B); double abl = AP.length();
 Point BA = (*B)-(*A);

 if (abl*apl == 0.0 || abl*bpl == 0.0) return apl;

 if (AB.getAngle(AP) > PI2) return apl;
 else if (BA.getAngle(BP) > PI2) return bpl;

 return distanceFromLine(A,B);
}

/////////////////// Distance from a segment ///////////////////////
//// 'cc' is initialized as the point of the segment whose     ////
//// distance from 'this' is minimum.                          ////

double Point::distanceFromEdge(const Point *A, const Point *B, Point *cc) const
{
 Point AP = (*A)-(*this); double apl = AP.length();
 Point BP = (*B)-(*this); double bpl = BP.length();

 if (apl == 0) {cc->setValue(A); return 0.0;}
 if (bpl == 0) {cc->setValue(B); return 0.0;}

 Point AB = (*A)-(*B); double abl = AP.length();
 Point BA = (*B)-(*A);

 if (abl*apl == 0.0 || abl*bpl == 0.0) {cc->setValue(A); return apl;}

 if (AB.getAngle(AP) > PI2) {cc->setValue(A); return apl;}
 else if (BA.getAngle(BP) > PI2) {cc->setValue(B); return bpl;}
 
 double t = (AB*AB);
 if (t == 0.0) {cc->setValue(A); return apl;}
 else t = (AP*AB)/(-t);
 cc->x = t*AB.x + A->x;
 cc->y = t*AB.y + A->y;
 cc->z = t*AB.z + A->z;
 return distanceFromLine(A,B);
}

///////////////// Angle between two vectors ///////////////

double Point::getAngle(const Point& p) const
{
 double ac, l = (length()*(p.length()));

 if (l == 0.0) {JMesh::warning("getAngle : One or both vectors are null !\n"); return -1.0;}

 ac = ((*this)*p)/l;
 if ((FABS((ac - 1.0))) < JMesh::acos_tolerance) ac = 0.0;
 else if ((FABS((1.0 + ac))) < JMesh::acos_tolerance) ac = M_PI;
 else ac = acos(ac);

 if (errno == EDOM)
 {
  JMesh::warning("Point::getAngle(): Warning. acos domain overflow !\n");
  errno = 0;
  return -1.0;
 }

 return ac;
}


//// Side: returns 1 if right, -1 if left, 0 aligned         ////
//// respect to the plane whose directional vector is "this" ////

int Point::side3D(const Point *p1, const Point *p2, const Point *p3) const
{
 Point A = ((*p1)-(*p2)), B = ((*p2)-(*p3));
 double l1 = A.length(), l2 = B.length();
 if (l1*l2 == 0.0) return 0;
 double ang = A.getAngle(B);
 if (ang == 0.0 || ang == M_PI) return 0;

 double db = ((A&B)*(*this));
 return (db > 0)?(1):(-1);
}
 
/////////// Distance of two straight lines ///////////////

double Point::distanceLineLine(const Point *A, const Point *A1, const Point *B1) const
{
 Point uu1 = ((*this)-(*A))&((*A1)-(*B1));
 double nom = ((*A)-(*A1))*(uu1);
 return FABS(nom)/(uu1.length());
}


/////////// Solution of a linear system 3 x 3    //////////
///// System Ax = d, where A = (a,b,c) rows, d = this /////

Point Point::linearSystem(const Point& a, const Point& b, const Point& c)
{
 Point ret;
 double det_A = a*(b&c);
 if (det_A == 0.0) return INFINITE_POINT;
 ret.x = (Point(x,a.y,a.z))*((Point(y,b.y,b.z))&(Point(z,c.y,c.z)));
 ret.y = (Point(a.x,x,a.z))*((Point(b.x,y,b.z))&(Point(c.x,z,c.z)));
 ret.z = (Point(a.x,a.y,x))*((Point(b.x,b.y,y))&(Point(c.x,c.y,z)));

 return (ret/det_A);
}


/// Initializes the point with the intersection between the edge AB and
/// the plane ax+by+cz+d=0. If there is no intersection returns 0. If the edge
/// lies entirely on the plane returns 2.

int Point::intersectionWithPlane(const Point *A, const Point *B, const double& a, const double& b, const double& c, const double& d)
{
 Point abc = Point(a,b,c);
 double pa = (abc*(*A))+d;
 double pb = (abc*(*B))+d;

 if (pa == 0) {x = A->x; y = A->y; z = A->z; return (pb==0)?(2):(1);}
 if (pb == 0) {x = B->x; y = B->y; z = B->z; return 1;}
 if (pa*pb > 0) return 0;
 pa = FABS(pa);
 pb = FABS(pb);

 Point o = ((*B)*pa + (*A)*pb)/(pa+pb);

 x = o.x; y = o.y; z = o.z; 
 return 1;
}

//// As the previous one, the plane is (direction, starting point) ////

int Point::intersectionWithPlane(const Point *A, const Point *B, const Point *nor, const Point *s)
{
 return intersectionWithPlane(A, B, nor->x, nor->y, nor->z, -((*nor)*(*s)));
}

//// Computes the closest points of the two lines 'this'-v1 and p1-p2  ////
//// Returns FALSE if the lines are parallel.                          ////

int Point::closestPoints(const Point *v1, const Point *p1, const Point *p2, Point *ptOnThis, Point *ptOnLine2) const
{
 Point pos1 = *this; Point dir1 = (*v1)-pos1;
 Point pos2 = *p1;   Point dir2 = (*p2)-pos2;
 double d1l = dir1.length(), d2l = dir2.length();

 if (d1l == 0.0 && d2l == 0.0)
  {ptOnThis->setValue(this); ptOnLine2->setValue(p1); return 1;}
 if (d1l*d2l == 0.0)
 {
  if (d1l <= d2l)
   {ptOnThis->setValue(this); distanceFromLine(p1, p2, ptOnLine2); return 1;}
  if (d2l <= d1l)
   {ptOnLine2->setValue(p1); p1->distanceFromLine(this, v1, ptOnThis); return 1;}
 }

 double ang = dir1.getAngle(dir2);
 if (ang == 0.0 || ang == M_PI) return 0;

 double s, t, A, B, C, D, E, F, denom;

 denom = ((dir1*dir2)/(d1l*d2l));
 denom = denom*denom - 1;

 dir1.normalize();
 dir2.normalize();

 A = E = dir1*dir2;
 B = dir1*dir1;
 C = (dir1*pos1) - (dir1*pos2);
 D = dir2*dir2;
 F = (dir2*pos1) - (dir2*pos2);

 s = ( C * D - A * F ) / denom;
 t = ( C * E - B * F ) / denom;
 *ptOnThis  = pos1 + (dir1*s);
 *ptOnLine2 = pos2 + (dir2*t);

// Uncomment the following to compute the distance between segments
// if (s < 0 || s > ((*v1)-(*this)).length() || t < 0 || t > ((*p2)-(*p1)).length())
//	return 0;	       // The points does not belong to the edges

 return 1;
}
