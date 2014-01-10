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

#include "triangle.h"
#include <stdlib.h>

//////////////////// Constructor //////////////////////

Triangle::Triangle(Edge *a, Edge *b, Edge *c)
{
 e1 = a;
 e2 = b;
 e3 = c;
 mask = 0;
 info = NULL;
}


//////////////////// Normal vector //////////////////////

Point Triangle::getNormal() const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 Point vd = (((*va)-(*vb))&((*vb)-(*vc)));
 double l = vd.length();

 if (l == 0) return Point(0,0,0);

 return vd/l;
}


////// Directional vector (more robust than the above one) ////////

Point Triangle::getVector() const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 Point vd1 = (((*va)-(*vb))&((*vb)-(*vc)));
 Point vd2 = (((*vb)-(*vc))&((*vc)-(*va)));
 Point vd3 = (((*vc)-(*va))&((*va)-(*vb)));
 Point vd = vd1+vd2+vd3;

 if (vd.isNull()) return Point(0,0,0);

 return vd;
}


/////////////////// Normal consistence check ////////////////////

bool Triangle::checkAdjNor(const Triangle *t) const
{
 Edge *e = commonEdge(t);
 if (e == NULL) return 1;

 Edge *ea = nextEdge(e);
 Edge *eb = t->nextEdge(e);
 if (ea->commonVertex(eb) == ea->commonVertex(e)) return 0;

 return 1;
}


//////////////////////// Triangle area /////////////////////////

double Triangle::area() const
{
 double a = e1->length(), b = e2->length(), c = e3->length();
 if (a==0.0 || b==0.0 || c==0.0) return 0.0;
 double p = (a+b+c)/2.0;
 p = p*(p-a)*(p-b)*(p-c); if (p<0) return 0.0;
 return sqrt(p);
}


/////////////////////// Triangle perimeter /////////////////////

double Triangle::perimeter() const
{
 return e1->length()+e2->length()+e3->length();
}


///////////// Barycenter ///////////////////////

Point Triangle::getCenter() const
{
 Point va = *v1(), vb = *v2(), vc = *v3();
 return (va+vb+vc)/3.0; 
}


///////////////////////// Circlecenter /////////////////////////

Point Triangle::getCircleCenter() const
{
 Point va = *v1(), vb = *v2(), vc = *v3();
 Point q1 = vb-va;
 Point q2 = vc-va;
 Point n = q2&q1;
 Point m1 = e2->getMidPoint();
 Point m2 = e1->getMidPoint();

 return Point(n*va,q1*m1,q2*m2).linearSystem(n,q1,q2);
}


/////// Check wether the point is inside the triangle's bounding ball /////

bool Triangle::inSphere(const Point *p) const
{
 Point c = getCircleCenter();
 double rad = c.distance(e1->v1);

 return (p->distance(c) < rad);
}


//////////////////// Angle at a vertex /////////////////////

double Triangle::getAngle(const Vertex *v) const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 if (v == va) return v->getAngle(vb, vc);
 if (v == vb) return v->getAngle(va, vc);
 if (v == vc) return v->getAngle(vb, va);

 return -1.0;
}


/////////// Angle between the two directional vectors /////////

double Triangle::getDAngle(const Triangle *t) const
{
 Point thisNormal = getVector();
 Point otherNormal = t->getVector();

 if (thisNormal.isNull() || otherNormal.isNull()) return -1.0;

 return thisNormal.getAngle(otherNormal);
}


///////////// Distance from the plane of the triangle //////////////

double Triangle::distanceFromPoint(const Point *p) const 
{
 Point n = getNormal();

 if (n.isNull()) return -1.0; 
 double d = (n*(*p))-(n*(*(e1->v1)));

 return FABS(d);
}

///////////// Squared distance from the plane of the triangle //////////////

double Triangle::squaredDistanceFromPoint(const Point *p) const 
{
 Point CA = e1->toVector()&e2->toVector();
 double CA2 = CA*CA;

 if (CA2 == 0) return -1.0; 
 double d = ((CA*(*p))-(CA*(*(e1->v1))));

 return (d*d)/CA2;
}


///////////// Distance of point from the triangle //////////////

double Triangle::pointTriangleDistance(const Point *p, Point *cp) const
{
 Point n = getNormal();
 if (n.isNull()) return -1.0;
 Vertex *va = v1(), *vb = v2(), *vc = v3();

 Point p1 = e1->toUnitVector();
 Point p2 = e2->toUnitVector();
 double d1 = (p1*(*p));
 double d2 = (p2*(*p));
 double d3 = (n*(*va));
 Point i = Point(d1,d2,d3).linearSystem(p1,p2,n);

 d1 = ((((*va)-(*vb))&((*vb)-i))*n);
 d2 = ((((*vb)-(*vc))&((*vc)-i))*n);
 d3 = ((((*vc)-(*va))&((*va)-i))*n);

 if (d1 >= 0 && d2 >= 0 && d3 >= 0) {if (cp) cp->setValue(i); return i.distance(p);}

 if (d2 < 0) {va=vb; vb=vc;}
 else if (d3 < 0) {vb=va; va=vc;}

 i = p->projection(va,vb); p1 = i-(*va); p2 = i-(*vb);

 if (p1*p2 <=0) {if (cp) cp->setValue(i); return i.distance(p);}

 d1=p1.squaredLength(); d2=p2.squaredLength();
 if (d1<d2) {if (cp) cp->setValue(va); return p->distance(va);}
 else {if (cp) cp->setValue(vb); return p->distance(vb);}
}


///////////// Distance of point from the triangle //////////////

double Triangle::pointTriangleSquaredDistance(const Point *p) const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 Point n(((*va)-(*vb))&((*vb)-(*vc)));
 if (n.x == 0 && n.y == 0 && n.z == 0) return -1.0;

 double d1 = ((((*va)-(*vb))&((*vb)-(*p)))*n);
 double d2 = ((((*vb)-(*vc))&((*vc)-(*p)))*n);
 double d3 = ((((*vc)-(*va))&((*va)-(*p)))*n);

 if (d1 >= 0 && d2 >= 0 && d3 >= 0) return squaredDistanceFromPoint(p);

 if (d2 < 0) {va=vb; vb=vc;}
 else if (d3 < 0) {vb=va; va=vc;}

 Point i(p->projection(va,vb));
 Point p1(i-(*va)); Point p2(i-(*vb));

 if (p1*p2 <=0) return i.squaredDistance(p);

 d1=p1.squaredLength(); d2=p2.squaredLength();
 if (d1<d2) return p->squaredDistance(va);
 else return p->squaredDistance(vb);
}


/////////// Projection of point 'p' on the plane of the triangle /////

Point Triangle::project(const Point *p) const
{
 Point n = getNormal();
 if (n.isNull()) return INFINITE_POINT;

 Point v = *(v1());
 Point p1 = e1->toUnitVector();
 double d1 = (p1*(*p));
 Point p2 = e2->toUnitVector();
 double d2 = (p2*(*p));
 double d3 = (n*v);
 Point i = Point(d1,d2,d3).linearSystem(p1,p2,n);
 return i;
}


///////// TRUE if the point prjection belongs to the triangle /////////

bool Triangle::isInside(const Point *p) const
{
 Vertex va = *v1(), vb = *v2(), vc = *v3();
 Point v = *p;
 if (v == va || v == vb || v == vc) return 1;
 Point n = getNormal();
 if (n.isNull()) return 0;

 Point pa1 = ((vb-va)&(v-va));
 if (pa1.isNull()) return (((va-v)*(vb-v)) <= 0);
 Point pa2 = ((vc-vb)&(v-vb));
 if (pa2.isNull()) return (((vb-v)*(vc-v)) <= 0);
 Point pa3 = ((va-vc)&(v-vc));
 if (pa3.isNull()) return (((vc-v)*(va-v)) <= 0);

 double a1 = n.getAngle(pa1);
 double a2 = n.getAngle(pa2);
 double a3 = n.getAngle(pa3);

 return ((a1 >= PI2 && a2 >= PI2 && a3 >= PI2) || (a1 <= PI2 && a2 <= PI2 && a3 <= PI2));
}


/// Checks wether the triangle intersects the edge 'e'. If it does, ///
/// the method returns 1 and 'out' is initialized with the point of ///
/// intersection. If there is no intersection or if the edge lies   ///
/// on the plane of the triangle, 0 is returned. 't' is a threshold ///
/// distance within which a point is considered to be on the face. ///

bool Triangle::intersectsEdge(const Edge *e, Point *out, const double t) const
{
 Point p, p1, p2, p3, nor=getNormal();
 p1.setValue(v1()); p2.setValue(v2()); p3.setValue(v3());
 double dv1 = ((*(e->v1))-p1)*nor;
 double dv2 = ((*(e->v2))-p1)*nor;
 double dv3;
 if (FABS(dv1) < t) dv1 = 0;
 if (FABS(dv2) < t) dv2 = 0;

 if (dv1*dv2 > 0) return 0;
 if (dv1 == 0 && dv2 == 0) return 0;

 dv1 = FABS(dv1); dv2 = FABS(dv2);
 if (dv1 == 0) p.setValue(e->v1);
 else if (dv2 == 0) p.setValue(e->v2);
 else p = (((*(e->v2))*dv1)+((*(e->v1))*dv2))/(dv1+dv2);

 Point r1 = (p2-p1)&nor; r1.normalize(); dv1 = (p-p1)*r1;
 Point r2 = (p3-p2)&nor; r2.normalize(); dv2 = (p-p2)*r2;
 Point r3 = (p1-p3)&nor; r3.normalize(); dv3 = (p-p3)*r3;

 if (FABS(dv1) < t) dv1 = 0;
 if (FABS(dv2) < t) dv2 = 0;
 if (FABS(dv3) < t) dv3 = 0;

 if (dv1 > 0 || dv2 > 0 || dv3 > 0) return 0;

 out->setValue(p);
 return 1;
}


////////// Cap-like degeneracy //////////

Edge *Triangle::isCap() const
{
 if (e1->length()==0 || e2->length()==0 || e3->length()==0) return NULL;

 if (getAngle(v1()) == M_PI) return e3;
 if (getAngle(v2()) == M_PI) return e1;
 if (getAngle(v3()) == M_PI) return e2;
 return 0;
}

////////// Needle-like degeneracy //////////

Edge *Triangle::isNeedle() const
{
 if (e1->length() == 0) return e1;
 if (e2->length() == 0) return e2;
 if (e3->length() == 0) return e3;

 double a1 = getAngle(v1());
 double a2 = getAngle(v2());
 double a3 = getAngle(v3());

 if (a1 == M_PI || a2 == M_PI || a3 == M_PI) return 0;
 if (a1 == 0) return e3;
 if (a2 == 0) return e1;
 if (a3 == 0) return e2;
 return 0;
}


//////////// Degenerate check //////////////////

bool Triangle::isDegenerate() const
{
 return (isCap() != NULL || isNeedle() != NULL);
}


///////////// Overlap check ////////////////////

bool Triangle::overlaps() const
{
 Triangle *tt1 = t1();
 Triangle *tt2 = t2();
 Triangle *tt3 = t3();

 if (tt1 && getDAngle(tt1) == M_PI) return 1;
 if (tt2 && getDAngle(tt2) == M_PI) return 1;
 if (tt3 && getDAngle(tt3) == M_PI) return 1;

 return 0;
}



/// Debug

void Triangle::printTriangle(FILE *fp) const
{
 v1()->printPoint(fp);
 v2()->printPoint(fp);
 v3()->printPoint(fp);
}
