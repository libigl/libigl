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

#include "vertex.h"
#include "edge.h"
#include "triangle.h"
#include <stdlib.h>
#include <errno.h>


//////////////////// Constructors ////////////////////////

Vertex::Vertex() : Point()
{
 e0 = NULL;
 mask = 0;
}


Vertex::Vertex(const coord& a, const coord& b, const coord& c) : Point(a,b,c)
{
 e0 = NULL;
 mask = 0;
}


Vertex::Vertex(const Point *p) : Point(p->x, p->y, p->z)
{
 e0 = NULL;
 mask = 0;
}


Vertex::Vertex(const Point& p) : Point(p.x, p.y, p.z)
{
 e0 = NULL;
 mask = 0;
}

///////////////////// Destructor ///////////////////////

Vertex::~Vertex()
{
}


/////////////// VE relation ////////////////////////////

List *Vertex::VE() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 List *ve = new List();

 if (e0 == NULL) return ve;

 e = e0;
 do
 {
  ve->appendTail(e);
  v = e->oppositeVertex(this);
  t = e->leftTriangle(this);
  if (t == NULL) break;
  e = t->oppositeEdge(v);
 } while (e != e0);

 if (e == e0 && ve->numels() > 1) return ve;

 ve->popHead();
 e = e0;

 do
 {
  ve->appendHead(e);
  v = e->oppositeVertex(this);
  t = e->rightTriangle(this);
  if (t == NULL) break;
  e = t->oppositeEdge(v);
 } while (e != e0);

 return ve;
}

/////////////// VV relation ////////////////////////////

List *Vertex::VV() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 List *vv = new List();

 if (e0 == NULL) return vv;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  vv->appendTail(v);
  t = e->leftTriangle(this);
  if (t == NULL) break;
  e = t->oppositeEdge(v);
 } while (e != e0);

 if (e == e0 && vv->numels() > 1) return vv;

 vv->popHead();
 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  vv->appendHead(v);
  t = e->rightTriangle(this);
  if (t == NULL) break;
  e = t->oppositeEdge(v);
 } while (e != e0);

 return vv;
}

/////////////// VT relation ////////////////////////////

List *Vertex::VT() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 List *vt = new List();

 if (e0 == NULL) return vt;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  t = e->leftTriangle(this);
  if (t == NULL) break;
  vt->appendTail(t);
  e = t->oppositeEdge(v);
 } while (e != e0);

 if (e == e0 && vt->numels() > 1) return vt;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  t = e->rightTriangle(this);
  if (t == NULL) break;
  vt->appendHead(t);
  e = t->oppositeEdge(v);
 } while (e != e0);

 return vt;
}

/////////////// Returns the edge (this,v2) ////////////////

Edge *Vertex::getEdge(const Vertex *v2) const
{
 List *ve = VE();
 Node *m;
 Edge *e;

 FOREACHVEEDGE(ve, e, m)
  if (e->oppositeVertex(this) == v2) {delete(ve); return e;}

 delete(ve);
 return NULL;
}


///////////// Returns the vertex valence ///////////////////////

int Vertex::valence() const
{
 List *ve = VE();
 int n = ve->numels();
 delete(ve);
 return n;
}

/////////////// Checks the boundary ////////////////////////////

int Vertex::isOnBoundary() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;

 if (e0 == NULL) return 0;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  t = e->leftTriangle(this);
  if (t == NULL) return 1;
  e = t->oppositeEdge(v);
 } while (e != e0);

 return 0;
}


/////////////// Next boundary edge ////////////////////////////

Edge *Vertex::nextBoundaryEdge() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;

 if (e0 == NULL) return NULL;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  t = e->leftTriangle(this);
  if (t == NULL) return e;
  e = t->oppositeEdge(v);
 } while (e != e0);

 return NULL;
}


/////////////// Next boundary vertex ////////////////////////////

Vertex *Vertex::nextOnBoundary() const
{
 Edge *e = nextBoundaryEdge();
 if (e != NULL) return e->oppositeVertex(this);

 return NULL;
}


/////////////// Previous boundary edge ////////////////////////////

Edge *Vertex::prevBoundaryEdge() const
{
 Triangle *t;
 Edge *e;
 Vertex *v;

 if (e0 == NULL) return NULL;

 e = e0;
 do
 {
  v = e->oppositeVertex(this);
  t = e->rightTriangle(this);
  if (t == NULL) return e;
  e = t->oppositeEdge(v);
 } while (e != e0);

 return NULL;
}


/////////////// Previous boundary vertex ////////////////////////////

Vertex *Vertex::prevOnBoundary() const
{
 Edge *e = prevBoundaryEdge();
 if (e != NULL) return e->oppositeVertex(this);

 return NULL;
}


///// Vertex normal as weighted average of the incident triangle normals ////
///// The weight is the incidence angle.				 ////
///// Possible degenerate triangles are not taken into account.          ////

Point Vertex::getNormal() const
{
 List *vt = VT();
 Node *n;
 Triangle *t;
 double pa;
 Point tnor, ttn;

 FOREACHVTTRIANGLE(vt, t, n)
 {
  pa = t->getAngle(this);
  ttn = t->getNormal();
  if (!ttn.isNull()) tnor = tnor+(ttn*pa);
 }
 delete(vt);

 if (tnor.isNull()) return Point(0,0,0);

 tnor.normalize();
 return tnor;
}


////////// Returns the angle between the two boundary edges ////////

double Vertex::getBoundaryAngle() const
{
 Edge *e1 = prevBoundaryEdge();
 Edge *e2 = nextBoundaryEdge();
 if (e1 == NULL || e2 == NULL) return -1.0;
 Vertex *v1 = e1->oppositeVertex(this);
 Vertex *v2 = e2->oppositeVertex(this);
 double ang = getAngle(v1,v2);

 return ang;
}


////////// Returns the discriminant for triangulation ////////

double Vertex::getAngleForTriangulation() const
{
 Edge *e1 = prevBoundaryEdge();
 Edge *e2 = nextBoundaryEdge();
 if (e1 == NULL || e2 == NULL) return DBL_MAX;
 Triangle *t1 = e1->getBoundaryTriangle();
 Triangle *t2 = e2->getBoundaryTriangle();
 Vertex *v1 = e1->oppositeVertex(this);
 Vertex *v2 = e2->oppositeVertex(this);
 if ((*v2)==(*v1)) return -2;
 if (distance(v1)*distance(v2) == 0.0) return -1;

 double ang = getAngle(v1,v2);
 if (ang == M_PI) return 3*M_PI;
 if (ang == 0)    return 0;

 Edge e3(v1,v2);
 Triangle t(e1,e2,&e3);
 double da1 = t.getDAngle(t1);
 double da2 = t.getDAngle(t2);

 if (da1==M_PI && da2==M_PI) return (DBL_MAX/2.0);
 if (da1==M_PI || da2==M_PI) return (DBL_MAX/4.0);

 return da1+da2+ang;
}


////////// Returns the AP discriminant for triangulation ////////

double Vertex::getAngleOnAveragePlane(Point *nor) const
{
 Edge *e1 = prevBoundaryEdge();
 Edge *e2 = nextBoundaryEdge();
 if (e1 == NULL || e2 == NULL) return DBL_MAX;
 Vertex *v1 = e1->oppositeVertex(this);
 Vertex *v2 = e2->oppositeVertex(this);
 Point p, p1, p2;
 p1.setValue(v1);
 p2.setValue(v2);
 p.setValue(this);
 p.project(nor);
 p1.project(nor);
 p2.project(nor);
 if (p.distance(p1)*p.distance(p2) == 0.0)
 {
  JMesh::warning("getAngleOnAveragePlane: coincident projections\n");
  return 0.0;
 }
 double ang = p.getAngle(&p1, &p2);
 if (nor->side3D(&p1, &p, &p2) < 0) ang = 2*M_PI-ang;

 return ang;
}


///// mean curvature at the vertex: sum of signed dihedral angles ////

double Vertex::totalDihedralAngle() const
{
 List *ve = VE();
 double mc = 0;
 Edge *e;
 Node *n;

 FOREACHVEEDGE(ve, e, n)
  if (e->isOnBoundary()) {delete(ve); return DBL_MAX;}
  else mc += (M_PI-e->dihedralAngle());
 mc /= ve->numels();

 delete(ve);

 return mc;
}


///// Sum of all the incident angles ////

double Vertex::totalAngle() const
{
 List *ve = VE();
 double ta = 0.0;
 Edge *e;
 Node *n;

 FOREACHVEEDGE(ve, e, n)
  if (e->isOnBoundary()) {delete(ve); return -1.0;}
  else ta += e->leftTriangle(this)->getAngle(this);

 delete(ve);

 return ta;
}


/////// Voronoi area around the vertex ///////////////////

double Vertex::voronoiArea() const
{
 List *vt = VT();
 Node *n;
 Triangle *t;
 double va=0.0;

 FOREACHVTTRIANGLE(vt, t, n) va += t->area();
 delete(vt);

 return va/3.0;
}


// Closes the gap starting from this vertex
// If 'check_geom' is true, the zipping stops
// whether the coordinates of the vertices to
// be zipped are not equal.

int Vertex::zip(const bool check_geom)
{
 Node *n;
 Edge *e;

 List *ve = VE();
 Edge *be1 = (Edge *)ve->head()->data;
 Edge *be2 = (Edge *)ve->tail()->data;
 delete(ve);
 if (!be1->isOnBoundary() || !be2->isOnBoundary()) return 0;
 Vertex *ov1 = be1->oppositeVertex(this);
 Vertex *ov2 = be2->oppositeVertex(this);

 if (check_geom && ((*ov1)!=(*ov2))) return 0;

 if (ov1 != ov2)
 {
  ve = ov2->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(ov2, ov1);
  delete(ve);
  ov2->e0 = NULL;
 }

 Triangle *t = (be2->t1!=NULL)?(be2->t1):(be2->t2);
 t->replaceEdge(be2, be1);
 be1->replaceTriangle(NULL, t);
 be2->v1=be2->v2=NULL;
 e0 = ov1->e0 = be1;
 return 1+ov1->zip(check_geom);
}


/////// Progressive Mesh: Vertex split ///////////////////

Edge *Vertex::inverseCollapse(Vertex *v2, Vertex *v3, Vertex *v4)
{
 Edge *e, *e1, *e2=NULL, *e3=NULL, *e4;
 Triangle *t1, *t2, *ta1, *ta4;
 Node *n;
 Vertex *tmp;

 List *ve = VE();
 FOREACHVEEDGE(ve, e, n)
 {
  tmp = e->oppositeVertex(this);
  if (tmp == v3) e2 = e;
  else if (tmp == v4) e3 = e;
 }

 if (!e2 || !e3) {delete(ve); return NULL;}

 ta1 = e2->rightTriangle(this);
 ta4 = e3->leftTriangle(this);

 FOREACHVEEDGE(ve, e, n) if (e == e3) break;
 FOREACHNODECIRCULAR((*ve), n, n)
 {
  e = ((Edge *)n->data);
  if (e == e2) break;
  else e->replaceVertex(this, v2);
 }
 delete(ve);

 e = new Edge(this, v2);
 e1 = new Edge(v2, v3);
 e4 = new Edge(v2, v4);
 t1 = new Triangle(e, e1, e2);
 t2 = new Triangle(e, e3, e4);

 e->t1 = t1; e->t2 = t2;
 e2->replaceTriangle(ta1, t1);
 e3->replaceTriangle(ta4, t2);
 if (ta1) ta1->replaceEdge(e2, e1);
 if (ta4) ta4->replaceEdge(e3, e4);
 e1->t1 = t1; e1->t2 = ta1;
 e4->t1 = ta4; e4->t2 = t2;
 v2->e0 = e0 = e;

// Point p = (*this)-((*v2)-(*this));
// x = p.x; y = p.y; z = p.z;

 return e;
}


/////// Progressive Mesh: Vertex split ///////////////////

Edge *Vertex::inverseCollapse(Vertex *v2, Edge *e, Edge *e1, Edge *e2, Edge *e3, Edge *e4, Triangle *t1, Triangle *t2)
{
 Triangle *ta1, *ta4;
 Node *n;
 Edge *f;

 ta1 = e2->rightTriangle(this);
 ta4 = e3->leftTriangle(this);

 List *ve = VE();
 FOREACHVEEDGE(ve, f, n) if (f == e3) break;
 FOREACHNODECIRCULAR((*ve), n, n)
 {
  f = ((Edge *)n->data);
  if (f == e2) break;
  else f->replaceVertex(this, v2);
 }
 delete(ve);

 e->v1 = this; e->v2 = v2;
 e1->v1 = v2; e1->v2 = e2->oppositeVertex(this);
 e4->v1 = v2; e4->v2 = e3->oppositeVertex(this);
 t1->e1 = e; t1->e2 = e1; t1->e3 = e2;
 t2->e1 = e; t2->e2 = e3; t2->e3 = e4;

 e->t1 = t1; e->t2 = t2;
 e2->replaceTriangle(ta1, t1);
 e3->replaceTriangle(ta4, t2);
 if (ta1) ta1->replaceEdge(e2, e1);
 if (ta4) ta4->replaceEdge(e3, e4);
 e1->t1 = t1; e1->t2 = ta1;
 e4->t1 = ta4; e4->t2 = t2;
 v2->e0 = e0 = e;

 return e;
}


//// Returns the number of boundaries of the spherical neighborhood ////
//// of radius 'r'.  						    ////

#define INCREASE_INFO_FIELD(v) ((v)->info = (void *)(((j_voidint)(v)->info)+1))
#define DECREASE_INFO_FIELD(v) ((v)->info = (void *)(((j_voidint)(v)->info)-2))

int Vertex::getTopology(const double& r) const
{
 int numBoundaries=0;
 double r2 = r*r;
 List *ve, triList, trisInSphere, edgsInSphere, edgsOnBoundary;
 Node *n;
 Triangle *t, *s;
 Edge *e, *f;
 Vertex *w;

 t = (e0->t1 != NULL)?(e0->t1):(e0->t2);
 triList.appendHead(t);
 MARK_BIT(t,3);

 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.popHead();
  trisInSphere.appendHead(t);

  if ((s = t->t1()) != NULL && !IS_BIT(s,3) && (squaredDistance(s->oppositeVertex(t->e1))<r2))
       {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t2()) != NULL && !IS_BIT(s,3) && (squaredDistance(s->oppositeVertex(t->e2))<r2))
       {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t3()) != NULL && !IS_BIT(s,3) && (squaredDistance(s->oppositeVertex(t->e3))<r2))
       {triList.appendHead(s); MARK_BIT(s,3);}
 }

 while ((s = (Triangle *)trisInSphere.popHead()) != NULL)
 {
  UNMARK_BIT(s, 3);
  e = s->e1; if (!IS_BIT(e, 3)) {MARK_BIT(e, 3); edgsInSphere.appendHead(e);} else UNMARK_BIT(e, 3);
  e = s->e2; if (!IS_BIT(e, 3)) {MARK_BIT(e, 3); edgsInSphere.appendHead(e);} else UNMARK_BIT(e, 3);
  e = s->e3; if (!IS_BIT(e, 3)) {MARK_BIT(e, 3); edgsInSphere.appendHead(e);} else UNMARK_BIT(e, 3);
 }

 while ((e = (Edge *)edgsInSphere.popHead()) != NULL) if (IS_BIT(e, 3)) edgsOnBoundary.appendHead(e);
 FOREACHVEEDGE((&edgsOnBoundary), e, n) e->v1->info = e->v2->info = 0;
 FOREACHVEEDGE((&edgsOnBoundary), e, n)
 {
  INCREASE_INFO_FIELD(e->v1);
  INCREASE_INFO_FIELD(e->v2);
 }

 while ((e = (Edge *)edgsOnBoundary.popHead()) != NULL) if (IS_BIT(e, 3))
 {
  numBoundaries++;
  DECREASE_INFO_FIELD(e->v1); UNMARK_BIT(e, 3); w=e->v2;
  while (w->info)
  {
   e = NULL;
   ve = w->VE();
   FOREACHVEEDGE(ve, f, n) if (IS_BIT(f, 3)) {e = f; UNMARK_BIT(e, 3); break;}
   delete(ve);
   DECREASE_INFO_FIELD(w);
   w = e->oppositeVertex(w);
  }
 }

 return numBoundaries;
}

#undef INCREASE_INFO_FIELD
#undef DECREASE_INFO_FIELD
