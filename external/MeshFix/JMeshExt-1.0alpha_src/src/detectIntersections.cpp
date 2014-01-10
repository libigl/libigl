/****************************************************************************
* JMeshExt                                                                  *
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

#include "detectIntersections.h"
#include <string.h>
#include <stdlib.h>
#include "jqsort.h"
#include "jrs_predicates.h"


inline double di_remeshOrient3D(Point *p1, Point *p2, Point *p3, Point *p4)
{
 return orient3d((double *)(p1), (double *)(p2), (double *)(p3), (double *)(p4));
}

/////////// Check for coplanar edge-triangles !!!!!

// Intersection point between two edges.
// Edges are assumed to intersect properly.

Point *di_cell::edgeEdgeIntersection(Edge *a, Edge *b)
{
 Edge *e1, *e2;

 if (lexEdgeCompare(a, b)>0) {e1=a; e2=b;}
 else {e1=b; e2=a;}

 double d1 = (((*(e1->v2))-(*(e2->v2)))&((*(e1->v1))-(*(e2->v2)))).length();
 double d2 = (((*(e1->v2))-(*(e2->v1)))&((*(e1->v1))-(*(e2->v1)))).length();

 static Point p;
 p = (((*(e2->v1))*d1)+((*(e2->v2))*d2))/(d1+d2);

 return &p;
}


// Intersection point between 'e' and 't'.
// Returns NULL if e does not intersect t or if e and t are coplanar.

Point *di_cell::edgeIntersectsTriangle(Edge *e, Triangle *t, Edge **te)
{
 if (t->hasEdge(e)) return NULL;

 Vertex *v1=t->v1(), *v2=t->v2(), *v3=t->v3(), *v4, *v0;
 Vertex *pv0, *nv0;
 static Point p, p1, p2;
 static bool jrs_toinit=1;
 if (jrs_toinit) {exactinit(); jrs_toinit=0;}

 if (t->hasVertex(e->v1) || t->hasVertex(e->v2))
 {
  if (t->hasVertex(e->v1)) {v0 = e->v1; v4 = e->v2;}
  else if (t->hasVertex(e->v2)) {v0 = e->v2; v4 = e->v1;}

  if (di_remeshOrient3D((v1), (v2), (v3), (v4))==0.0)
  {
   pv0 = (t->prevVertex(v0));
   nv0 = (t->nextVertex(v0));
   p = (*v4)-(*v0);
   p1 = ((*pv0)-(*v4));
   p2 = ((*nv0)-(*v4));
   if (((p1&p)*(p2&p))>0.0) return NULL;
   p = (*pv0)-(*nv0);
   p2 = (*pv0)-(*v0);
   if (((p1&p)*(p2&p))>0.0) return NULL;
   return edgeEdgeIntersection(e, t->oppositeEdge(v0));
  }
  else return NULL;
 }

 (*te)=NULL;

 double d1 = di_remeshOrient3D((v1), (v2), (v3), ((e->v1)));
 double d2 = di_remeshOrient3D((v1), (v2), (v3), ((e->v2)));

 if (d1 == 0 && d2 == 0) return NULL;
 if ((d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)) return NULL;

 double e2 = di_remeshOrient3D((v1), (v2), ((e->v1)), ((e->v2)));
 double e3 = di_remeshOrient3D((v2), (v3), ((e->v1)), ((e->v2)));
 double e1 = di_remeshOrient3D((v3), (v1), ((e->v1)), ((e->v2)));

 if (e1==0 && e2==0 && e3==0) return NULL;

 if ((e1 >= 0 && e2 >= 0 && e3 >= 0) || (e1 <= 0 && e2 <= 0 && e3 <= 0))
 {
  if (d1==0) return e->v1;
  if (d2==0) return e->v2;
  if (e1==0 && e2==0) return v1;
  if (e2==0 && e3==0) return v2;
  if (e3==0 && e1==0) return v3;
// Edge-edge intersection must be exactly commutative (maybe sorting edges lexycographically)
  if (e1==0) {(*te)=t->e1; return edgeEdgeIntersection(e, t->e1);}
  if (e2==0) {(*te)=t->e2; return edgeEdgeIntersection(e, t->e2);}
  if (e3==0) {(*te)=t->e3; return edgeEdgeIntersection(e, t->e3);}
  d1 = FABS(d1); d2 = FABS(d2);
  p = (((*(e->v2))*d1)+((*(e->v1))*d2))/(d1+d2);
  return &p;
 }

 return NULL;
}



di_cell::di_cell(Triangulation *tin, bool useAll)
{
 Node *n;
 Vertex *v;
 Triangle *t;
 Mp.x = -DBL_MAX, mp.x = DBL_MAX;
 Mp.y = -DBL_MAX, mp.y = DBL_MAX;
 Mp.z = -DBL_MAX, mp.z = DBL_MAX;
 FOREACHVVVERTEX((&(tin->V)), v, n) if (useAll || IS_VISITED(v))
 {
  if (v->x < mp.x) mp.x = v->x;
  if (v->x > Mp.x) Mp.x = v->x;
  if (v->y < mp.y) mp.y = v->y;
  if (v->y > Mp.y) Mp.y = v->y;
  if (v->z < mp.z) mp.z = v->z;
  if (v->z > Mp.z) Mp.z = v->z;
 }

 mp -= DI_EPSILON_POINT;
 Mp += DI_EPSILON_POINT;

 FOREACHVTTRIANGLE((&(tin->T)), t, n) if (useAll || IS_VISITED(t)) triangles.appendTail(t);
}

bool di_cell::is_triangleBB_in_cell(Triangle *t)
{
 Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
 double mx = MIN(v1->x, MIN(v2->x, v3->x));
 double Mx = MAX(v1->x, MAX(v2->x, v3->x));
 double my = MIN(v1->y, MIN(v2->y, v3->y));
 double My = MAX(v1->y, MAX(v2->y, v3->y));
 double mz = MIN(v1->z, MIN(v2->z, v3->z));
 double Mz = MAX(v1->z, MAX(v2->z, v3->z));

 return (mx <= Mp.x && Mx >= mp.x && my <= Mp.y && My >= mp.y && mz <= Mp.z && Mz >= mp.z);
}

di_cell *di_cell::fork()
{
 Triangle *t;
 Point e = Mp-mp;
 di_cell *nc = new di_cell;
 List tmp;

 if (e.x > e.y && e.x > e.z)
  {nc->mp = mp; nc->Mp = Mp; nc->Mp.x -= (e.x/2); mp.x += (e.x/2);}
 else if (e.y > e.x && e.y > e.z)
  {nc->mp = mp; nc->Mp = Mp; nc->Mp.y -= (e.y/2); mp.y += (e.y/2);}
 else
  {nc->mp = mp; nc->Mp = Mp; nc->Mp.z -= (e.z/2); mp.z += (e.z/2);}
 
 while ((t=(Triangle *)triangles.popHead()) != NULL)
 {
  if (is_triangleBB_in_cell(t)) tmp.appendHead(t);
  if (nc->is_triangleBB_in_cell(t)) nc->triangles.appendHead(t);
 }
 triangles.joinTailList(&tmp);

 return nc;
}


// Returns TRUE if the cell does not contain self intersections for sure.
// Sufficient conditions to state this are (TO VERIFY):
// 1) All the triangles within the cell belong to the same connected component
// 2) There are no boundary edges
// 3) All the normals dot-multiply positively with the average normal

bool di_cell::doesNotIntersectForSure()
{
 Triangle *t, *s;
 Node *n;
 Point anor;
 int ns = 0;

 List tec;
 FOREACHVTTRIANGLE((&triangles), t, n) MARK_BIT(t, 3);
 FOREACHVTTRIANGLE((&triangles), t, n) if (!IS_VISITED2(t))
 {
  ns++; if (ns>1) break;
  tec.appendHead(t); MARK_VISIT2(t);
  while ((t=(Triangle *)tec.popHead()) != NULL)
  {
   s = t->t1(); if (s != NULL && IS_BIT(s, 3) && !IS_VISITED2(s)) {tec.appendTail(s); MARK_VISIT2(s);}
   s = t->t2(); if (s != NULL && IS_BIT(s, 3) && !IS_VISITED2(s)) {tec.appendTail(s); MARK_VISIT2(s);}
   s = t->t3(); if (s != NULL && IS_BIT(s, 3) && !IS_VISITED2(s)) {tec.appendTail(s); MARK_VISIT2(s);}
  }
 }
 FOREACHVTTRIANGLE((&triangles), t, n) {UNMARK_VISIT2(t); UNMARK_BIT(t, 3);}
 if (ns>1) return 0;

 FOREACHVTTRIANGLE((&triangles), t, n)
 {
  if (t->e1->isOnBoundary() || t->e2->isOnBoundary() || t->e3->isOnBoundary()) return 0;
  anor += DI_STORED_NORMAL(t);
 }

 FOREACHVTTRIANGLE((&triangles), t, n) if (anor*DI_STORED_NORMAL(t) <= 0) return 0;

 return 1;
}

// Brute force all-with-all intersection test of the triangles in 'triangles'.
void di_cell::di_selectIntersections()
{
 Triangle *t, *y;
 Edge *e, *te;
 Node *n, *m;
 List edges;

 FOREACHVTTRIANGLE((&triangles), t, n) MARK_VISIT2(t);
 FOREACHVTTRIANGLE((&triangles), t, n)
 {
  e = t->e1; y = t->t1(); if (!IS_VISITED2(e) && ((y!=NULL && IS_VISITED2(y)) || y==NULL)) {MARK_VISIT2(e); edges.appendHead(e);}
  e = t->e2; y = t->t2(); if (!IS_VISITED2(e) && ((y!=NULL && IS_VISITED2(y)) || y==NULL)) {MARK_VISIT2(e); edges.appendHead(e);}
  e = t->e3; y = t->t3(); if (!IS_VISITED2(e) && ((y!=NULL && IS_VISITED2(y)) || y==NULL)) {MARK_VISIT2(e); edges.appendHead(e);}
 }
 FOREACHVTTRIANGLE((&triangles), t, n) UNMARK_VISIT2(t);
 FOREACHVEEDGE((&edges), e, n) UNMARK_VISIT2(e);

 FOREACHVTTRIANGLE((&triangles), t, n)
  FOREACHVEEDGE((&edges), e, m)
   if (edgeIntersectsTriangle(e, t, &te))
   {
    MARK_VISIT(t);
    if (e->t1 != NULL) MARK_VISIT(e->t1);
    if (e->t2 != NULL) MARK_VISIT(e->t2);
   }
}


/////////////////////////////////////////////////////////////////////////
//                                                                     ||
////////////////////// Select   Intersections ///////////////////////////
//                                                                     ||
/////////////////////////////////////////////////////////////////////////

int ExtTriMesh::selectIntersectingTriangles(UINT16 tris_per_cell)
{
 Triangle *t;
 Vertex *v;
 Node *n;
 bool isSelection=0;
 List *selT = new List, *selV = new List;

 JMesh::begin_progress();
 JMesh::report_progress(NULL);

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  isSelection=1;
  selT->appendTail(t);
  v=t->v1(); if (!IS_VISITED(v)) {MARK_VISIT(v); selV->appendTail(v);}
  v=t->v2(); if (!IS_VISITED(v)) {MARK_VISIT(v); selV->appendTail(v);}
  v=t->v3(); if (!IS_VISITED(v)) {MARK_VISIT(v); selV->appendTail(v);}
 }
 JMesh::report_progress(NULL);

 if (!isSelection) {delete(selT); delete(selV); selT=&T; selV=&V;}

 // Store triangle normal in each triangle's info field
 FOREACHVTTRIANGLE(selT, t, n) DI_STORED_PANORMAL(t) = new Point(t->getNormal());

 di_cell *c2, *c = new di_cell(this, !isSelection);
 List cells, todo(c);
 int i=0;

 while ((c = (di_cell *)todo.popHead()) != NULL)
 {
  if (i>DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= tris_per_cell) cells.appendHead(c);
  else
  {
   i++;
   JMesh::report_progress(NULL);
   c2 = c->fork();
   if (c->doesNotIntersectForSure()) delete(c); else todo.appendTail(c);
   if (c2->doesNotIntersectForSure()) delete(c2); else todo.appendTail(c2);
  }
 }

 // Deselect everything and select only intersecting triangles
 deselectTriangles();
 i=0; FOREACHNODE(cells, n)
 {
  (((di_cell *)n->data)->di_selectIntersections());
  JMesh::report_progress("%d %% done   ",((i++)*100)/cells.numels());
 }
 JMesh::end_progress();

 // Dispose memory allocated for cells
 while (cells.numels()) delete((di_cell *)cells.popHead());

 // Count selected triangles for final report and delete stored normals
 int its=0;
 FOREACHVTTRIANGLE(selT, t, n) {delete(DI_STORED_PNORMAL(t)); t->info = NULL; if (IS_VISITED(t)) its++;}

 if (its) JMesh::info("%d intersecting triangles have been selected.\n",its);
 else JMesh::info("No intersections detected.\n");

 FOREACHVVVERTEX(selV, v, n) UNMARK_VISIT(v);
 if (isSelection) {delete(selT); delete(selV);}

 return its;
}
