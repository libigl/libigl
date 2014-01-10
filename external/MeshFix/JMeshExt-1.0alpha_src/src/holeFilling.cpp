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

#include "holeFilling.h"


////////// Generic method for patching holes. Heuristic. /////////////
////////// Small angles are patched first, if possible.  /////////////

int ExtTriMesh::TriangulateHole(Edge *e)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0;
 Triangle *t;
 v = e->v1;

 t = (e->t1!=NULL)?(e->t1):(e->t2);
 if (t->nextEdge(e)->isOnBoundary() && t->prevEdge(e)->isOnBoundary()) return 0;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);

 while (bvs.numels() > 2)
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&bvs), v, n)
   if (!IS_VISITED(v) && v->e0 && (ang = v->getAngleForTriangulation()) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   JMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&bvs), v, n) UNMARK_VISIT(v);
   int i=0; FOREACHTRIANGLE(t, n) if (i++==nt) break; else unlinkTriangle(t);
   removeUnlinkedElements();
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  if ((t=EulerEdgeTriangle(e1,e2))==NULL) MARK_VISIT(v);
  else {bvs.removeCell(gn); UNMARK_VISIT(v1); UNMARK_VISIT(v2); MARK_VISIT(t); nt++;}
 }

 return nt;
}


// Fills the hole identified by 'e' and leaves the new triangle selected.
// 'refine_and_smooth' is for bi-laplacian fairing.

void ExtTriMesh::FillHole(Edge *e, bool refine_and_smooth)
{
 int i, nt;
 Node *n;
 Triangle *t;
 Vertex *v;

 deselectTriangles();
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
 nt = TriangulateHole(e);
 if (!nt) return;

 i=0; FOREACHTRIANGLE(t, n) if (i++==nt) break; else MARK_VISIT(t);

 if (refine_and_smooth)
 {
  t = (Triangle *)T.head()->data;
  refineSelectedHolePatches(t);
  fairSelection(t);
 }
}

//// Triangulate Small Boundaries (with less than 'nbe' edges) /////

int ExtTriMesh::fillSmallBoundaries(int nbe, bool refine_patches, bool smooth_patches)
{
 Vertex *v,*w;
 Triangle *t;
 Node *n;
 int grd, is_selection=0, tbds = 0, pct = 100;
 List bdrs;

 JMesh::begin_progress();
 JMesh::report_progress("0%% done ");

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {is_selection=1; break;}

 if (is_selection) FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
  {MARK_VISIT2(t->v1()); MARK_VISIT2(t->v2()); MARK_VISIT2(t->v3());}

 FOREACHVERTEX(v, n)
 {
  grd = 0;
  if (!IS_VISITED2(v) && v->isOnBoundary())
  {
   tbds++;
   w = v;
   do
   {
    if (IS_VISITED2(w)) grd=nbe+1;
    MARK_VISIT2(w);
    grd++;
    w = w->nextOnBoundary();
   } while (w != v);
   if (grd <= nbe) bdrs.appendHead(w->nextBoundaryEdge());
  }
 }
 FOREACHVERTEX(v, n) {UNMARK_VISIT(v); UNMARK_VISIT2(v);}

 deselectTriangles();

 pct=0; FOREACHNODE(bdrs, n)
 {
  if (TriangulateHole((Edge *)n->data) && refine_patches)
  {
   t = (Triangle *)T.head()->data;
   if (!refineSelectedHolePatches(t) && smooth_patches) fairSelection(t);
  }
  JMesh::report_progress("%d%% done ",((++pct)*100)/bdrs.numels());
 }

 grd = bdrs.numels();

 JMesh::end_progress();

 return grd;
}


// Inserts new vertices in the current selection so as
// to reflect the density of the surrounding mesh.
// This method assumes that the selection has no internal vertices.

int ExtTriMesh::refineSelectedHolePatches(Triangle *t0)
{
 Node *n, *m;
 Triangle *t, *t1, *t2;
 Edge *e, *f;
 Vertex *v;
 List *ve, toswap, reg, all_edges, interior_edges, boundary_edges, boundary_vertices, interior_vertices;
 double sigma, l, sv1, sv2, sv3, dv1, dv2, dv3;
 int swaps, totits, nee, ntb, nnt=-1, pnnt, gits=0;
 const double alpha = 1.4142136;
 Point vc;

 if (t0 != NULL)
 {
  UNMARK_VISIT(t0); toswap.appendHead(t0);
  while ((t=(Triangle *)toswap.popHead()) != NULL)
  {
   reg.appendHead(t);
   t1=t->t1(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t2(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t3(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
  }
  FOREACHVTTRIANGLE((&reg), t, n) MARK_VISIT(t);
 }
 else FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) reg.appendHead(t);

 FOREACHVTTRIANGLE((&reg), t, n)
 {
  e = t->e1; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);
  e = t->e2; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);
  e = t->e3; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);
 }

 while (all_edges.numels())
 {
  e = (Edge *)all_edges.popHead();
  if (IS_VISITED2(e)) {boundary_edges.appendHead(e); UNMARK_VISIT2(e);}
  else {interior_edges.appendHead(e); MARK_VISIT2(e);}
 }

 FOREACHVEEDGE((&boundary_edges), e, n)
 {
  v = e->v1; if (!IS_VISITED2(v)) {MARK_VISIT2(v); boundary_vertices.appendHead(v);}
  v = e->v2; if (!IS_VISITED2(v)) {MARK_VISIT2(v); boundary_vertices.appendHead(v);}
 }

 FOREACHVVVERTEX((&boundary_vertices), v, n) UNMARK_VISIT2(v);

 // Due to the above definitions, interior edges are VISITED2

 FOREACHVVVERTEX((&boundary_vertices), v, n)
 {
  ve = v->VE();
  sigma=0; nee=0; FOREACHVEEDGE(ve, e, m) if (!IS_VISITED2(e)) {nee++; sigma += e->length();}
  sigma /= nee; v->info = new double(sigma);
  delete(ve);
 }

 FOREACHVEEDGE((&interior_edges), e, n) UNMARK_VISIT2(e);
 FOREACHVEEDGE((&boundary_edges), e, n) MARK_BIT(e, 3);

 do
 {
  pnnt=nnt;
  nnt=0;
  FOREACHVTTRIANGLE((&reg), t, n)
  {
   vc = t->getCenter();
   sv1 = (*(double *)t->v1()->info);
   sv2 = (*(double *)t->v2()->info);
   sv3 = (*(double *)t->v3()->info); 
   sigma = (sv1+sv2+sv3)/3.0;
   dv1 = alpha*(t->v1()->distance(&vc));
   dv2 = alpha*(t->v2()->distance(&vc));
   dv3 = alpha*(t->v3()->distance(&vc));
   if (dv1>sigma && dv1>sv1 && dv2>sigma && dv2>sv2 && dv3>sigma && dv3>sv3)   
   {
    ntb = T.numels();
    v = splitTriangle(t,&vc,1);
    nnt += (T.numels()-ntb);
    if (T.numels() == ntb+2)
    {
     v->info = new double(sigma);
     interior_vertices.appendHead(v);
     interior_edges.appendHead(v->e0);
     interior_edges.appendHead(v->e0->leftTriangle(v)->prevEdge(v->e0));
     interior_edges.appendHead(v->e0->rightTriangle(v)->nextEdge(v->e0));
     t1 = ((Triangle *)T.head()->data);
     t2 = ((Triangle *)T.head()->next()->data);
     t1->mask = t2->mask = t->mask;
     reg.appendHead(t1); reg.appendHead(t2);
    }
   }
  }

  FOREACHVEEDGE((&interior_edges), e, n) {MARK_VISIT2(e); toswap.appendHead(e);}
  totits=0; swaps=1;
  while (swaps && totits++ < 10)
  {
   swaps = 0; 
   while ((e=(Edge *)toswap.popHead())!=NULL)
   {
    UNMARK_VISIT2(e);
    l = e->squaredLength();
    if (e->swap())
    {
     if (e->squaredLength() >= l*0.999999) e->swap(1);
     else
     {
      swaps++;
      toswap.appendTail(e);
      f = e->t1->nextEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
      f = e->t1->prevEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
      f = e->t2->nextEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
      f = e->t2->prevEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
     }
    }
   }
  }

  if (pnnt==nnt) gits++;
 } while (nnt && gits<10);

 FOREACHVEEDGE((&boundary_edges), e, n) UNMARK_BIT(e, 3);
 FOREACHVVVERTEX((&boundary_vertices), v, n) {delete((double *)v->info); v->info=NULL;}
 FOREACHVVVERTEX((&interior_vertices), v, n) {delete((double *)v->info); v->info=NULL;}

 if (gits>=10) {JMesh::warning("Fill holes: Refinement stage failed to converge. Breaking.\n"); return 1;}

 return 0;
}


// Fairs the inner vertices of the selection using a second-order umbrella operator
// similar to a boundary constrained bi-laplacian smoothing

void fs_sparseSystem::solve(List *vl)
{
 if (kterm_size != 3) JMesh::error("fs_sparseSystem::solve(List *): Known term size is not 3!\n");
 Node *n;
 int i, nv = vl->numels();
 if (nv != num_variables) JMesh::error("fs_sparseSystem::solve(List *): Vertex list size does not match system size!\n");
 double *x = new double[nv];
 sparseSystem::solve(x, 0);
 for(i=0, n=vl->head(); i<nv; i++, n=n->next()) ((Vertex *)n->data)->x = -x[i];
 sparseSystem::solve(x, 1);
 for(i=0, n=vl->head(); i<nv; i++, n=n->next()) ((Vertex *)n->data)->y = -x[i];
 sparseSystem::solve(x, 2);
 for(i=0, n=vl->head(); i<nv; i++, n=n->next()) ((Vertex *)n->data)->z = -x[i];
 delete(x);
}

void ExtTriMesh::fairSelection(Triangle *t0)
{
 Node *n, *m, *o;
 Triangle *t;
 Edge *e1, *e2;
 Vertex *v, *vicino, *vicino2;
 List *vt, interior_vertices, all_vertices, *vicini, *vicini2;
 int i, j, isb, niv;
 double peso, W_j, W_vicino, W_j_vicino, W_vicino_vicino2;

// FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 if (t0 != NULL)
 {
  List toswap, reg;
  Triangle *t1;
  UNMARK_VISIT(t0); toswap.appendHead(t0);
  while ((t=(Triangle *)toswap.popHead()) != NULL)
  {
   reg.appendHead(t);
   v = t->v1(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
   v = t->v2(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
   v = t->v3(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
   t1=t->t1(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t2(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t3(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
  }
  FOREACHVTTRIANGLE((&reg), t, n) MARK_VISIT(t);
 }
 else
 {
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
  {
   v = t->v1(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
   v = t->v2(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
   v = t->v3(); if (!IS_VISITED(v)) {MARK_VISIT(v); all_vertices.appendHead(v);}
  }
 }

 FOREACHVVVERTEX((&all_vertices), v, n)
 {
  isb=0;
  vt = v->VT();
  FOREACHVTTRIANGLE(vt, t, m) if (!IS_VISITED(t)) {isb=1; break;}
  delete(vt);
  if (isb) UNMARK_VISIT(v);
  else interior_vertices.appendHead(v);
 }

 if (!interior_vertices.numels()) return;

 niv = interior_vertices.numels();
 fs_sparseSystem sps(niv);

 for(i=0, n=interior_vertices.head(); i<niv; i++, n=n->next()) 
 { 
  sps.sumCoefficient(1, i, i);
  ((Vertex *)n->data)->info = (void *)i;
 }

 for(j=0, o=interior_vertices.head(); j<niv; j++, o=o->next())
 {
  v = (Vertex *)o->data;
  vicini = v->VE();
  W_j=0; FOREACHVEEDGE(vicini, e1, n) W_j += e1->length();
  FOREACHVEEDGE(vicini, e1, n)
  {
   vicino = e1->oppositeVertex(v);
   vicini2 = vicino->VE();
   W_vicino=0; FOREACHVEEDGE(vicini2, e2, m) W_vicino += e2->length();
   W_j_vicino = e1->length();
   
   // Alec: replaced "int" with "j_voidint"
   if (IS_VISITED(vicino)) sps.sumCoefficient(-2*W_j_vicino/W_j, j, (j_voidint)vicino->info);
   else
   {
    sps.sumKnownTerm(2*W_j_vicino/W_j*vicino->x, j, 0);
    sps.sumKnownTerm(2*W_j_vicino/W_j*vicino->y, j, 1);
    sps.sumKnownTerm(2*W_j_vicino/W_j*vicino->z, j, 2);
   }
   
   FOREACHVEEDGE(vicini2, e2, m)
   {
    vicino2 = e2->oppositeVertex(vicino);
    W_vicino_vicino2 = e2->length();
    peso = W_j_vicino * W_vicino_vicino2 / ( W_j * W_vicino);
    // Alec: replaced "int" with "j_voidint"
    if (IS_VISITED(vicino2)) sps.sumCoefficient(peso, j, (j_voidint)vicino2->info);
    else
    {
     sps.sumKnownTerm(-peso*vicino2->x, j, 0);
     sps.sumKnownTerm(-peso*vicino2->y, j, 1);
     sps.sumKnownTerm(-peso*vicino2->z, j, 2);
    }
   }
   delete vicini2;
  }
  delete vicini;
 }

 sps.solve(&interior_vertices);

 FOREACHVVVERTEX((&all_vertices), v, n) UNMARK_VISIT(v);
}


// This method looks if exactly two vertices belonging to two different boundary loops are
// selected. If so, such a pair is joined through an edge and a pair of triangles is added
// to change the topology of the mesh (2 boundaries -> 1 boundary).
// On success, the return value is the new edge connecting the two vertices selected.
// Returns NULL if it could not detect such a pair of vertices.
// If 'justconnect', the remaining hole is not patched.
// If 'refine', the patching triangles are refined to reproduce neighboring density
// If 'fair', the refined vertices are moved so as get tangential continuity

Edge *ExtTriMesh::joinBoundaryLoops(bool justconnect, bool refine, bool fair)
{
 Vertex *v, *gv=NULL, *gw=NULL;
 Node *n;

 FOREACHVERTEX(v, n) if (IS_VISITED(v))
  {if (gv==NULL) gv=v; else if (gw==NULL) gw=v; else return NULL;}

 return joinBoundaryLoops(gv, gw, justconnect, refine, fair);
}

Edge *ExtTriMesh::joinBoundaryLoops(Vertex *gv, Vertex *gw, bool justconnect, bool refine, bool fair)
{
 Vertex *v, *gvn, *gwn;
 Edge *e, *gve, *gwe;
 Triangle *t;
 Node *n;
 double tl1=0.0, tl2=0.0, pl1, pl2;

 if (gv==NULL || gw==NULL || !gv->isOnBoundary() || !gw->isOnBoundary()) return NULL;

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
 deselectTriangles();

 v=gv;
 if (!justconnect)
 {
  do {v=v->nextOnBoundary(); if (v==gw) return NULL;} while (v!=gv);
 }
 else
 {
  gvn=gv->nextOnBoundary(); gwn=gv->prevOnBoundary();
  if (gw==gvn || gw==gwn) return NULL;
  if (gw == gvn->nextOnBoundary())
   {t=EulerEdgeTriangle(gvn->prevBoundaryEdge(), gvn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gvn);}
  if (gw == gwn->prevOnBoundary())
   {t=EulerEdgeTriangle(gwn->prevBoundaryEdge(), gwn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gwn);}
 }

 gve = gv->prevBoundaryEdge();
 gvn = gve->oppositeVertex(gv);
 gwe = gw->nextBoundaryEdge();
 gwn = gwe->oppositeVertex(gw);

 Edge *je = CreateEdge(gv, gw);
 Edge *je1 = CreateEdge(gv, gwn);
 Edge *je2 = CreateEdge(gwn, gvn);

 t = CreateTriangle(je, gwe, je1); MARK_VISIT(t);
 t = CreateTriangle(je1, je2, gve); MARK_VISIT(t);

 if (justconnect) return je;

 v=gv; do {e=v->nextBoundaryEdge(); v=e->oppositeVertex(v); tl1+=e->length();} while (v!=gv);
 v=gw; do {e=v->nextBoundaryEdge(); v=e->oppositeVertex(v); tl2+=e->length();} while (v!=gw);
 pl1=tl1; pl2=tl2;

 double c1, c2;

 e=je;
 while (e->isOnBoundary())
 {
  gv = (e->t2 != NULL)?(e->v2):(e->v1); gve=gv->nextBoundaryEdge();
  gw = (e->t1 != NULL)?(e->v2):(e->v1); gwe=gw->prevBoundaryEdge();
  c1 = fabs((pl1-gve->length())*tl2-pl2*tl1);
  c2 = fabs((pl2-gwe->length())*tl1-pl1*tl2);
  if (c1<c2)
  {
   t=EulerEdgeTriangle(e, gve); MARK_VISIT(t);
   pl1 -= gve->length();
   e = t->nextEdge(gve);
  }
  else
  {
   t=EulerEdgeTriangle(gwe, e); MARK_VISIT(t);
   pl2 -= gwe->length();
   e = t->prevEdge(gwe);
  }
 }

 if (refine)
 {
  refineSelectedHolePatches();
  if (fair) fairSelection();
 }

 return je;
}

/*
Edge *ExtTriMesh::joinBoundaryLoops(bool justconnect, bool refine, bool fair)
{
 Vertex *v, *gv=NULL, *gw=NULL, *gvn, *gwn;
 Edge *e, *gve, *gwe;
 Triangle *t;
 Node *n;

 FOREACHVERTEX(v, n) if (IS_VISITED(v))
  {if (gv==NULL) gv=v; else if (gw==NULL) gw=v; else return NULL;}

 if (gv==NULL || gw==NULL || !gv->isOnBoundary() || !gw->isOnBoundary()) return NULL;

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 v=gv;
 if (!justconnect)
 {
  do {v=v->nextOnBoundary(); if (v==gw) return NULL;} while (v!=gv);
 }
 else
 {
  gvn=gv->nextOnBoundary(); gwn=gv->prevOnBoundary();
  if (gw==gvn || gw==gwn) return NULL;
  if (gw == gvn->nextOnBoundary())
   {t=EulerEdgeTriangle(gvn->prevBoundaryEdge(), gvn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gvn);}
  if (gw == gwn->prevOnBoundary())
   {t=EulerEdgeTriangle(gwn->prevBoundaryEdge(), gwn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gwn);}
 }

 gve = gv->prevBoundaryEdge();
 gvn = gve->oppositeVertex(gv);
 gwe = gw->nextBoundaryEdge();
 gwn = gwe->oppositeVertex(gw);

 Edge *je = CreateEdge(gv, gw);
 Edge *je1 = CreateEdge(gv, gwn);
 Edge *je2 = CreateEdge(gwn, gvn);

 t = CreateTriangle(je, gwe, je1); MARK_VISIT(t);
 t = CreateTriangle(je1, je2, gve); MARK_VISIT(t);

 if (justconnect) return je;

 double lb1=0.0, lb2=0.0;
 int nb1=0, nb2=0;
 Edge *ne;

 v=gv;
 while ((e=v->nextBoundaryEdge()) != je2)
  { lb1 += e->length(); nb1++; v=e->oppositeVertex(v); }
 v=e->oppositeVertex(v);
 while ((e=v->nextBoundaryEdge()) != je)
  { lb2 += e->length(); nb2++; v=e->oppositeVertex(v); }

 e=je;
 while (nb1+nb2)
 {
  if (lb1>lb2)
  {
   v = (e->t2 != NULL)?(e->v2):(e->v1);
   ne=v->nextBoundaryEdge(); t=EulerEdgeTriangle(e, ne); MARK_VISIT(t);
   lb1 -= ne->length(); nb1--;
   e = t->nextEdge(ne);
  }
  else
  {
   v = (e->t1 != NULL)?(e->v2):(e->v1);
   ne=v->prevBoundaryEdge(); t=EulerEdgeTriangle(ne, e); MARK_VISIT(t);
   lb2 -= ne->length(); nb2--;
   e = t->prevEdge(ne);
  }
 }

 if (refine)
 {
  refineSelectedHolePatches();
  if (fair) fairSelection();
 }

 return je;
}
*/
