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

#include "tin.h"
#include <stdlib.h>
#include <string.h>


//////////////////////////////////////////////////////////////////
//                                                              //
//    C L A S S   C O N S T R U C T O R S                       //
//                                                              //
//////////////////////////////////////////////////////////////////


///////////////////// Constructor (Empty) ////////////////////

Triangulation::Triangulation()
{
 n_boundaries = n_handles = n_shells = 0;
 d_boundaries = d_handles = d_shells = 0;
}


//////////////////// Constructor (Pre-defined) ///////////////

Triangulation::Triangulation(const char *tin_definition)
{
 if (!strcmp(tin_definition, "tetrahedron"))
 {
  Vertex *v1 = new Vertex(-1,-1.4142136,0);
  Vertex *v2 = new Vertex(-1,1.4142136,0);
  Vertex *v3 = new Vertex(1,0,-1.4142136);
  Vertex *v4 = new Vertex(1,0,1.4142136);
  Edge   *e1 = new Edge(v1,v2); v1->e0 = e1;
  Edge   *e2 = new Edge(v2,v3); v2->e0 = e2;
  Edge   *e3 = new Edge(v3,v1); v3->e0 = e3;
  Edge   *e4 = new Edge(v1,v4); v4->e0 = e4;
  Edge   *e5 = new Edge(v2,v4);
  Edge   *e6 = new Edge(v3,v4);
  Triangle *t1 = new Triangle(e1,e2,e3);
  Triangle *t2 = new Triangle(e1,e4,e5);
  Triangle *t3 = new Triangle(e2,e5,e6);
  Triangle *t4 = new Triangle(e3,e6,e4);
  e1->t1 = t1; e1->t2 = t2;
  e2->t1 = t1; e2->t2 = t3;
  e3->t1 = t1; e3->t2 = t4;
  e4->t1 = t2; e4->t2 = t4;
  e5->t1 = t3; e5->t2 = t2;
  e6->t1 = t4; e6->t2 = t3;
  V.appendHead(v1); V.appendHead(v2); V.appendHead(v3); V.appendHead(v4);
  T.appendHead(t1); T.appendHead(t2); T.appendHead(t3); T.appendHead(t4);
  E.appendHead(e1); E.appendHead(e2); E.appendHead(e3);
  E.appendHead(e4); E.appendHead(e5); E.appendHead(e6);

  n_boundaries = 0;
  n_handles = 0;
  n_shells = 1;
  d_boundaries = d_handles = d_shells = 0;
 }
 else JMesh::error("Unknown triangulation type '%s'\n",tin_definition);
}


///////////////////// Cloning TIN ///////////////////////////

Triangulation::Triangulation(const Triangulation *tin, const bool clone_info)
{
 Node *n;
 Vertex *v, *nv;
 Edge *e, *ne;
 Triangle *t, *nt;

 int i;
 void **t_info = new void *[tin->T.numels()];
 i=0; FOREACHVTTRIANGLE((&(tin->T)), t, n) t_info[i++]=t->info;
 void **e_info = new void *[tin->E.numels()];
 i=0; FOREACHVEEDGE((&(tin->E)), e, n) e_info[i++]=e->info;
 void **v_info = new void *[tin->V.numels()];
 i=0; FOREACHVVVERTEX((&(tin->V)), v, n) v_info[i++]=v->info;

 FOREACHVVVERTEX((&(tin->V)), v, n)
  {nv=new Vertex(v); V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&(tin->E)), e, n)
  {ne=new Edge((Vertex *)e->v1->info, (Vertex *)e->v2->info); E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&(tin->T)), t, n)
  {nt=new Triangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); T.appendTail(nt); t->info = nt;}

 FOREACHVVVERTEX((&(tin->V)), v, n) {((Vertex *)v->info)->e0 = (Edge *)v->e0->info; v->info = NULL;}

 FOREACHVEEDGE((&(tin->E)), e, n)
  {((Edge *)e->info)->t1 = (e->t1)?((Triangle *)e->t1->info):(NULL); ((Edge *)e->info)->t2 = (e->t2)?((Triangle *)e->t2->info):(NULL); e->info = NULL;}

 i=0; FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info=t_info[i++];
 i=0; FOREACHVEEDGE((&(tin->E)), e, n) e->info=e_info[i++];
 i=0; FOREACHVVVERTEX((&(tin->V)), v, n) v->info=v_info[i++];

 if (clone_info)
 {
  i=0; FOREACHTRIANGLE(t, n) t->info=t_info[i++];
  i=0; FOREACHEDGE(e, n) e->info=e_info[i++];
  i=0; FOREACHVERTEX(v, n) v->info=v_info[i++];  
 }
 delete [] t_info; delete [] e_info; delete [] v_info;

 d_boundaries = d_handles = d_shells = 1;
}


//// Creates a new Triangulation out of a connected component of an existing Triangulation.
//// If 'keep_reference' is set to 'true', each element of the existing mesh keeps a
//// pointer to the corresponding new element in the 'info' field.

Triangulation::Triangulation(const Triangle *t0, const bool keep_reference)
{
 List todo(t0), st, sv, se;
 Node *n;
 Triangle *t, *nt;
 Vertex *v, *nv;
 Edge *e, *ne;

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  st.appendHead(t);
  nt=t->t1(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
  nt=t->t2(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
  nt=t->t3(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
 }

 FOREACHVTTRIANGLE((&st), t, n)
 {
  UNMARK_VISIT2(t);
  e = t->e1; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  e = t->e2; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  e = t->e3; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  v = t->v1(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
  v = t->v2(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
  v = t->v3(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
 }

 FOREACHVVVERTEX((&sv), v, n)
  {UNMARK_VISIT2(v); nv=new Vertex(v); V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&se), e, n)
  {UNMARK_VISIT2(e); ne=new Edge((Vertex *)e->v1->info, (Vertex *)e->v2->info); E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&st), t, n)
  {nt=new Triangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); T.appendTail(nt); t->info = nt;}

 FOREACHVVVERTEX((&sv), v, n) ((Vertex *)v->info)->e0 = (Edge *)v->e0->info;

 FOREACHVEEDGE((&se), e, n)
  {((Edge *)e->info)->t1 = (e->t1)?((Triangle *)e->t1->info):(NULL); ((Edge *)e->info)->t2 = (e->t2)?((Triangle *)e->t2->info):(NULL);}

 if (!keep_reference)
 {
  FOREACHVVVERTEX((&sv), v, n) v->info = NULL;
  FOREACHVEEDGE((&se), e, n) e->info = NULL;
  FOREACHVTTRIANGLE((&st), t, n) t->info = NULL;
 }

 eulerUpdate();
}


///////////////////// Destructor ///////////////////////////

Triangulation::~Triangulation()
{
 T.freeNodes();
 V.freeNodes();
 E.freeNodes();
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    P R I M I T I V E   C O N S T R U C T I O N               //
//                                                              //
//////////////////////////////////////////////////////////////////


//////////////////// Creates an edge ////////////////////////////

Edge *Triangulation::CreateEdge(Vertex *v1, Vertex *v2)
{
 Edge *e;
 
 if ((e = v1->getEdge(v2)) != NULL) return e;
 
 e = new Edge(v1,v2);
 v1->e0 = e;
 v2->e0 = e;
 E.appendHead(e);

 return e;
}


//////////////////// Creates an edge ////////////////////////////

Edge *Triangulation::CreateEdge(ExtVertex *v1, ExtVertex *v2, const bool check)
{
 Edge *e;
 Node *n;

 if (check)
  FOREACHVEEDGE((&(v1->VE)), e, n)
   if (e->oppositeVertex(v1->v) == v2->v) return e;

 e = new Edge(v1->v,v2->v);
 if (v1->v->e0 == NULL) v1->v->e0 = e;
 if (v2->v->e0 == NULL) v2->v->e0 = e;
 v1->VE.appendHead(e);
 v2->VE.appendHead(e);
 E.appendHead(e);

 return e;
}


///////////////////// Creates a triangle //////////////////////////

Triangle *Triangulation::CreateTriangle(Edge *e1, Edge *e2, Edge *e3)
{
 Triangle *tt, **at1, **at2, **at3;
 
 if (e1->commonVertex(e2) == e1->v2 && e1->t1 == NULL) at1 = &(e1->t1);
 else if (e1->commonVertex(e2) == e1->v1 && e1->t2 == NULL) at1 = &(e1->t2);
 else return NULL;
 if (e2->commonVertex(e3) == e2->v2 && e2->t1 == NULL) at2 = &(e2->t1);
 else if (e2->commonVertex(e3) == e2->v1 && e2->t2 == NULL) at2 = &(e2->t2);
 else return NULL;
 if (e3->commonVertex(e1) == e3->v2 && e3->t1 == NULL) at3 = &(e3->t1);
 else if (e3->commonVertex(e1) == e3->v1 && e3->t2 == NULL) at3 = &(e3->t2);
 else return NULL;

 tt = new Triangle(e1,e2,e3);
 *at1 = *at2 = *at3 = tt;
 T.appendHead(tt);

 d_boundaries = d_handles = d_shells = 1;

 return tt;
}


///////////////////// Creates an unoriented triangle //////////////////////////

Triangle *Triangulation::CreateUnorientedTriangle(Edge *e1, Edge *e2, Edge *e3)
{
 Triangle *tt, **at1, **at2, **at3;
 
 if (e1->t1 == NULL) at1 = &(e1->t1);
 else if (e1->t2 == NULL) at1 = &(e1->t2);
 else return NULL;
 if (e2->t1 == NULL) at2 = &(e2->t1);
 else if (e2->t2 == NULL) at2 = &(e2->t2);
 else return NULL;
 if (e3->t1 == NULL) at3 = &(e3->t1);
 else if (e3->t2 == NULL) at3 = &(e3->t2);
 else return NULL;

 tt = new Triangle(e1,e2,e3);
 *at1 = *at2 = *at3 = tt;
 T.appendHead(tt);

 return tt;
}


////////////// Euler operatior: Create edge and triangle //////////////////////

Triangle *Triangulation::EulerEdgeTriangle(Edge *e2, Edge *e3)
{
 Vertex *cv = e2->commonVertex(e3);
 Triangle *adj = (e2->t1 == NULL)?(e2->t2):(e2->t1);
 if (cv == NULL || !e2->isOnBoundary() || !e3->isOnBoundary()) return NULL;

 Edge *e1 = CreateEdge(e2->oppositeVertex(cv), e3->oppositeVertex(cv));
 if (adj->nextEdge(e2)->hasVertex(cv)) return CreateTriangle(e1,e3,e2);
 return CreateTriangle(e1,e2,e3);
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    P R I M I T I V E   D E S T R U C T I O N                 //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Unlinks the triangle (elements are not removed from the lists) //////////

void Triangulation::unlinkTriangle(Triangle *t)
{
 Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
 Edge *e1 = t->e1, *e2 = t->e2, *e3 = t->e3;

 int v1nm = (v1->isOnBoundary() && !e1->isOnBoundary() && !e2->isOnBoundary());
 int v2nm = (v2->isOnBoundary() && !e2->isOnBoundary() && !e3->isOnBoundary());
 int v3nm = (v3->isOnBoundary() && !e3->isOnBoundary() && !e1->isOnBoundary());

 v1->e0 = ((e2->isOnBoundary())?(e1):(e2));
 v2->e0 = ((e3->isOnBoundary())?(e2):(e3));
 v3->e0 = ((e1->isOnBoundary())?(e3):(e1));

 e1->replaceTriangle(t, NULL);
 e2->replaceTriangle(t, NULL);
 e3->replaceTriangle(t, NULL);

 if (e1->isIsolated() && e2->isIsolated()) v1->e0 = NULL;
 if (e2->isIsolated() && e3->isIsolated()) v2->e0 = NULL;
 if (e3->isIsolated() && e1->isIsolated()) v3->e0 = NULL;
 if (e1->isIsolated()) e1->v1 = e1->v2 = NULL;
 if (e2->isIsolated()) e2->v1 = e2->v2 = NULL;
 if (e3->isIsolated()) e3->v1 = e3->v2 = NULL;
 t->e1 = t->e2 = t->e3 = NULL;

 Vertex *nv;
 Edge *e;
 List *ve;
 Node *n;

 if (v1nm)
 {
  nv = new Vertex(v1->x, v1->y, v1->z);
  nv->e0 = v1->e0;
  ve = v1->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v1, nv);
  delete(ve);
  v1->e0 = e1;
  V.appendHead(nv);
 }
 if (v2nm)
 {
  nv = new Vertex(v2->x, v2->y, v2->z);
  nv->e0 = v2->e0;
  ve = v2->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v2, nv);
  delete(ve);
  v2->e0 = e2;
  V.appendHead(nv);
 }
 if (v3nm)
 {
  nv = new Vertex(v3->x, v3->y, v3->z);
  nv->e0 = v3->e0;
  ve = v3->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v3, nv);
  delete(ve);
  v3->e0 = e3;
  V.appendHead(nv);
 }
}


//// Unlinks the triangle (elements are not removed from the lists)              ////
//// Differently from the above method, non-manifold vertices are not duplicated ////

void Triangulation::unlinkTriangleNoManifold(Triangle *t)
{
 Edge *e1 = t->e1, *e2 = t->e2, *e3 = t->e3;

 e1->replaceTriangle(t, NULL);
 e2->replaceTriangle(t, NULL);
 e3->replaceTriangle(t, NULL);

 if (e1->isIsolated()) e1->v1 = e1->v2 = NULL;
 if (e2->isIsolated()) e2->v1 = e2->v2 = NULL;
 if (e3->isIsolated()) e3->v1 = e3->v2 = NULL;
 t->e1 = t->e2 = t->e3 = NULL;
}


///// Removes all the triangles with NULL edges /////

int Triangulation::removeTriangles()
{
 Node *n;
 Triangle *t;
 int r = 0;

 n = T.head();
 do
 {
  t = (Triangle *)n->data;
  n = n->next();
  if (t->e1 == NULL || t->e2 == NULL || t->e3 == NULL)
  {
   r++;
   T.removeCell((n!=NULL)?(n->prev()):T.tail());
   delete t;
  }
 } while (n != NULL);

 d_boundaries = d_handles = d_shells = 1;

 return r;
}


///// Removes all the edges with NULL vertices /////

int Triangulation::removeEdges()
{
 Node *n;
 Edge *e;
 int r = 0;

 n = E.head();
 do
 {
  e = (Edge *)n->data;
  n = n->next();
  if (e->v1 == NULL || e->v2 == NULL)
  {
   r++;
   E.removeCell((n!=NULL)?(n->prev()):E.tail());
   delete e;
  }
 } while (n != NULL);

 d_boundaries = d_handles = d_shells = 1;

 return r;
}


/////////// Removes all the vertices with e0 field = NULL ////////////

int Triangulation::removeVertices()
{
 Node *n;
 Vertex *v;
 int r = 0;

 n = V.head();
 do
 {
  v = (Vertex *)n->data;
  n = n->next();
  if (v->e0 == NULL)
  {
   r++;
   V.removeCell((n!=NULL)?(n->prev()):V.tail());
   delete v;
  }
 } while (n != NULL);

 d_boundaries = d_handles = d_shells = 1;

 return r;
}



//////////////////////////////////////////////////////////////////
//                                                              //
//    S E L E C T I O N   M A N A G E M E N T                   //
//                                                              //
//////////////////////////////////////////////////////////////////


//////////////// Deselect all the triangles /////////////

void Triangulation::deselectTriangles()
{
 Triangle *t;
 Node *n;
 FOREACHTRIANGLE(t, n) UNMARK_VISIT(t);
}


//////// Removes all the selected (IS_VISITED) triangles /////////////

void Triangulation::removeSelectedTriangles()
{
 Node *n;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) unlinkTriangle(t);
 removeUnlinkedElements();
}


// Mark all the triangles having at least one boundary vertex as 'selected'

int Triangulation::selectBoundaryTriangles()
{
 Node *n;
 Edge *e;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;
 int ns=0;

 FOREACHEDGE(e, n) if (e->isOnBoundary()) {MARK_VISIT(e->v1); MARK_VISIT(e->v2);}

 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) {MARK_VISIT(t); ns++;}
 }
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 return ns;
}


// Grows the current selection (1 triangle width)

int Triangulation::growSelection()
{
 Node *n;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;
 int ns=0;

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  MARK_VISIT(v1); MARK_VISIT(v2); MARK_VISIT(v3);
 }
 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) {MARK_VISIT(t); ns++;}
 }

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 return ns;
}


// Shrinks the current selection (1 triangle width)

void Triangulation::shrinkSelection()
{
 Node *n;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  MARK_VISIT(v1); MARK_VISIT(v2); MARK_VISIT(v3);
 }
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) UNMARK_VISIT(t);
 }

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
}


// Toggles the selection status of the triangles

void Triangulation::invertSelection(Triangle *t0)
{
 Node *n;
 Triangle *t;

 if (t0 != NULL)
 {
  List totoggle(t0);
  Triangle *s;
  bool unmark = IS_VISITED(t0);
  if (unmark) UNMARK_VISIT(t0); else MARK_VISIT(t0);
  while ((t = (Triangle *)totoggle.popHead()) != NULL)
  {
   if ((s = t->t1()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
   if ((s = t->t2()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
   if ((s = t->t3()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
  }
 }
 else
 {
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) UNMARK_VISIT(t); else MARK_VISIT(t);
 }
}


void Triangulation::reselectSelection(Triangle *t0)
{
 if (!IS_VISITED(t0)) return;

 Node *n;
 Triangle *t, *s;
 List triList(t0);
 MARK_VISIT2(t0);

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  if ((s = t->t1()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
 }

 FOREACHTRIANGLE(t, n) if (!IS_VISITED2(t)) UNMARK_VISIT(t); else UNMARK_VISIT2(t);
}

// Creates a new mesh out of a selection.

Triangulation *Triangulation::createSubMeshFromSelection(Triangle *t0, bool keep_ref)
{
 Triangle *t,*s, *nt;
 Node *n;

 if (t0 != NULL && !IS_VISITED(t0)) return NULL;

 Triangulation *tin = new Triangulation;
 Vertex *v,*nv;
 Edge *e, *ne;
 List triList, sT, sE, sV;

 if (t0 != NULL)
 {
  triList.appendHead(t0); MARK_BIT(t0,3);
  while(triList.numels())
  {
   t = (Triangle *)triList.popHead();
   sT.appendHead(t);
   if (!IS_BIT(t->e1, 3)) {sE.appendHead(t->e1); MARK_BIT(t->e1, 3);}
   if (!IS_BIT(t->e2, 3)) {sE.appendHead(t->e2); MARK_BIT(t->e2, 3);}
   if (!IS_BIT(t->e3, 3)) {sE.appendHead(t->e3); MARK_BIT(t->e3, 3);}
   if ((v=t->v1()) && !IS_BIT(v, 3)) {sV.appendHead(v); MARK_BIT(v, 3);}
   if ((v=t->v2()) && !IS_BIT(v, 3)) {sV.appendHead(v); MARK_BIT(v, 3);}
   if ((v=t->v3()) && !IS_BIT(v, 3)) {sV.appendHead(v); MARK_BIT(v, 3);}
   if ((s = t->t1()) != NULL && !IS_BIT(s,3) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,3);}
   if ((s = t->t2()) != NULL && !IS_BIT(s,3) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,3);}
   if ((s = t->t3()) != NULL && !IS_BIT(s,3) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,3);}
  }
 }
 else
 {
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
   {sT.appendHead(t); MARK_BIT(t->e1, 3); MARK_BIT(t->e2, 3); MARK_BIT(t->e3, 3);}
  FOREACHEDGE(e, n) if (IS_BIT(e,3))
   {sE.appendHead(e); MARK_BIT(e->v1, 3); MARK_BIT(e->v2, 3);}
  FOREACHVERTEX(v, n) if (IS_BIT(v,3)) sV.appendHead(v);
 }

 FOREACHVEEDGE((&sE), e, n) e->v1->e0 = e->v2->e0 = e;

 int i;
 void **v_info = NULL, **e_info = NULL, **t_info = NULL;
 if (!keep_ref)
 {
  v_info = new void *[sV.numels()];
  i=0; FOREACHVVVERTEX((&sV), v, n) v_info[i++] = v->info;
  e_info = new void *[sE.numels()];
  i=0; FOREACHVEEDGE((&sE), e, n) e_info[i++] = e->info;
  t_info = new void *[sT.numels()];
  i=0; FOREACHVTTRIANGLE((&sT), t, n) t_info[i++] = t->info;
 }

 FOREACHVVVERTEX((&sV), v, n)
  {nv=new Vertex(v); tin->V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&sE), e, n)
  {ne=new Edge((Vertex *)e->v1->info, (Vertex *)e->v2->info); tin->E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&sT), t, n)
  {nt=new Triangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); tin->T.appendTail(nt); t->info = nt; nt->info = t;}

 FOREACHVVVERTEX((&sV), v, n) ((Vertex *)v->info)->e0 = (Edge *)v->e0->info;

 FOREACHVEEDGE((&sE), e, n)
 {
  ((Edge *)e->info)->t1 = (e->t1 && IS_VISITED(e->t1))?((Triangle *)e->t1->info):(NULL);
  ((Edge *)e->info)->t2 = (e->t2 && IS_VISITED(e->t2))?((Triangle *)e->t2->info):(NULL);
 }

 i=0; if (!keep_ref) FOREACHVVVERTEX((&sV), v, n) v->info = v_info[i++];
 i=0; if (!keep_ref) FOREACHVEEDGE((&sE), e, n) e->info = e_info[i++];
 i=0; if (!keep_ref) FOREACHVTTRIANGLE((&sT), t, n) t->info = t_info[i++];

 FOREACHVTTRIANGLE((&sT), t, n) UNMARK_BIT(t, 3);
 FOREACHVEEDGE((&sE), e, n) UNMARK_BIT(e, 3);
 FOREACHVVVERTEX((&sV), v, n) UNMARK_BIT(v, 3);

 if (!sT.numels()) {delete(tin); return NULL;}

 tin->duplicateNonManifoldVertices();
 tin->eulerUpdate();

 return tin;
}


///// Marks all the triangles within distance L as selected //////

int Triangulation::selectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;
 int nt=0;

 FOREACHVTTRIANGLE(reg, s, n) {MARK_VISIT(s); nt++;}
 delete(reg);

 return nt;
}


///// Deselects all the triangles within distance L //////

int Triangulation::deselectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;
 int nt=0;

 FOREACHVTTRIANGLE(reg, s, n) {UNMARK_VISIT(s); nt++;}
 delete(reg);

 return nt;
}


///// Selects all the triangles within distance L which were already //////
///// selected. Deselects the others.                                //////

void Triangulation::reselectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;

 FOREACHVTTRIANGLE(reg, s, n) MARK_VISIT2(s);
 FOREACHTRIANGLE(s, n) if (IS_VISITED(s) && !IS_VISITED2(s)) UNMARK_VISIT(s);
 FOREACHVTTRIANGLE(reg, s, n) UNMARK_VISIT2(s);
 delete(reg);
}


//// Remove all the selected triangles and re-triangulate /////

bool Triangulation::retriangulateSelectedRegion()
{
 List ttbr;
 Node *n;
 Triangle *u;
 Point nor;
 FOREACHTRIANGLE(u, n) if (IS_VISITED(u))
  {ttbr.appendHead(u); nor = nor+(u->getNormal()*u->area());}

 if (ttbr.numels() < 2)
 {
  JMesh::warning("retriangulateRegion: Nothing to retriangulate.\n");
  return 0;
 }

 FOREACHVTTRIANGLE((&(ttbr)), u, n)
  if (u->getNormal()*nor <= 0.0)
  {
   JMesh::warning("retriangulateRegion: Too complex geometry. Can't retriangulate.\n");
   return 0;
  }

 if (!isSelectionSimple(&ttbr))
 {
  JMesh::warning("retriangulateRegion: Non-simple region. Can't retriangulate.\n");
  return 0;
 }

 List *ms = getRegionInternalVertices(&ttbr);

 FOREACHVTTRIANGLE((&(ttbr)), u, n) unlinkTriangle(u);
 Edge *e = ((Edge *)ms->head()->data);
 List *vl = ((List *)ms->head()->next()->data);
 TriangulateHole(e, vl);
 delete(vl);
 delete(ms);
 removeUnlinkedElements();

 return 1;
}


////////// Check wether 's' represents a simple selection ////////

bool Triangulation::isSelectionSimple(List *s)
{
 if (!s->numels()) return 0; // Empty region is not simple

 Node *n;
 Triangle *ta, *t = (Triangle *)s->head()->data;
 List bdr, top(t);
 MARK_VISIT2(t);
 int nv=0;

 while (top.numels())
 {
  t = (Triangle *)top.popHead();
  nv++;
  ta=t->t1(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e1);
  ta=t->t2(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e2);
  ta=t->t3(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e3);
 }

 FOREACHVTTRIANGLE(s, t, n) UNMARK_VISIT2(t);
 if (top.numels()) return 0; // Mesh-boundary in selection
 if (nv != s->numels()) return 0; // Disconnected selection

 Edge *e, *f, *ge=NULL, *e0;
 List *ve;
 FOREACHVEEDGE((&(bdr)), e, n) MARK_VISIT(e);
 int nae;

 nv = 0;
 e = e0 = (Edge *)bdr.head()->data;
 Vertex *v = e->v1;

 do
 {
  nv++;
  v = e->oppositeVertex(v);
  ve = v->VE();
  nae=0; FOREACHVEEDGE(ve, f, n) if (f!=e && IS_VISITED(f)) {ge=f; nae++;}
  delete(ve);
  if (nae > 1) break;
  e=ge;
 } while (e != e0);

 FOREACHVEEDGE((&(bdr)), e, n) UNMARK_VISIT(e);
 if (nv != bdr.numels()) return 0; // Non-simple selection

 return 1;
}


//// Unmarks all the elements of the triangulation ////
//// but leaves selected triangles marked.         ////

void Triangulation::unmarkEverythingButSelections()
{
 Vertex *v;
 Edge *e;
 Triangle *t;
 Node *n;
 FOREACHVERTEX(v, n) v->mask = 0;
 FOREACHEDGE(e, n) e->mask = 0;
 FOREACHTRIANGLE(t, n) t->mask &= (unsigned char)1;
}


//// Selects all the triangles of the shell containing 't' ////

int Triangulation::selectConnectedComponent(Triangle *t0, bool sos)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;
 int ns = 0;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (!IS_VISITED(t))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && !IS_VISITED(t1) && (!(sos && IS_SHARPEDGE(t->e1)))) todo.appendHead(t1);
   if (t2 != NULL && !IS_VISITED(t2) && (!(sos && IS_SHARPEDGE(t->e2)))) todo.appendHead(t2);
   if (t3 != NULL && !IS_VISITED(t3) && (!(sos && IS_SHARPEDGE(t->e3)))) todo.appendHead(t3);

   MARK_VISIT(t); ns++;
  }
 }

 return ns;
}


//// Deselects all the triangles of the shell containing 't' ////

int Triangulation::deselectConnectedComponent(Triangle *t0, bool sos)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;
 int ns = 0;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (IS_VISITED(t))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && IS_VISITED(t1) && (!(sos && IS_SHARPEDGE(t->e1)))) todo.appendHead(t1);
   if (t2 != NULL && IS_VISITED(t2) && (!(sos && IS_SHARPEDGE(t->e2)))) todo.appendHead(t2);
   if (t3 != NULL && IS_VISITED(t3) && (!(sos && IS_SHARPEDGE(t->e3)))) todo.appendHead(t3);

   UNMARK_VISIT(t); ns++;
  }
 }

 return ns;
}


// Append to the current mesh a copy of all the elements of 'src'.
// The newly created elements form a new selection.

void Triangulation::append(Triangulation *src)
{
 deselectTriangles();
 Triangulation cb(src);
 cb.invertSelection();
 V.joinTailList(&(cb.V));
 E.joinTailList(&(cb.E));
 T.joinTailList(&(cb.T));
 d_boundaries = d_handles = d_shells = 1;
}

//////////////////////////////////////////////////////////////////
//                                                              //
//    R E G I O N   M A N I P U L A T I O N                     //
//                                                              //
//////////////////////////////////////////////////////////////////

///// Make a list with all the triangles within distance L //////

List *Triangulation::getRegion(Triangle *t, const double L, const Point *center)
{
 List triList, *toRemove = new List;
 if (t->v1()->distance(center) > L) return toRemove;
 if (t->v2()->distance(center) > L) return toRemove;
 if (t->v3()->distance(center) > L) return toRemove;

 Triangle *s;
 Node *n;

 triList.appendHead(t);
 MARK_BIT(t,3);

 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.head()->data;
  triList.removeCell(triList.head());
  toRemove->appendHead(t);

  if ((s = t->t1()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e1)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t2()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e2)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t3()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e3)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
 }

 FOREACHVTTRIANGLE(toRemove, s, n) UNMARK_BIT(s, 3);

 return toRemove;
}


///// Unlink all the triangles within distance L //////

void Triangulation::removeRegion(Triangle *t, const double L, const Point *center)
{
 List triList, toRemove;
 Node *n;
 Triangle *s;

 triList.appendHead(t);
 MARK_VISIT(t);

 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.head()->data;
  triList.removeCell(triList.head());
  toRemove.appendHead(t);

  if ((s = t->t1()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e1)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e2)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e3)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
 }

 for (n = toRemove.tail(); n != NULL; n=n->prev())
 {
  s = ((Triangle *)n->data);
  unlinkTriangle(s);
 }

 removeUnlinkedElements();
}


//////// Next region's boundary vertex /////////////

Vertex *Triangulation::nextVertexOnRegionBoundary(Vertex *sv) const
{
 Triangle *lt, *rt;
 Edge *e;
 List *ve = sv->VE();
 Node *n;

 FOREACHVEEDGE(ve, e, n)
 {
  lt = e->leftTriangle(sv);
  rt = e->rightTriangle(sv);
  if (lt != NULL && IS_VISITED(lt) && (rt == NULL || !IS_VISITED(rt)))
   {delete(ve); return e->oppositeVertex(sv);}
 }
 delete(ve);

 return NULL;
}


//// This method returns a list containing an edge of the region's boundary ////
//// as its first element, and all the internal vertices as the remaining   ////

List *Triangulation::getRegionInternalVertices(List *reg)
{
 List *iVertices = new List;
 List *outList = new List;
 Edge *bEdge = NULL;
 Triangle *s, *t;
 Node *n;
 Vertex *v1, *v2, *v3;

 FOREACHVTTRIANGLE(reg, t, n) {MARK_VISIT(t); MARK_BIT(t, 3);}

 FOREACHVTTRIANGLE(reg, t, n)
 {
  if (IS_BIT(t,3))
  {
   UNMARK_BIT(t,3);
   if ((s = t->t1()) != NULL && !IS_VISITED(s)) {bEdge = t->e1; MARK_BIT(t->e1->v1, 3); MARK_BIT(t->e1->v2, 3);}
   if ((s = t->t2()) != NULL && !IS_VISITED(s)) {bEdge = t->e2; MARK_BIT(t->e2->v1, 3); MARK_BIT(t->e2->v2, 3);}
   if ((s = t->t3()) != NULL && !IS_VISITED(s)) {bEdge = t->e3; MARK_BIT(t->e3->v1, 3); MARK_BIT(t->e3->v2, 3);}
  }
 }

 FOREACHVTTRIANGLE(reg, s, n)
 {
  v1 = s->v1(); v2 = s->v2(); v3 = s->v3();
  if (!IS_BIT(v1, 3)) {iVertices->appendHead(v1); MARK_BIT(v1, 3);}
  if (!IS_BIT(v2, 3)) {iVertices->appendHead(v2); MARK_BIT(v2, 3);}
  if (!IS_BIT(v3, 3)) {iVertices->appendHead(v3); MARK_BIT(v3, 3);}
 }
 FOREACHVTTRIANGLE(reg, s, n)
 {
  v1 = s->v1(); v2 = s->v2(); v3 = s->v3();
  UNMARK_BIT(v1, 3); UNMARK_BIT(v2, 3); UNMARK_BIT(v3, 3); 
 }

 outList->appendHead(iVertices);
 outList->appendHead(bEdge);

 return outList;
}


// Transforms only the shell indicated by 't0'

void Triangulation::transformShell(Triangle *t0, const Matrix4x4& m)
{
 List todo(t0), st, sv;
 Triangle *t, *nt;
 Vertex *v;
 double x, y, z, w;

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  st.appendHead(t);
  nt=t->t1(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
  nt=t->t2(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
  nt=t->t3(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
 }

 while (st.numels())
 {
  t = (Triangle *)st.popHead();
  UNMARK_VISIT(t);
  v = t->v1(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
  v = t->v2(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
  v = t->v3(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
 }

 while (sv.numels())
 {
  v = (Vertex *)sv.popHead();
  UNMARK_VISIT(v);
  x = ((*v)*Point(m.matrix[0][0],m.matrix[1][0],m.matrix[2][0]))+m.matrix[3][0];
  y = ((*v)*Point(m.matrix[0][1],m.matrix[1][1],m.matrix[2][1]))+m.matrix[3][1];
  z = ((*v)*Point(m.matrix[0][2],m.matrix[1][2],m.matrix[2][2]))+m.matrix[3][2];
  w = ((*v)*Point(m.matrix[0][3],m.matrix[1][3],m.matrix[2][3]))+m.matrix[3][3];
  v->x = x/w; v->y = y/w; v->z = z/w; 
 }
}


//// Removes all the triangles of the shell containing 't0' ////

void Triangulation::removeShell(Triangle *t0)
{
 List todo(t0);
 Triangle *t, *t1, *t2, *t3;

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

  if (t1 != NULL && !IS_VISITED2(t1)) {MARK_VISIT2(t1); todo.appendHead(t1);}
  if (t2 != NULL && !IS_VISITED2(t2)) {MARK_VISIT2(t2); todo.appendHead(t2);}
  if (t3 != NULL && !IS_VISITED2(t3)) {MARK_VISIT2(t3); todo.appendHead(t3);}

  unlinkTriangle(t);
 }

 removeUnlinkedElements();
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    G L O B A L   O P E R A T I O N S                         //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Tags as sharp all the edges exceeding the given curvature ///////////

void Triangulation::sharpEdgeTagging(const double ta)
{
 Node *n;
 Edge *e;
 FOREACHEDGE(e, n)
  if (e->curvature() > ta) TAG_SHARPEDGE(e);
  else UNTAG_SHARPEDGE(e);
}


//// Unmarks all the elements of the triangulation ////

void Triangulation::unmarkEverything()
{
 Vertex *v;
 Edge *e;
 Triangle *t;
 Node *n;
 FOREACHVERTEX(v, n) v->mask = 0;
 FOREACHEDGE(e, n) e->mask = 0;
 FOREACHTRIANGLE(t, n) t->mask = 0;
}


///// Compute the bounding box and return its max edge /////

double Triangulation::getBoundingBox(Point& mp, Point& Mp) const
{
 Vertex *v; Node *n;
 Mp.x = -DBL_MAX, mp.x = DBL_MAX;
 Mp.y = -DBL_MAX, mp.y = DBL_MAX;
 Mp.z = -DBL_MAX, mp.z = DBL_MAX;
 FOREACHVERTEX(v, n)
 {
  if (v->x < mp.x) mp.x = v->x;
  if (v->x > Mp.x) Mp.x = v->x;
  if (v->y < mp.y) mp.y = v->y;
  if (v->y > Mp.y) Mp.y = v->y;
  if (v->z < mp.z) mp.z = v->z;
  if (v->z > Mp.z) Mp.z = v->z;
 }

 return MAX(Mp.x-mp.x,MAX(Mp.y-mp.y,Mp.z-mp.z));
}


///// Compute the approximate bounding ball radius /////

double Triangulation::getBoundingBallRadius() const
{
 Vertex *v; Node *n;
 Point tc, mp, Mp;
 double tb, bsr = getBoundingBox(mp, Mp)/2;
 Point bsc = (Mp+mp)/2;

 FOREACHVERTEX(v, n)
  if ((tb = ((*v)-bsc).length()) > bsr)
  {
   tc = ((*v)-bsc); tc.normalize();
   tb = ((tb-bsr)/2);
   bsc = bsc+(tc*tb);
   bsr += tb;
  }

 return bsr;
}


////// Returns the surface area of the mesh //////

double Triangulation::area() const
{
 Triangle *t;
 Node *n;
 double a=0.0;
 FOREACHTRIANGLE(t, n) a += t->area();

 return a;
}


////// Returns the volume of the mesh //////

double Triangulation::volume() const
{
 Triangle *t;
 Node *n;
 double v=0.0;
 FOREACHTRIANGLE(t, n)
  v += (t->getCenter()*t->getNormal())*t->area();

 return v/3;
}


///// Places the mesh into the unit cube by translating and resizing /////
///// so that all the coordinates are between 0 and mc (default =1). /////

void Triangulation::normalize(double mc)
{
 Vertex *v;
 Node *n;
 Point mp, Mp;
 double mel = getBoundingBox(mp, Mp)/mc;
 FOREACHVERTEX(v, n) v->setValue(((*v)-mp)/mel);	// Shift and normalize
}


/////// Transforms all the vertices using the 4x4 matrix 'm' ///////////

void Triangulation::transform(const Matrix4x4& m)
{
 Node *n;
 Vertex *v;
 double x,y,z,w;

 FOREACHVERTEX(v, n)
 {
  x = ((*v)*Point(m.matrix[0][0],m.matrix[1][0],m.matrix[2][0]))+m.matrix[3][0];
  y = ((*v)*Point(m.matrix[0][1],m.matrix[1][1],m.matrix[2][1]))+m.matrix[3][1];
  z = ((*v)*Point(m.matrix[0][2],m.matrix[1][2],m.matrix[2][2]))+m.matrix[3][2];
  w = ((*v)*Point(m.matrix[0][3],m.matrix[1][3],m.matrix[2][3]))+m.matrix[3][3];
  v->x = x/w; v->y = y/w; v->z = z/w; 
 }
}


// Add noise in the normal direction. Normal displacement is
// bounded by ns% of the bounding ball radius.

void Triangulation::addNormalNoise(double ns)
{
 Vertex *v;
 Node *n;
 Point np;
 int i;
 double noise;
 coord *xyz = (coord *)malloc(sizeof(coord)*V.numels()*3);
 ns *= (getBoundingBallRadius()/100.0);

 i=0; FOREACHVERTEX(v, n)
 {
  noise = ns*(((((double)rand()))-(((double)RAND_MAX)/2.0))/((double)RAND_MAX));
  np = (*v)+((v->getNormal())*noise);
  xyz[i++]=np.x; xyz[i++]=np.y; xyz[i++]=np.z;
 }
 i=0; FOREACHVERTEX(v, n) {v->x=xyz[i++]; v->y=xyz[i++]; v->z=xyz[i++];}

 free(xyz);
}


// Iteratively swap edges to maximize the minimum angle.
// Checks and avoids normal inversion.

bool Triangulation::iterativeEdgeSwaps()
{
 Node *n;
 Edge *e, *f;
 double l;
 int swaps=1, totits=1;
 Point n1, n2, nor;
 List toswap;

 bool selection=0;
 Triangle *t;
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {selection=1; break;}

 FOREACHEDGE(e, n) if (!IS_SHARPEDGE(e) && !e->isOnBoundary())
 {
  MARK_VISIT(e); if ((!selection || (IS_VISITED(e->t1) && IS_VISITED(e->t2)))) toswap.appendTail(e);
 }

 JMesh::begin_progress();
 while (swaps && totits++ < 10)
 {
  swaps = 0; for (n=toswap.head(); n!=NULL; )
  {
   e = (Edge *)n->data;
   if (n==toswap.tail()) {toswap.removeCell(toswap.tail()); n=NULL;}
   else {n=n->next(); toswap.removeCell(n->prev());}
   UNMARK_VISIT(e);
   if (!e->t1->isNeedle() && !e->t2->isNeedle())
   {
    n1 = e->t1->getNormal();
    n2 = e->t2->getNormal();
    nor = n1+n2;
    l = e->delaunayMinAngle();
    if (e->swap())
    {
     if (e->t1->isNeedle() || e->t2->isNeedle() || e->delaunayMinAngle() <= l*1.000001 || nor*e->t1->getNormal() <= 0 || nor*e->t2->getNormal() <= 0) e->swap(1);
     else
     {
      swaps++;
      f = e->t1->nextEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t1->prevEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t2->nextEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t2->prevEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
     }
    }
   }
  }
  JMesh::report_progress("Swaps: %d      ", swaps);
 }
 JMesh::end_progress();

 FOREACHEDGE(e, n) UNMARK_VISIT(e);

 if (totits >= 10)
 {
  JMesh::warning("Optimization did not converge after 10 iterations! Stopping.\n");
  JMesh::warning("You may try to run the method again.\n");
  return 0;
 }

 return 1;
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    T O P O L O G Y   M A N I P U L A T I O N                 //
//                                                              //
//////////////////////////////////////////////////////////////////


////// Invertes all the triangle normals and edge orientations ///////

void Triangulation::flipNormals()
{
 Node *n;
 Edge *e;
 Triangle *t;

 FOREACHTRIANGLE(t, n) t->invert();
 FOREACHEDGE(e, n) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
}


//// Marks all the triangles of the shell containing 't' ////

void Triangulation::flipNormals(Triangle *t0)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (!IS_BIT(t,2))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && !IS_BIT(t1,2)) todo.appendHead(t1);
   if (t2 != NULL && !IS_BIT(t2,2)) todo.appendHead(t2);
   if (t3 != NULL && !IS_BIT(t3,2)) todo.appendHead(t3);

   t->invert();
   if (!IS_BIT(t->e1,2)) p_swap((void **)(&(t->e1->v1)), (void **)(&(t->e1->v2)));
   if (!IS_BIT(t->e2,2)) p_swap((void **)(&(t->e2->v1)), (void **)(&(t->e2->v2)));
   if (!IS_BIT(t->e3,2)) p_swap((void **)(&(t->e3->v1)), (void **)(&(t->e3->v2)));
   MARK_BIT(t->e1,2); MARK_BIT(t->e2,2); MARK_BIT(t->e3,2);
   MARK_BIT(t,2);
  }
 }

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (IS_BIT(t,2))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && IS_BIT(t1,2)) todo.appendHead(t1);
   if (t2 != NULL && IS_BIT(t2,2)) todo.appendHead(t2);
   if (t3 != NULL && IS_BIT(t3,2)) todo.appendHead(t3);

   UNMARK_BIT(t->e1,2); UNMARK_BIT(t->e2,2); UNMARK_BIT(t->e3,2);
   UNMARK_BIT(t,2);
  }
 }
}


//// Returns the top triangle of the mesh (max. z) ////

Triangle *Triangulation::topTriangle(Triangle *t0)
{
 Node *n;
 Vertex *v, *hv = NULL, *v1, *v2, *v3;
 Edge *e, *fe = NULL;
 coord az, Mz = -COORD_MAX;
 Triangle *t, *t1, *t2, *t3;
 List *ve, todo, tlist, elist, vlist;

 todo.appendHead(t0); MARK_BIT(t0,2);

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead(); tlist.appendHead(t);
  t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (!IS_VISITED(v1)) {MARK_VISIT(v1); vlist.appendHead(v1);}
  if (!IS_VISITED(v2)) {MARK_VISIT(v2); vlist.appendHead(v2);}
  if (!IS_VISITED(v3)) {MARK_VISIT(v3); vlist.appendHead(v3);}

  if (!IS_VISITED(t->e1)) {MARK_VISIT(t->e1); elist.appendHead(t->e1);}
  if (!IS_VISITED(t->e2)) {MARK_VISIT(t->e2); elist.appendHead(t->e2);}
  if (!IS_VISITED(t->e3)) {MARK_VISIT(t->e3); elist.appendHead(t->e3);}

  if (t1 != NULL && !IS_BIT(t1,2)) {MARK_BIT(t1,2); todo.appendHead(t1);}
  if (t2 != NULL && !IS_BIT(t2,2)) {MARK_BIT(t2,2); todo.appendHead(t2);}
  if (t3 != NULL && !IS_BIT(t3,2)) {MARK_BIT(t3,2); todo.appendHead(t3);}
 }

 ve = new List;

 FOREACHVVVERTEX((&(vlist)), v, n) {UNMARK_VISIT(v); if ((az = v->z) > Mz) {Mz=az; hv = v;}}
 Mz = COORD_MAX;
 FOREACHVEEDGE((&(elist)), e, n) {UNMARK_VISIT(e); if (e->hasVertex(hv) && e->length() != 0) ve->appendHead(e);}
 FOREACHVTTRIANGLE((&(tlist)), t, n) UNMARK_BIT(t, 2);

 FOREACHVEEDGE(ve, e, n)
  if ((az = (hv->z - e->oppositeVertex(hv)->z)/e->length()) < Mz) {Mz=az; fe = e;}
 delete(ve);

 if (fe == NULL) fe = hv->e0;
 if (fe->t1 == NULL || fe->t2 == NULL) return NULL;

 return (fabs(fe->t1->getNormal().z) > fabs(fe->t2->getNormal().z))?(fe->t1):(fe->t2);
}


/////// Computes boundaries and handles ///////////

void Triangulation::eulerUpdate()
{
 Vertex *v,*w;
 Edge *e;
 Triangle *t,*s;
 List triList;
 Node *n;
 n_boundaries = n_shells = n_handles = 0;

 FOREACHTRIANGLE(t, n) if (!IS_BIT(t, 3))
 {
  n_shells++;
  triList.appendHead(t);
  MARK_BIT(t,3);

  while(triList.numels())
  {
   t = (Triangle *)triList.popHead();
   if ((s = t->t1()) != NULL && !IS_BIT(s,3)) {triList.appendHead(s); MARK_BIT(s,3);}
   if ((s = t->t2()) != NULL && !IS_BIT(s,3)) {triList.appendHead(s); MARK_BIT(s,3);}
   if ((s = t->t3()) != NULL && !IS_BIT(s,3)) {triList.appendHead(s); MARK_BIT(s,3);}
  }
 }
 FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 3);

 FOREACHEDGE(e, n) if (e->isOnBoundary())
  {MARK_BIT(e->v1, 3); MARK_BIT(e->v2, 3);}

 FOREACHVERTEX(v, n) if (IS_BIT(v,3))
 {
  n_boundaries++;
  for (w=v; IS_BIT(w, 3); w = w->nextOnBoundary()) UNMARK_BIT(w, 3);
 }

 n_handles = (E.numels() - V.numels() - T.numels() + 2*n_shells - n_boundaries)/2;
 d_boundaries = d_handles = d_shells = 0;
}


// Makes the mesh equivalent to a topological disk.
// Edges and vertices are duplicated when necessary.

void Triangulation::openToDisk()
{
 Triangle *t = (Triangle *)T.head()->data;
 Triangle *s;
 List triList, *ve;
 Vertex *v, *w;
 Edge *e, *ne;
 Node *n;
 triList.appendHead(t);
 MARK_BIT(t,3);

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  if ((s = t->t1()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e1,3);}
  if ((s = t->t2()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e2,3);}
  if ((s = t->t3()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e3,3);}
 }
 FOREACHTRIANGLE (t, n) UNMARK_BIT(t, 3);

 FOREACHVERTEX(v, n) v->info = new List;

 FOREACHEDGE(e, n) if (!IS_BIT(e, 3))
 {
  ((List *)e->v1->info)->appendHead(e);
  ((List *)e->v2->info)->appendHead(e);
 }

 FOREACHVERTEX(v, n) if (((List *)v->info)->numels()==1) triList.appendHead(v);
 if (!triList.numels()) JMesh::error("Triangulation::openToDisk: Couldn't find a root.\n");

 while(triList.numels())
 {
  v = (Vertex *)triList.popHead();
  ve = ((List *)v->info);
  if (ve->numels())
  {
   e = (Edge *)(ve->head()->data);
   MARK_BIT(e, 3);
   ve->popHead();
   w = e->oppositeVertex(v);
   ve = ((List *)w->info);
   ve->removeNode(e);
   if (ve->numels() == 1) triList.appendHead(w);
  }
  else
  {
   ve = v->VE();
   e = (Edge *)ve->head()->data; UNMARK_BIT(e, 3); ((List *)v->info)->appendHead(e);
   e = (Edge *)ve->head()->next()->data; UNMARK_BIT(e, 3); ((List *)v->info)->appendHead(e);
   delete(ve);
  }
 }

 FOREACHEDGE(e, n) if (!IS_BIT(e, 3) && !e->isOnBoundary())
 {
  ne = new Edge(e->v1, e->v2);
  ne->t1 = e->t1; e->t1 = NULL; E.appendHead(ne);
  ne->t1->replaceEdge(e, ne);
 }

 FOREACHEDGE(e, n) UNMARK_BIT(e, 3);

 FOREACHVERTEX(v, n) if (v->info) {delete(((List *)v->info)); v->info = NULL;}
 duplicateNonManifoldVertices();
 d_boundaries = d_handles = d_shells = 1;
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    T R I A N G U L A T I O N   M E T H O D S                 //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Computes the center of mass of the hole and star-patches ////////

int Triangulation::StarTriangulateHole(Edge *e)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n;
 Edge *e1, *e2, *e3;
 Point np;
 Vertex *v, *nv, *v1, *v2;
 int nt=0;

 v = e->v1;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);

 FOREACHVVVERTEX((&(bvs)), v, n) np = np+(*v);

 np = np/bvs.numels();
 nv = new Vertex(&np);
 V.appendHead(nv);

 v1 = ((Vertex *)bvs.head()->data);
 e1 = CreateEdge(nv, v1);
 for (n=bvs.head()->next(); n!=NULL; n=n->next())
 {
  v2 = ((Vertex *)n->data);
  e2 = CreateEdge(nv, v2);
  e3 = v1->getEdge(v2);
  CreateTriangle(e1, e2, e3);
  nt++;
  v1 = v2;
  e1 = e2;
 }
 v2 = ((Vertex *)bvs.head()->data);
 e2 = nv->getEdge(v2);
 e3 = v1->getEdge(v2);
 CreateTriangle(e1, e2, e3);
 nt++;

 return nt;
}


///// Patch holes using 2D Delaunay triangulation on the plane 'nor' /////

int Triangulation::TriangulateHole(Edge *e, Point *nor)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0;

 v = e->v1;
 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);

 while (bvs.numels() > 2)
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&(bvs)), v, n)
   if (!IS_VISITED(v) && v->e0 && (ang = v->getAngleOnAveragePlane(nor)) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   JMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&(bvs)), v, n) UNMARK_VISIT(v);
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  if (!EulerEdgeTriangle(e1,e2)) MARK_VISIT(v);
  else {bvs.removeCell(gn); UNMARK_VISIT(v1); UNMARK_VISIT(v2); nt++;}
 }

 int i, skips;
 do
 {
  skips = 0;
  for (n=E.head(), i=2*nt*nt; i<nt; n=n->next(), i--)
  {
   e = ((Edge *)n->data);
   ang = e->delaunayMinAngle();
   if (e->swap())
    {if (e->delaunayMinAngle() <= ang) e->swap(1); else skips++;}
  }
  if (i < 0) {JMesh::warning("Optimization is taking too long. I give up.\n"); break;}
 } while (skips);

 return nt;
}


///// Triangulates a hole using the additional vertices in 'vl' /////

int Triangulation::TriangulateHole(Edge *e, List *vl)
{
 if (!e->isOnBoundary()) return 0;

 List bvs, ovbs, nedg;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0, neb;

 v = e->v1;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);
 ovbs.appendList(&bvs);

 while (bvs.numels() > 2) // While there are more than two boundary vertices
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&(bvs)), v, n)
   if (!IS_VISITED(v) && v->e0 && (ang = v->getAngleForTriangulation()) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   JMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&(bvs)), v, n) UNMARK_VISIT(v);
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  neb = E.numels();
  if (!EulerEdgeTriangle(e1,e2)) MARK_VISIT(v);
  else
  {
   bvs.removeCell(gn);
   UNMARK_VISIT(v1);
   UNMARK_VISIT(v2);
   nt++;
   if (E.numels() > neb) nedg.appendHead(E.head()->data);
  }
 }

// if (nt < 2) return nt;

 // Calcolo una normale per il buco come media dei nuovi triangoli
 int i;
 Point nor;

 for (i=0, n=T.head(); i<nt; i++, n=n->next()) nor = nor+((Triangle *)n->data)->getNormal();
 if (nor.isNull())
 {
  JMesh::warning("TriangulateHole: Unable to compute an average normal. Can't optimize.\n");
  return nt;
 }
 nor.normalize();
 
 // Memorizzo da qualche parte la posizione originale dei vertici
 // e Proietto il boundary sul piano con normale quella appena calcolata
 coord *ovps = (coord *)malloc((ovbs.numels())*3*sizeof(coord));
 int j = 0;
 FOREACHVVVERTEX((&(ovbs)), v, n)
 {
  ovps[j++] = v->x; ovps[j++] = v->y; ovps[j++] = v->z;
  v->project(&nor);
 }

 // Proietto i punti interni sul piano d'appoggio
 Point *p;
 coord *ovpsi = (coord *)malloc((vl->numels())*3*sizeof(coord));
 j = 0;
 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  ovpsi[j++] = p->x; ovpsi[j++] = p->y; ovpsi[j++] = p->z;
 }
 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  p->project(&nor);
 }

 // Ottimizzo secondo Delaunay vincolato al boundary la nuova regione

 int sw;
 do
 {
  sw = 0;
  FOREACHVEEDGE((&(nedg)), e1, n)
  {
   ang = e1->delaunayMinAngle();
   if (e1->swap())
    {if (e1->delaunayMinAngle() <= ang) e1->swap(1); else sw++;}
  }
 } while (sw);

 // Inserisco i punti interni
 int ntt = T.numels()-nt;
 List ivs;

 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  ivs.appendTail(watsonInsert(p, &T, T.numels()-ntt));
 }
 nt = (T.numels() - ntt);

 // Riporto i vertici interni al loro posto
 j=0; FOREACHVVVERTEX((&(ivs)), v, n)
 {
  if (v != NULL)
  {
   v->x = ovpsi[j++]; v->y = ovpsi[j++]; v->z = ovpsi[j++];
  } else j+=3;
 }
 free(ovpsi);

 // Riporto i vertici del boundary al loro posto
 j=0; FOREACHVVVERTEX((&(ovbs)), v, n)
 {
  v->x = ovps[j++]; v->y = ovps[j++]; v->z = ovps[j++];
 }
 free(ovps);

 return nt;
}


////// Inserts the point 'p' in the Delaunay triangulation 'tR' with 'nt' triangles ////

Vertex *Triangulation::watsonInsert(Point *p, List *tR, int nt)
{
 Node *n, *m;
 Edge *e;
 Triangle *t;
 List bdr, bdrs, todo, *ve;
 Vertex *v1, *v2, *v3;
 int i;

 for (i=0, n = T.head(); i<nt; n=n->next(), i++)
 {
  t = ((Triangle *)n->data);
  if (t->e1 != NULL && t->inSphere(p))
  {
   v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
   if (!IS_VISITED(v1)) bdr.appendHead(v1);
   if (!IS_VISITED(v2)) bdr.appendHead(v2);
   if (!IS_VISITED(v3)) bdr.appendHead(v3);
   MARK_VISIT(v1); MARK_VISIT(v2); MARK_VISIT(v3);
   MARK_BIT(t, 3);
   todo.appendHead(t);
  }
 }
 if (bdr.numels() == 0) return NULL;

 FOREACHVVVERTEX((&(bdr)), v1, n)
 {
  ve = v1->VE();
  FOREACHVEEDGE(ve, e, m) if (!IS_BIT(e->t1, 3) || !IS_BIT(e->t2, 3)) v1->e0 = e;
  delete(ve);
 }

 while (todo.numels())
 {
  t = ((Triangle *)todo.head()->data);
  todo.removeCell(todo.head());
  unlinkTriangleNoManifold(t);
 }

 Node *tmp;
 for (i=0, n = T.head(); i<nt; i++)
 {
  t = ((Triangle *)n->data);
  if (t->e1 == NULL) {tmp = n; n=n->next(); T.freeCell(tmp);}
  else n=n->next();
 }

 for (n = bdr.head(); n!=NULL;)
 {
  v1 = ((Vertex *)n->data);
  if (v1->e0 == NULL) {tmp = n; n=n->next(); bdr.removeCell(tmp);}
  else n=n->next();
 }

 v1 = v2 = ((Vertex *)bdr.head()->data);
 do
 {
  bdrs.appendHead(v1);
  v1 = v1->nextOnBoundary();
 } while (v1 != v2);

 Vertex *v = new Vertex(p->x, p->y, p->z);
 V.appendHead(v);

 v1 = ((Vertex *)bdrs.head()->data);
 v->e0 = e = new Edge(v, v1);
 UNMARK_VISIT(v1);
 E.appendHead(e);

 for (n = bdrs.head()->next(); n!=NULL; n=n->next())
 {
  v1 = ((Vertex *)n->data);
  UNMARK_VISIT(v1);
  v2 = ((Vertex *)n->prev()->data);
  e = new Edge(v, v1);
  CreateTriangle(e, v1->getEdge(v2), (Edge *)E.head()->data);
  E.appendHead(e);
 }
 EulerEdgeTriangle(v->e0, (Edge *)E.head()->data);

 return v;
}


//// Removes a vertex and retriangulates its VT ////

int Triangulation::retriangulateVT(Vertex *v)
{
 Point nor;
 Edge *e, *e0 = v->e0->t1->oppositeEdge(v);
 List *vt = v->VT();
 List oe;
 Triangle *t;
 Node *m, *n;
 int i, nt;

 FOREACHVTTRIANGLE(vt, t, m)
 {
  e = t->oppositeEdge(v);
  oe.appendTail(t->prevEdge(e));
  oe.appendTail(e);
  oe.appendTail(t->nextEdge(e));
  nor = nor+t->getNormal();
  unlinkTriangle(t);
 }
 nor.normalize();
 nt = TriangulateHole(e0, &nor);

 for (m=T.head(), i=0; i<nt; i++, m=m->next())
 {
  t = ((Triangle *)m->data);
  if (t->overlaps() || t->isDegenerate()) break;
 }
 if (i<nt)
 {
  JMesh::warning("Re-triangulation failed. Restoring..\n");
  for (m=T.head(), i=0; i<nt; i++, m=m->next())
   unlinkTriangle(((Triangle *)m->data));
  n = oe.head();
  FOREACHVTTRIANGLE(vt, t, m)
  {
   t->e1 = ((Edge *)n->data); n=n->next();
   t->e2 = ((Edge *)n->data); n=n->next();
   t->e3 = ((Edge *)n->data); n=n->next();
   t->e1->v1 = v; t->e1->v2 = (t->e2->t1 == NULL)?(t->e2->v1):(t->e2->v2);
   t->e3->v1 = v; t->e3->v2 = (t->e2->t1 == NULL)?(t->e2->v2):(t->e2->v1);
   ((t->e2->t1 == NULL)?(t->e2->t1):(t->e2->t2)) = t;
   t->e1->t1 = t;
   t->e3->t2 = t;
  }
  v->e0 = ((Triangle *)vt->head()->data)->e1;
 }
 delete(vt);

 return 1;
}


//// Splits the edge 'e' by the point 'p'. ////

Vertex *Triangulation::splitEdge(Edge *e, Point *p, bool copy_mask)
{
 if ((*p)==(*(e->v1))) return e->v1;
 if ((*p)==(*(e->v2))) return e->v2;
 Vertex *v3 = (e->t1 != NULL)?(e->t1->oppositeVertex(e)):(NULL);
 Vertex *v4 = (e->t2 != NULL)?(e->t2->oppositeVertex(e)):(NULL);
 Edge *be1 = (e->t1 != NULL)?(e->t1->nextEdge(e)):(NULL);
 Edge *be4 = (e->t2 != NULL)?(e->t2->prevEdge(e)):(NULL);
 Vertex *v = new Vertex(p->x, p->y, p->z);
 Edge *ne = new Edge(v, e->v2);
 Edge *ne1 = (e->t1 != NULL)?(new Edge(v, v3)):(NULL);
 Edge *ne2 = (e->t2 != NULL)?(new Edge(v, v4)):(NULL);
 Triangle *nt1 = (e->t1 != NULL)?(new Triangle(ne1, ne,be1)):(NULL);
 Triangle *nt2 = (e->t2 != NULL)?(new Triangle(ne, ne2,be4)):(NULL);

 ne->t1 = nt1; ne->t2 = nt2;
 if (ne1 != NULL) {ne1->t1 = e->t1; ne1->t2 = nt1;}
 if (ne2 != NULL) {ne2->t1 = nt2; ne2->t2 = e->t2;}
 if (be1 != NULL) be1->replaceTriangle(e->t1, nt1);
 if (be4 != NULL) be4->replaceTriangle(e->t2, nt2);
 e->v2->e0 = (be1 != NULL)?(be1):(be4);
 e->v2 = v;
 v->e0 = e;
 if (e->t1 != NULL) e->t1->replaceEdge(be1, ne1);
 if (e->t2 != NULL) e->t2->replaceEdge(be4, ne2);

 if (copy_mask)
 {
  ne->mask = e->mask;
  if (nt1 != NULL) nt1->mask = e->t1->mask;
  if (nt2 != NULL) nt2->mask = e->t2->mask;
 }

 V.appendHead(v);
 E.appendHead(ne);
 if (ne1 != NULL) E.appendHead(ne1);
 if (ne2 != NULL) E.appendHead(ne2);
 if (nt1 != NULL) T.appendHead(nt1);
 if (nt2 != NULL) T.appendHead(nt2);

 return v;
}


//// Splits the trianlge 't' by the point 'p'. If 'corr' is TRUE ////
//// the method does not perform the split if this causes the   ////
//// creation of a triangle with a vertex angle < MIN_ANG_TRI.  ////

Vertex *Triangulation::splitTriangle(Triangle *t, Point *p, bool only_proper)
{
 Vertex *v1 = t->v1();
 Vertex *v2 = t->v2();
 Vertex *v3 = t->v3();

 if ((*p) == (*v1)) return (only_proper)?(NULL):v1;
 if ((*p) == (*v2)) return (only_proper)?(NULL):v2;
 if ((*p) == (*v3)) return (only_proper)?(NULL):v3;
 if (!only_proper)
 {
  if (p->getAngle(v1,v2) == M_PI) return splitEdge(t->e2, p);
  if (p->getAngle(v2,v3) == M_PI) return splitEdge(t->e3, p);
  if (p->getAngle(v3,v1) == M_PI) return splitEdge(t->e1, p);
 }

 Vertex *v = new Vertex(p->x, p->y, p->z);
 Edge *ne1 = new Edge(v, v1);
 Edge *ne2 = new Edge(v, v2);
 Edge *ne3 = new Edge(v, v3);
 Triangle *nt1 = new Triangle(ne2, t->e3,ne3);
 Triangle *nt2 = new Triangle(ne3, t->e1,ne1);
 t->e3->replaceTriangle(t, nt1);
 t->e1->replaceTriangle(t, nt2);
 t->replaceEdge(t->e3, ne2);
 t->replaceEdge(t->e1, ne1);
 ne1->t1 = t; ne1->t2 = nt2;
 ne2->t1 = nt1; ne2->t2 = t;
 ne3->t1 = nt2; ne3->t2 = nt1;
 v->e0 = ne1;

 V.appendHead(v);
 E.appendHead(ne1);
 E.appendHead(ne2);
 E.appendHead(ne3);
 T.appendHead(nt1);
 T.appendHead(nt2);

 return v;
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    D E B U G   AND WORK-IN-PROGRESS                          //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Prints general information ///////////

void Triangulation::printReport()
{
 eulerUpdate();

 JMesh::info("*** Triangulation Report ***\n");
 JMesh::info("V: %d\n",V.numels());
 JMesh::info("E: %d\n",E.numels());
 JMesh::info("T: %d\n",T.numels());

 JMesh::info("Boundary: %d components.\n",boundaries());
 JMesh::info("Handles: %d.\n",handles());
 JMesh::info("Shells: %d.\n",shells());
}
