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

#include "jmesh.h"
#include "jqsort.h"
#include <stdlib.h>
#include <string.h>


////////////// Performs some checks and attempts to fix possible errors or degeneracies //////////////

int Triangulation::checkAndRepair()
{
 if (mergeCoincidentEdges()) JMesh::warning("Some coincident edges have been merged.\n");
 if (removeDegenerateTriangles()) JMesh::warning("Some edges have been swapped or collapsed to eliminate degenerate triangles.\n");

 printReport();

 return 1;
}

////////////// Checks the triangulation's connectivity //////////////
//
// This method should be used when implementing new algorithms to
// check the consistency of the connectivity graph. Because an
// inconsistent graph is not assumed by all the other methods, such
// a flaw is considered critical and the program should terminate.
// If connectivity is ok, NULL is returned, otherwise a string
// describing the error is returned.
//
/////////////////////////////////////////////////////////////////////

const char *Triangulation::checkConnectivity()
{
 Vertex *v;
 Edge *e,*e2;
 Triangle *t;
 Node *n,*m;
 List *ve;

 FOREACHVERTEX(v, n)
 {
  if (v == NULL) return "checkConnectivity: detected NULL element in V list!";
  if (v->e0 == NULL) return "checkConnectivity: detected NULL e0 pointer for a vertex!";
  if (!v->e0->hasVertex(v)) return "checkConnectivity: detected wrong e0 pointer for a vertex!";
 }

 FOREACHEDGE(e, n)
 {
  if (e == NULL) return "checkConnectivity: detected NULL element in E list!";
  if (e->v1 == NULL || e->v2 == NULL) return "checkConnectivity: detected edge with one or two NULL end-points!";
  if (e->v1 == e->v2) return "checkConnectivity: detected edge with two coincident end-points!";
  if (e->t1 == NULL && e->t2 == NULL) return "checkConnectivity: detected edge with no incident triangles!";
  if (e->t1 != NULL)
  {
   if (!e->t1->hasEdge(e)) return "checkConnectivity: detected wrong t1 triangle at an edge";
   if (e->commonVertex(e->t1->nextEdge(e)) == e->v1)
	return "checkConnectivity: Edge orientation does not match t1 normal";
  }
  if (e->t2 != NULL)
  {
   if (!e->t2->hasEdge(e)) return "checkConnectivity: detected wrong t2 triangle at an edge";
   if (e->commonVertex(e->t2->nextEdge(e)) == e->v2)
	return "checkConnectivity: Edge orientation does not match t2 normal";
  }
 }

 FOREACHTRIANGLE(t, n)
 {
  if (t == NULL) return "checkConnectivity: detected NULL element in T list!";
  if (t->e1 == NULL || t->e2 == NULL || t->e3 == NULL) return "checkConnectivity: detected NULL as a triangle edge!";
  if (t->e1 == t->e2 || t->e1 == t->e3 || t->e2 == t->e3) return "checkConnectivity: detected triangle with two coincident edges!";
  if (t->v1() == NULL || t->v2() == NULL || t->v3() == NULL) return "checkConnectivity: triangle edges do not share vertices!";
  if (t->e1->t1 != t && t->e1->t2 != t) return "checkConnectivity: detected triangle with 1st edge not pointing to the triangle itself!";
  if (t->e2->t1 != t && t->e2->t2 != t) return "checkConnectivity: detected triangle with 2nd edge not pointing to the triangle itself!";
  if (t->e3->t1 != t && t->e3->t2 != t) return "checkConnectivity: detected triangle with 3rd edge not pointing to the triangle itself!";
 }

 FOREACHEDGE(e, n)
 {
  ve = e->v1->VE();
  FOREACHVEEDGE(ve, e2, m)
  {
   if (e2 != e && e2->oppositeVertex(e->v1) == e->v2) return "checkConnectivity: detected duplicate edge!";
  }
  if (ve->containsNode(e) == NULL) return "checkConnectivity: detected non manifold vertex!";
  delete(ve);
  ve = e->v2->VE();
  FOREACHVEEDGE(ve, e2, m)
  {
   if (e2 != e && e2->oppositeVertex(e->v2) == e->v1) return "checkConnectivity: detected duplicate edge!";
  }
  if (ve->containsNode(e) == NULL) return "checkConnectivity: detected non manifold vertex!";
  delete(ve);
 }

 return NULL;
}


////////////// Duplicate non-manifold vertices /////////////////////
//
// If a vertex is topologically non-manifold, this data structure
// does not guarantee its functionality. Therefore, in order to use
// the same triangle mesh, this method allows to duplicate such
// vertices. Notice that the data-structure cannot code non-manifold
// edges.
//
////////////////////////////////////////////////////////////////////

int Triangulation::duplicateNonManifoldVertices()
{
 Vertex *v;
 Edge *e, *f;
 Node *n, *m;
 List *ve;
 int dv = 0;

 FOREACHEDGE(e, n)
 {
  ve = e->v1->VE();
  if (ve->containsNode(e) == NULL)
  {
   V.appendHead(v = new Vertex(e->v1));
   FOREACHVEEDGE(ve, f, m) f->replaceVertex(e->v1, v);
   v->e0 = e->v1->e0;
   e->v1->e0 = e;
   dv++;
  }
  delete(ve);
 }
 FOREACHEDGE(e, n)
 {
  ve = e->v2->VE();
  if (ve->containsNode(e) == NULL)
  {
   V.appendHead(v = new Vertex(e->v2));
   FOREACHVEEDGE(ve, f, m) f->replaceVertex(e->v2, v);
   v->e0 = e->v2->e0;
   e->v2->e0 = e;
   dv++;
  }
  delete(ve);
 }

 if (dv) d_boundaries = d_handles = d_shells = 1;

 return dv;
}

////////////// Checks the triangulation geometry //////////////
//// 							 //////
//// Looks for coincident vertices, degenerate triangles //////
//// and overlapping triangles.				 //////
//// If something is wrong returns the closest vertex.   //////
//// 							 //////
///////////////////////////////////////////////////////////////

Vertex *Triangulation::checkGeometry()
{
 int i;
 Vertex *ret = NULL;
 double ang, minda = 0;
 Triangle *t;
 Edge *e;
 Vertex **varr = (Vertex **)V.toArray();
 Edge **evarr;
 Vertex *v1, *v2;
 Node *n;

 if (varr == NULL) JMesh::warning("checkGeometry: Not enough memory. Can't check for coincident vertices.\n");
 else
 {
  jqsort((void **)varr, V.numels(), xyzCompare);
  for (i=0; i<(V.numels()-1); i++)
  {
   v1 = ((Vertex *)varr[i]);
   v2 = ((Vertex *)varr[i+1]);
   if ((*v1)==(*v2))
   {
    ret = v1;
    JMesh::warning("checkGeometry: detected coincident vertices.\n");
    if (v1->getEdge(v2))
    {
     JMesh::warning("               and there is an edge connecting them!\n");
     free(varr);
     return v1;
    }
   }
  }
  free(varr);
 }

 evarr = (Edge **)E.toArray();
 if (evarr == NULL) JMesh::warning("checkGeometry: Not enough memory. Can't check for coincident edges.\n");
 else
 {
  jqsort((void **)evarr, E.numels(), lexEdgeCompare);
  for (i=0; i<(E.numels()-1); i++)
  {
   if (!lexEdgeCompare(evarr[i], evarr[i+1]))
   {
    ret = ((Edge *)evarr[i])->v1;
    JMesh::warning("checkGeometry: detected coincident edges.\n");
   }
  }
  free(evarr);
 }

 FOREACHTRIANGLE(t, n)
 {
  ang = t->getAngle(t->v1());
  if (ang == 0 || ang == M_PI) {JMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v1();}
  ang = t->getAngle(t->v2());
  if (ang == 0 || ang == M_PI) {JMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v2();}
  ang = t->getAngle(t->v3());
  if (ang == 0 || ang == M_PI) {JMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v3();}
 }

 ang = minda = 0;
 FOREACHEDGE(e, n)
  if (e->t1 != NULL && e->t2 != NULL && (ang = e->t1->getDAngle(e->t2)) == M_PI)
   {JMesh::warning("checkGeometry: overlapping triangles detected.\n"); return e->v1;}
  else minda = MAX(minda,ang);

 JMesh::info("checkGeometry: minimum dihedral angle = %f (%f DEGs)\n", M_PI-minda, ((M_PI-minda)*360)/(2*M_PI));

 return ret;
}


///// Merges possible coincident edges //////////

int Triangulation::mergeCoincidentEdges()
{
 List bes;
 Node *n;
 Edge *e1, *e2;

 FOREACHEDGE(e1, n) if (e1->isOnBoundary()) bes.appendHead(e1);
 if (!bes.numels()) return 0;

 Edge **evarr = (Edge **)bes.toArray();
 int i, ret = 0, wrn = 0;

 jqsort((void **)evarr, bes.numels(), lexEdgeCompare);

  for (i=0; i<(bes.numels()-1); i++)
  {
   e1 = evarr[i];
   e2 = evarr[i+1];
   if (!e2->isLinked()) evarr[i+1]=e1;
   else if (!lexEdgeCompare(e1, e2) &&
            (
             (((*(e1->v1))==(*(e2->v1))) && ((e1->t1 && e2->t2) || (e1->t2 && e2->t1))) ||
             (((*(e1->v1))==(*(e2->v2))) && ((e1->t1 && e2->t1) || (e1->t2 && e2->t2)))
            )
           )
           {
            if (!e2->merge(e1)) wrn++;
            else
            {
             ret += e2->v1->zip();
             ret += e2->v2->zip();
            }
           }
  }

 free(evarr);
 removeEdges();
 removeVertices();

 if (wrn) JMesh::warning("mergeCoincidentEdges: Couldn't merge unconsistently oriented edges.\n");

 if (ret) d_boundaries = d_handles = d_shells = 1;

 return ret;
}


//////// Eliminates duplicated triangles (i.e. having the same vertices) /////////

int Triangulation::removeDuplicatedTriangles()
{
 Edge *e;
 Node *n;
 Point p;
 int i=0;

 FOREACHEDGE(e, n)
  if (!e->isOnBoundary() && e->t1->oppositeVertex(e) == e->t2->oppositeVertex(e))
  {
   p = e->t2->getCenter();
   splitTriangle(e->t2, &p, 1);
   i++;
  }

 if (i)  d_boundaries = d_handles = d_shells = 1;

 return i;
}


//////// Flip or collapse edges to eliminate degenerate triangles /////////

int Triangulation::removeDegenerateTriangles()
{
 Triangle *t;
 Edge *e;
 int collapses, swaps, degn, tcs = 0;
 List todo, *vt;
 Node *n;
int i=0;
 do
 {
  collapses = swaps = 0;

  FOREACHTRIANGLE(t, n) if (t->isDegenerate())
   {MARK_VISIT2(t); todo.appendHead(t);}

  while (todo.numels())
  {
i++;
   t = (Triangle *)todo.popHead();
   UNMARK_VISIT2(t);
   if (t->isLinked())
   {
    if ((e = t->isCap()) != NULL)
    {
     if (e->isOnBoundary()) {unlinkTriangle(t); collapses++;}
     else if (e->swap())
     {
      if (e->t1->overlaps() || e->t2->overlaps() || e->t1->isCap() || e->t2->isCap()) e->swap(1);
      else
      {
       swaps++;
       if (!IS_VISITED2(e->t1)) {MARK_VISIT2(e->t1); todo.appendTail(e->t1);}
       if (!IS_VISITED2(e->t2)) {MARK_VISIT2(e->t2); todo.appendTail(e->t2);}
      }
     }
    }
    else if ((e = t->isNeedle()) != NULL)
    {
     vt = e->v2->VT();
     if (e->collapse())
     {
      collapses++;
      FOREACHVTTRIANGLE(vt, t, n)
       if (t->isLinked() && !IS_VISITED2(t)) {MARK_VISIT2(t); todo.appendTail(t);}
     }
     else if (!e->isOnBoundary() && e->oppositeTriangle(t)->oppositeVertex(e)->valence() == 3)
     {
      if (e->oppositeTriangle(t)->nextEdge(e)->collapse())
      {
       MARK_VISIT2(t); todo.appendHead(t);
       t = e->oppositeTriangle(t);
       if (t && !IS_VISITED2(t)) {MARK_VISIT2(t); todo.appendTail(t);}
       collapses++;
      }
     }
     delete(vt);
    }
   }
  }

  if (collapses) removeUnlinkedElements();
  tcs += collapses; tcs += swaps;
 } while (collapses+swaps);

 degn=0;
 FOREACHTRIANGLE(t, n) if (t->isDegenerate()) degn++;

 if (degn)
 {
  FOREACHTRIANGLE(t, n) if (t->isDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
  JMesh::warning("removeDegenerateTriangles: %d degeneracies couldn't be removed.\n",degn);
  JMesh::warning("removeDegenerateTriangles: and have been selected.\n");
  return -tcs;
 }

 return tcs;
}


//// If the mesh is made of more than one connected component ////
//// keep only the biggest one and remove all the others.     ////

int Triangulation::removeSmallestComponents()
{
 Node *n,*m;
 List todo;
 List components;
 List *component, *biggest = NULL;
 Triangle *t, *t1, *t2, *t3;
 int nt = 0, gnt = 0;

 t = ((Triangle *)T.head()->data);
 n = T.head();
 do
 {
  component = new List;
  components.appendHead(component);
  todo.appendHead(t);
  while (todo.numels())
  {
   t = (Triangle *)todo.head()->data;
   todo.removeCell(todo.head());
   if (!IS_VISITED2(t))
   {
    t1 = t->t1();
    t2 = t->t2();
    t3 = t->t3();

    if (t1 != NULL && !IS_VISITED2(t1)) todo.appendHead(t1);
    if (t2 != NULL && !IS_VISITED2(t2)) todo.appendHead(t2);
    if (t3 != NULL && !IS_VISITED2(t3)) todo.appendHead(t3);

    MARK_VISIT2(t);
    component->appendTail(t);
   }
  }
  todo.removeNodes();
  for (; n != NULL; n=n->next()) {t = ((Triangle *)n->data); if (!IS_VISITED2(t)) break;}
 }
 while (n != NULL);

 FOREACHNODE(components, n)
  if ((nt = ((List *)n->data)->numels()) > gnt) {gnt=nt; biggest = (List *)n->data;}

 FOREACHTRIANGLE(t, n) UNMARK_VISIT2(t);

 nt = 0;
 FOREACHNODE(components, n)
  if (((List *)n->data) != biggest)
   FOREACHVTTRIANGLE(((List *)n->data), t, m)
   {
    if (t->e1->v1 != NULL) t->e1->v1->e0 = NULL;
    if (t->e1->v2 != NULL) t->e1->v2->e0 = NULL;
    if (t->e2->v1 != NULL) t->e2->v1->e0 = NULL;
    if (t->e2->v2 != NULL) t->e2->v2->e0 = NULL;
    if (t->e3->v1 != NULL) t->e3->v1->e0 = NULL;
    if (t->e3->v2 != NULL) t->e3->v2->e0 = NULL;
    t->e1->v1 = t->e1->v2 = t->e2->v1 = t->e2->v2 = t->e3->v1 = t->e3->v2 = NULL;
    t->e1 = t->e2 = t->e3 = NULL;
    nt++;
   }

 FOREACHNODE(components, n) delete((List *)n->data);

 if (nt)
 {
  d_boundaries = d_handles = d_shells = 1;
  removeUnlinkedElements();
  return 1;
 }

 return 0;
}


//// Traverses the triangulation and inverts normals in order ////
//// to make the adjacences consistent.			      ////

int Triangulation::forceNormalConsistence()
{
 int ret = 0;
 Node *n;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (!IS_VISITED2(t))
  ret |= forceNormalConsistence(t);
 FOREACHTRIANGLE(t, n) UNMARK_VISIT2(t);
 return ret;
}

int Triangulation::forceNormalConsistence(Triangle *t0)
{
 Node *n;
 Edge *e;
 List todo, elist;
 Triangle *t, *t1, *t2, *t3;
 int tmp1, tmp2, r=0, wrn = 0, isclosed = 1;

 todo.appendHead(t0);

 while (todo.numels())
 {
  t = (Triangle *)todo.head()->data;
  todo.removeCell(todo.head());
  if (!IS_VISITED2(t))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
   if (!IS_VISITED2(t->e1)) {MARK_VISIT2(t->e1); elist.appendHead(t->e1);}
   if (!IS_VISITED2(t->e2)) {MARK_VISIT2(t->e2); elist.appendHead(t->e2);}
   if (!IS_VISITED2(t->e3)) {MARK_VISIT2(t->e3); elist.appendHead(t->e3);}

   if (t1 != NULL && !IS_VISITED2(t1)) {todo.appendHead(t1); if (!t->checkAdjNor(t1)) {t1->invert(); r=1;}}
   if (t2 != NULL && !IS_VISITED2(t2)) {todo.appendHead(t2); if (!t->checkAdjNor(t2)) {t2->invert(); r=1;}}
   if (t3 != NULL && !IS_VISITED2(t3)) {todo.appendHead(t3); if (!t->checkAdjNor(t3)) {t3->invert(); r=1;}}

   MARK_VISIT2(t);
  }
 }

 FOREACHVEEDGE((&(elist)), e, n)
 {
  UNMARK_VISIT2(e);
  if (isclosed && e->isOnBoundary()) isclosed = 0;
  tmp1 = (e->t1 != NULL)?((e->commonVertex(e->t1->nextEdge(e)) == e->v1)?(-1):(1)):(0);
  tmp2 = (e->t2 != NULL)?((e->commonVertex(e->t2->nextEdge(e)) == e->v2)?(-1):(1)):(0);

  if (tmp1*tmp2 < 0)
  {
   wrn++;
   if (tmp1 == -1) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
   Edge *ne = new Edge(e->v2, e->v1);
   E.appendHead(ne);
   e->t2->replaceEdge(e, ne);
   ne->t2 = e->t2; e->t2 = NULL;
  } else if (tmp1 == -1 || tmp2 == -1) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
 }

 if (wrn)
 {
  d_boundaries = d_handles = d_shells = 1;
  JMesh::warning("forceNormalConsistence: Triangulation was not orientable. Cut performed.\n");
 }

 if (isclosed)
 {
  t = topTriangle(t0);
  if (t->getNormal().z < 0) {flipNormals(t0); r=1;}
 }

 return r;
}


// If possible, swap edges to remove overlaps. When it is not
// enough, remove the overlapping triangles from the mesh.

int Triangulation::removeOverlappingTriangles()
{
 Node *n;
 Edge *e, *ea;
 List oved;

 FOREACHEDGE(e, n)
  if (!e->isOnBoundary() && e->t1->getDAngle(e->t2) == M_PI)
   oved.appendHead(e);

 oved.sort(edgeCompare);

 for (n=oved.tail(); n!=NULL; n=n->prev())
 {
  e = (Edge *)n->data;
  if (e->t1->getDAngle(e->t2) == M_PI && e->swap())
  {
   if (e->t1->isDegenerate() || e->t2->isDegenerate()) {e->swap(1); continue;}
   ea = e->t1->nextEdge(e);
   if (!ea->isOnBoundary() && ea->t1->getDAngle(ea->t2) == M_PI) {e->swap(1); continue;}
   ea = e->t1->prevEdge(e);
   if (!ea->isOnBoundary() && ea->t1->getDAngle(ea->t2) == M_PI) {e->swap(1); continue;}
   ea = e->t2->nextEdge(e);
   if (!ea->isOnBoundary() && ea->t1->getDAngle(ea->t2) == M_PI) {e->swap(1); continue;}
   ea = e->t2->prevEdge(e);
   if (!ea->isOnBoundary() && ea->t1->getDAngle(ea->t2) == M_PI) {e->swap(1); continue;}
  }
 }

 for (n=oved.tail(); n!=NULL; n=n->prev())
 {
  e = (Edge *)n->data;
  if (!e->isOnBoundary() && e->t1->getDAngle(e->t2) == M_PI) {unlinkTriangle(e->t1); unlinkTriangle(e->t2);}
 }
 removeUnlinkedElements();
 d_boundaries = d_handles = d_shells = 1;

 return 0;
}


//// Selects all the handles within a radius L ////

int Triangulation::selectTinyHandles(double L)
{
 if (shells() > 1) {JMesh::warning("Sorry. This feature was not implemented for multi-component meshes.\n"); return 1;}

 int bdrr, nh = 0;
 Triangle *t, *lt, *rt;
 Edge *e, *se, *le, *re, *gate = NULL;
 Vertex *v, *ov;
 Node *n, *m;

 FOREACHTRIANGLE(t, n)
  if (!t->v1()->isOnBoundary() && !t->v2()->isOnBoundary() && !t->v3()->isOnBoundary()) break;
 if (n==NULL) return 1;

 unmarkEverything();
 t->mask = t->v1()->mask = t->v2()->mask = t->v3()->mask = 1;
 t->e1->mask = t->e2->mask = t->e3->mask = 1;
 List todo((gate = t->e1));

 while ((gate = (Edge *)todo.popHead()) != NULL)
 {
  while(1)
  {
   bdrr = 0;
   t = (gate->t1 == NULL || gate->t1->mask)?((gate->t2 == NULL || gate->t2->mask)?(NULL):(gate->t2)):(gate->t1);
   if (t == NULL)
   {
    if (gate->t1) MARK_VISIT2(gate->t1); else MARK_VISIT2(gate->t2);
    nh++;
    break;
   }
   t->mask = 1; t->e1->mask = t->e2->mask = t->e3->mask = 1;
   le = t->prevEdge(gate);
   re = t->nextEdge(gate);
   lt = t->leftTriangle(gate);
   rt = t->rightTriangle(gate);
   ov = t->oppositeVertex(gate);
   if (lt && rt && !lt->mask && !rt->mask && !ov->mask && ov->isOnBoundary())
   {
    bdrr = 1;
    e = se = ov->nextBoundaryEdge();
    v = ov;
    do
    {
     e->mask = e->v1->mask = e->v2->mask = 1;
     v = e->oppositeVertex(v);
     e = v->nextBoundaryEdge();
    } while (e != se);
   }

   if (ov->mask)
   {
    if ((!lt || lt->mask) && (!rt || rt->mask)) break;	// E
    else if ((!lt || lt->mask)) gate = re;		// L
    else if ((!rt || rt->mask)) gate = le;		// R
    else {if (!bdrr) todo.appendHead(le); gate = re;}	// S
   }
   else gate = re;					// C

   t->v1()->mask = t->v2()->mask = t->v3()->mask = 1;
  }
 }

 nh /= 2;

 FOREACHTRIANGLE(t, n) UNMARK_VISIT(t);
 FOREACHEDGE(e, n) UNMARK_VISIT(e);
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 Point center;
 Triangulation *stin;
 int sh=0;

 FOREACHTRIANGLE(t, n) if (IS_VISITED2(t))
 {
  center = t->getCenter();
  selectSphericalRegion(t, L, &center);
 }

 FOREACHTRIANGLE(t, n) if (IS_VISITED2(t))
 {
  UNMARK_VISIT2(t);
  if ((stin = createSubMeshFromSelection(t)) != NULL)
  {
   if (stin->handles() < 1) invertSelection(t);
   else
   {
    sh += stin->handles();
    FOREACHVTTRIANGLE((&(stin->T)), lt, m) {rt = (Triangle *)lt->info; UNMARK_VISIT2(rt);}
   }
   delete(stin);
  }
 }

 JMesh::info("Selected %d out of %d handles\n",sh,nh);
 JMesh::info("Radius: %f\n",L);

 return 0;
}
