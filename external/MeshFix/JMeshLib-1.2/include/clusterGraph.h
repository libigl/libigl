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

#ifndef _CLUSTER_GRAPH_H
#define _CLUSTER_GRAPH_H

#include "heap.h"
#include "graph.h"
#include "float.h"

class clusterEdge : public graphEdge
{
 public:
 
 int index;
 double cost;
 
 clusterEdge(graphNode *n1, graphNode *n2, int i) : graphEdge(n1,n2) {index = i; cost = 0.0;}
};


class clusterHeap : abstractHeap
{
 public:
 clusterEdge **edges;

 clusterHeap(int n, clusterEdge **edg) : abstractHeap(n)
  {positions = new int[n+1]; edges = edg; for (int i=0; i<=n; i++) positions[i]=0;}

 ~clusterHeap() {delete(positions);}

 void push(clusterEdge *e) {insert((void *)(e->index));}
 clusterEdge *getFirst() {return (numels)?(edges[(j_voidint)getHead()]):(NULL);}
 clusterEdge *popHead() {return (numels)?(edges[(j_voidint)removeHead()]):(NULL);}
 int isEmpty() {return (numels==0);}

 void remove(clusterEdge *e)
  {e->cost = -1; if (positions[e->index]) {upheap(positions[e->index]); removeHead();}}

 void update(clusterEdge *e)
  {if (!positions[e->index]) push(e); else downheap(upheap(positions[e->index]));}

 int compare(const void *, const void *);
};


class clusterGraph : public Graph
{
 int curEdgeIndex, maxNumEdges;
 clusterEdge **ces;
 clusterHeap *ch;
 double (*costFunction)(const void *, const void *);
 
 public:

 clusterGraph(int n, double (*cf)(const void *, const void *))
 {
  ces = new clusterEdge *[n];
  maxNumEdges = n;
  curEdgeIndex = 0;
  ch = new clusterHeap(n, ces);
  costFunction = cf;
 }
 ~clusterGraph() {delete ces; delete ch;}

 clusterEdge *createEdge(graphNode *n1, graphNode *n2);
 
 clusterEdge *getFirstEdge();
 double getLowestCost() {clusterEdge *e = getFirstEdge(); return (e!=NULL)?(e->cost):(DBL_MAX);}
 int collapseFirstEdge(void (*)(const void *, const void *) = NULL);
};

#endif // _CLUSTER_GRAPH_H
