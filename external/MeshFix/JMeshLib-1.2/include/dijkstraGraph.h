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

#ifndef _DIJKSTRA_GRAPH_H
#define _DIJKSTRA_GRAPH_H

#include "heap.h"
#include "graph.h"
#include "float.h"

class dijkstraNode : public graphNode
{
 public:

 int index;
 double dist;

 dijkstraNode() : graphNode() {dist = DBL_MAX;}
};


class dijkstraEdge : public graphEdge
{
 public:
 
 double cost;
 
 dijkstraEdge(graphNode *n1, graphNode *n2, double c) : graphEdge(n1,n2) {cost = c;}
};


class dijkstraHeap : abstractHeap
{
 public:
 dijkstraNode **nodes;

 dijkstraHeap(int n, dijkstraNode **nod) : abstractHeap(n)
  {positions = new int[n+1]; nodes = nod; for (int i=0; i<=n; i++) positions[i]=0;}

 ~dijkstraHeap() {delete(positions);}

 void push(dijkstraNode *n) {insert((void *)(n->index));}
 dijkstraNode *getHead() {return (dijkstraNode *)getHead();}
 dijkstraNode *popHead() {return (numels)?(nodes[(j_voidint)removeHead()]):(NULL);}
 int isEmpty() {return (numels==0);}

 void remove(dijkstraNode *n)
  {n->dist = -1; if (positions[n->index]) {upheap(positions[n->index]); removeHead();}}

 void update(dijkstraNode *n)
  {if (!positions[n->index]) push(n); else downheap(upheap(positions[n->index]));}

 int compare(const void *, const void *);
};


class dijkstraGraph : public Graph
{
 int curNodeIndex, maxNumNodes;
 dijkstraNode **nds;
 dijkstraHeap *ch;
 
 public:

 dijkstraGraph(int n)
 {
  nds = new dijkstraNode *[n];
  maxNumNodes = n;
  curNodeIndex = 0;
  ch = new dijkstraHeap(n, nds);
 }
 ~dijkstraGraph() {delete nds; delete ch;}

 dijkstraNode *addNode(dijkstraNode *n);
 dijkstraEdge *createEdge(dijkstraNode *n1, dijkstraNode *n2, double c);
 
 dijkstraNode *popFirstNode() {return ch->popHead();}
 dijkstraNode **getNodes() const {return nds;}
 int getNumNodes() const {return curNodeIndex;}

 // If 'use_distances' is set to true, current distances associated to the nodes
 // are not set to infinity before running the algorithm.
 void runDijkstra(dijkstraNode *s, bool use_distances = false);

 double computeDistance(dijkstraNode *source, dijkstraNode *destination);
};

#endif // _DIJKSTRA_GRAPH_H
