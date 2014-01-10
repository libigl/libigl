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

#include "dijkstraGraph.h"


int dijkstraHeap::compare(const void *n1, const void *n2)
{
 dijkstraNode *a = nodes[(j_voidint)n1];
 dijkstraNode *b = nodes[(j_voidint)n2];
 double l1 = a->dist;
 double l2 = b->dist;

 if (l1 < l2) return -1;
 if (l1 > l2) return 1;
 return 0;
}


dijkstraNode *dijkstraGraph::addNode(dijkstraNode *n)
{
 nodes.appendHead(n);
 if (curNodeIndex >= maxNumNodes) return NULL;
 nds[curNodeIndex] = n;
 n->index = curNodeIndex;
 curNodeIndex++;
 return n;
}

dijkstraEdge *dijkstraGraph::createEdge(dijkstraNode *n1, dijkstraNode *n2, double cost)
{
 Node *n;
 FOREACHNODE(n1->edges, n)
  if (((dijkstraEdge *)n->data)->hasNode(n2))
   return (dijkstraEdge *)n->data;

 dijkstraEdge *ne = new dijkstraEdge(n1, n2, cost);
 edges.appendHead(ne);
 
 return ne;
}

void dijkstraGraph::runDijkstra(dijkstraNode *n0, bool use_distances)
{
 Node *n;
 dijkstraEdge *de;
 dijkstraNode *dn, *dd;
 double d;
 FOREACHNODE(nodes, n) {dn=((dijkstraNode *)n->data); dn->mask=0; if (!use_distances) dn->dist = DBL_MAX;}

 n0->dist = 0.0;
 ch->push(n0);

 while ((dn=ch->popHead())!=NULL)
 {
  dn->mask=1;
  FOREACHNODE(dn->edges, n)
  {
   de = ((dijkstraEdge *)n->data);
   dd = (dijkstraNode *)de->oppositeNode(dn);
   if (dd->mask == 0)
   {
    d = dn->dist + de->cost;
    if (d < dd->dist)
    {
     dd->dist = d;
     ch->update(dd);
    }
   }
  }
 }
}


double dijkstraGraph::computeDistance(dijkstraNode *n0, dijkstraNode *n1)
{
 Node *n;
 dijkstraEdge *de;
 dijkstraNode *dn, *dd;
 double d;
 FOREACHNODE(nodes, n) {dn=((dijkstraNode *)n->data); dn->mask=0; dn->dist = DBL_MAX;}

 n0->dist = 0.0;
 ch->push(n0);

 while ((dn=ch->popHead())!=NULL)
 {
  if (dn == n1) {while (ch->popHead()!=NULL); return n1->dist;}
  dn->mask=1;
  FOREACHNODE(dn->edges, n)
  {
   de = ((dijkstraEdge *)n->data);
   dd = (dijkstraNode *)de->oppositeNode(dn);
   if (dd->mask == 0)
   {
    d = dn->dist + de->cost;
    if (d < dd->dist)
    {
     dd->dist = d;
     ch->update(dd);
    }
   }
  }
 }

 return DBL_MAX;
}
