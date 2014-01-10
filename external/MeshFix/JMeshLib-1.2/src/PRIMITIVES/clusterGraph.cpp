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

#include "clusterGraph.h"


int clusterHeap::compare(const void *e1, const void *e2)
{
 clusterEdge *a = edges[(j_voidint)e1];
 clusterEdge *b = edges[(j_voidint)e2];
 double l1 = a->cost;
 double l2 = b->cost;

 if (l1 < l2) return -1;
 if (l1 > l2) return 1;
 return 0;
}


clusterEdge *clusterGraph::createEdge(graphNode *n1, graphNode *n2)
{
 Node *n;
 FOREACHNODE(n1->edges, n)
  if (((clusterEdge *)n->data)->hasNode(n2))
   return (clusterEdge *)n->data;

 if (curEdgeIndex >= maxNumEdges) return NULL;

 clusterEdge *ne = new clusterEdge(n1, n2, curEdgeIndex);
 edges.appendHead(ne);
 ces[curEdgeIndex] = ne;
 ne->cost = costFunction(ne->n1, ne->n2);
 ch->update(ne);
 curEdgeIndex++;
 
 return ne;
}


clusterEdge *clusterGraph::getFirstEdge()
{
 while (!ch->isEmpty() && ch->getFirst()->isUnlinked()) ch->popHead();
 return ch->getFirst();
}

int clusterGraph::collapseFirstEdge(void (*mergenodes)(const void *, const void *))
{
 clusterEdge *e;

 while ((e=ch->popHead())!=NULL) if (!e->isUnlinked()) break;
 if (e==NULL) return 0;

 if (mergenodes != NULL) mergenodes(e->n1, e->n2);
  
 graphNode *gn = e->n1;
 e->collapse();
 Node *n;
 clusterEdge *ne;
 FOREACHNODE(gn->edges, n)
 {
  ne = (clusterEdge *)n->data;
  ne->cost = costFunction(ne->n1, ne->n2);
  ch->update(ne);
 }
 return 1;
}
