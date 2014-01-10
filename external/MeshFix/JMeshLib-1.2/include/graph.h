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

#ifndef _GRAPH_H
#define _GRAPH_H

#include "list.h"


//! Base class type for nodes of non-oriented graphs

class graphNode
{
 public:

 //! List of incident edges.
 List edges;

 //! Generic 8-bit mask for marking purposes.
 unsigned char mask;

 graphNode() {mask=0;}
 virtual ~graphNode() {}

 //! Returns TRUE if the node is isolated. O(1).
 bool isIsolated() {return (edges.numels()==0);}

 //! Returns the edge connecting this with 'n'. NULL if not connected. O(degree).
 class graphEdge *getEdge(graphNode *n);
};


//! Base class type for edges of non-oriented graphs

class graphEdge
{
 public:

 //! Edge's end-points
 graphNode *n1, *n2;

 //! Generic 8-bit mask for marking purposes.
 unsigned char mask;

 //! Constructor.
 graphEdge(graphNode *, graphNode *);
 virtual ~graphEdge() {}

 //! Returns the node oppsite to 'n'. O(1).
 graphNode *oppositeNode(graphNode *n) {return (n1==n)?(n2):((n2==n)?(n1):(NULL));}

 //! Returns TRUE if this edge does not connect points. O(1).
 bool isUnlinked() {return (n1==NULL);}

 //! Returns TRUE if 'n' is a node of this edge. O(1).
 bool hasNode(graphNode *n) {return (n1==n || n2==n);}

 //! Makes this edge as 'unlinked' from the graph. O(1).
 void makeUnlinked() {n1=NULL; n2=NULL;}

 //! Edge collapse. O(degree of neighbors).

 //! After one (or a series of) collapse, remember to call Graph::deleteUnlinkedElements()
 //! to make the graph coherent with its node and edge lists.
 void collapse();
};


//! Base class type for non oriented graphs

class Graph
{
 public:

 //! Nodes and edges of the graph.
 List nodes, edges;

 ~Graph();

 //! Adds an existing isolated node to the graph. O(1).
 graphNode *addNode(graphNode *n) {nodes.appendHead(n); return n;}

 //! Creates a new edge out of a pair of nodes. O(degree of nodes).

 //! If the edges already exists, no new edge is created and the
 //! existing one is returned. Otherwise the newly created edge is returned.
 graphEdge *createEdge(graphNode *n1, graphNode *n2);

 //! Eliminates isolated nodes and unlinked edges from the lists. O(N).
 //! The eliminated elements are deleted too.
 void deleteUnlinkedElements();
};


#endif // _GRAPH_H
