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

#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "jqsort.h"


// Create a new node containing 'd', and link it to       //
// 'p' on the left (prev) and to 'n' on the right (next). //

Node::Node(const Node *p, const void *d, const Node *n)
{
 data=(void *)d;
 if ((n_prev=(Node *)p) != NULL) n_prev->n_next = this;
 if ((n_next=(Node *)n) != NULL) n_next->n_prev = this;
}


// Destroy and unlink the node

Node::~Node()
{
 if (n_prev != NULL) n_prev->n_next = n_next;
 if (n_next != NULL) n_next->n_prev = n_prev;
}


/////////// Constructor from list ///////////////////

List::List(const void **d, int n)
{
 l_head = l_tail = NULL; l_numels = 0;
 for (int i=0; i<n; i++) appendTail(d[i]);
}

///////////////////////// Destructor //////////////////////////

List::~List()
{
 while (l_head != NULL) removeCell(l_head);
}

////////////////// Append an element //////////////////

void List::appendHead(const void *d)
{
 Node *nn = new Node(NULL, d, l_head);
 l_head = nn;
 if (l_tail == NULL) l_tail = l_head;
 l_numels++;
}

void List::appendTail(const void *d)
{
 Node *nn = new Node(l_tail, d, NULL);
 l_tail = nn;
 if (l_head == NULL) l_head = l_tail;
 l_numels++;
}

////////////////// Appends a list //////////////////

void List::appendList(const List *l)
{
 Node *n = l->l_tail;
 
 while (n != NULL)
 {
  appendHead(n->data);
  n=n->prev();
 }
}

////////////////// Joins a list to the l_tail //////////////////

void List::joinTailList(List *l)
{
 if (l->l_numels == 0) return;
 if (l_tail != NULL)
 {
  l_tail->n_next = l->l_head; l->l_head->n_prev = l_tail; l_tail = l->l_tail;
  l_numels += l->l_numels;
 }
 else
 {
  l_head = l->l_head; l_tail = l->l_tail; l_numels = l->l_numels;
 }
 l->l_head = l->l_tail = NULL; l->l_numels = 0;
}

//// Removes the first node and returns the corresponding data /////

void *List::popHead()
{
 void *data = (l_head != NULL)?(l_head->data):(NULL);
 if (l_head != NULL) removeCell(l_head);
 return data;
}

//// Removes the last node and returns the corresponding data /////

void *List::popTail()
{
 void *data = (l_tail != NULL)?(l_tail->data):(NULL);
 if (l_tail != NULL) removeCell(l_tail);
 return data;
}

//////////////////// Removes an element //////////////////

int List::removeNode(const void *d)
{
 Node *tmp = l_head;
 int i=1;

 while (tmp != NULL)
  if (tmp->data == d)
  {
   removeCell(tmp);
   return i;
  }
  else {tmp=tmp->n_next; i++;}

 return 0;
}


//////////////////// Removes an element //////////////////

int List::removeNode(int i)
{
 Node *tmp = l_head;

 while (tmp!=NULL && i--) tmp=tmp->n_next;
 if (tmp==NULL) return 0;

 removeCell(tmp);
 return 1;
}


//////////////////// Gets a node //////////////////

Node *List::getNode(int i) const
{
 Node *tmp = l_head;

 while (tmp!=NULL && i--) tmp=tmp->n_next;
 return tmp;
}


//////////////////// Removes a node //////////////////

void List::removeCell(Node *n)
{
 if (n==l_head) l_head = n->n_next;
 if (n==l_tail) l_tail = n->n_prev;
 delete(n);
 l_numels--;
}


////////////////// Garbage collection //////////////

void List::freeCell(Node *n)
{
 free(n->data);
 removeCell(n);
}

void List::freeNode(void *d)
{
 free(d);
 removeNode(d);
}

//////////////////// Belonging check /////////////////

Node *List::containsNode(const void *d) const
{
 Node *tmp = l_head;
 
 while (tmp != NULL)
  if (tmp->data == d) return tmp;
  else tmp=tmp->n_next;
 
 return NULL;
}

//////////////////// Replaces a node /////////////////

Node *List::replaceNode(const void *od, const void *nd)
{
 Node *tmp = containsNode(od);
 if (tmp != NULL) {tmp->data = (void *)nd; return tmp;}
 appendTail(nd);
 return l_tail;
}

//////////////////////// Garbage collector /////////////////////

void List::freeNodes()
{
 while (l_head != NULL) freeCell(l_head);
}

//////////////////////// Garbage collector /////////////////////

void List::removeNodes()
{
 while (l_head != NULL) removeCell(l_head);
}


///// Conversion to array ///////

void **List::toArray() const
{
 Node *n = l_head;
 int i;
 void **array;

 if (l_numels == 0) return NULL;
 array = (void **)malloc(sizeof(void *)*l_numels);
 if (array == NULL) return NULL;
 for (i=0; i<l_numels; i++, n=n->n_next) array[i] = n->data;

 return array;
}

///// Sorts the list /////////

int List::sort(int (*comp)(const void *, const void *))
{
 void **array;
 int ne = l_numels-1;

 if (l_numels < 2) return 0;
 if ((array = toArray()) == NULL) return 1;

 jqsort(array, l_numels, comp);
 removeNodes();
 for (; ne >= 0; ne--) appendHead(array[ne]);
 free(array);

 return 0;
}

