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

#ifndef _JLIST_H
#define _JLIST_H

#include <stdio.h>

/////////////////////////////////////////////////////////////////////////////////////////////

//! Generic node of a doubly likned list.


class Node
{
 friend class List; // This is to make methods in 'List' able to modify n_prev and n_next

 public :
 void *data;			//!< Actual data stored in the node

 //! Creates an isolated node storing 'd'
 Node(const void *d) {data=(void *)d; n_prev=n_next=NULL;}

 //! Creates a new node storing 'd' and links it to a previous node 'p' and to a next one 'n'.
 Node(const Node *p, const void *d, const Node *n);
 ~Node();			//!< Standard destructor

 Node *prev() {return n_prev;}	//!< Returns the previous node in the list, possibly NULL
 Node *next() {return n_next;}	//!< Returns the next node in the list, possibly NULL

 protected:
 Node *n_prev,*n_next;		//!< Previous and next node pointers
};


/////////////////////////////////////////////////////////////////////////////////////////////

//! Doubly linked list.


class List
{ 
 protected :

 Node *l_head;			//!< First node pointer
 Node *l_tail;			//!< Last node pointer
 int l_numels;			//!< Number of elements in the list
 
 public :

 //! Creates an empty list
 List() {l_head = l_tail = NULL; l_numels = 0;}

 //! Creates a list containing an element 'd' (singleton)
 List(const void *d) {l_head = l_tail = new Node(d); l_numels = 1;}

 //! Creates a list out of an array 'd' made of 'n' elements.
 List(const void **d, int n);

 //! Creates a duplicated list.
 List(List& l) {l_head = l_tail = NULL; l_numels = 0; appendList(&l);}

 //! Creates a duplicated list.
 List(List* l) {l_head = l_tail = NULL; l_numels = 0; appendList(l);}

 //! Destructor
 ~List();

 Node *head() const {return l_head;}	//!< Gets the first node, NULL if empty. \n O(1).
 Node *tail() const {return l_tail;}	//!< Gets the last node, NULL if empty. \n O(1).
 int numels() const {return l_numels;}	//!< Gets the number of elements. \n O(1).

 void appendHead(const void *d);	//!< Appends a new node storing 'd' to the head. \n O(1).
 void appendTail(const void *d);	//!< Appends a new node storing 'd' to the tail. \n O(1).

 //! Deletes and removes the node containing 'd'. Returns its position, 0 if 'd' is not in the list. \n O(numels()).
 int  removeNode(const void *d);

 //! Deletes and i'th node (starting from 0). Returns 0 if the list has less than i+1 nodes. \n O(numels()).
 int  removeNode(int i);

 //! Returns the node at position 'i' (starting from 0). Returns NULL if the list has less than i+1 nodes. \n O(numels()).
 Node *getNode(int i) const;

 //! Deletes and removes the node 'n' from the list. \n O(1).
 void removeCell(Node *n);

 //! Appends a list 'l' to the head by duplicating nodes in 'l'. \n O(l->numels()).
 void appendList(const List *l);

 //! Appends a list 'l' to the tail by linking the first node of 'l' to the last one of this list. 'l' becomes empty. \n O(1).
 void joinTailList(List *l);

 void *popHead();		//!< Deletes and removes the first node. Returns its data. \n O(1).
 void *popTail();		//!< Deletes and removes the last node. Returns its data. \n O(1).

 //! Deletes and removes the node 'n' from the list and frees data memory. \n O(1).

 //! Warning. This method uses the free() function to to dispose the memory space
 //! used by the data stored in the node. This means that such data should have
 //! been allocated through malloc(), calloc() or realloc(), and not through the
 //! 'new' operator. On some systems, however, the 'delete' operator simply calls 'free()'
 //! right after the execution of the proper object destructor so, if the object
 //! does not need to free internally allocated memory, it is safe to dispose the
 //! memory trhough free() although the object was allocated by 'new'. This works
 //! on Linux Fedora Core 2 distributions.
 void freeCell(Node *n);

 //! Deletes and removes the node storing 'd' and frees the memory occupied by 'd' itself. \n O(numels()).

 //! Warning. Read the comment for the method 'freeCell()'
 void freeNode(void *d);

 //! Returns the node storing 'd'. NULL if not found. \n O(numels()).
 Node *containsNode(const void *d) const;

 //! Replaces old_n with new_n. The Node containing new_n is returned. \n O(numels()).
 Node *replaceNode(const void *old_n, const void *new_n);

 //! Deletes and removes all the nodes and frees data memory. \n O(numels()).

 //! Warning. Read the comment for the method 'freeCell()'
 void freeNodes();

 void removeNodes();		//!< Deletes and removes all the nodes. \n O(numels()).

 void **toArray() const;		//!< Creates an array out of the list. \n O(numels()).

 //! Sorts the list using 'comp' as comparison function for two elements. \n O(numels()^2).

 //! This method uses the QuickSort algorithm for sorting, thus the complexity is N^2 in the
 //! worst case, but it is actually much faster in the average case. If, however, there is
 //! the need to have a guaranteed O(NlogN) complexity, it is possible to implement a heap
 //! based on the 'abstractHeap' class. See the documentation of the standard 'qsort' library
 //! function for details on the prototype of the comparison function 'comp'.
 int sort(int (*comp)(const void *, const void *));
};

//! Convenience macro to scan the nodes of a list.
#define FOREACHNODE(l, n) for ((n) = (l).head(); (n) != NULL; (n)=(n)->next())

//! Convenience macro to circulate around the nodes of a list 'l' starting from node 'm'. Must exit with break or return.
#define FOREACHNODECIRCULAR(l, m, n) for ((n) = (m); ; (n)=((n)!=(l).tail())?((n)->next()):((l).head()))

#endif // _JLIST_H

