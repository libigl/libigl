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

#ifndef _HEAP_H
#define _HEAP_H

#include "j_mesh.h"

//! Heap base class type.

//! abstractHeap is the base class for implementing heaps.
//! Each implementation (class extension) must define the method 
//! compare to be used for sorting the heap. If the objects being 
//! sorted are non-negative numbers, a special implementation may
//! use the field positions to record the index of each element
//! within the heap. This feature is useful when there is a need to 
//! re-sort an element whose cost changes after its insertion into 
//! the heap. The array 'positions' must be allocated by 
//! the extended class constructor, and must be able to contain NMAX+1 
//! integer numbers, where NMAX is the maximum value that can be 
//! assumed by an object.


class abstractHeap
{
 protected:
 void **heap;		//!< Heap data is stored here
 int numels;		//!< Current number of elements
 int maxels;		//!< Maximum number of elements
 int *positions;	//!< Optional pointer to an array of positions

 int upheap(int i);	//!< Moves the i'th object up on the heap
 int downheap(int i);	//!< Moves the i'th object down on the heap

 //! Comparison of two heap elements

 //! This function must be implemented in the extended class.
 //! The return value must be <0 if a<b, >0 if a>b or 0 if a=b.

 virtual int compare(const void *a, const void *b) = 0;

 public :
 
 abstractHeap(int n);	//!< Creates a heap which can contain up to 'n' elements
 virtual ~abstractHeap() = 0;		//!< Default destructor
 
 //! Inserts 'e' into the heap

 //! Inserts an element 'e' into the heap in the correct position, according to the
 //! method compare. If the insertion fails because the heap is full, -1 is
 //! returned, otherwise the index position of the newly inserted element is
 //! returned.

 int insert(void *e);

 int isEmpty() const {return (numels == 0);}	//!< Returns TRUE if the heap is empty
 void *getHead() const {return heap[1];}	//!< Returns the first element of the heap
 void *removeHead();				//!< Removes and returns the first element after rearranging the heap
 void flush() {numels=0;}			//!< Removes all the elements
};

#endif // _HEAP_H
