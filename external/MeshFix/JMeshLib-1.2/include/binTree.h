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

#ifndef _BINTREE_H
#define _BINTREE_H

//! Generic binary tree.

//! This class defines a node of a generic
//! binary tree.

class binTree
{
 protected:
 void *data;				//!< Data within the node
 binTree *parent, *left, *right;	//!< Connectivity

 public:

 binTree(void *v);			//!< Constructor (singleton)

 //! Constructor (join two sub-trees)

 //! Creates a new node for 'v' and makes it the parent
 //! of the existing binary trees 'l' and 'r'.

 binTree(void *v, binTree *l, binTree *r);
 ~binTree();				//!< Destructor

 void setValue(void *v) {data=v;}			//!< Sets the value of this node
 void *getValue() const {return data;} 			//!< Returns the value of this node
 binTree *getLeftChild() const {return left;}		//!< Returns the left child of this node
 binTree *getRightChild() const {return right;}		//!< Returns the right child of this node
 binTree *getParentNode() const {return parent;}	//!< Returns the parent of this node
};

#endif // _BINTREE_H
