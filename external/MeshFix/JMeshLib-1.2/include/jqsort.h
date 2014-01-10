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

//! \file
//! \brief Declaration of a generic QuickSort function.
//!
//! The  jqsort()  function sorts an array with numels elements.
//! The v argument points to the start of the array of elements casted to void *.
//! The contents of the array are sorted in ascending order according to  a
//! comparison  function  pointed  to  by  comp, which is called with two
//! arguments that point to the objects being compared.
//! The comparison function must return an integer less than, equal to,  or
//! greater  than  zero  if  the first argument is considered to be respectively
//! less than, equal to, or greater than the second.  If two members
//! compare as equal, their order in the sorted array is undefined.
//! See the manpage of the standard library qsort() function for further information.

extern void jqsort(void *v[], int numels, int (*comp)(const void *, const void *));
