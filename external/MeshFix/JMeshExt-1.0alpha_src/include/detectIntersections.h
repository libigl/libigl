#ifndef DETECT_INTERSECTIONS_H
/****************************************************************************
* JMeshExt                                                                  *
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

#define DETECT_INTERSECTIONS_H

#include "exttrimesh.h"

#define DI_STORED_PANORMAL(t) (((t)->info))
#define DI_STORED_PNORMAL(t) (((Point *)DI_STORED_PANORMAL(t)))
#define DI_STORED_NORMAL(t) (*DI_STORED_PNORMAL(t))

#define DI_MAX_NUMBER_OF_CELLS	10000

#define DI_TEINT_EPS 1.0e-15

#define DI_EPSILON_POINT Point(1.0e-9, 1.0e-9, 1.0e-9)


class di_cell
{
 public:
 Point mp, Mp;
 List triangles;

 di_cell() {}
 di_cell(Triangulation *tin, bool useAll=1);

 inline bool is_Point_in_cell(Point *p)
  {return (p->x >= mp.x && p->x <= Mp.x && p->y >= mp.y && p->y <= Mp.y && p->z >= mp.z && p->z <= Mp.z);}

 bool is_triangleBB_in_cell(Triangle *t);

 di_cell *fork();

 bool doesNotIntersectForSure();
 void di_selectIntersections();

 bool collinearPoints(Point *, Point *, Point *);
 Point *edgeIntersectsTriangle(Edge *, Triangle *, Edge **);
 Point *edgeEdgeIntersection(Edge *, Edge *);
};

#endif // DETECT_INTERSECTIONS_H
