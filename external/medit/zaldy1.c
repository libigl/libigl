/*
 *  memory allocation for data structures
 *
 *  Written by Pascal J. Frey, Inria-Rocquencourt
 *  Copyright (c) Inria, 1999.  All rights reserved. 
*/

#ifndef IGL
#ifdef __cplusplus
extern "C" {
#endif
#endif

#include <assert.h>
#include <stdlib.h>
#include "medit.h"
#include "extern.h"
#include "sproto.h"


int zaldy1(pMesh mesh) {

  /* default */
  if ( ddebug ) printf("allocate %d points\n",mesh->np);
  mesh->point = (pPoint)M_calloc(mesh->np+1,sizeof(Point),"zaldy1.point");
  assert(mesh->point);

  if ( ddebug ) printf("allocate %d tria\n",mesh->nt);
  mesh->tria = (pTriangle)M_calloc(mesh->nt+1,sizeof(Triangle),"zaldy1.tria");
  assert(mesh->tria);

  if ( ddebug ) printf("allocate %d quad\n",mesh->nq);
  if ( mesh->nq ) {
    mesh->quad = (pQuad)M_calloc(mesh->nq+1,sizeof(Quad),"zaldy1.quad");
	assert(mesh->quad);
  }

  if ( mesh->ntet ) {
    if ( ddebug ) printf("allocate %d tetra\n",mesh->ntet);
    mesh->tetra = (pTetra)M_calloc(mesh->ntet+1,sizeof(Tetra),"zaldy1.tetra");
	assert(mesh->tetra);
  }

  if ( mesh->nhex ) {
    if ( ddebug ) printf("allocate %d hexa\n",mesh->nhex);
    mesh->hexa = (pHexa)M_calloc(mesh->nhex+1,sizeof(Hexa),"zaldy1.hexa");
	assert(mesh->hexa);
  }

  if ( mesh->na ) {
    if ( ddebug ) printf("allocate %d edges\n",mesh->na);
    mesh->edge = (pEdge)M_calloc(mesh->na+1,sizeof(Edge),"zaldy1.edge");
	assert(mesh->edge);
  }

  if ( mesh->nvn || mesh->ntg ) {
    mesh->extra = (pExtra)M_calloc(1,sizeof(Extra),"zaldy1.extra");
	assert(mesh->extra);
  
	if ( mesh->nvn ) {
  	  mesh->extra->n = (float*)M_calloc(3*mesh->nvn+1,sizeof(float),"inmesh");
	  assert(mesh->extra->n);
    }

    if ( mesh->ntg ) {
      mesh->extra->t = (float*)M_calloc(3*mesh->ntg+1,sizeof(float),"inmesh");
	  assert(mesh->extra->n);
    }
  }

  return(1);
}
