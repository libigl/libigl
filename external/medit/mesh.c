#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define FLOAT_MAX  1.e20
ubyte   dosurf;


double volTet(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];

  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];

  vol = (c2[0]-c1[0]) * (ay*bz - az*by) \
      + (c2[1]-c1[1]) * (az*bx - ax*bz) \
      + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return(vol / 6.0);
}


/* dump mesh info */
void meshInfo(pMesh mesh) {
  fprintf(stdout,"    Vertices   %8d",mesh->np);
  if ( mesh->nc ) fprintf(stdout,"  Corners  %d",mesh->nc);
  if ( mesh->nr ) fprintf(stdout,"  Required %d",mesh->nr);
  fprintf(stdout,"\n");
  if ( mesh->nt )   fprintf(stdout,"    Triangles  %8d\n",mesh->nt);
  if ( mesh->nq )   fprintf(stdout,"    Quads      %8d\n",mesh->nq);
  if ( mesh->ntet ) fprintf(stdout,"    Tetrahedra %8d\n",mesh->ntet);
  if ( mesh->nhex ) fprintf(stdout,"    Hexahedra  %8d\n",mesh->nhex);
  if ( mesh->na ) {
    fprintf(stdout,"    Edges      %8d",mesh->na);
    if ( mesh->nri ) fprintf(stdout,"  Ridges   %d",mesh->nri);
    if ( mesh->nre ) fprintf(stdout,"  Required %d",mesh->nre);
    fprintf(stdout,"\n");
  }
  if ( mesh->nvn || mesh->ntg) {
    fprintf(stdout,"    Normals    %8d",mesh->nvn);
    fprintf(stdout,"  Tangents   %d\n",mesh->ntg);
  }
  fprintf(stdout,"    Bounding box:  x:[%g  %g]  y:[%g  %g]  z:[%g  %g]\n",
          mesh->xmin,mesh->xmax,mesh->ymin,mesh->ymax,mesh->zmin,mesh->zmax);
}


void meshCoord(pMesh mesh,int displ) {
  pTetra     pt1;
  pPoint     ppt;
  pSolution  ps;
  double    *c1,*c2,*c3,*c4,mul,vold,volf;
  int        k;

  /* default */
  if ( ddebug ) printf("change mesh coord\n");
  mul = 2.0*displ - 1.;
  if ( mesh->dim == 2 ) {   
   for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
	  ps  = &mesh->sol[k];
	  ppt->c[0] = ppt->c[0] + mul*ps->m[0];
	  ppt->c[1] = ppt->c[1] + mul*ps->m[1];
    }
  }
  else {
	vold = 0.0;
	for (k=1; k<=mesh->ntet; k++) {
	  pt1 = &mesh->tetra[k];
	  if ( !pt1->v[0] )  continue;
	  c1 = &mesh->point[pt1->v[0]].c[0];
	  c2 = &mesh->point[pt1->v[1]].c[0];
	  c3 = &mesh->point[pt1->v[2]].c[0];
	  c4 = &mesh->point[pt1->v[3]].c[0];
	  vold += volTet(c1,c2,c3,c4);
	}
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
	  ps  = &mesh->sol[k];
	  ppt->c[0] = mesh->xtra + ppt->c[0] + mul*ps->m[0];
	  ppt->c[1] = mesh->ytra + ppt->c[1] + mul*ps->m[1];
	  ppt->c[2] = mesh->ztra + ppt->c[2] + mul*ps->m[2];
    }
	volf = 0.0;
	for (k=1; k<=mesh->ntet; k++) {
	  pt1 = &mesh->tetra[k];
	  if ( !pt1->v[0] )  continue;
	  c1 = &mesh->point[pt1->v[0]].c[0];
	  c2 = &mesh->point[pt1->v[1]].c[0];
	  c3 = &mesh->point[pt1->v[2]].c[0];
	  c4 = &mesh->point[pt1->v[3]].c[0];
	  volf += volTet(c1,c2,c3,c4);
	}
	fprintf(stdout,"  Volume: initial %E    final %E\n",vold,volf);
  }
}


void meshBox(pMesh mesh,int bb) {
  pPoint     ppt;
  int        k;

  /* default */
  if ( ddebug ) printf("compute mesh box\n");
  if ( bb ) {
    mesh->xmin = mesh->ymin = mesh->zmin =  FLOAT_MAX;
    mesh->xmax = mesh->ymax = mesh->zmax = -FLOAT_MAX;

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( ppt->tag == M_UNUSED && mesh->ne )  continue;
      if ( ppt->c[0] < mesh->xmin ) mesh->xmin = ppt->c[0];
      if ( ppt->c[0] > mesh->xmax ) mesh->xmax = ppt->c[0];
      if ( ppt->c[1] < mesh->ymin ) mesh->ymin = ppt->c[1];
      if ( ppt->c[1] > mesh->ymax ) mesh->ymax = ppt->c[1];
      if ( ppt->c[2] < mesh->zmin ) mesh->zmin = ppt->c[2];
      if ( ppt->c[2] > mesh->zmax ) mesh->zmax = ppt->c[2];
    }
  }

  /* translate mesh at center */
  mesh->xtra = 0.5 * (mesh->xmin+mesh->xmax);
  mesh->ytra = 0.5 * (mesh->ymin+mesh->ymax);
  mesh->ztra = 0.5 * (mesh->zmin+mesh->zmax);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    /*if ( ppt->tag == M_UNUSED && mesh->ne )  continue;*/
    ppt->c[0] -= mesh->xtra; 
    ppt->c[1] -= mesh->ytra; 
    ppt->c[2] -= mesh->ztra;
  }
}


int meshSurf(pMesh mesh) {
  fprintf(stderr,"****meshSurf****\n");
  pTetra     ptt,pt2;
  pHexa      ph;
  pTriangle  pt;
  pQuad      pq;
  int       *adj,i,k,iadr;
  ubyte      i1,i2,i3;
  static int idirt[7] = {0,1,2,3,0,1,2};
  static int ch[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, 
			  {1,2,6,5}, {2,3,7,6}, {0,3,7,4} };

  if ( !dosurf )  return(1);

  /* extract surface */
  if ( mesh->ntet > 0 && !mesh->nt ) {
    if ( !hashTetra(mesh) )  return(0);

    for (k=1; k<=mesh->ntet; k++) {
      ptt = &mesh->tetra[k];
      if ( !ptt->v[0] )  continue;
      iadr = 4*(k-1)+1;
      adj  = &mesh->adja[iadr];
      for (i=0; i<4; i++)
        if ( !adj[i] ) ++mesh->nt;
	else {
	  pt2 = &mesh->tetra[adj[i]];
	  if ( ptt->ref != pt2->ref && k < adj[i] ) ++mesh->nt;
	}
    }

    /* memory alloc */
    if ( ddebug ) printf("allocate %d tria\n",mesh->nt);
    mesh->tria = (pTriangle)calloc(mesh->nt+1,sizeof(struct striangle));
    if ( !mesh->tria ) {
      fprintf(stdout,"  ## No triangle\n");
      return(0);
    }

    /* find triangles */
    mesh->nt = 0;
    for (k=1; k<=mesh->ntet; k++) {
      ptt = &mesh->tetra[k];
      if ( !ptt->v[0] )  continue;
      iadr = 4*(k-1)+1;
      adj  = &mesh->adja[iadr];
      for (i=0; i<4; i++) {
        if ( !adj[i] ) {
          pt = &mesh->tria[++mesh->nt];
          i1 = idirt[i+1];
          i2 = idirt[i+2];
          i3 = idirt[i+3];
          pt->v[0] = ptt->v[i1];
          pt->v[1] = ptt->v[i2];
          pt->v[2] = ptt->v[i3];
	      pt->ref  = 0;
        }
	    else {
	      pt2 = &mesh->tetra[adj[i]];
	      if ( (ptt->ref != pt2->ref) && (k < adj[i]) ) {
            pt = &mesh->tria[++mesh->nt];
            i1 = idirt[i+1];
            i2 = idirt[i+2];
            i3 = idirt[i+3];
            pt->v[0] = ptt->v[i1];
            pt->v[1] = ptt->v[i2];
            pt->v[2] = ptt->v[i3];
	        pt->ref  = MEDIT_MAX(ptt->ref,pt2->ref);
	  } 
	}
      }	
    } 
  }

  if ( mesh->nhex > 0 && !mesh->nq ) {
    if ( mesh->adja ) {
      free(mesh->adja);
      free(mesh->voy);
      mesh->adja = 0;
      mesh->voy  = 0;
    }
    if ( !hashHexa(mesh) )  return(0);

    for (k=1; k<=mesh->nhex; k++) {
      ph = &mesh->hexa[k];
      if ( !ph->v[0] )  continue;
      iadr = 6*(k-1)+1;
      adj  = &mesh->adja[iadr];
      for (i=0; i<6; i++)
        if ( !adj[i] ) ++mesh->nq;
    }

    /* memory alloc */
    if ( ddebug ) printf("allocate %d quads\n",mesh->nq);
    mesh->quad = (pQuad)calloc(mesh->nq+1,sizeof(struct squad));
    if ( !mesh->quad ) {
      fprintf(stdout,"  ## No quad\n");
      return(0);
    }

    /* find quadrilaterals */
    mesh->nq = 0;
    for (k=1; k<=mesh->nhex; k++) {
      ph = &mesh->hexa[k];
      if ( !ph->v[0] )  continue;
      iadr = 6*(k-1)+1;
      adj  = &mesh->adja[iadr];
      for (i=0; i<6; i++)
        if ( !adj[i] ) {
          pq = &mesh->quad[++mesh->nq];
          pq->v[0] = ph->v[ch[i][0]];
          pq->v[1] = ph->v[ch[i][1]];
          pq->v[2] = ph->v[ch[i][2]];
          pq->v[3] = ph->v[ch[i][3]];
        }
    } 
  }

  return(1);
}

void meshRef(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pTriangle  pt;
  pQuad      pq;
  pTetra     pte;
  pHexa      ph;
  pPoint     ppt;
  int       *old,i,k,m,nmat;

  /* default */
  if ( ddebug ) printf("    assign %d references\n",sc->par.nbmat-1);

  /* sort elements by references */
  if ( !quiet ) fprintf(stdout,"   Identifying sub-domains\n");
  old = (int*)calloc(sc->par.nbmat+1,sizeof(int));
  if ( !old )  exit(2);

  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    pm->ext[0] = mesh->xmax-mesh->xtra;
    pm->ext[1] = mesh->ymax-mesh->ytra;
    pm->ext[2] = mesh->zmax-mesh->ztra;
    pm->ext[3] = mesh->xmin-mesh->xtra;
    pm->ext[4] = mesh->ymin-mesh->ytra;
    pm->ext[5] = mesh->zmin-mesh->ztra;
  }

  for (k=mesh->nt; k>0; k--) {
    pt   = &mesh->tria[k];
	if ( !pt->v[0] )  continue;
    nmat = matRef(sc,pt->ref);
    /*nmat = !pt->ref ? DEFAULT_MAT : 1+(pt->ref-1)%(sc->par.nbmat-1);*/
    pt->nxt   = old[nmat];
    old[nmat] = k;
    pm = &sc->material[nmat];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      pm->ext[0] = MEDIT_MIN(pm->ext[0],ppt->c[0]);
      pm->ext[1] = MEDIT_MIN(pm->ext[1],ppt->c[1]);
      pm->ext[2] = MEDIT_MIN(pm->ext[2],ppt->c[2]);
      pm->ext[3] = MEDIT_MAX(pm->ext[3],ppt->c[0]);
      pm->ext[4] = MEDIT_MAX(pm->ext[4],ppt->c[1]);
      pm->ext[5] = MEDIT_MAX(pm->ext[5],ppt->c[2]);
    }
  }

  for (m=0;m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    pm->depmat[LTria] = old[m];
    old[m] = 0;
  }

  for (k=mesh->nq; k>0; k--) {
    pq   = &mesh->quad[k];
	if ( !pq->v[0] )  continue;
    nmat = matRef(sc,pq->ref);
    /*nmat = !pq->ref ? DEFAULT_MAT : 1+(pq->ref-1)%(sc->par.nbmat-1);*/
    pq->nxt   = old[nmat];
    old[nmat] = k;
    pm = &sc->material[nmat];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pq->v[i]];
      pm->ext[0] = MEDIT_MIN(pm->ext[0],ppt->c[0]);
      pm->ext[1] = MEDIT_MIN(pm->ext[1],ppt->c[1]);
      pm->ext[2] = MEDIT_MIN(pm->ext[2],ppt->c[2]);
      pm->ext[3] = MEDIT_MAX(pm->ext[3],ppt->c[0]);
      pm->ext[4] = MEDIT_MAX(pm->ext[4],ppt->c[1]);
      pm->ext[5] = MEDIT_MAX(pm->ext[5],ppt->c[2]);
    }
  }

  for (m=0;m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    pm->depmat[LQuad] = old[m];
    old[m] = 0;
  }

  for (k=mesh->ntet; k>0; k--) {
    pte  = &mesh->tetra[k];
	if ( !pte->v[0] )  continue;
    nmat = matRef(sc,pte->ref);
    /*nmat = !pte->ref ? DEFAULT_MAT : 1+(pte->ref-1)%(sc->par.nbmat-1);*/
    pte->nxt  = old[nmat];
    old[nmat] = k;
    pm = &sc->material[nmat];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pte->v[i]];
      pm->ext[0] = MEDIT_MIN(pm->ext[0],ppt->c[0]);
      pm->ext[1] = MEDIT_MIN(pm->ext[1],ppt->c[1]);
      pm->ext[2] = MEDIT_MIN(pm->ext[2],ppt->c[2]);
      pm->ext[3] = MEDIT_MAX(pm->ext[3],ppt->c[0]);
      pm->ext[4] = MEDIT_MAX(pm->ext[4],ppt->c[1]);
      pm->ext[5] = MEDIT_MAX(pm->ext[5],ppt->c[2]);
    }
  }

  for (m=0;m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    pm->depmat[LTets] = old[m];
    old[m] = 0;
  }

  for (k=mesh->nhex; k>0; k--) {
    ph   = &mesh->hexa[k];
    nmat = matRef(sc,ph->ref);
    /*nmat = !ph->ref ? DEFAULT_MAT : 1+(ph->ref-1)%(sc->par.nbmat-1);*/
    ph->nxt   = old[nmat];
    old[nmat] = k;
    pm = &sc->material[nmat];
    for (i=0; i<8; i++) {
      ppt = &mesh->point[ph->v[i]];
      pm->ext[0] = MEDIT_MIN(pm->ext[0],ppt->c[0]);
      pm->ext[1] = MEDIT_MIN(pm->ext[1],ppt->c[1]);
      pm->ext[2] = MEDIT_MIN(pm->ext[2],ppt->c[2]);
      pm->ext[3] = MEDIT_MAX(pm->ext[3],ppt->c[0]);
      pm->ext[4] = MEDIT_MAX(pm->ext[4],ppt->c[1]);
      pm->ext[5] = MEDIT_MAX(pm->ext[5],ppt->c[2]);
    }
  }
  for (m=0;m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    pm->depmat[LHexa] = old[m];
    old[m] = 0;
  }

  free(old);

  /* remove unused materials */
/*
  for (m=1; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    if ( !pm->depmat[LTria] && !pm->depmat[LQuad] &&
         !pm->depmat[LTets] && !pm->depmat[LHexa] ) {
      pm1 = &sc->material[sc->par.nbmat-1];
      memcpy(pm,pm1,sizeof(struct material));
      sc->par.nbmat--;
      m--;
    }
  }
*/
  if ( ddebug )
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      if ( pm->depmat[LTria] || pm->depmat[LQuad] || pm->depmat[LTets] ||
           pm->depmat[LHexa] )
        fprintf(stdout,"    depart[%d], ref %d = %d %d %d %d\n",
                m,pm->ref,pm->depmat[LTria],pm->depmat[LQuad],
                pm->depmat[LTets],pm->depmat[LHexa]);
    }
}

int meshUpdate(pScene sc,pMesh mesh) {
  int     k,ret;
  clock_t  ct;

  /* release mesh structure */
  if ( ddebug ) fprintf(stdout,"loadMesh: update mesh\n");
  if ( mesh->tria )  { M_free(mesh->tria);  mesh->tria  = 0; }
  if ( mesh->quad )  { M_free(mesh->quad);  mesh->quad  = 0; }
  if ( mesh->edge )  { M_free(mesh->edge);  mesh->edge  = 0; }
  if ( mesh->tetra ) { M_free(mesh->tetra); mesh->tetra = 0; }
  if ( mesh->hexa )  { M_free(mesh->hexa);  mesh->hexa  = 0; }
  if ( mesh->adja )  { M_free(mesh->adja);  mesh->adja  = 0; }
  if ( mesh->voy  )  { M_free(mesh->voy);   mesh->voy   = 0; }
  if ( mesh->point ) { M_free(mesh->point); mesh->point = 0; }
  if ( mesh->extra ) {
    if ( mesh->extra->iv )  M_free(mesh->extra->nv);
    if ( mesh->extra->it )  M_free(mesh->extra->nt);
    if ( mesh->extra->iq )  M_free(mesh->extra->nq);
    if ( mesh->extra->n )   M_free(mesh->extra->n);
    M_free(mesh->extra);
    mesh->extra = 0;
  }
  if ( mesh->nbb && mesh->sol ) {
    if ( (mesh->dim==2 && mesh->nfield==3) ||
         (mesh->dim==3 && mesh->nfield==6) )
      for (k=1; k<=mesh->nbb; k++)
        free(mesh->sol[k].m);
    M_free(mesh->sol);
    mesh->sol = 0;
  }

  /* read mesh */
  fprintf(stdout," Loading data file(s)\n");
  ct = clock();
  ret = 0;
  switch ( mesh->typ ) {
  case 0:
    ret = loadMesh(mesh);
    break;
  case 1:
    ret = inmsh2(mesh);
    break;
  case 2:
    ret = loadGIS(mesh);
    break;
  }
  if ( !ret )  return(0);

  /* compute mesh box */
  if ( (mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq) )  
    meshSurf(mesh);
  meshBox(mesh,1);
  /*parsop(sc,mesh);*/
  if ( !quiet )  meshInfo(mesh);

  /* read metric */
  if ( !loadSol(mesh,mesh->name,1) )
    bbfile(mesh);
  if ( !quiet && mesh->nbb )
    fprintf(stdout,"    Solutions  %8d\n",mesh->nbb);

  ct = difftime(clock(),ct);
  fprintf(stdout,"  Input seconds:     %.2f\n",
          (double)ct/(double)CLOCKS_PER_SEC);

  return(1);
}
