#include "medit.h"
#ifdef IGL
extern "C"{
#endif
#include "libmesh5.h"
#ifdef IGL
}
#endif
#include "extern.h"
#ifdef IGL
#include "eigenv.h"
#endif


int loadMesh(pMesh mesh) {
  pPoint      ppt;
  pEdge       pr;
  pTriangle   pt;
  pQuad       pq;
  pTetra      ptet;
  pHexa       ph;
  double      d,dp1,dp2,dp3,dn[3];
  float      *n,fp1,fp2,fp3;
  int         i,ia,ib,inm,ref,is,k,disc,nn,nt,nq;
  char       *ptr,data[256];

  printf("use loadMesh\n");
  /* default */
  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(-1);
      }
    }
  }
  else if (!(inm = GmfOpenMesh(data,GmfRead,&mesh->ver,&mesh->dim)) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return(-1);
  }
  if ( !quiet )  fprintf(stdout,"  Reading %s\n",data);

  /* parse keywords */
  mesh->np   = GmfStatKwd(inm,GmfVertices);
  printf("mesh->np=%i\n",mesh->np);
  mesh->nt   = GmfStatKwd(inm,GmfTriangles);
  printf("mesh->nt=%i\n",mesh->nt);
  mesh->nq   = GmfStatKwd(inm,GmfQuadrilaterals);
  mesh->ntet = GmfStatKwd(inm,GmfTetrahedra);
  mesh->nhex = GmfStatKwd(inm,GmfHexahedra);
  mesh->nc   = GmfStatKwd(inm,GmfCorners);
  mesh->nr   = GmfStatKwd(inm,GmfRequiredVertices);
  mesh->na   = GmfStatKwd(inm,GmfEdges);
  mesh->nri  = GmfStatKwd(inm,GmfRidges);
  mesh->nre  = GmfStatKwd(inm,GmfRequiredEdges);
  mesh->nvn  = GmfStatKwd(inm,GmfNormals);
  mesh->ntg  = GmfStatKwd(inm,GmfTangents);
  mesh->ne = mesh->nt + mesh->nq + mesh->ntet + mesh->nhex;

  /* check space dimension */
  if ( mesh->dim < 2 || mesh->dim > 3 ) {
	fprintf(stdout,"  ## Wrong dim\n");
    GmfCloseMesh(inm);
    return(-1);
  }
  /* check if vertices and elements found */
  if ( !mesh->np ) {
    fprintf(stdout,"  ## No vertex\n");
    GmfCloseMesh(inm);
    return(-1);
  }

  /* memory allocation for mesh */
  if ( !zaldy1(mesh) ) {
    GmfCloseMesh(inm);
    return(-1);
  }

  /* read mesh vertices */
  GmfGotoKwd(inm,GmfVertices);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( mesh->ver == GmfFloat ) {
      if ( mesh->dim == 2 ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
  	    ppt->c[0] = fp1;
  	    ppt->c[1] = fp2;
	    printf("vertices %f %f %i\n",fp1,fp2,ref);
      }
      else {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
  	    ppt->c[0] = fp1;
  	    ppt->c[1] = fp2;
	    ppt->c[2] = fp3;
      }
    }
    else {
      if ( mesh->dim == 2 ) {
        GmfGetLin(inm,GmfVertices,&dp1,&dp2,&ref);
        ppt->c[0] = dp1;
        ppt->c[1] = dp2;
      }
      else {
        GmfGetLin(inm,GmfVertices,&dp1,&dp2,&dp3,&ref);
        ppt->c[0] = dp1;
        ppt->c[1] = dp2;
        ppt->c[2] = dp3;
      }
    }
    ppt->ref = ref & 0x7fff;
    ppt->tag = M_UNUSED;
  }

  /* read mesh triangles */
  disc = 0 ;
  GmfGotoKwd(inm,GmfTriangles);
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    GmfGetLin(inm,GmfTriangles,&pt->v[0],&pt->v[1],&pt->v[2],&ref);
    pt->ref  = ref & 0x7fff;
    for (i=0; i<3; i++) {    
      if ( pt->v[i] < 1 || pt->v[i] > mesh->np ) {
        disc++;
    pt->v[0] = 0;
        break;
      }
      else {
        ppt = &mesh->point[pt->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  /* mesh quadrilaterals */
  GmfGotoKwd(inm,GmfQuadrilaterals);
  for (k=1; k<=mesh->nq; k++) {
    pq = &mesh->quad[k];
    GmfGetLin(inm,GmfQuadrilaterals,&pq->v[0],&pq->v[1],&pq->v[2],&pq->v[3],&ref);
    pq->ref  = ref & 0x7fff;
    for (i=0; i<4; i++) {    
      if ( pq->v[i] < 1 || pq->v[i] > mesh->np ) {
        disc++;
        pq->v[0] = 0;
        break;
      }
      else {
        ppt = &mesh->point[pq->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  /* mesh tetrahedra */
  GmfGotoKwd(inm,GmfTetrahedra);
  for (k=1; k<=mesh->ntet; k++) {
    ptet = &mesh->tetra[k];
    GmfGetLin(inm,GmfTetrahedra,&ptet->v[0],&ptet->v[1],&ptet->v[2],&ptet->v[3],&ref);
    ptet->ref  = ref & 0x7fff;
    for (i=0; i<4; i++) {    
      if ( ptet->v[i] < 1 || ptet->v[i] > mesh->np ) {
        disc++;
        ptet->v[0] = 0;
        break;
      }
      else {
        ppt = &mesh->point[ptet->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  /* mesh hexahedra */
  GmfGotoKwd(inm,GmfHexahedra);
  for (k=1; k<=mesh->nhex; k++) {
    ph = &mesh->hexa[k];
    GmfGetLin(inm,GmfHexahedra,&ph->v[0],&ph->v[1],&ph->v[2],&ph->v[3],
              &ph->v[4],&ph->v[5],&ph->v[6],&ph->v[7],&ref);
    ph->ref  = ref & 0x7fff;
    for (i=0; i<8; i++) {    
      if ( ph->v[i] < 1 || ph->v[i] > mesh->np ) {
        disc++;
        ph->v[0] = 0;
        break;
      }
      else {
        ppt = &mesh->point[ph->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  /* mesh corners */
  GmfGotoKwd(inm,GmfCorners);
  for (k=1; k<=mesh->nc; k++) {
    GmfGetLin(inm,GmfCorners,&is);
    if ( is < 1 || is > mesh->np )
      disc++;
    else {
      ppt = &mesh->point[is];
      ppt->tag |= M_CORNER;
      ppt->tag &= ~M_UNUSED;
    }
  }

  /* required vertices */
  GmfGotoKwd(inm,GmfRequiredVertices);
  for (k=1; k<=mesh->nr; k++) {
    GmfGetLin(inm,GmfRequiredVertices,&is);
    if ( is < 1 || is > mesh->np )
      disc++;
    else {
      ppt = &mesh->point[is];
      ppt->tag |= M_REQUIRED;
      ppt->tag &= ~M_UNUSED;
    }
  }

  /*mesh edges */
  GmfGotoKwd(inm,GmfEdges);
  for (k=1; k<=mesh->na; k++) {
    GmfGetLin(inm,GmfEdges,&ia,&ib,&ref);
    if ( ia < 1 || ia > mesh->np || ib < 1 || ib > mesh->np )
      disc++;
    else {
      pr = &mesh->edge[k];
      pr->v[0] = ia;
      pr->v[1] = ib;
      pr->ref  = ref & 0x7fff;
      pr->tag  = !pr->ref ? M_NOTAG : M_TAG;
      ppt = &mesh->point[ia];
      ppt->tag &= ~M_UNUSED;
      ppt = &mesh->point[ib];
      ppt->tag &= ~M_UNUSED;
    }
  }

  /* mesh ridges */
  GmfGotoKwd(inm,GmfRidges);
  for (k=1; k<=mesh->nri; k++) {
    GmfGetLin(inm,GmfRidges,&is);
    if ( is < 1 || is > mesh->na )
      disc++;
    else {
      pr = &mesh->edge[is];
      pr->tag |= M_RIDGE;
    }
  }

  /* required edges */
  GmfGotoKwd(inm,GmfRequiredEdges);
  for (k=1; k<=mesh->nre; k++) {
    GmfGetLin(inm,GmfRequiredEdges,&is);
    if ( is < 1 || is > mesh->na )
      disc++;
    else {
      pr = &mesh->edge[is];
      pr->tag |= M_REQUIRED;
    }
  }


  /* mesh normals */
  GmfGotoKwd(inm,GmfNormals);
  for (k=1; k<=mesh->nvn; k++) {
    n = &mesh->extra->n[3*(k-1)+1];
    if ( mesh->ver == GmfFloat )
      GmfGetLin(inm,GmfNormals,&n[0],&n[1],&n[2]);
    else {
      GmfGetLin(inm,GmfNormals,&dn[0],&dn[1],&dn[2]);
      n[0] = dn[0];
      n[1] = dn[1];
      n[2] = dn[2];
    }
    d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if ( d > 0.0 ) {
      d = 1.0 / sqrt(d);
      n[0] *= d;
      n[1] *= d;
      n[2] *= d;
    }
  }

  if ( mesh->nvn ) {
    /* normals at vertices */
    mesh->extra->iv = GmfStatKwd(inm,GmfNormalAtVertices);
    mesh->extra->nv = (int*)M_calloc(mesh->np+1,sizeof(int),"inmesh");
    assert(mesh->extra->nv);
    GmfGotoKwd(inm,GmfNormalAtVertices);
    for (k=1; k<=mesh->extra->iv; k++) {
      GmfGetLin(inm,GmfNormalAtVertices,&nn,&is);
      if ( nn < 1 || nn > mesh->np )
        disc++;
      else
    mesh->extra->nv[nn] = is;
    }

    /* normals at triangle vertices */
    mesh->extra->it = GmfStatKwd(inm,GmfNormalAtTriangleVertices);
    mesh->extra->nt = (int*)M_calloc(3*mesh->nt+1,sizeof(int),"inmesh");
    assert(mesh->extra->nt);
    GmfGotoKwd(inm,GmfNormalAtTriangleVertices);
    for (k=1; k<=mesh->extra->it; k++) {
      GmfGetLin(inm,GmfNormalAtTriangleVertices,&nt,&is,&nn);
      if ( nt < 1 || nt > mesh->nt || is < 1 || is > 3 || nn < 1 || nn > mesh->nvn )
        disc++;
      else
    mesh->extra->nt[3*(nt-1)+is] = nn;
    }

    /*normals at quadrilateral vertices */
    mesh->extra->iq = GmfStatKwd(inm,GmfNormalAtQuadrilateralVertices);
    mesh->extra->nq = (int*)M_calloc(4*mesh->nq+1,sizeof(int),"inmesh");
    assert(mesh->extra->nq);
    GmfGotoKwd(inm,GmfNormalAtQuadrilateralVertices);
    for (k=1; k<=mesh->extra->iq; k++) {
      GmfGetLin(inm,GmfNormalAtQuadrilateralVertices,&nq,&is,&nn);
      if ( nq < 1 || nq > mesh->nq || is < 1 || is > 4 || nn < 1 || nn > mesh->nvn )
        disc++;
      else
    mesh->extra->nq[3*(nq-1)+is] = nn;
    }
  }

  /*mesh tangents */
  GmfGotoKwd(inm,GmfTangents);
  for (k=1; k<=mesh->ntg; k++) {
    n = &mesh->extra->t[3*(k-1)+1];
    if ( mesh->ver == GmfFloat )
      GmfGetLin(inm,GmfTangents,&n[0],&n[1],&n[2]);
    else {
      GmfGetLin(inm,GmfTangents,&dn[0],&dn[1],&dn[2]);
      n[0] = dn[0];
      n[1] = dn[1];
      n[2] = dn[2];
    }
    d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if ( d > 0.0 ) {
      d = 1.0 / sqrt(d);
      n[0] *= d;
      n[1] *= d;
      n[2] *= d;
    }
  }

  if (mesh->ntg ) {
    /*tangents at vertices */
    mesh->extra->jv = GmfStatKwd(inm,GmfTangentAtVertices);
    mesh->extra->tv = (int*)M_calloc(mesh->np+1,sizeof(int),"inmesh");
    assert(mesh->extra->tv);
    GmfGotoKwd(inm,GmfTangentAtVertices);
    for (k=1; k<=mesh->extra->jv; k++) {
      GmfGetLin(inm,GmfTangentAtVertices,&nn,&is);
      if ( nn < 1 || nn > mesh->np )
        disc++;
      else
    mesh->extra->tv[nn] = is;
    }
    
    /* tangent at edge vertices */
    mesh->extra->je = GmfStatKwd(inm,GmfTangentAtEdgeVertices);
    mesh->extra->te = (int*)M_calloc(2*mesh->na+1,sizeof(int),"inmesh");
    assert(mesh->extra->te);
    GmfGotoKwd(inm,GmfTangentAtEdgeVertices);
    for (k=1; k<=mesh->extra->je; k++) {
      GmfGetLin(inm,GmfTangentAtEdgeVertices,&nt,&is,&nn);
      if ( nt < 1 || nt > mesh->np || is < 1 || is > 2 || nn < 1 || nn > mesh->ntg )
        disc++;
      else
    mesh->extra->te[3*(nt-1)+is] = is;
    }
  }

  GmfCloseMesh(inm);
  if ( disc > 0 ) {
    fprintf(stdout,"  ## %d entities discarded\n",disc);
  }
  return(1);
}


/*mark clipped elements */
static int markPt(pMesh mesh) {
  pTriangle  pt;
  pQuad      pq;
  pPoint     ppt;
  int        i,k,pos,neg,nul;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_UNUSED )  continue;
    ppt->tmp = 0;
  }

  for (k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    pos = neg = nul = 0;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->clip == 2 )       pos++;
      else if ( ppt->clip == 1 )  neg++;
      else                        nul++;
    }
    if ( pos && pos+nul < 4 ) {
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt->v[i]];
        ppt->tmp = 1;
      }
    }
  }

  for (k=1; k<=mesh->nq; k++ ) {
    pq = &mesh->quad[k];
    if ( !pq->v[0] )  continue;
    pos = neg = nul = 0;
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pq->v[i]];
      if ( ppt->clip == 2 )       pos++;
      else if ( ppt->clip == 1 )  neg++;
      else                        nul++;
    }
    if ( pos && pos+nul < 5 ) {
      for (i=0; i<4; i++) {
        ppt = &mesh->point[pq->v[i]];
        ppt->tmp = 1;
      }
    }
  }

  return(1);
}


/* save (part of) mesh to disk */
int saveMesh(pScene sc,pMesh mesh,char *fileout,ubyte clipon) {
  pPoint     ppt;
  pTriangle  pt;
  pQuad      pq;
  pTetra     ptt;
  pMaterial  pm;
  float      fp1,fp2,fp3;
  int        outm,i,k,m,ver,np,nt,nq,ref;

  ver = GmfFloat;
  if ( !(outm = GmfOpenMesh(fileout,GmfWrite,ver,mesh->dim)) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",fileout);
    return(0);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",fileout);

  /*compact vertices */
  if ( clipon )  
    markPt(mesh);
  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_UNUSED )  continue;
    ppt->tmp = ++np;
  }

  GmfSetKwd(outm,GmfVertices,np);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & M_UNUSED )  continue;
    ref = ppt->ref;
    if ( mesh->dim == 2 ) {
      fp1 = ppt->c[0] + mesh->xtra;
      fp2 = ppt->c[1] + mesh->ytra;
      ref = ppt->ref;
      GmfSetLin(outm,GmfVertices,fp1,fp2,ref);
    }
    else {
      fp1 = ppt->c[0] + mesh->xtra;
      fp2 = ppt->c[1] + mesh->ytra;
      fp3 = ppt->c[2] + mesh->ztra;
      ref = ppt->ref;
      GmfSetLin(outm,GmfVertices,fp1,fp2,fp3,ref);
    }
  }

  /* write triangles */
  nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    m  = matRef(sc,pt->ref);
    pm = &sc->material[m];
    if ( pm->flag )  continue;
    for (i=0; i<3; i++)
      if ( !mesh->point[pt->v[i]].tmp )  break;
    if ( i == 3 )  nt++;
  }
  GmfSetKwd(outm,GmfTriangles,nt);
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    m  = matRef(sc,pt->ref);
    pm = &sc->material[m];
    if ( pm->flag )  continue;
    for (i=0; i<3; i++)
      if ( !mesh->point[pt->v[i]].tmp )  break;
    if ( i < 3 )  continue;
    ref = pt->ref;
    GmfSetLin(outm,GmfTriangles,mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp,
              mesh->point[pt->v[2]].tmp,ref);
  }

  /* write quads */
  nq = 0;
  for (k=1; k<=mesh->nq; k++) {
    pq = &mesh->quad[k];
    if ( !pq->v[0] )  continue;
    m  = matRef(sc,pq->ref);
    pm = &sc->material[m];
    if ( pm->flag )  continue;
    for (i=0; i<4; i++)
      if ( !mesh->point[pq->v[i]].tmp )  break;
    if ( i == 4 )  nq++;
  }
  GmfSetKwd(outm,GmfQuadrilaterals,nq);
  for (k=1; k<=mesh->nq; k++) {
    pq = &mesh->quad[k];
    if ( !pq->v[0] )  continue;
    m  = matRef(sc,pq->ref);
    pm = &sc->material[m];
    if ( pm->flag )  continue;
    for (i=0; i<4; i++)
      if ( !mesh->point[pq->v[i]].tmp )  break;
    if ( i < 4 )  continue;
    ref = pq->ref;
    GmfSetLin(outm,GmfQuadrilaterals,mesh->point[pq->v[0]].tmp,mesh->point[pq->v[1]].tmp,
              mesh->point[pq->v[2]].tmp,mesh->point[pq->v[3]].tmp,ref);
  }

  /* write tetrahedra */
  GmfSetKwd(outm,GmfTetrahedra,mesh->ntet);
  for (k=1; k<=mesh->ntet; k++) {
    ptt = &mesh->tetra[k];
    if ( !ptt->v[0] )  continue;
    m  = matRef(sc,ptt->ref);
    pm = &sc->material[m];
    if ( pm->flag )  continue;
    for (i=0; i<4; i++)
      if ( !mesh->point[ptt->v[i]].tmp )  break;
    if ( i < 4 )  continue;
    ref = ptt->ref;
    GmfSetLin(outm,GmfTetrahedra,mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp,
              mesh->point[ptt->v[2]].tmp,mesh->point[ptt->v[3]].tmp,ref);
  }

  /* write hexahedra */

  if ( !quiet ) {
    fprintf(stdout,"     TOTAL NUMBER OF VERTICES   %8d\n",np);
    fprintf(stdout,"     TOTAL NUMBER OF TRIANGLES  %8d\n",nt);
    fprintf(stdout,"     TOTAL NUMBER OF QUADS      %8d\n",nq);
    fprintf(stdout,"     TOTAL NUMBER OF TETRA      %8d\n",mesh->ntet);
  }

  GmfCloseMesh(outm);
  return(1);
}


/* load solution (metric) */
int loadSol(pMesh mesh,char *filename,int numsol) {
  pSolution    sol;
  double       dbuf[ GmfMaxTyp ];
  float        fbuf[ GmfMaxTyp ];
  double       m[6],lambda[3],eigv[3][3],vp[2][2];
  int          inm,k,i,key,nel,size,type,iord,off,typtab[GmfMaxTyp],ver,dim;
  char        *ptr,data[128];

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".solb");
  if (!(inm = GmfOpenMesh(data,GmfRead,&ver,&dim)) ) {
    ptr  = strstr(data,".sol");
    *ptr = '\0';
    strcat(data,".sol");
    if (!(inm = GmfOpenMesh(data,GmfRead,&ver,&dim)) )
      return(0);
  }
  if ( !quiet )  fprintf(stdout,"  Reading %s\n",data);

  if ( dim != mesh->dim ) {
    fprintf(stderr,"  %%%% Wrong dimension %d.\n",dim);
    GmfCloseMesh(inm);
    return(0);
  }

  nel = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,typtab);
  if ( nel ) {
    if ( nel > mesh->np ) {
      fprintf(stderr,"  %%%% Wrong number: %d Solutions discarded\n",nel-mesh->np);
    }
    mesh->typage = 2;
    key = GmfSolAtVertices;
  }
  else {
    mesh->typage = 1;
    if ( mesh->dim == 2 && mesh->nt ) {
      //if( mesh->nt){
      nel = GmfStatKwd(inm,GmfSolAtTriangles,&type,&size,typtab);
      if ( nel && nel != mesh->nt ) {
        fprintf(stderr,"  %%%% Wrong number %d.\n",nel);
        GmfCloseMesh(inm);
        return(0);
      }
      key = GmfSolAtTriangles;
    }
    else if ( mesh->dim == 2 && mesh->nq ) {
      nel = GmfStatKwd(inm,GmfSolAtQuadrilaterals,&type,&size,typtab);
      if ( nel && nel != mesh->nq ) {
        fprintf(stderr,"  %%%% Wrong number %d.\n",nel);
        GmfCloseMesh(inm);
        return(0);
      }
      key = GmfSolAtQuadrilaterals;
    }
    else if ( mesh->ntet ) {
      nel = GmfStatKwd(inm,GmfSolAtTetrahedra,&type,&size,typtab);
      if ( nel && nel != mesh->ntet ) {
        fprintf(stderr,"  %%%% Wrong number %d.\n",nel);
        GmfCloseMesh(inm);
        return(0);
      }
      key = GmfSolAtTetrahedra;
    }
    else if ( mesh->nhex ) {
      nel = GmfStatKwd(inm,GmfSolAtHexahedra,&type,&size,typtab);
      if ( nel && nel != mesh->nhex ) {
        fprintf(stderr,"  %%%% Wrong number %d.\n",nel);
        GmfCloseMesh(inm);
        return(0);
      }
      key = GmfSolAtHexahedra;
    }
  }
  if ( !nel )  return(1);
  if(ddebug) printf("numsol=%i, type=%i, size=%i\n",numsol,type,size); 
  if ( numsol > type )  numsol = 1;
  numsol--;
  mesh->nbb    = nel;
  mesh->bbmin  =  1.0e20;
  mesh->bbmax  = -1.0e20;
  if ( !zaldy2(mesh) ) {
    GmfCloseMesh(inm);
    return(0);
  }
  sol = mesh->sol;
  sol->dim = dim;
  sol->ver = ver;

  off = 0;
  for (i=0; i<numsol; i++) {
    if(ddebug) printf("typtab[%i]=%i",i,typtab[i]); 
    switch(typtab[i]) {
      case GmfSca:  
        off++; break;
      case GmfVec:
        off += sol->dim; break;
      case GmfSymMat:
        off += sol->dim*(sol->dim+1)/2;  break;
    }
  }
  if(ddebug) printf("typtab[%i]=%i, off%i",i,typtab[i],off); 

  GmfGotoKwd(inm,key);
  if(ddebug) printf("numsol=%i,typtab[i]=%i\n",numsol,typtab[i]); 
  switch(typtab[numsol]) {
    case GmfSca:
      if ( ddebug )  printf("   scalar field\n");
      mesh->nfield = 1;
      for (k=1; k<=nel; k++) {
        if ( sol->ver == GmfFloat )
          GmfGetLin(inm,key,fbuf);
        else {
          GmfGetLin(inm,key,dbuf);
          for (i=0; i<GmfMaxTyp; i++)
            fbuf[i] = dbuf[off+i];
        }
        mesh->sol[k].bb = fbuf[off];
	if(ddebug) printf("valeur donneer %i %f\n",k,fbuf[off]);
	if ( mesh->sol[k].bb < mesh->bbmin )  mesh->bbmin = mesh->sol[k].bb;
        if ( mesh->sol[k].bb > mesh->bbmax )  mesh->bbmax = mesh->sol[k].bb;
      }
      break;

    case GmfVec:
      if ( ddebug )  printf("   vector field\n");
      mesh->nfield = sol->dim;
      for (k=1; k<=nel; k++) {
	mesh->sol[k].bb = 0.0;
        if ( sol->ver == GmfFloat )
          GmfGetLin(inm,key,fbuf);
        else {
          GmfGetLin(inm,key,dbuf);
          for (i=0; i<GmfMaxTyp; i++)
            fbuf[i] = dbuf[off+i];
        }
        for (i=0; i<sol->dim; i++) {
          mesh->sol[k].m[i] = fbuf[off+i];
	  if(ddebug) printf("valeur donner %i composante %i %f\n",k,i,fbuf[off+i]);
          mesh->sol[k].bb  += fbuf[off+i]*fbuf[off+i];
        }
        mesh->sol[k].bb = sqrt(mesh->sol[k].bb);
        if ( mesh->sol[k].bb < mesh->bbmin )  mesh->bbmin = mesh->sol[k].bb;
        if ( mesh->sol[k].bb > mesh->bbmax )  mesh->bbmax = mesh->sol[k].bb;
      }
      break;

    case GmfSymMat:
      if ( ddebug )  printf("   metric field\n");
      mesh->nfield = sol->dim*(sol->dim+1) / 2;
      for (k=1; k<=nel; k++) {
        if ( sol->ver == GmfFloat )
          GmfGetLin(inm,key,fbuf);
        else {
          GmfGetLin(inm,key,dbuf);
          for (i=0; i<GmfMaxTyp; i++)
            fbuf[i] = dbuf[off+i];
        }
        
        if ( sol->dim == 2 ) {
          for (i=0; i<3; i++)  
            mesh->sol[k].m[i] = m[i] = fbuf[off+i];
          iord = eigen2(m,lambda,vp);
          mesh->sol[k].bb = MEDIT_MIN(lambda[0],lambda[1]);
          if ( mesh->sol[k].bb < mesh->bbmin )  mesh->bbmin = mesh->sol[k].bb;
          if ( mesh->sol[k].bb > mesh->bbmax )  mesh->bbmax = mesh->sol[k].bb;
        }
        else {
          for (i=0; i<6; i++)  
            mesh->sol[k].m[i] = fbuf[off+i];
          mesh->sol[k].m[2] = fbuf[off+3];
          mesh->sol[k].m[3] = fbuf[off+2];
          for (i=0; i<6; i++)  m[i] = mesh->sol[k].m[i];
          iord = eigenv(1,m,lambda,eigv);
          if ( iord ) {
            mesh->sol[k].bb = lambda[0];
            mesh->sol[k].bb = MEDIT_MAX(mesh->sol[k].bb,lambda[1]);
            mesh->sol[k].bb = MEDIT_MAX(mesh->sol[k].bb,lambda[2]);
            if ( mesh->sol[k].bb < mesh->bbmin )  mesh->bbmin = mesh->sol[k].bb;
            if ( mesh->sol[k].bb > mesh->bbmax )  mesh->bbmax = mesh->sol[k].bb;
          }
        }
    }
    break;
  }

  GmfCloseMesh(inm);
  return(1);  
}


