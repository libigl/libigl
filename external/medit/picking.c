#include "medit.h"
#include "extern.h"
#include "sproto.h"


typedef struct color {
  GLuint  rMask,gMask,bMask,aMask;
  int     rShift,gShift,bShift;
  int     rBits,gBits,bBits,aBits;
} Color;


int    refmat=-1,reftype=-1,refitem=0,numel=0,refval=0;

static int ch[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, 
			{1,2,6,5}, {2,3,7,6}, {0,3,7,4} };
static int ct[4][3] = { {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3} };
extern int refpick;


static void drawTria(pScene sc,pMesh mesh,int k) {
  pMaterial    pm;
  pTriangle    pt;
  pPoint       p0,p1,p2;
  double       ax,ay,az,bx,by,bz,dd;
  float        shrink,cx,cy,cz,n[3];

  /* default */
  if ( ddebug ) printf("draw triangle %d\n",k);
  if ( k < 1 || k > mesh->nt )  return;

  pt = &mesh->tria[k];
  if ( refpick > 0 )  pt->ref = refpick;
  refmat = matRef(sc,pt->ref);
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  pm = &sc->material[refmat];
  cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.;
  cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.;
  cz = (p0->c[2] + p1->c[2] + p2->c[2]) / 3.;
  shrink = 0.95 * sc->shrink;

  glBegin(GL_TRIANGLES);
  glColor3f(1.0-pm->dif[0],1.0-pm->dif[1],1.0-pm->dif[2]);

  /* compute normal */
  ax = p1->c[0] - p0->c[0];
  ay = p1->c[1] - p0->c[1];
  az = p1->c[2] - p0->c[2];
  bx = p2->c[0] - p0->c[0];
  by = p2->c[1] - p0->c[1];
  bz = p2->c[2] - p0->c[2];
  n[0] = ay*bz - az*by;
  n[1] = az*bx - ax*bz;
  n[2] = ax*by - ay*bx;
  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd > 0.0f ) {
    dd = 1.0f / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  }
  glNormal3fv(n);
  glVertex3f(shrink*(p0->c[0]-cx)+cx,
             shrink*(p0->c[1]-cy)+cy,
             shrink*(p0->c[2]-cz)+cz);
  glNormal3fv(n);
  glVertex3f(shrink*(p1->c[0]-cx)+cx,
             shrink*(p1->c[1]-cy)+cy,
             shrink*(p1->c[2]-cz)+cz);
  glNormal3fv(n);
  glVertex3f(shrink*(p2->c[0]-cx)+cx,
             shrink*(p2->c[1]-cy)+cy,
             shrink*(p2->c[2]-cz)+cz);
  glEnd();

  glColor3f(1.0-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
  output3(p0->c[0],p0->c[1],p0->c[2],"%d",pt->v[0]);
  output3(p1->c[0],p1->c[1],p1->c[2],"%d",pt->v[1]);
  output3(p2->c[0],p2->c[1],p2->c[2],"%d",pt->v[2]);
}

static void drawQuad(pScene sc,pMesh mesh,int k) {
  pMaterial  pm;
  pQuad      pq;
  pPoint     p0,p1,p2,p3;
  double     ax,ay,az,bx,by,bz,dd;
  float      shrink,cx,cy,cz,n[3];

  /* default */
  if ( ddebug ) printf("draw quad %d\n",k);
  if ( k < 1 || k > mesh->nq )  return;
  pq = &mesh->quad[k];
  if ( refpick > 0 )  pq->ref = refpick;
  refmat = matRef(sc,pq->ref);

  p0 = &mesh->point[pq->v[0]];
  p1 = &mesh->point[pq->v[1]];
  p2 = &mesh->point[pq->v[2]];
  p3 = &mesh->point[pq->v[3]];
  pm = &sc->material[refmat];
  cx = (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]) / 4.;
  cy = (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]) / 4.;
  cz = (p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]) / 4.;
  shrink = 0.95 * sc->shrink;

  glBegin(GL_QUADS);
  glColor3f(1.0-pm->dif[0],1.0-pm->dif[1],1.0-pm->dif[2]);

  /* compute normal */
  ax = p1->c[0] - p0->c[0];
  ay = p1->c[1] - p0->c[1];
  az = p1->c[2] - p0->c[2];
  bx = p2->c[0] - p0->c[0];
  by = p2->c[1] - p0->c[1];
  bz = p2->c[2] - p0->c[2];
  n[0] = ay*bz - az*by;
  n[1] = az*bx - ax*bz;
  n[2] = ax*by - ay*bx;
  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd > 0.0f ) {
    dd = 1.0f / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  }
  glNormal3fv(n);
  glVertex3f(shrink*(p0->c[0]-cx)+cx,
             shrink*(p0->c[1]-cy)+cy,
             shrink*(p0->c[2]-cz)+cz);
  glNormal3fv(n);
  glVertex3f(shrink*(p1->c[0]-cx)+cx,
             shrink*(p1->c[1]-cy)+cy,
             shrink*(p1->c[2]-cz)+cz);
  glNormal3fv(n);
  glVertex3f(shrink*(p2->c[0]-cx)+cx,
             shrink*(p2->c[1]-cy)+cy,
             shrink*(p2->c[2]-cz)+cz);
  glNormal3fv(n);
  glVertex3f(shrink*(p3->c[0]-cx)+cx,
             shrink*(p3->c[1]-cy)+cy,
             shrink*(p3->c[2]-cz)+cz);
  glEnd();

  /* display vertex number */
  glColor3f(1.0-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
  output3(p0->c[0],p0->c[1],p0->c[2],"%d",pq->v[0]);
  output3(p1->c[0],p1->c[1],p1->c[2],"%d",pq->v[1]);
  output3(p2->c[0],p2->c[1],p2->c[2],"%d",pq->v[2]);
  output3(p3->c[0],p3->c[1],p3->c[2],"%d",pq->v[3]);
}

static void drawTets(pScene sc,pMesh mesh,int k) {
  pMaterial  pm;
  pTetra     pt;
  pPoint     p0,p1,p2,p3;
  float      ax,ay,az,bx,by,bz,d,n[3];
  float      shrink,cx,cy,cz;
  int        l;

  /* default */
  if ( ddebug ) printf("draw tetra %d\n",k);
  if ( k < 1 || k > mesh->ntet )  return;
  pt = &mesh->tetra[k];
  if ( refpick > 0 )  pt->ref = refpick;
  refmat = matRef(sc,pt->ref);

  pm = &sc->material[refmat];
  shrink = 0.95*sc->shrink;

  glBegin(GL_TRIANGLES);
  glColor3f(1.0-pm->dif[0],1.0-pm->dif[1],1.0-pm->dif[2]);
  for (l=0; l<4; l++) {
    p0 = &mesh->point[pt->v[ct[l][0]]];
    p1 = &mesh->point[pt->v[ct[l][1]]];
    p2 = &mesh->point[pt->v[ct[l][2]]];
    cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.;
    cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.;
    cz = (p0->c[2] + p1->c[2] + p2->c[2]) / 3.;

    /* compute face normal */
    ax = p1->c[0] - p0->c[0]; ay = p1->c[1] - p0->c[1]; az = p1->c[2] - p0->c[2];
    bx = p2->c[0] - p0->c[0]; by = p2->c[1] - p0->c[1]; bz = p2->c[2] - p0->c[2];
    n[0] = ay*bz - az*by;
    n[1] = az*bx - ax*bz;
    n[2] = ax*by - ay*bx;
    d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if ( d > 0.0f ) {
      d = 1.0 / sqrt(d);
      n[0] *= d;
      n[1] *= d;
      n[2] *= d;
    }
    glNormal3fv(n);
    glVertex3f(shrink*(p0->c[0]-cx)+cx,
               shrink*(p0->c[1]-cy)+cy,
               shrink*(p0->c[2]-cz)+cz);
    glNormal3fv(n);
    glVertex3f(shrink*(p1->c[0]-cx)+cx,
               shrink*(p1->c[1]-cy)+cy,
               shrink*(p1->c[2]-cz)+cz);
    glNormal3fv(n);
    glVertex3f(shrink*(p2->c[0]-cx)+cx,
               shrink*(p2->c[1]-cy)+cy,
               shrink*(p2->c[2]-cz)+cz);
  }
  glEnd();

  /* display vertex number */
  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  p3 = &mesh->point[pt->v[3]];

  glColor3f(1.0-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
  output3(p0->c[0],p0->c[1],p0->c[2],"%d",pt->v[0]);
  output3(p1->c[0],p1->c[1],p1->c[2],"%d",pt->v[1]);
  output3(p2->c[0],p2->c[1],p2->c[2],"%d",pt->v[2]);
  output3(p3->c[0],p3->c[1],p3->c[2],"%d",pt->v[3]);

  /*if ( mesh->nfield == 6 )  drawEllipse(sc,mesh,LTets,k);*/
  if ( !mesh->nbb )
    circumSphere(sc,mesh,LTets,k);
}

static void drawHexa(pScene sc,pMesh mesh,int k) {
  pMaterial  pm;
  pHexa      ph;
  pPoint     p0,p1,p2,p3;
  float      ax,ay,az,bx,by,bz,d,n[3];
  float      shrink,cx,cy,cz;
  int        l;

  /* default */
  if ( ddebug ) printf("draw hexa %d\n",k);
  if ( k < 1 || k > mesh->nhex )  return;
  ph = &mesh->hexa[k];
  if ( refpick > 0 )  ph->ref = refpick;
  refmat = matRef(sc,ph->ref);
  pm = &sc->material[refmat];
  shrink = 0.95*sc->shrink;
  
  glBegin(GL_QUADS);
  glColor3f(1.0-pm->dif[0],1.0-pm->dif[1],1.0-pm->dif[2]);
  for (l=0; l<6; l++) {
    p0 = &mesh->point[ph->v[ch[l][0]]];
    p1 = &mesh->point[ph->v[ch[l][1]]];
    p2 = &mesh->point[ph->v[ch[l][2]]];
    p3 = &mesh->point[ph->v[ch[l][3]]];
    cx = (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]) / 4.;
    cy = (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]) / 4.;
    cz = (p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]) / 4.;

    /* compute face normal */
    ax = p1->c[0] - p0->c[0]; ay = p1->c[1] - p0->c[1]; az = p1->c[2] - p0->c[2];
    bx = p2->c[0] - p0->c[0]; by = p2->c[1] - p0->c[1]; bz = p2->c[2] - p0->c[2];
    n[0] = ay*bz - az*by;
    n[1] = az*bx - ax*bz;
    n[2] = ax*by - ay*bx;
    d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if ( d > 0.0f ) {
      d = 1.0 / sqrt(d);
      n[0] *= d;  
      n[1] *= d;  
      n[2] *= d;
    }
    glNormal3fv(n);
    glVertex3f(sc->shrink*(p0->c[0]-cx)+cx,
               sc->shrink*(p0->c[1]-cy)+cy,
               sc->shrink*(p0->c[2]-cz)+cz);
    glNormal3fv(n);
    glVertex3f(sc->shrink*(p1->c[0]-cx)+cx,
               sc->shrink*(p1->c[1]-cy)+cy,
               sc->shrink*(p1->c[2]-cz)+cz);
    glNormal3fv(n);
    glVertex3f(sc->shrink*(p2->c[0]-cx)+cx,
               sc->shrink*(p2->c[1]-cy)+cy,
               sc->shrink*(p2->c[2]-cz)+cz);
    glNormal3fv(n);
    glVertex3f(sc->shrink*(p3->c[0]-cx)+cx,
               sc->shrink*(p3->c[1]-cy)+cy,
               sc->shrink*(p3->c[2]-cz)+cz);
  }
  glEnd();

  /* display vertex number */
  glColor3f(1.0-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
  for (l=0; l<8; l++) {
    p0 = &mesh->point[ph->v[l]];
    output3(p0->c[0],p0->c[1],p0->c[2],"%d",ph->v[l]);
  }
}

static void drawPoint(pScene sc,pMesh mesh,int k) {
  pPoint      pt;

  pt = &mesh->point[k];

  /*glDisable(GL_DEPTH_TEST);*/
  glDisable(GL_LIGHTING);
  glPointSize(6.0);
  glColor3f(1.0,0.,0.);
  glBegin(GL_POINTS);
    glVertex3f(pt->c[0],pt->c[1],pt->c[2]);
  glEnd();
  output3(pt->c[0],pt->c[1],pt->c[2],"%d",refitem);
  glEnable(GL_LIGHTING);
  /*glEnable(GL_DEPTH_TEST);*/
}


static void infoData(pScene sc,pMesh mesh,int k,int typel) {
  pSolution  ps;

  if ( !mesh->nbb )  return;
  ps = &mesh->sol[k];
  if ( mesh->nfield == 1 )
    fprintf(stdout,"  Data (scalar): %f\n",ps->bb);
  else if ( mesh->nfield == mesh->dim ) {
    fprintf(stdout,"  Data (vector): %f %f",ps->m[0],ps->m[1]);
    if (mesh->dim == 3 )  fprintf(stdout," %f",ps->m[2]); 
    fprintf(stdout,"\n");
  }
  else if ( mesh->dim == 2 && mesh->nfield == 3 ) {
    fprintf(stdout,"  Data (tensor): %f %f %f\n",
            ps->m[0],ps->m[1],ps->m[2]);
    drawEllipse(sc,mesh,typel,k);
  }
  else if ( mesh->dim == 3 && mesh->nfield == 6 ) {
    if ( mesh->ne)
      fprintf(stdout,"  Data (tensor): %f %f %f %f %f %f\n",
              ps->m[0],ps->m[1],ps->m[2],ps->m[3],ps->m[4],ps->m[5]);
    drawEllipsoid(sc,mesh,typel,k);
  }
  fflush(stdout); /* add J. Morice 12/2008 */
}


static void infoEntity(pScene sc,pMesh mesh,int k,int type) {
  pMaterial  pm;
  pTriangle  pt;
  pTetra     ptt;
  pHexa      ph;
  pQuad      pq;
  pPoint     p0;
  int        i;

  if ( mesh->ne)  fprintf(stdout,"\n Picking result :\n");
  pm = &sc->material[refmat];
  switch(type) {
  case LPoint:
    p0 = &mesh->point[k];
    if ( mesh->ne)
      fprintf(stdout,"  Vertex %5d : %f, %f, %f    ref : %d\n",
	      k,p0->c[0]+mesh->xtra,p0->c[1]+mesh->ytra,p0->c[2]+mesh->ztra,p0->ref);
    if ( mesh->nbb && mesh->typage == 2 )  infoData(sc,mesh,k,LPoint);
    break;
 
  case LTria:
    pt = &mesh->tria[k];
    fprintf(stdout,"  Triangle %5d : %d, %d, %d    ref : %d [%s]\n",
	    k,pt->v[0],pt->v[1],pt->v[2],pt->ref,pm->name);
    if ( mesh->nbb && mesh->typage == 1 )  infoData(sc,mesh,k,LTria);
    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      fprintf(stdout,"  vertex   %5d : %f %f %f   ref %d\n",
	    pt->v[i],p0->c[0]+mesh->xtra,p0->c[1]+mesh->ytra,p0->c[2]+mesh->ztra,p0->ref);
      if ( mesh->nbb && mesh->typage == 2 )  infoData(sc,mesh,pt->v[i],LPoint);
    }
    break;
  
  case LQuad:
    pq = &mesh->quad[k];
    fprintf(stdout,"  Quad  %5d : %d, %d, %d, %d    ref : %d [%s]\n",
	    k,pq->v[0],pq->v[1],pq->v[2],pq->v[3],pq->ref,pm->name);
    if ( mesh->nbb && mesh->typage == 1 )  infoData(sc,mesh,k,LQuad);
    for (i=0; i<4; i++) {
      p0 = &mesh->point[pq->v[i]];
      fprintf(stdout,"  vertex   %5d : %f %f %f   ref %d\n",
	    pq->v[i],p0->c[0]+mesh->xtra,p0->c[1]+mesh->ytra,p0->c[2]+mesh->ztra,p0->ref);
      if ( mesh->nbb && mesh->typage == 2 )  infoData(sc,mesh,pq->v[i],LPoint);
    }
    break;
  
  case LTets:
    ptt = &mesh->tetra[k];
    fprintf(stdout,"  Tetra  %5d : %d, %d, %d, %d    ref : %d [%s]\n",
	    k,ptt->v[0],ptt->v[1],ptt->v[2],ptt->v[3],ptt->ref,pm->name);
    if ( mesh->nbb && mesh->typage == 1 )  infoData(sc,mesh,k,LTets);
    for (i=0; i<4; i++) {
      p0 = &mesh->point[ptt->v[i]];
      fprintf(stdout,"  vertex   %5d : %f %f %f   ref %d\n",
	      ptt->v[i],p0->c[0]+mesh->xtra,p0->c[1]+mesh->ytra,p0->c[2]+mesh->ztra,p0->ref);
      if ( mesh->nbb && mesh->typage == 2 )  infoData(sc,mesh,ptt->v[i],LPoint);
    }
    break;  
  
  case LHexa:
    ph = &mesh->hexa[k];
    fprintf(stdout,"  Hexa   %5d : %d, %d, %d, %d, %d, %d, %d, %d    ref : %d [%s]\n",
	    k,ph->v[0],ph->v[1],ph->v[2],ph->v[3],ph->v[4],ph->v[5],
            ph->v[6],ph->v[7],ph->ref,pm->name);
    if ( mesh->nbb && mesh->typage == 1 )  infoData(sc,mesh,k,LHexa);
    for (i=0; i<8; i++) {
      p0 = &mesh->point[ph->v[i]];
      fprintf(stdout,"  vertex   %5d : %f %f %f   ref %d\n",
	      ph->v[i],p0->c[0]+mesh->xtra,p0->c[1]+mesh->ytra,p0->c[2]+mesh->ztra,p0->ref);
      if ( mesh->nbb && mesh->typage == 2 )  infoData(sc,mesh,ph->v[i],LPoint);
    }
    break;  
  }
  fflush(stdout); /* add J. Morice 12/2008 */
}


static int getColorRange(Color *c,pMesh mesh) {
  GLubyte    mask;
  GLint      rBits,gBits,bBits,aBits;
  long       nbmax;
  int        i,nbits;

  glGetIntegerv(GL_RED_BITS, &rBits);
  glGetIntegerv(GL_GREEN_BITS, &gBits);
  glGetIntegerv(GL_BLUE_BITS, &bBits);
  glGetIntegerv(GL_ALPHA_BITS, &aBits);
  nbits = rBits+gBits+bBits+aBits;

  if ( nbits < 32 ) {
    nbmax = 2 << (nbits-1);
    if ( nbmax < mesh->nt+mesh->nq+mesh->ntet+mesh->nhex ) {
      fprintf(stderr,"  Sorry. Picking disabled. (%ld,%d)\n",nbmax,mesh->nt+mesh->nq);
      return(0);
    }
    else if ( nbmax < 0.1*(mesh->ntet+mesh->nhex) ) {
      fprintf(stderr,"  Sorry. Picking disabled. (%ld,%d)\n",nbmax,mesh->nt+mesh->nq);
      return(0);
    }
  }

  mask  = 0;
  nbits = nbits - rBits;
  for (i=0; i<rBits; i++)  mask |= 1 << i;
  c->rMask   = mask << nbits;
  c->rShift  = nbits;
  c->rBits   = 8 - rBits;

  mask  = 0;
  nbits = nbits - gBits;
  for (i=0; i<gBits; i++)  mask |= 1 << i;
  c->gMask  = mask << nbits;
  c->gShift = nbits;
  c->gBits  = 8 - gBits;

  mask  = 0;
  nbits = nbits - bBits;
  for (i=0; i<bBits; i++)  mask |= 1 << i;
  c->bMask   = mask << aBits;
  c->bShift  = nbits;
  c->bBits   = 8 - bBits;

  mask = 0;
  for (i=0; i<aBits; i++)  mask |= 1 << i;
  c->aMask  = mask;
  c->aBits  = 8 - aBits;
}


static void displayPoint(pScene sc,pMesh mesh,Color *c) {
  pPoint       ppt;
  unsigned int k,kk;

  glPointSize(5);
  glBegin(GL_POINTS);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];

    kk = 2*k+1;
    glColor4ub((kk & c->rMask) >> c->rShift << c->rBits,
               (kk & c->gMask) >> c->gShift << c->gBits,
               (kk & c->bMask) >> c->bShift << c->bBits,
               (kk & c->aMask)              << c->aBits);

    if ( mesh->dim == 2 )
      glVertex2dv(ppt->c);
    else
      glVertex3dv(ppt->c);
  }
  glEnd();
  glPointSize(1);
}


static void displayTria(pScene sc,pMesh mesh,Color *c) {
  pTriangle    pt;
  pPoint       p0,p1,p2;
  pMaterial    pm;
  int          m;
  unsigned int k,kk;

  glBegin(GL_TRIANGLES);
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTria];
    if ( !k || pm->flag )  continue;

    while ( k != 0 ) {
      pt = &mesh->tria[k];
      if ( pt->v[0] ) {
        p0 = &mesh->point[pt->v[0]];
        p1 = &mesh->point[pt->v[1]];
        p2 = &mesh->point[pt->v[2]];

        kk = 2*k+1;
	    glColor4ub((kk & c->rMask) >> c->rShift << c->rBits,
 	               (kk & c->gMask) >> c->gShift << c->gBits,
	               (kk & c->bMask) >> c->bShift << c->bBits,
		           (kk & c->aMask)              << c->aBits);

	    glVertex3dv(p0->c);
        glVertex3dv(p1->c);
        glVertex3dv(p2->c);
      }
      k = pt->nxt;
    }
  }
  glEnd();
}

static int displayQuad(pScene sc,pMesh mesh,Color *c) {
  pQuad        pq;
  pPoint       p0,p1,p2,p3;
  pMaterial    pm;
  int          k,m,base;
  unsigned int kk;

  base = mesh->nt;

  glBegin(GL_QUADS);
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LQuad];
    if ( !k || pm->flag )  continue;

    while ( k != 0 ) {
      pq = &mesh->quad[k];
      if ( pq->v[0] ) {
        p0 = &mesh->point[pq->v[0]];
        p1 = &mesh->point[pq->v[1]];
        p2 = &mesh->point[pq->v[2]];
        p3 = &mesh->point[pq->v[3]];

        kk = base + k;
        kk = 2*kk + 1;
	glColor4ub((kk & c->rMask) >> c->rShift << c->rBits,
 	           (kk & c->gMask) >> c->gShift << c->gBits,
	           (kk & c->bMask) >> c->bShift << c->bBits,
		   (kk & c->aMask)              << c->aBits);

	glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
        glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
        glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
        glVertex3f(p3->c[0],p3->c[1],p3->c[2]);
      }
      k = pq->nxt;
    }
  }
  glEnd();
}

static int displayTets(pScene sc,pMesh mesh,Color *c) {
  pTetra       ptt;
  pPoint       p0,p1,p2;
  pMaterial    pm;
  pClip        clip;
  int          k,l,m,base;
  unsigned int kk;

  clip = sc->clip;
  base = mesh->nt + mesh->nq;

  glBegin(GL_TRIANGLES);
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;

    while ( k != 0 ) {
      ptt = &mesh->tetra[k];
      if ( ptt->v[0] ) {
        if ( (clip->active & C_ON && ptt->clip) || !(clip->active & C_ON) ) {
          kk = base + k;
          kk = 2*kk + 1;
  	  for (l=0; l<4; l++) {
	    p0 = &mesh->point[ptt->v[ct[l][0]]];
	    p1 = &mesh->point[ptt->v[ct[l][1]]];
	    p2 = &mesh->point[ptt->v[ct[l][2]]];

	    glColor4ub((kk & c->rMask) >> c->rShift << c->rBits,
 	               (kk & c->gMask) >> c->gShift << c->gBits,
	               (kk & c->bMask) >> c->bShift << c->bBits,
		       (kk & c->aMask)              << c->aBits);

	    glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
            glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
            glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
	  }
	}
      }
      k = ptt->nxt;
    }
  }
  glEnd();
}

static int displayHexa(pScene sc,pMesh mesh,Color *c) {
  pHexa        ph;
  pPoint       p0,p1,p2,p3;
  pMaterial    pm;
  pClip        clip;
  int          k,l,m,base;
  unsigned int kk;

  clip = sc->clip;
  base = mesh->nt + mesh->nq + mesh->ntet;

  glBegin(GL_QUADS);
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LHexa];
    if ( !k || pm->flag )  continue;

    while ( k != 0 ) {
      ph = &mesh->hexa[k];
      if ( ph->v[0] ) {
        if ( (clip->active & C_ON && ph->clip) || !(clip->active & C_ON) ) {
          kk = base + k;
          kk = 2*kk + 1;
  	  for (l=0; l<6; l++) {
	    p0 = &mesh->point[ph->v[ch[l][0]]];
	    p1 = &mesh->point[ph->v[ch[l][1]]];
	    p2 = &mesh->point[ph->v[ch[l][2]]];
	    p3 = &mesh->point[ph->v[ch[l][3]]];

	    glColor4ub((kk & c->rMask) >> c->rShift << c->rBits,
 	               (kk & c->gMask) >> c->gShift << c->gBits,
	               (kk & c->bMask) >> c->bShift << c->bBits,
		       (kk & c->aMask)              << c->aBits);

	    glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
            glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
            glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
            glVertex3f(p3->c[0],p3->c[1],p3->c[2]);
	  }
	}
      }
      k = ph->nxt;
    }
  }
  glEnd();
}


void drawModelSimple(pScene sc,pMesh mesh,Color *c) {

  glDisable(GL_BLEND);
  glDisable(GL_LIGHTING);
  glShadeModel(GL_FLAT);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  if ( sc->clip->active & C_ON ) {
    glClipPlane(GL_CLIP_PLANE0,sc->clip->eqn);
    glEnable(GL_CLIP_PLANE0);
    if ( mesh->nt )    displayTria(sc,mesh,c);
    if ( mesh->nq )    displayQuad(sc,mesh,c);
    glDisable(GL_CLIP_PLANE0);
    if ( sc->clip->active & C_VOL ) {
      if ( mesh->ntet )  displayTets(sc,mesh,c);
      if ( mesh->nhex )  displayHexa(sc,mesh,c);
    }
  }
  else {
    if ( mesh->nt )  displayTria(sc,mesh,c);
    if ( mesh->nq )  displayQuad(sc,mesh,c);
    if ( mesh->dim == 3 ) {
      if ( !mesh->nt ) displayTets(sc,mesh,c);
      if ( !mesh->nq ) displayHexa(sc,mesh,c);
    }
  }

  if ( !mesh->ne )  displayPoint(sc,mesh,c);
    
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
}


static int closestPoint(pScene sc,pMesh mesh,int x,int y,int item,int type) {
  pTriangle   pt;
  pQuad       pq;
  pTetra      pte;
  pHexa       ph;
  pPoint      ppt;
  GLdouble    projmat[16],winx,winy,winz,matrix[16];
  GLint       viewport[4];
  double      dd,dmin;
  int         i,id;

  glGetDoublev(GL_PROJECTION_MATRIX,projmat);
  glGetIntegerv(GL_VIEWPORT,viewport);
  glGetDoublev(GL_MODELVIEW_MATRIX,matrix);
  y = viewport[3] - y;

  id   = -1;
  dmin = 1.e20;
  switch(type) {
  case LTria:
    pt   = &mesh->tria[item];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      gluProject(ppt->c[0],ppt->c[1],ppt->c[2],
                 matrix,projmat,viewport,
                 &winx,&winy,&winz);
      dd = (winx - x)*(winx - x) + (winy - y)*(winy - y);
      if ( dd < dmin ) {
        dmin = dd;
	id   = i;
      }
    }
    return(pt->v[id]);
  
  case LQuad:
    pq   = &mesh->quad[item];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pq->v[i]];
      gluProject(ppt->c[0]-sc->cx,ppt->c[1]-sc->cy,ppt->c[2]-sc->cz,
                 matrix,projmat,viewport,
                 &winx,&winy,&winz);
      dd = (winx - x)*(winx - x) + (winy - y)*(winy - y);
      if ( dd < dmin ) {
        dmin = dd;
	id   = i;
      }
    }
    return(pq->v[id]);
  
  case LTets:
    pte  = &mesh->tetra[item];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pte->v[i]];
      gluProject(ppt->c[0]-sc->cx,ppt->c[1]-sc->cy,ppt->c[2]-sc->cz,
                 matrix,projmat,viewport,
                 &winx,&winy,&winz);
      dd = (winx - x)*(winx - x) + (winy - y)*(winy - y);
      if ( dd < dmin ) {
        dmin = dd;
	id   = i;
      }
    }
    return(pte->v[id]);

  case LHexa:
    ph   = &mesh->hexa[item];
    for (i=0; i<8; i++) {
      ppt = &mesh->point[ph->v[i]];
      gluProject(ppt->c[0]-sc->cx,ppt->c[1]-sc->cy,ppt->c[2]-sc->cz,
                 matrix,projmat,viewport,
                 &winx,&winy,&winz);
      dd = (winx - x)*(winx - x) + (winy - y)*(winy - y);
      if ( dd < dmin ) {
        dmin = dd;
	id   = i;
      }
    }
    return(ph->v[id]);
  }

  return(0);
}


GLuint pickingScene(pScene sc,int x,int y,int ident) {
  pMesh         mesh;
  pClip         clip;
  GLint         viewport[4];
  GLubyte       pixel[4];
  GLuint        dlist;
  Color         c;
  int           k;
  unsigned int  item;

  dlist   = 0;  
  refitem = 0;
  refmat  = -1;
  mesh    = cv.mesh[sc->idmesh];
  clip    = sc->clip;

  if ( !getColorRange(&c,mesh) )  return(dlist);

  glGetIntegerv(GL_VIEWPORT,viewport);
  glDrawBuffer(GL_BACK_LEFT);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.,0.,-sc->persp->depth, 0.,0.,0., 0.0,1.0,0.0);

  setupView(sc);
  glMultMatrixf(sc->view->matrix);
  glTranslatef(sc->cx,sc->cy,sc->cz);
  drawModelSimple(sc,mesh,&c);
  glFlush();

  /* process color */
  glReadBuffer(GL_BACK_LEFT);
  glReadPixels(x,viewport[3]-y,1,1,
               GL_RGBA,GL_UNSIGNED_BYTE,(void *)pixel);
  glDrawBuffer(GL_BACK_LEFT);

  if ( ddebug )  printf("pixel %d %d %d %d\n",pixel[0],pixel[1],pixel[2],pixel[3]);

  item  = pixel[0] >> c.rBits << c.rShift & c.rMask;
  item += pixel[1] >> c.gBits << c.gShift & c.gMask;
  item += pixel[2] >> c.bBits << c.bShift & c.bMask;
  item += pixel[3] >> c.aBits             & c.aMask;
  item /= 2;
  if ( !item )  return(dlist);

  if ( ddebug )  printf("item %d\n",item);

  if ( !mesh->ne ) {
    if ( item <= mesh->np ) {
      refitem = item;
      reftype = LPoint;
    }
  }
  else if ( item <= mesh->nt ) {
    refitem = item;
    reftype = LTria;
  }
  else if ( item <= mesh->nt+mesh->nq ) {
    refitem = item - mesh->nt;
    reftype = LQuad;
  }
  else if ( item <= mesh->nt+mesh->nq+mesh->ntet ) {
    refitem = item - (mesh->nt+mesh->nq);
    reftype = LTets;
  }
  else if ( item <= mesh->nt+mesh->nq+mesh->ntet+mesh->nhex ) {
    refitem = item - (mesh->nt+mesh->nq+mesh->ntet);
    reftype = LHexa;
  }

  /* peculiar case: vertex */
  if ( refitem > 0 && ident == LPoint ) {
    if ( mesh->ne ) {
      refitem = closestPoint(sc,mesh,x,y,refitem,reftype);
      reftype = LPoint;
    }
  }

  if ( refitem > 0 ) {
    dlist = glGenLists(1);
    glNewList(dlist,GL_COMPILE);
    if ( glGetError() )  return(0);

    switch(reftype) {
    case LPoint:
      if ( refitem <= mesh->np && !mesh->ne ) {
        for (k=1; k<=mesh->np; k++) {
          /*drawPoint(sc,mesh,k);*/
          infoEntity(sc,mesh,k,LPoint);
        }
      }
      else {
        drawPoint(sc,mesh,refitem);
        infoEntity(sc,mesh,refitem,LPoint);
      }
      break;
    case LTria:
      drawTria(sc,mesh,refitem);
      infoEntity(sc,mesh,refitem,LTria);
      break;
    case LQuad:
      drawQuad(sc,mesh,refitem);
      infoEntity(sc,mesh,refitem,LQuad);
      break;
    case LTets:
      drawTets(sc,mesh,refitem);
      infoEntity(sc,mesh,refitem,LTets);
      break;
    case LHexa:
      drawHexa(sc,mesh,refitem);
      infoEntity(sc,mesh,refitem,LHexa);
      break;
    }

    glEndList();
    return(dlist);
  }
  else
    refmat = -1;

  return(dlist); 
}


GLuint pickItem(pMesh mesh,pScene sc,int numit) {
  pMaterial  pm;
  pPoint     ppt;
  GLuint     dlist;
  int        nm;

  /* default */
  if ( ddebug ) printf("create pickitem list\n");
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  if ( numit <= mesh->np && mesh->ne ) {
    ppt = &mesh->point[numit];
    pm  = &sc->material[DEFAULT_MAT];
    if ( ppt->tag != M_UNUSED ) {
      if ( ppt->ref ) {
        nm = matRef(sc,ppt->ref);
        /*nm = 1+(ppt->ref-1) % (sc->par.nbmat-1);*/
        pm = &sc->material[nm];
      }
      glPointSize(3.0);
      glColor3f(pm->dif[0],pm->dif[1],pm->dif[2]);
      output3(ppt->c[0],ppt->c[1],ppt->c[2],"%d",numit);
    }
  }

  if ( numit <= mesh->nt )   drawTria(sc,mesh,numit);
  if ( numit <= mesh->nq )   drawQuad(sc,mesh,numit);
  if ( numit <= mesh->ntet ) drawTets(sc,mesh,numit);
  if ( numit <= mesh->nhex ) drawHexa(sc,mesh,numit);

  glEndList(); 
  return(dlist);
}
