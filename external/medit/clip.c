#include "medit.h"
#include "extern.h"
#include "sproto.h"

static pClip     cclip   = 0;
static ubyte     curclip = 0;
static GLfloat   plane[4]  = {-1.0, 0.0, 0.0, 0.0};


static void drawCap(pScene sc,pClip clip,GLboolean docap) {

  if ( !docap ) {  
    if ( clip->active & C_EDIT )        glColor3f(1.0,0.0,1.0);
    else if ( clip->active & C_FREEZE ) glColor3f(0.0,0.6,0.9);
    else                                glColor3f(0.0,1.0,0.0);

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(3.);
    glBegin(GL_QUADS);
      glVertex3f(0.,-1.,-1.);
      glVertex3f(0., 1.,-1.);
      glVertex3f(0., 1.,1.);
      glVertex3f(0.,-1.,1.);
    glEnd();

    glLineWidth(1.);
    glColor3f(1.,0.7,0.);
    glBegin(GL_LINES);
      glVertex3f(0.,0.,-1.);
      glVertex3f(0.,0.,1.);
      glColor3f(1.,0.7,0.);
      glVertex3f(0.,-1.,0.);
      glVertex3f(0.,1.,0.);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  }
  else {
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);    
    glColor4fv(sc->material[DEFAULT_MAT].dif);
    glRecti(-100,-100,100,100);
  }
}


/* display clipping plane */
void updateClip(pClip clip,pMesh mesh) {
  pScene      sc;
  pTransform  cliptr,view;
  GLfloat     dd,dmax,inv[16],axis[4],trans[4],matrix[16];
  int         idw;

  /* default */
  if ( ddebug ) printf("updateClip\n");

  /* retrieve context */
  idw    = currentScene();
  sc     = cv.scene[idw];
  view   = sc->view;
  cliptr = clip->cliptr;
#ifdef IGL
  if(sc->igl_params->render_on_C_UPDATE)
  {
    sc->igl_params->render_on_next = true;
  }
#endif

  /* compute clip transform */
  if ( clip->active & C_EDIT ) {
    invertMatrix(view->matrix,inv);
    inv[3] = inv[7] = inv[11] = 0.0;  inv[15] = 1.0;
    
    /* rotation cumulative */
    if ( cliptr->angle != 0.0 ) {
      transformVector(trans,cliptr->axis,inv);
      glPushMatrix();
      glLoadIdentity();
      glRotatef(cliptr->angle,trans[0],trans[1],trans[2]);
      glMultMatrixf(cliptr->rot);
      glGetFloatv(GL_MODELVIEW_MATRIX,cliptr->rot);
      glPopMatrix();
    }

    /* translation cumulative */
    axis[0] = cliptr->panx;
    axis[1] = cliptr->pany;
    axis[2] = 0.0;
    axis[3] = 1.0;
    if ( cliptr->manim == GL_FALSE )
      cliptr->panx = cliptr->pany = 0.0;
    transformVector(trans,axis,inv);

    dd = trans[0]*clip->eqn[0]+trans[1]*clip->eqn[1]+trans[2]*clip->eqn[2];
    trans[0] = dd*clip->eqn[0];
    trans[1] = dd*clip->eqn[1];
    trans[2] = dd*clip->eqn[2];

    cliptr->tra[12] += trans[0];
    cliptr->tra[13] += trans[1];
    cliptr->tra[14] += trans[2];

    /* truncation */
    dmax = mesh->xmax - mesh->xmin;
    dmax = MEDIT_MAX(dmax,mesh->ymax - mesh->ymin);
    dmax = MEDIT_MAX(dmax,mesh->zmax - mesh->zmin) / 1.8;
    if ( fabs(cliptr->tra[12]) > dmax || fabs(cliptr->tra[13]) > dmax ||
         fabs(cliptr->tra[14]) > dmax ) {
      if ( cliptr->manim == GL_TRUE ) {
        cliptr->panx = -cliptr->panx;
        cliptr->pany = -cliptr->pany;
      }
      else {
        cliptr->tra[12] = MEDIT_MAX(-dmax,MEDIT_MIN(dmax,cliptr->tra[12]));
        cliptr->tra[13] = MEDIT_MAX(-dmax,MEDIT_MIN(dmax,cliptr->tra[13]));
        cliptr->tra[14] = MEDIT_MAX(-dmax,MEDIT_MIN(dmax,cliptr->tra[14]));
      }
    }

    /* final transformation */
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixf(cliptr->tra);
    glMultMatrixf(cliptr->rot);
    glGetFloatv(GL_MODELVIEW_MATRIX,cliptr->matrix);
    glPopMatrix();

    /* compute plane equation */
    invertMatrix(cliptr->matrix,inv);
    transformPoint(clip->eqn,plane,inv);

    if ( clip->active & C_REDO )
      clipVertices(mesh,sc,clip);
   /* clip->active |= C_REDO;*/
  }

  else if ( clip->active & C_FREEZE ) {
    glPushMatrix();
    glLoadMatrixf(view->matrix);
    glTranslatef(sc->cx,sc->cy,sc->cz);
    glGetFloatv(GL_MODELVIEW_MATRIX,matrix);
    glPopMatrix();
    invertMatrix(matrix,inv);
    inv[3] = inv[7] = inv[11] = 0.0;  inv[15] = 1.0;

    glPushMatrix();
    glLoadMatrixf(inv);
    glMultMatrixf(view->oldmat);
    glTranslatef(sc->cx,sc->cy,sc->cz);
    glMultMatrixf(cliptr->matrix);
    glGetFloatv(GL_MODELVIEW_MATRIX,cliptr->matrix);
    glPopMatrix();

    /* compute plane equation */
    invertMatrix(cliptr->matrix,inv);
    transformPoint(clip->eqn,plane,inv);
    clip->active |= C_REDO;
  }
  if ( !cliptr->manim ) {
    cliptr->angle = 0.0;
    clip->active ^= C_UPDATE;
  }
}


void clipVertices(pMesh mesh,pScene sc,pClip clip) {
  pTetra     pt;
  pHexa      ph;
  pPoint     p0;
  double     dd1,zero;
  int        k,l,nbpos,nbneg,nbnul;

#ifdef IGL
  if(sc->igl_params->hot_dog_view)
  {
    /* check points in plane */
    zero = sc->dmax*1.e-13;
    const double width = sc->igl_params->width(mesh);
    const double hot_dog_ratio = sc->igl_params->hot_dog_ratio;
    for (k=1; k<=mesh->np; k++) {
      for(int h = 0;h<sc->igl_params->num_hot_dog_slices;h++)
      {
        p0 = &mesh->point[k];
        p0->clip = 0;
        if ( p0->tag & M_UNUSED )  continue;
        p0->hd_dd1[h] = p0->c[0]*clip->eqn[0] + p0->c[1]*clip->eqn[1] \
          + p0->c[2]*clip->eqn[2] + clip->eqn[3] + 
          (h%2==0 ?  width*h : width*(h-1) + width*hot_dog_ratio*2);
        if ( p0->hd_dd1[h] > zero )      p0->hd_clip[h] = 2;
        else if ( p0->hd_dd1[h] < zero ) p0->hd_clip[h] = 1;
        else                   p0->hd_clip[h] = 0;
      }
    }

    /* update tetrahedra */
    for (k=1; k<=mesh->ntet; k++) {
      pt = &mesh->tetra[k];
      pt->clip = 0;

      for(int h = 0;h<sc->igl_params->num_hot_dog_slices;h++)
      {
        nbpos = nbneg = nbnul = 0;
        for (l=0; l<4; l++) {
          p0  = &mesh->point[pt->v[l]];
          if ( p0->hd_clip[h]== 2 )       nbpos++;
          else if ( p0->hd_clip[h] == 1 )  nbneg++;
          else                       nbnul++;
        }
        if ( nbpos && nbpos+nbnul < 4 )  pt->clip = 1;
      }
    }
  }else{
#endif

  /* check points in plane */
  zero = sc->dmax*1.e-13;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( p0->tag & M_UNUSED )  continue;
    dd1 = p0->c[0]*clip->eqn[0] + p0->c[1]*clip->eqn[1] \
        + p0->c[2]*clip->eqn[2] + clip->eqn[3];
    if ( dd1 > zero )      p0->clip = 2;
    else if ( dd1 < zero ) p0->clip = 1;
    else                   p0->clip = 0;
  }

  /* update tetrahedra */
  for (k=1; k<=mesh->ntet; k++) {
    pt = &mesh->tetra[k];
    pt->clip = 0;
    nbpos = nbneg = nbnul = 0;
    for (l=0; l<4; l++) {
      p0  = &mesh->point[pt->v[l]];
      if ( p0->clip == 2 )       nbpos++;
      else if ( p0->clip == 1 )  nbneg++;
      else                       nbnul++;
    }
    if ( nbpos && nbpos+nbnul < 4 )  pt->clip = 1;
  }
#ifdef IGL
  }
#endif

  /* update hexahedra */
  for (k=1; k<=mesh->nhex; k++) {
    ph = &mesh->hexa[k];
    ph->clip = 0;
    nbpos = nbneg = nbnul = 0;
    for (l=0; l<8; l++) {
      p0  = &mesh->point[ph->v[l]];
      if ( p0->clip == 2 )       nbpos++;
      else if ( p0->clip == 1 )  nbneg++;
      else                       nbnul++;
    }
    if ( nbpos && nbpos+nbnul < 8 )  ph->clip = 1;
  }
}


void drawClip(pScene sc,pClip clip,pMesh mesh,GLboolean docap) {
  pTransform  cliptr,view;
  GLfloat     scale;

  /* default */
  if ( ddebug ) printf("drawClip\n");
  view   = sc->view;
  cliptr = clip->cliptr;

  if ( clip->active & C_UPDATE )  updateClip(clip,mesh);
  if ( clip->active & C_REDO ) {
    if ( !animate )  glutSetCursor(GLUT_CURSOR_WAIT);
    /* update clip plane */
    clipVertices(mesh,sc,clip);

    /* build display lists */
    if ( clip->active & C_VOL ) {
      if ( sc->clist[LTets] )  glDeleteLists(sc->clist[LTets],1);
      if ( sc->clist[LHexa] )  glDeleteLists(sc->clist[LHexa],1);
      if ( sc->cmlist[LTets] ) glDeleteLists(sc->cmlist[LTets],1);
      if ( sc->cmlist[LHexa] ) glDeleteLists(sc->cmlist[LHexa],1);
      sc->clist[LTets] = sc->cmlist[LTets] = (GLuint)0;
      sc->clist[LHexa] = sc->cmlist[LHexa] = (GLuint)0;

      if ( clip->active & C_CAP ) {
        sc->clist[LTets]  = capTetra(mesh);
        if ( sc->mode & S_MAP )
          sc->cmlist[LTets] = capTetraMap(mesh);
        else if ( sc->isotyp & S_ISOLINE )
          sc->ilist[LTets] = capTetraIso(mesh);
      }
      else {
        sc->clist[LTets] = listTetra(sc,mesh,1);
        sc->clist[LHexa] = listHexa(sc,mesh,1);
        if ( sc->mode & S_MAP ) {
          sc->cmlist[LTets] = listTetraMap(sc,mesh,1);
          sc->cmlist[LHexa] = listHexaMap(sc,mesh,1);
        }
      }
      if ( !animate )  clip->active ^= C_REDO;
    }
    else clip->active ^= C_REDO;
    
    if ( sc->isotyp & S_VECTOR ) {
      if ( sc->vlist[LTets] )  glDeleteLists(sc->vlist[LTets],1);
        sc->vlist[LTets] = listClipTetraVector(mesh);
      if ( sc->vlist[LHexa] )  glDeleteLists(sc->vlist[LHexa],1);
        sc->vlist[LHexa] = listClipHexaVector(mesh);
    }
    if ( !animate )  glutSetCursor(GLUT_CURSOR_INHERIT);
  }

  /* display plane frame */
  if ( clip->active & C_HIDE )  return;
  glPushMatrix();
  glMultMatrixf(cliptr->matrix);
  scale = 0.3*(sc->dmax+sc->dmin);
  glScalef(scale,scale,scale);

  drawCap(sc,clip,docap);

  glPopMatrix();
}


void copyClip(pClip clip) {
  if ( !cclip ) {
    cclip = (pClip)M_calloc(1,sizeof(struct clip),"clip");
    if ( !clip )  exit(2);
  }
  cclip = (pClip)memcpy(cclip,clip,sizeof(struct clip));
  if ( !cclip )  exit(2);
  curclip = 1;
}

int pasteClip(pClip clip) {
  if ( !curclip )  return(0);
  clip = (pClip)memcpy(clip,cclip,sizeof(struct clip));
  if ( !clip )  exit(2);
  clip->active = 0;
  curclip = 0;
  return(1);
}

void tiltClip(pScene sc,pClip clip) { 
  float    axis[4];
  
  axis[0] = clip->cliptr->axis[0];
  axis[1] = clip->cliptr->axis[1];
  axis[2] = clip->cliptr->axis[2];

  transformVector(clip->cliptr->axis,axis,sc->view->matrix);
  /*clip->cliptr->angle = 90.0;*/

  clip->active |= C_REDO;
}

/* change clip orientation */
void invertClip(pScene sc,pClip clip) {
  clip->eqn[0] = -clip->eqn[0];
  clip->eqn[1] = -clip->eqn[1];
  clip->eqn[2] = -clip->eqn[2];
  clip->eqn[3] = -clip->eqn[3];
  plane[0]     = -plane[0];
}

void resetClip(pScene sc,pClip clip,pMesh mesh) {
  double    dd;

  resetTransform(clip->cliptr);
  dd = sc->par.clip[0]*sc->par.clip[3] + sc->par.clip[1]*sc->par.clip[4] \
     + sc->par.clip[2]*sc->par.clip[5];

  clip->active |= C_REDO;
  clip->eqn[0] = sc->par.clip[3];
  clip->eqn[1] = sc->par.clip[4];
  clip->eqn[2] = sc->par.clip[5];
  clip->eqn[3] = -dd;
  /*clip->active = C_ON + C_VOL;*/
}


/* create a clipping plane */
pClip createClip(pScene sc,pMesh mesh) {
  pClip    clip;

  /* default */
  clip = (pClip)M_calloc(1,sizeof(struct clip),"clip");
  assert(clip);
  clip->cliptr = (pTransform)M_calloc(1,sizeof(struct transform),"clip");
  assert(clip->cliptr);

  resetClip(sc,clip,mesh);
  return(clip);
}


