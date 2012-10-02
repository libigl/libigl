#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define SCALV      1.0e-04
#define EPST       1.e-14
#define COSPI15    .9781476007338057
#define SINPI15    .2079116908177593


/* draw a vector in 3D: u = unit vector */
void drawVector3D(float p[3],double u[3],double scale) {
  double   c[3],cc[3],m[9],dd,scal5;

  /* local frame */
  m[0] = u[0];  m[1] = u[1];  m[2] = u[2];

  if ( fabs(u[0]) > EPS ) {
    m[3] = -(u[1]+u[2]) / u[0];
    m[4] = m[5] = 1.0f;
  }
  else if ( fabs(u[1]) > EPS ) {
    m[4] = -(u[0]+u[2]) / u[1];
    m[3] = m[5] = 1.0f;
  }
  else {
    m[5] = -(u[0]+u[1]) / u[2];
    m[3] = m[4] = 1.0;
  }
  dd = 1.0 / sqrt(m[3]*m[3] + m[4]*m[4] + m[5]*m[5]);
  m[3] *= dd;
  m[4] *= dd;
  m[5] *= dd;
  
  m[6] = m[1]*m[5] - m[2]*m[4];
  m[7] = m[2]*m[3] - m[0]*m[5];
  m[8] = m[0]*m[4] - m[3]*m[1];

  u[0] *= scale;
  u[1] *= scale;
  u[2] *= scale;
  c[0] = p[0] + u[0];
  c[1] = p[1] + u[1];
  c[2] = p[2] + u[2];
  glVertex3fv(p);
  glVertex3dv(c);

  /* draw 4 arrows */
  scal5 = scale / 10.0;
  cc[0] = -COSPI15*scal5;
  cc[1] =  0.0;
  cc[2] =  SINPI15*scal5;

  /* M^-1 . X */
  glVertex3dv(c);
  glVertex3d(c[0]+(m[0]*cc[0]+m[6]*cc[2]),
             c[1]+(m[1]*cc[0]+m[7]*cc[2]),
             c[2]+(m[2]*cc[0]+m[8]*cc[2]));

  cc[0] = -COSPI15*scal5;
  cc[2] = -cc[2];
  glVertex3dv(c);
  glVertex3d(c[0]+(m[0]*cc[0]+m[6]*cc[2]),
             c[1]+(m[1]*cc[0]+m[7]*cc[2]),
             c[2]+(m[2]*cc[0]+m[8]*cc[2]));

  cc[0] = -COSPI15*scal5;
  cc[1] =  SINPI15*scal5;
  cc[2] =  0.0;
  glVertex3dv(c);
  glVertex3d(c[0]+(m[0]*cc[0]+m[3]*cc[1]),
             c[1]+(m[1]*cc[0]+m[4]*cc[1]),
             c[2]+(m[2]*cc[0]+m[5]*cc[1]));

  cc[1] = -cc[1];
  glVertex3dv(c);
  glVertex3d(c[0]+(m[0]*cc[0]+m[3]*cc[1]),
             c[1]+(m[1]*cc[0]+m[4]*cc[1]),
             c[2]+(m[2]*cc[0]+m[5]*cc[1]));
}

void drawVector2D(float p[2],double u[2],double scale) {
  double   c[2],dx,dy;

  u[0] *= scale;
  u[1] *= scale;
  c[0] = p[0] + u[0];
  c[1] = p[1] + u[1];
  glVertex2fv(p);
  glVertex2dv(c);

  dx = ( COSPI15*u[0] + SINPI15*u[1]) / 3.0;
  dy = (-SINPI15*u[0] + COSPI15*u[1]) / 3.0;
  glVertex2dv(c);
  glVertex2d(c[0]-dx,c[1]-dy);

  dx = ( COSPI15*u[0] - SINPI15*u[1]) / 3.0;
  dy = ( SINPI15*u[0] + COSPI15*u[1]) / 3.0;
  glVertex2dv(c);
  glVertex2d(c[0]-dx,c[1]-dy);
}


GLuint listClipTetraVector(pMesh mesh) {
  pMaterial   pm;
  pTetra      pt;
  pPoint      ppt;
  pSolution   ps0;
  pScene      sc;
  pClip       clip;
  double      rgb[3],u[3],epsra,iso,kc,dd,scal,scalemin,scalemax;
  float       cp[3];
  GLuint      dlist = 0;
  int         ia,k,l,m;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* default */
  if ( !mesh->ntet || !mesh->nbb )  return(0);
  if ( ddebug ) printf("create vector list for clip\n");

  sc   = cv.scene[currentScene()];
  clip = sc->clip;
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* create display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  scalemin = sc->dmax * SCALV;
  scalemax = 10.0*scalemin;
  mesh->mark++;
  glLineWidth(2.0);

  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;
 
    glBegin(GL_LINES);
    while ( k != 0 ) {
      pt = &mesh->tetra[k];
      if ( !pt->v[0] || !pt->clip ) {
        k = pt->nxt;
        continue;
      }
 
      /* element size */
      scal  = 0.5*sizeTetra(mesh,k);
      epsra = EPST * scal;

      /* linear interpol. */
      if ( mesh->typage == 2 )
        for (l=0; l<4; l++) {
          ppt = &mesh->point[pt->v[l]];
          if ( ppt->mark == mesh->mark )  continue;
          ppt->mark = mesh->mark;

          ps0 = &mesh->sol[pt->v[l]];
          dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
                         + ps0->m[2]*ps0->m[2]);
          if ( dd < epsra )  continue;

          /* color =  norm */
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          /* 3D vectors */
          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          u[2]  = ps0->m[2] * dd;
          cp[0] = ppt->c[0];
          cp[1] = ppt->c[1];
          cp[2] = ppt->c[2];
          drawVector3D(cp,u,scal);
        }
      else {
        cp[0] = cp[1] = cp[2] = 0.0;
        for (l=0; l<4; l++) {
          ppt = &mesh->point[pt->v[l]];
          cp[0] += ppt->c[0];
          cp[1] += ppt->c[1];
          cp[2] += ppt->c[2];
        }
        cp[0] *= 0.25;
        cp[1] *= 0.25;
        cp[2] *= 0.25;
        
        /* color = norm */
        ps0 = &mesh->sol[k];
        dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
                        +ps0->m[2]*ps0->m[2]);
        if ( dd > epsra ) {
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          u[2]  = ps0->m[2] * dd;
          drawVector3D(cp,u,scal);
        }
      }
      k = pt->nxt;
    }
    glEnd();
  }
  glLineWidth(1.0);
  glEndList();

  return(dlist);
}

GLuint listClipHexaVector(pMesh mesh) {
  pMaterial   pm;
  pHexa       ph;
  pPoint      ppt;
  pSolution   ps0;
  pScene      sc;
  pClip       clip;
  double      rgb[3],u[3],epsra,iso,kc,dd,scal,scalemin,scalemax;
  float       cp[3];
  GLuint      dlist = 0;
  int         ia,k,l,m;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* default */
  if ( !mesh->nhex || !mesh->nbb )  return(0);
  if ( ddebug ) printf("create vector list for clip\n");

  sc   = cv.scene[currentScene()];
  clip = sc->clip;
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* create display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  scalemin = sc->dmax * SCALV;
  scalemax = 10.0*scalemin;
  mesh->mark++;
  glLineWidth(2.0);

  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LHexa];
    if ( !k || pm->flag )  continue;
 
    glBegin(GL_LINES);
    while ( k != 0 ) {
      ph = &mesh->hexa[k];
      if ( !ph->v[0] || !ph->clip ) {
        k = ph->nxt;
        continue;
      }
 
      /* element size */
      scal  = 0.5*sizeHexa(mesh,k);
      epsra = EPST * scal;

      /* linear interpol. */
      if ( mesh->typage == 2 )
        for (l=0; l<8; l++) {
          ppt = &mesh->point[ph->v[l]];
          if ( ppt->mark == mesh->mark )  continue;
          ppt->mark = mesh->mark;

          ps0 = &mesh->sol[ph->v[l]];
          dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
                         + ps0->m[2]*ps0->m[2]);
          if ( dd < epsra )  continue;

          /* color =  norm */
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          /* 3D vectors */
          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          u[2]  = ps0->m[2] * dd;
          cp[0] = ppt->c[0];
          cp[1] = ppt->c[1];
          cp[2] = ppt->c[2];
          drawVector3D(cp,u,scal);
        }
      else {
        cp[0] = cp[1] = cp[2] = 0.0;
        for (l=0; l<8; l++) {
          ppt = &mesh->point[ph->v[l]];
          cp[0] += ppt->c[0];
          cp[1] += ppt->c[1];
          cp[2] += ppt->c[2];
        }
        cp[0] *= 0.125;
        cp[1] *= 0.125;
        cp[2] *= 0.125;
        
        /* color = norm */
        ps0 = &mesh->sol[k];
        dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
                        +ps0->m[2]*ps0->m[2]);
        if ( dd > epsra ) {
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          u[2]  = ps0->m[2] * dd;
          drawVector3D(cp,u,scal);
        }
      }
      k = ph->nxt;
    }
    glEnd();
  }
  glLineWidth(1.0);
  glEndList();

  return(dlist);
}


GLuint listTria2dVector(pMesh mesh) {
  pMaterial   pm;
  pTriangle   pt;
  pPoint      ppt;
  pSolution   ps0;
  pScene      sc;
  double      rgb[3],u[2],epsra,iso,kc,dd,scalemin,scalemax,scal;
  float       cp[2];
  GLuint      dlist = 0;
  int         ia,i,k,l,m;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* default */
  if ( !mesh->nt || !mesh->nbb )  return(0);

  sc    = cv.scene[currentScene()];
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("create vector list\n");

  /* create display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  scalemin = sc->dmax * SCALV;
  scalemax = 15*scalemin;
  mesh->mark++;
  glLineWidth(3.0);

  if ( mesh->typage == 2 ) {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTria];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        scal  = 0.5*sizeTria(mesh,k);
        epsra = EPST * scal;

        for (i=0; i<3; i++) {
          ppt = &mesh->point[pt->v[i]];
          if ( ppt->mark == mesh->mark )  continue;
          ppt->mark = mesh->mark;

          ps0 = &mesh->sol[pt->v[i]];
          dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]);
          if ( dd < epsra )  continue;

          /* color = norm */
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          cp[0] = ppt->c[0];
          cp[1] = ppt->c[1];
          drawVector2D(cp,u,scal);
        }
        k = pt->nxt; 
      }
    }
    glEnd();
  }
  else {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTria];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        scal  = 0.5*sizeTria(mesh,k);
        epsra = EPST * scal;
        cp[0] = cp[1] = 0.0;
        for (l=0; l<3; l++) {
          ppt = &mesh->point[pt->v[l]];
          cp[0] += ppt->c[0];
          cp[1] += ppt->c[1];
        }
        cp[0] /= 3.0;
        cp[1] /= 3.0;

        /* color = norm */
        ps0 = &mesh->sol[k];
        dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]);
        if ( dd < epsra ) {
          k = pt->nxt; 
          continue;
        }
        if ( iso < sc->iso.val[0] ) 
          iso = sc->iso.val[0];  
        else if ( iso > sc->iso.val[MAXISO-1] )
          iso = sc->iso.val[MAXISO-1];
        for (ia=0; ia<MAXISO-1; ia++)
          if ( iso < sc->iso.val[ia] )  break;
        kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
        hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
        hsvrgb(hsv,rgb);
        glColor3dv(rgb);
        
        dd   = 1.0 / dd;
        u[0] = ps0->m[0] * dd;
        u[1] = ps0->m[1] * dd;
        drawVector2D(cp,u,scal);
        k = pt->nxt;
      }
    }
    glEnd();
  }

  glLineWidth(1.0);
  glEndList();
  return(dlist);
}


GLuint listQuad2dVector(pMesh mesh) {
  pMaterial   pm;
  pQuad       pq;
  pPoint      ppt;
  pSolution   ps0;
  pScene      sc;
  double      rgb[3],u[2],epsra,iso,kc,dd,scalemin,scalemax,scal;
  float       cp[2];
  GLuint      dlist = 0;
  int         ia,i,k,l,m;
  static double hsv[3] = { 0.0, 1.0, 0.80 };

  /* default */
  if ( !mesh->nq || !mesh->nbb )  return(0);

  sc    = cv.scene[currentScene()];
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("create vector list\n");

  /* create display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  scalemin = sc->dmax * SCALV;
  scalemax = 15*scalemin;
  mesh->mark++;
  glLineWidth(3.0);

  if ( mesh->typage == 2 ) {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LQuad];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pq = &mesh->quad[k];
        if ( !pq->v[0] ) {
          k = pq->nxt;
          continue;
        }

        scal  = 0.5*sizeQuad(mesh,k);
        epsra = EPST * scal;

        for (i=0; i<4; i++) {
          ppt = &mesh->point[pq->v[i]];
          if ( ppt->mark == mesh->mark )  continue;
          ppt->mark = mesh->mark;

          ps0 = &mesh->sol[pq->v[i]];
          dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]);
          if ( dd < epsra )  continue;

          /* color = norm */
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          cp[0] = ppt->c[0];
          cp[1] = ppt->c[1];
          drawVector2D(cp,u,scal);
        }
        k = pq->nxt; 
      }
    }
    glEnd();
  }
  else {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LQuad];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pq = &mesh->quad[k];
        if ( !pq->v[0] ) {
          k = pq->nxt;
          continue;
        }

        scal  = 0.5*sizeQuad(mesh,k);
        epsra = EPST * scal;
        cp[0] = cp[1] = 0.0;
        for (l=0; l<4; l++) {
          ppt = &mesh->point[pq->v[l]];
          cp[0] += ppt->c[0];
          cp[1] += ppt->c[1];
        }
        cp[0] *= 0.25;
        cp[1] *= 0.25;

        /* color = norm */
        ps0 = &mesh->sol[k];
        dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]);
        if ( dd < epsra ) {
          k = pq->nxt; 
          continue;
        }
        if ( iso < sc->iso.val[0] ) 
          iso = sc->iso.val[0];  
        else if ( iso > sc->iso.val[MAXISO-1] )
          iso = sc->iso.val[MAXISO-1];
        for (ia=0; ia<MAXISO-1; ia++)
          if ( iso < sc->iso.val[ia] )  break;
        kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
        hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
        hsvrgb(hsv,rgb);
        glColor3dv(rgb);
        
        dd   = 1.0 / dd;
        u[0] = ps0->m[0] * dd;
        u[1] = ps0->m[1] * dd;
        drawVector2D(cp,u,scal);
        k = pq->nxt;
      }
    }
    glEnd();
  }

  glLineWidth(1.0);
  glEndList();
  return(dlist);
}


GLuint listTria3dVector(pMesh mesh) {
  pMaterial   pm;
  pTriangle   pt;
  pPoint      ppt;
  pSolution   ps0;
  pScene      sc;
  double      rgb[3],u[3],epsra,iso,kc,dd,scalemin,scalemax,scal;
  float       cp[3];
  GLuint      dlist = 0;
  int         ia,i,k,l,m;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* default */
  if ( !mesh->nbb )  return(0);

  sc    = cv.scene[currentScene()];
  /*if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);*/
  if ( ddebug ) printf("create vector list\n");

  /* create display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  scalemin = sc->dmax * SCALV;
  scalemax = 15*scalemin;
  mesh->mark++;
  glLineWidth(2.0);

  if ( mesh->typage == 2 ) {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTria];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        scal  = sizeTria(mesh,k);
        epsra = EPST * scal;

        for (i=0; i<3; i++) {
          ppt = &mesh->point[pt->v[i]];
          if ( ppt->mark == mesh->mark )  continue;
          ppt->mark = mesh->mark;

          ps0 = &mesh->sol[pt->v[i]];
          dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
              + ps0->m[2]*ps0->m[2]);
          if ( dd < epsra )  continue;

          /* color = norm */
          if ( iso < sc->iso.val[0] ) 
            iso = sc->iso.val[0];  
          else if ( iso > sc->iso.val[MAXISO-1] )
            iso = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( iso < sc->iso.val[ia] )  break;
          kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);
          glColor3dv(rgb);

          dd    = 1.0 / dd;
          u[0]  = ps0->m[0] * dd;
          u[1]  = ps0->m[1] * dd;
          u[2]  = ps0->m[2] * dd;
          cp[0] = ppt->c[0];
          cp[1] = ppt->c[1];
          cp[2] = ppt->c[2];
          drawVector3D(cp,u,scal);
        }
        k = pt->nxt; 
      }
    }
    glEnd();
  }

  else {
    glBegin(GL_LINES);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTria];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        scal  = 0.5*sizeTria(mesh,k);
        epsra = EPST * scal;
        cp[0] = cp[1] = 0.0;
        for (l=0; l<3; l++) {
          ppt = &mesh->point[pt->v[l]];
          cp[0] += ppt->c[0];
          cp[1] += ppt->c[1];
          cp[2] += ppt->c[2];
        }
        cp[0] /= 3.0;
        cp[1] /= 3.0;
        cp[2] /= 3.0;

        /* color = norm */
        ps0 = &mesh->sol[k];
        dd  = iso = sqrt(ps0->m[0]*ps0->m[0] + ps0->m[1]*ps0->m[1]\
            + ps0->m[2]*ps0->m[2]);
        if ( dd < epsra ) {
          k = pt->nxt; 
          continue;
        }
        if ( iso < sc->iso.val[0] ) 
          iso = sc->iso.val[0];  
        else if ( iso > sc->iso.val[MAXISO-1] )
          iso = sc->iso.val[MAXISO-1];
        for (ia=0; ia<MAXISO-1; ia++)
          if ( iso < sc->iso.val[ia] )  break;
        kc = (iso-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
        hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
        hsvrgb(hsv,rgb);
        glColor3dv(rgb);
        
        dd   = 1.0 / dd;
        u[0] = ps0->m[0]/dd;
        u[1] = ps0->m[1]/dd;
        u[2] = ps0->m[2]/dd;
        drawVector3D(cp,u,scal);
        k = pt->nxt;
      }
    }
    glEnd();
  }

  glLineWidth(1.0);
  glEndList();
  return(dlist);
}
