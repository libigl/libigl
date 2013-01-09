#include "medit.h"
#include "extern.h"
#include "sproto.h"
#include "eigenv.h"

#ifndef M_PI2
#define M_PI2  (2.0*M_PI)
#endif

extern int refmat;
extern int eigen2(double m[3],double lambda[2],double vp[2][2]);

static GLfloat IdMatrix[16] = {
   1.0, 0.0, 0.0, 0.0,
   0.0, 1.0, 0.0, 0.0,
   0.0, 0.0, 1.0, 0.0,
   0.0, 0.0, 0.0, 1.0
};

void drawEllipsoid(pScene sc,pMesh mesh,int typel,int k) {
  pMaterial    pm;
  pSolution    ps;
  pTriangle    pt;
  pTetra       pt1;
  pPoint       p0;
  GLfloat      mat[16],cx,cy,cz;
  double       m[6],lambda[3],v[3][3];
  int          i,j,l,iord;

  /* compute average size */
  if ( mesh->nfield != 6 || !mesh->nbb )  return;

  /* draw average ellipse at element */
  if ( typel == LPoint ) {
    p0 = &mesh->point[k];
    /*pm = &sc->material[refmat];*/
    pm = &sc->material[p0->ref];
    ps = &mesh->sol[k];
    for (j=0; j<6; j++) m[j] = ps->m[j];

    iord = eigenv(1,m,lambda,v);
    if ( !iord )  return;

    if ( mesh->ne ) {    
      fprintf(stdout,"  Eigenvectors :\n   vp1 : %f  %f  %f\n",
              v[0][0],v[0][1],v[0][2]);
      fprintf(stdout,"   vp2 : %f  %f  %f\n",v[1][0],v[1][1],v[1][2]);
      fprintf(stdout,"   vp3 : %f  %f  %f\n",v[2][0],v[2][1],v[2][2]);
      fprintf(stdout,"  Eigenvalues  : %f  %f  %f    %d\n",
              lambda[0],lambda[1],lambda[2],iord);
      if ( lambda[0] <= 0.0 || lambda[1] <= 0.0 || lambda[2] <= 0.0)
        return;
      fprintf(stdout,"  Sizes        : %f  %f  %f\n",
              1.0/sqrt(lambda[0]),1.0/sqrt(lambda[1]),1.0/sqrt(lambda[2]));
    }
    else if ( lambda[0] <= 0.0 || lambda[1] <= 0.0 || lambda[2] <= 0.0)
      return;

    lambda[0] = MEDIT_MAX(EPS,0.5/sqrt(lambda[0]));
    lambda[1] = MEDIT_MAX(EPS,0.5/sqrt(lambda[1]));
    lambda[2] = MEDIT_MAX(EPS,0.5/sqrt(lambda[2]));
    memcpy(mat,IdMatrix,16*sizeof(GLfloat));
    for (j=0; j<3; j++) 
      for (l=0; l<3; l++)
        mat[j*4+l] = v[j][l];

    glDisable(GL_LIGHTING);
    glPushMatrix();
    glColor4fv(pm->dif);
    glTranslatef(p0->c[0],p0->c[1],p0->c[2]);
    glMultMatrixf(mat);
    glScalef(lambda[0],lambda[1],lambda[2]);
    glutWireSphere(1.,30,30);
    glPopMatrix();
    glEnable(GL_LIGHTING);
  }

  else if ( typel == LTria ) {
    pt = &mesh->tria[k];
    if ( mesh->nbb == mesh->np )
      pm = &sc->material[refmat];
    else if ( mesh->nbb == mesh->nt )
	  pm = &sc->material[k];
	else
	  return;
    glColor4fv(pm->dif);
    for (j=0; j<6; j++) m[j] = 0.;
    cx = cy = cz = 0.0;
    for (i=0; i<3; i++) {
      ps = &mesh->sol[pt->v[i]];
      p0 = &mesh->point[pt->v[i]];
      cx += p0->c[0];
      cy += p0->c[1];
      cz += p0->c[2];
      for (j=0; j<6; j++)
        m[j] += ps->m[j];
    }
    cx /= 3.;  cy /= 3.;  cz /= 3.;
    for (j=0; j<6; j++)  m[j] /= 6.;

    if ( !eigenv(1,m,lambda,v) )  return;
    lambda[0] = MEDIT_MAX(EPS,0.5/sqrt(lambda[0]));
    lambda[1] = MEDIT_MAX(EPS,0.5/sqrt(lambda[1]));
    lambda[2] = MEDIT_MAX(EPS,0.5/sqrt(lambda[2]));

    memcpy(mat,IdMatrix,16*sizeof(GLfloat));
    for (j=0; j<3; j++) 
      for (l=0; l<3; l++)
        mat[j*4+l] = v[j][l];

    glDisable(GL_LIGHTING);
    glPushMatrix();
    glColor4fv(pm->dif);
    glTranslatef(cx,cy,cz);
    glMultMatrixf(mat);
    glScalef(lambda[0],lambda[1],lambda[2]);
    glutWireSphere(1.,30,30);
    glPopMatrix();
    glEnable(GL_LIGHTING);
  }

  else if ( typel == LTets ) {
    pt1 = &mesh->tetra[k];
    pm = &sc->material[refmat];
    glColor4fv(pm->dif);
    for (j=0; j<6; j++)  m[j] = 0.;
    cx = cy = cz = 0.0;
    for (i=0; i<4; i++) {
      if ( mesh->nbb == mesh->np )
        ps = &mesh->sol[pt1->v[i]];
	  else if ( mesh->nbb == mesh->ntet )
        ps = &mesh->sol[k];
      else
		return;
      p0 = &mesh->point[pt1->v[i]];
      cx += p0->c[0];
      cy += p0->c[1];
      cz += p0->c[2];
      for (j=0; j<6; j++)
        m[j] += ps->m[j];
    }
    cx /= 4.;  cy /= 4.;  cz /= 4.;
    for (j=0; j<6; j++)  m[j] /= 6.;

    if ( !eigenv(1,m,lambda,v) )  return;
    lambda[0] = MEDIT_MAX(EPS,0.5/sqrt(lambda[0]));
    lambda[1] = MEDIT_MAX(EPS,0.5/sqrt(lambda[1]));
    lambda[2] = MEDIT_MAX(EPS,0.5/sqrt(lambda[2]));

    memcpy(mat,IdMatrix,16*sizeof(GLfloat));
    for (j=0; j<3; j++) 
      for (l=0; l<3; l++)
        mat[j*4+l] = v[j][l];
    glDisable(GL_LIGHTING);
    glPushMatrix();
    glColor4fv(pm->dif);
    glTranslatef(cx,cy,cz);
    glMultMatrixf(mat);
    glScalef(lambda[0],lambda[1],lambda[2]);
    glutWireSphere(1.,30,30);
    glPopMatrix();
    glEnable(GL_LIGHTING);
	
  }
  else return;
}


void glCircle(float radius) {
  float   ang,ux,uy;

  ux = 0.0f;
  uy = radius;

  glBegin(GL_LINE_STRIP);
  for (ang=0.0f; ang<=2*M_PI+0.2; ang+=0.2) {
    ux = radius*(float)sin((double)ang);
    uy = radius*(float)cos((double)ang);
    glVertex2f(ux,uy);
  }
  glEnd();
}


void drawEllipse(pScene sc,pMesh mesh,int typel,int k) {
  pMaterial    pm;
  pSolution    ps;
  pPoint       p0;
  double       m[3],vp[2][2],lambda[2],dd1,dd2;
  float        theta;

  /* draw ellipse at vertex */
  if ( typel == LPoint ) {
    ps = &mesh->sol[k];
    p0 = &mesh->point[k];
    pm = &sc->material[refmat];

    m[0] = ps->m[0];
    m[1] = ps->m[1];
    m[2] = ps->m[2];
    if ( !eigen2(m,lambda,vp) ) return;

    /* consider eigenvalues as sizes */
    dd1 = 1.0 / sqrt(fabs(lambda[0]));
    dd2 = 1.0 / sqrt(fabs(lambda[1]));

    glDisable(GL_LIGHTING);
    glPushMatrix();
    glLineWidth(1.0);
    glBegin(GL_LINES);
      glColor3fv(pm->dif);
      glVertex3f(p0->c[0],p0->c[1],0.0);
      glColor3fv(pm->dif);
      glVertex3f(p0->c[0]+dd1*vp[0][0],p0->c[1]+dd1*vp[0][1],0.0);

      glColor3fv(pm->dif);
      glVertex3f(p0->c[0],p0->c[1],0.0);
      glColor3fv(pm->dif);
      glVertex3f(p0->c[0]+dd2*vp[1][0],p0->c[1]+dd2*vp[1][1],0.0);
    glEnd();

    theta = atan2(vp[0][1],vp[0][0])*RTOD;
    glTranslatef(p0->c[0],p0->c[1],0.0);
    glRotatef(theta,0.0,0.0,1.0);
    glScaled(dd1,dd2,0.0);

    glColor3fv(pm->dif);
    glCircle(1.0);
    glLineWidth(1.0);
    glPopMatrix();
    glEnable(GL_LIGHTING);

    /* print out info */
    fprintf(stdout,"  Eigenvectors :\n");
    fprintf(stdout,"   vp1 : %f  %f\n",vp[0][0],vp[0][1]);
    fprintf(stdout,"   vp2 : %f  %f\n",vp[1][0],vp[1][1]);
    fprintf(stdout,"  Eigenvalues  : %f  %f\n",lambda[0],lambda[1]);
    fprintf(stdout,"  Sizes        : %f  %f\n",dd1,dd2);
  }
}


GLuint drawAllEllipse(pScene sc,pMesh mesh) {
  GLuint       dlist;
  pSolution    ps;
  pMaterial    pm;
  pTriangle    pt;
  pPoint       p0;
  double       m[3],vp[2][2],lambda[2],dd1,dd2;
  float        theta,cx,cy;
  int          k,i,ref;

  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* draw ellipse at vertex */
  glDisable(GL_LIGHTING);
  glLineWidth(1.0);
  
  if ( mesh->typage == 1 ) {
    for (k=1; k<=mesh->ne; k++) {
      ps = &mesh->sol[k];
      pt = &mesh->tria[k];
      if ( !pt->v[0] )  continue;

      ref = matRef(sc,pt->ref);
      pm  = &sc->material[ref];
      if ( pm->flag )   continue;

      cx = cy = 0.0;
      for (i=0; i<3; i++) {
        p0 = &mesh->point[pt->v[i]];
        cx += p0->c[0];
        cy += p0->c[1];
      }
      cx *= 1. / 3.;
      cy *= 1. / 3.;

      m[0] = ps->m[0];
      m[1] = ps->m[1];
      m[2] = ps->m[2];
      if ( !eigen2(m,lambda,vp) ) return (0) ;

      /* consider eigenvalues as sizes */
      dd1 = 1.0 / sqrt(fabs(lambda[0]));
      dd2 = 1.0 / sqrt(fabs(lambda[1]));

      glPushMatrix();
      theta = atan2(vp[0][1],vp[0][0])*RTOD;
      glTranslatef(cx,cy,0.0);
      glRotatef(theta,0.0,0.0,1.0);
      glScaled(dd1,dd2,0.0);
      glColor3fv(pm->dif);
      glCircle(1.0);
      glPopMatrix();
    }
  }
  else {
    for (k=1; k<=mesh->np; k++) {
      ps = &mesh->sol[k];
      p0 = &mesh->point[k];

      ref = matRef(sc,p0->ref);
      pm  = &sc->material[ref];
      if ( pm->flag )   continue;

      m[0] = ps->m[0];
      m[1] = ps->m[1];
      m[2] = ps->m[2];
      if ( !eigen2(m,lambda,vp) ) return (0) ;

      /* consider eigenvalues as sizes */
      dd1 = 1.0 / sqrt(fabs(lambda[0]));
      dd2 = 1.0 / sqrt(fabs(lambda[1]));

      glPushMatrix();
      theta = atan2(vp[0][1],vp[0][0])*RTOD;
      glTranslatef(p0->c[0],p0->c[1],0.0);
      glRotatef(theta,0.0,0.0,1.0);
      glScaled(dd1,dd2,0.0);
      glColor3fv(pm->dif);
      glCircle(1.0);
      glPopMatrix();
    }
  }
  glEnable(GL_LIGHTING);

  glEndList();
  return(dlist);
}


void circumSphere(pScene sc,pMesh mesh,int typel,int k) {
  pMaterial    pm;
  double       c[3],rad;

  cenrad(mesh,k,c,&rad);
  rad = sqrt(rad);
  pm = &sc->material[refmat];

  glDisable(GL_LIGHTING);
  glPushMatrix();
  glColor4fv(pm->dif);
  glTranslated(c[0],c[1],c[2]);
  glScalef(rad,rad,rad);
  glutWireSphere(1.,30,30);
  glPopMatrix();
  glEnable(GL_LIGHTING);
}

