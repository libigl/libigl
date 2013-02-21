#include "medit.h"
#include "extern.h"
#include "sproto.h"

GLuint geomList(pScene sc,pMesh mesh) {
  GLuint     list = 0;
  pMaterial  pm;
  pEdge      pr;
  pPoint     ppt,pp0,pp1;
  double     dd;
  float      n[3];
  int        k,it = 0,nm;
  static float green[4] = {0.0, 1.0, 0.0, 1.0};
  static float rouge[4] = {1.0, 0.0, 0.0, 1.0};
  static float jaune[4] = {1.0, 1.0, 0.0, 1.0};

  /* default */
  if ( mesh->na+mesh->nc+mesh->np == 0 )  return(0);

  /* create display list */
  list = glGenLists(1);
  if ( !list )  return(0);
  glNewList(list,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* draw corners, ridges and required items */
  if ( ddebug ) printf("construct point list\n");
  if ( mesh->ne ) {
    /*glPointSize(3);*/
    glPointSize(sc->par.pointsize);  /* pour Herve LeDret */
    glBegin(GL_POINTS);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( ppt->tag & M_UNUSED && !ppt->ref )  continue;
      if ( ppt->tag == M_CORNER )
        glColor3fv(rouge);
      else if ( ppt->tag == M_REQUIRED )
        glColor3fv(green);
      else continue;
      it++;
      if ( sc->par.linc == 1 )  glColor3fv(sc->par.edge);
      glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
    }
    glEnd();
    glPointSize(1);
  }
  else {
    pm = &sc->material[DEFAULT_MAT];
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);
    glBegin(GL_POINTS);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      n[0] = ppt->c[0] - sc->cx;
      n[1] = ppt->c[1] - sc->cy;
      n[2] = ppt->c[2] - sc->cz;
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd > 0.0f ) {
        dd = 1.0 / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }
      glNormal3fv(n);
      glVertex3f(ppt->c[0],ppt->c[1],ppt->c[2]);
    }
    glEnd(); 
    it = mesh->np;
  }

  /* draw edges */
#ifdef IGL
  for(int pass = 0;pass<2;pass++)
  {
#endif
#ifdef IGL
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
#endif
  if ( ddebug )  printf("construct edge list\n");
  glLineWidth(sc->par.linewidth);
#ifdef IGL
  if(pass == 0)
  {
    glPointSize(sc->igl_params->dot_size);
    glEnable( GL_POINT_SMOOTH );
    glBegin(GL_POINTS);
  }else
  {
    glBegin(GL_LINES);
  }
#endif
  //glBegin(GL_LINES);
  for (k=1; k<=mesh->na; k++) {
    pr = &mesh->edge[k];
    if ( pr->v[0] > mesh->np || pr->v[1] > mesh->np )
      continue;

    if ( pr->tag & M_RIDGE ) {
      if ( pr->tag & M_TAG )
	    glColor3fv(jaune);  /* ridge + ref en jaune */
      else
	    glColor3fv(rouge);  /* ridges en rouge */
    }
    else if ( !pr->ref ) {
      glColor3fv(sc->par.edge);
    }
    else {
      nm = matRef(sc,pr->ref);
      pm = &sc->material[nm];
      glColor3fv(pm->dif);
    }
    if ( sc->par.linc == 1 )  glColor3fv(sc->par.edge);
#ifdef IGL
      if(pr->ref == 1)
      {
        // Non-manifold
        glColor3fv(sc->igl_params->nme_color);
      }else if(pr->ref == 2)
      {
        // Boundary
        glColor3fv(sc->igl_params->open_color);
      }
#endif
    pp0 = &mesh->point[pr->v[0]];
    pp1 = &mesh->point[pr->v[1]];
    glVertex3f(pp0->c[0],pp0->c[1],pp0->c[2]);
    glVertex3f(pp1->c[0],pp1->c[1],pp1->c[2]);
    it++;
  }
  glEnd();
#ifdef IGL
  }
#endif
  glLineWidth(1.0);
  glEndList();

  if ( it == 0 ) {
    glDeleteLists(list,1);
    return(0);
  }
  else
    return(list);
}

