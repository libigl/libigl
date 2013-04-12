#include "medit.h"
#include "extern.h"
#include "sproto.h"

/* list of element faces */
static int ct[4][3] = { {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3} };
static int ch[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, 
			{1,2,6,5}, {2,3,7,6}, {0,3,7,4} };

static float redcol[4]   = {1.0, 0.0, 0.0, 1.0};
static float greencol[4] = {0.0, 0.6, 0.0, 1.0}; 


/* build list of triangles */
GLuint listTria(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pTriangle  pt;
  pPoint     p0,p1,p2;
  GLuint     dlist;
  double     dd,ax,ay,az,bx,by,bz;
  float      cx,cy,cz,n[3],pp0[3],pp1[3],pp2[3];
  int        k,m,mm,is0,is1,is2,transp;

  /* default */
  if ( !mesh->nt ) return(0);
  if ( ddebug ) printf("create display list / TRIA\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    mm = sc->matsort[m];
    pm = &sc->material[mm];
    k  = pm->depmat[LTria];
    if ( !k || pm->flag )  continue;
    transp = 0;

    if ( !(sc->mode & S_MATERIAL) )
      pm = &sc->material[DEFAULT_MAT];
    transp = pm->amb[3] < 0.999 || pm->dif[3] < 0.999 || pm->spe[3] < 0.999;

#ifdef IGL
    int old_depth_func =0;
    glGetIntegerv(GL_DEPTH_FUNC,&old_depth_func);
#endif
    if ( transp ) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
#ifdef IGL
      glDepthFunc(GL_ALWAYS);
#else
      glDepthMask(GL_FALSE);
#endif
    }

    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);

    glBegin(GL_TRIANGLES);
    if ( sc->type & S_FLAT ) {
      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }
        p0 = &mesh->point[pt->v[0]];
        p1 = &mesh->point[pt->v[1]];
        p2 = &mesh->point[pt->v[2]];

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
        if ( dd > 0.0 ) {
          dd = 1.0 / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }

        if ( sc->shrink < 1.0 ) {
          cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.0;
          cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.0;
          cz = (p0->c[2] + p1->c[2] + p2->c[2]) / 3.0;
          pp0[0] = sc->shrink*(p0->c[0]-cx)+cx;
          pp0[1] = sc->shrink*(p0->c[1]-cy)+cy;
          pp0[2] = sc->shrink*(p0->c[2]-cz)+cz;
          pp1[0] = sc->shrink*(p1->c[0]-cx)+cx;
          pp1[1] = sc->shrink*(p1->c[1]-cy)+cy;
          pp1[2] = sc->shrink*(p1->c[2]-cz)+cz;
          pp2[0] = sc->shrink*(p2->c[0]-cx)+cx;
          pp2[1] = sc->shrink*(p2->c[1]-cy)+cy;
          pp2[2] = sc->shrink*(p2->c[2]-cz)+cz;
          glNormal3fv(n);  glVertex3fv(pp0);
          glNormal3fv(n);  glVertex3fv(pp1);
          glNormal3fv(n);  glVertex3fv(pp2);
        }
        else {
          glNormal3fv(n);  glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
          glNormal3fv(n);  glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
          glNormal3fv(n);  glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
        }
        k = pt->nxt;
      }
    }
    else {
      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }
        p0 = &mesh->point[pt->v[0]];
        p1 = &mesh->point[pt->v[1]];
        p2 = &mesh->point[pt->v[2]];

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
        if ( dd > 0.0 ) {
          dd = 1.0 / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }

        is0 = is1 = is2 = 0;
        if ( mesh->extra->iv ) {
          if ( pt->v[0] <= mesh->nvn )
            is0 = mesh->extra->nv[pt->v[0]];
          if ( pt->v[1] <= mesh->nvn )
            is1 = mesh->extra->nv[pt->v[1]];
          if ( pt->v[2] <= mesh->nvn )
            is2 = mesh->extra->nv[pt->v[2]];
        }
        if ( !is0 && pt->v[0] <= mesh->extra->it )
          is0 = mesh->extra->nt[3*(k-1)+1];
        if ( !is1 && pt->v[1] <= mesh->extra->it )
          is1 = mesh->extra->nt[3*(k-1)+2];
        if ( !is2 && pt->v[2] <= mesh->extra->it )
          is2 = mesh->extra->nt[3*(k-1)+3];

        if ( sc->shrink < 1.0 ) {
          cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.;
          cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.;
          cz = (p0->c[2] + p1->c[2] + p2->c[2]) / 3.;
          pp0[0] = sc->shrink*(p0->c[0]-cx)+cx;
          pp0[1] = sc->shrink*(p0->c[1]-cy)+cy;
          pp0[2] = sc->shrink*(p0->c[2]-cz)+cz;
          pp1[0] = sc->shrink*(p1->c[0]-cx)+cx;
          pp1[1] = sc->shrink*(p1->c[1]-cy)+cy;
          pp1[2] = sc->shrink*(p1->c[2]-cz)+cz;
          pp2[0] = sc->shrink*(p2->c[0]-cx)+cx;
          pp2[1] = sc->shrink*(p2->c[1]-cy)+cy;
          pp2[2] = sc->shrink*(p2->c[2]-cz)+cz;
          
          if ( !is0 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is0-1)+1]);
          glVertex3fv(pp0);
          if ( !is1 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is1-1)+1]);
          glVertex3fv(pp1);
          if ( !is2 )
            glNormal3fv(n);
          else
            glNormal3fv(&mesh->extra->n[3*(is2-1)+1]);
          glVertex3fv(pp2);
        }
        else {
          if ( !is0 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is0-1)+1]);
          glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
          if ( !is1 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is1-1)+1]);
          glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
          if ( !is2 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is2-1)+1]);
          glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
        }
        k = pt->nxt;
      }
    }
    glEnd();
    if ( transp ) {
#ifdef IGL
      glDepthFunc(old_depth_func);
#else
      glDepthMask(GL_TRUE);
#endif
      glDisable(GL_BLEND);
    }
  }

  glEndList();
  return(dlist);
}


/* build list of quadrilaterals */
GLuint listQuad(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pQuad      pq;
  pPoint     p0,p1,p2,p3;
  GLuint     dlist = 0;
  double     ax,ay,az,bx,by,bz,dd;
  float      cx,cy,cz,n[3],pp0[3],pp1[3],pp2[3],pp3[3];
  int        k,m,mm,is0,is1,is2,is3,transp;

  /* default */
  if ( !mesh->nq ) return(0);
  if ( ddebug ) printf("create display list / QUADS\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    mm = sc->matsort[m];
    pm = &sc->material[mm];
    k  = pm->depmat[LQuad];
    if ( !k || pm->flag )  continue;
    transp = 0;

    if ( !(sc->mode & S_MATERIAL) )
      pm = &sc->material[DEFAULT_MAT];
    transp = pm->amb[3] < 0.999 || pm->dif[3] < 0.999 || pm->spe[3] < 0.999;
    if ( transp ) {
      glEnable(GL_BLEND);
      glDepthMask(GL_FALSE);
      glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    }
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);      
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);

    glBegin(GL_QUADS);
    if ( sc->type & S_FLAT ) {
      while ( k != 0 ) {
        pq = &mesh->quad[k];
        if ( pq->v[0] == 0 ) {
          k = pq->nxt;
          continue;
        }
        p0 = &mesh->point[pq->v[0]];
        p1 = &mesh->point[pq->v[1]];
        p2 = &mesh->point[pq->v[2]];
        p3 = &mesh->point[pq->v[3]];

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
        if ( dd > 0.0 ) {
          dd = 1.0 / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }

        if ( sc->shrink < 1.0 ) {
          cx = (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]) / 4.;
          cy = (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]) / 4.;
          cz = (p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]) / 4.;
          pp0[0] = sc->shrink*(p0->c[0]-cx)+cx;
          pp0[1] = sc->shrink*(p0->c[1]-cy)+cy;
          pp0[2] = sc->shrink*(p0->c[2]-cz)+cz;
          pp1[0] = sc->shrink*(p1->c[0]-cx)+cx;
          pp1[1] = sc->shrink*(p1->c[1]-cy)+cy;
          pp1[2] = sc->shrink*(p1->c[2]-cz)+cz;
          pp2[0] = sc->shrink*(p2->c[0]-cx)+cx;
          pp2[1] = sc->shrink*(p2->c[1]-cy)+cy;
          pp2[2] = sc->shrink*(p2->c[2]-cz)+cz;
          pp3[0] = sc->shrink*(p3->c[0]-cx)+cx;
          pp3[1] = sc->shrink*(p3->c[1]-cy)+cy;
          pp3[2] = sc->shrink*(p3->c[2]-cz)+cz;
          glNormal3fv(n);  glVertex3fv(pp0);
          glNormal3fv(n);  glVertex3fv(pp1);
          glNormal3fv(n);  glVertex3fv(pp2);
          glNormal3fv(n);  glVertex3fv(pp3);
        }
        else {
          glNormal3fv(n);  glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
          glNormal3fv(n);  glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
          glNormal3fv(n);  glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
          glNormal3fv(n);  glVertex3f(p3->c[0],p3->c[1],p3->c[2]);
        }
        k = pq->nxt;
      }
    }
    else {
      while ( k != 0 ) {
        pq = &mesh->quad[k];
        if ( !pq->v[0] ) {
          k = pq->nxt;
          continue;
        }
        p0 = &mesh->point[pq->v[0]];
        p1 = &mesh->point[pq->v[1]];
        p2 = &mesh->point[pq->v[2]];
        p3 = &mesh->point[pq->v[3]];

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
        if ( dd > 0.0 ) {
          dd = 1.0 / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }
        is0 = is1 = is2 = is3 = 0;
        if ( mesh->extra->iv ) {
          if ( pq->v[0] <= mesh->nvn )
            is0 = mesh->extra->nv[pq->v[0]];
          if ( pq->v[1] <= mesh->nvn )
            is1 = mesh->extra->nv[pq->v[1]];
          if ( pq->v[2] <= mesh->nvn )
            is2 = mesh->extra->nv[pq->v[2]];
          if ( pq->v[3] <= mesh->nvn )
            is3 = mesh->extra->nv[pq->v[3]];
        }
        if ( !is0 && pq->v[0] <= mesh->extra->iq )
          is0 = mesh->extra->nq[4*(k-1)+1];
        if ( !is1 && pq->v[1] <= mesh->extra->iq )
          is1 = mesh->extra->nq[4*(k-1)+2];
        if ( !is2 && pq->v[2] <= mesh->extra->iq )
          is2 = mesh->extra->nq[4*(k-1)+3];
        if ( !is3 && pq->v[3] <= mesh->extra->iq )
          is3 = mesh->extra->nq[4*(k-1)+4];

        if ( sc->shrink < 1.0 ) {
          cx = (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]) / 4.;
          cy = (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]) / 4.;
          cz = (p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]) / 4.;
          pp0[0] = sc->shrink*(p0->c[0]-cx)+cx;
          pp0[1] = sc->shrink*(p0->c[1]-cy)+cy;
          pp0[2] = sc->shrink*(p0->c[2]-cz)+cz;
          pp1[0] = sc->shrink*(p1->c[0]-cx)+cx;
          pp1[1] = sc->shrink*(p1->c[1]-cy)+cy;
          pp1[2] = sc->shrink*(p1->c[2]-cz)+cz;
          pp2[0] = sc->shrink*(p2->c[0]-cx)+cx;
          pp2[1] = sc->shrink*(p2->c[1]-cy)+cy;
          pp2[2] = sc->shrink*(p2->c[2]-cz)+cz;
          pp3[0] = sc->shrink*(p3->c[0]-cx)+cx;
          pp3[1] = sc->shrink*(p3->c[1]-cy)+cy;
          pp3[2] = sc->shrink*(p3->c[2]-cz)+cz;

          if ( !is0 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is0-1)+1]);
          glVertex3fv(pp0);
          if ( !is1 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is1-1)+1]);
          glVertex3fv(pp1);
          if ( !is2 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is2-1)+1]);
          glVertex3fv(pp2);
          if ( !is3 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is3-1)+1]);
          glVertex3fv(pp3);
        }
        else {

          if ( !is0 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is0-1)+1]);
          glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
          if ( !is1 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is1-1)+1]);
          glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
          if ( !is2 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is2-1)+1]);
          glVertex3f(p2->c[0],p2->c[1],p2->c[2]);
          if ( !is3 )
            glNormal3fv(n);  
          else
            glNormal3fv(&mesh->extra->n[3*(is3-1)+1]);
          glVertex3f(p3->c[0],p3->c[1],p3->c[2]);
        }
        k = pq->nxt;
      }
    }
    glEnd();
    if ( transp ) {
      glDepthMask(GL_TRUE);
      glDisable(GL_BLEND);
    }
  }

  glEndList();
  return(dlist);
}


/* build list of tetrahedra */
GLuint listTetra(pScene sc,pMesh mesh,ubyte clip) {
  pMaterial  pm;
  pTetra     pt;
  pPoint     p0,p1,p2;
  GLuint     dlist;
  double     ax,ay,az,bx,by,bz,d;
  float      n[3],shrink,cx,cy,cz;
  int        k,l,m,mm;

  if ( !mesh->ntet )  return(0);
  if ( ddebug ) printf("create 3d display list w. materials\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
  shrink = MEDIT_MIN(sc->shrink,0.995f);

  /* scan tetra */
  for (m=0; m<sc->par.nbmat; m++) {
    mm = sc->matsort[m];
    pm = &sc->material[mm];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;

    bool transp = 0;
#ifdef IGL
    int old_depth_func =0;
    glGetIntegerv(GL_DEPTH_FUNC,&old_depth_func);
#endif

    transp = sc->igl_params->tet_color[3]<0.999;
    if ( transp ) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
#ifdef IGL
      glDepthFunc(GL_ALWAYS);
#else
      glDepthMask(GL_FALSE);
#endif
    }
    if ( sc->mode & S_MATERIAL ) {

      //if ( pm->dif[3] < 0.999 ) {
      //  glDepthMask(GL_FALSE);
      //  glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      //  glEnable(GL_BLEND);
      //}
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
      glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);
    }
    else
#ifdef IGL
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,sc->igl_params->tet_color);
#else
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,redcol);
#endif
    
    /* display triangular faces */
    glBegin(GL_TRIANGLES);
    while ( k != 0 ) {
      pt = &mesh->tetra[k];
      if ( !pt->v[0] || (clip && !pt->clip) ) {
        k = pt->nxt;
        continue;
      }
      /* build 4 faces */
      cx = cy = cz = 0.;
      for (l=0; l<4; l++) {
        p0  = &mesh->point[pt->v[l]];
        cx += p0->c[0];
        cy += p0->c[1];
        cz += p0->c[2];
      }
      cx *= 0.25; 
      cy *= 0.25;
      cz *= 0.25;

      for (l=0; l<4; l++) {
	p0 = &mesh->point[pt->v[ct[l][0]]];
	p1 = &mesh->point[pt->v[ct[l][1]]];
	p2 = &mesh->point[pt->v[ct[l][2]]];

	/* compute face normal */
	ax = p1->c[0] - p0->c[0]; ay = p1->c[1] - p0->c[1]; az = p1->c[2] - p0->c[2];
	bx = p2->c[0] - p0->c[0]; by = p2->c[1] - p0->c[1]; bz = p2->c[2] - p0->c[2];
	n[0] = ay*bz - az*by;
	n[1] = az*bx - ax*bz;
	n[2] = ax*by - ay*bx;
	d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
	if ( d > 0.0f ) {
	  d = 1.0f / sqrt(d);
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
      k = pt->nxt;
    }
    glEnd();
    if ( transp ) {
#ifdef IGL
      glDepthFunc(old_depth_func);
#else
      glDepthMask(GL_TRUE);
#endif
      glDisable(GL_BLEND);
    }
  }

  glEndList();
  return(dlist);
}


/* build list of hexahedra */
GLuint listHexa(pScene sc,pMesh mesh,ubyte clip) {
  pMaterial  pm;
  pHexa      ph;
  pPoint     p0,p1,p2,p3;
  GLuint     dlist;
  double     ax,ay,az,bx,by,bz,d;
  float      n[3],cx,cy,cz,shrink;
  int        k,l,m,mm;

  if ( !mesh->nhex )  return(0);
  if ( ddebug ) printf("create 3d display list w. materials\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
  shrink = MEDIT_MIN(sc->shrink,0.995f);

  /* scan hexa */
  for (m=0; m<sc->par.nbmat; m++) {
    mm = sc->matsort[m];
    pm = &sc->material[mm];
    k  = pm->depmat[LHexa];
    if ( !k || pm->flag )  continue;

    if ( sc->mode & S_MATERIAL ) {
      if ( pm->dif[3] < 0.999) {
        glDepthMask(GL_FALSE);
        glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
      }
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
      glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);
    }
    else
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,greencol);

    /* display quadrilateral faces */
    glBegin(GL_QUADS);
    while ( k != 0 ) {
      ph = &mesh->hexa[k];
      if ( !ph->v[0] || (clip && !ph->clip) ) {
        k = ph->nxt;
        continue;
      }

      /* build 6 faces */
      cx = cy = cz = 0.;
      for (l=0; l<8; l++) {
        p0  = &mesh->point[ph->v[l]];
        cx += p0->c[0];
        cy += p0->c[1];
        cz += p0->c[2];
      }
      cx /= 8.;
      cy /= 8.;
      cz /= 8.;
      for (l=0; l<6; l++) {
        p0 = &mesh->point[ph->v[ch[l][0]]];
        p1 = &mesh->point[ph->v[ch[l][1]]];
        p2 = &mesh->point[ph->v[ch[l][2]]];
        p3 = &mesh->point[ph->v[ch[l][3]]];

        /* compute face normal */
        ax = p1->c[0] - p0->c[0]; ay = p1->c[1] - p0->c[1]; az = p1->c[2] - p0->c[2];
        bx = p2->c[0] - p0->c[0]; by = p2->c[1] - p0->c[1]; bz = p2->c[2] - p0->c[2];
        n[0] = ay*bz - az*by;
        n[1] = az*bx - ax*bz;
        n[2] = ax*by - ay*bx;
        d = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if ( d > 0.0f ) {
          d = 1.0f / sqrt(d);
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
        glNormal3fv(n);
        glVertex3f(shrink*(p3->c[0]-cx)+cx,
                   shrink*(p3->c[1]-cy)+cy,
                   shrink*(p3->c[2]-cz)+cz);
      }
      k = ph->nxt;
    }
    glEnd();
    if ( sc->mode & S_MATERIAL && pm->dif[3] < 0.999 ) {
      glDepthMask(GL_TRUE);
      glDisable(GL_BLEND);
    }
  }

  glEndList();
  return(dlist);
}
