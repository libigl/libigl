#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define NBCOL    9

/* build list for capping */
GLuint capTetra(pMesh mesh) {
  pScene     sc;
  pClip      clip;
  GLuint     dlist = 0;
  pTetra     pt;
  pPoint     p0,p1;
  pMaterial  pm;
  double     dd1[6],d,ax,ay,az,bx,by,bz;
  double     cx[4],cy[4],cz[4],cc;
  float      n[3];
  int        m,k,k1,k2,l,transp,pos[6],neg[6],nbpos,nbneg,nbnul;
  static int tn[4] = {0,0,1,1};
  static int tp[4] = {0,1,1,0};

  /* default */
  if (!mesh->ntet )  return(0);
  if ( ddebug ) printf("create capping list / TETRA\n");
 
  /* build display list */
  sc   = cv.scene[currentScene()];
  clip = sc->clip;
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;

    transp = 0;
    if ( !(sc->mode & S_MATERIAL) )
      pm = &sc->material[DEFAULT_MAT+1];

    /* check transparency */
    transp = pm->amb[3] < 0.999 || pm->dif[3] < 0.999 || pm->spe[3] < 0.999;
    if ( transp ) {
      glDepthMask(GL_FALSE);
      glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
    }
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,pm->dif);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,pm->amb);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,pm->spe);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,pm->emi);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&pm->shininess);
 
    while ( k != 0 ) {
      pt = &mesh->tetra[k];
      if ( !pt->v[0] || !pt->clip ) {
        k = pt->nxt;
        continue;
      }
      nbpos = nbneg = nbnul = 0;
      for (l=0; l<4; l++) {
        p0 = &mesh->point[pt->v[l]];
        if ( p0->clip == 2 )      pos[nbpos++] = l;
        else if (p0->clip == 1 )  neg[nbneg++] = l;
        else                      nbnul++;
        dd1[l] = p0->c[0]*clip->eqn[0] + p0->c[1]*clip->eqn[1] \
               + p0->c[2]*clip->eqn[2] + clip->eqn[3];
      }

      if ( nbneg == 2 && nbpos == 2 ) {
        /* display quadrilateral */
        for (l=0; l<4; l++) {
          k1 = neg[tn[l]];
          k2 = pos[tp[l]];
          p0 = &mesh->point[pt->v[k1]];
          p1 = &mesh->point[pt->v[k2]];
          cc = 1.0f;
          if ( dd1[k2]-dd1[k1] != 0.0f )
            cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
          cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
          cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
          cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
        }

        /* compute face normal */
        ax = cx[1]-cx[0]; ay = cy[1]-cy[0]; az = cz[1]-cz[0];
        bx = cx[2]-cx[0]; by = cy[2]-cy[0]; bz = cz[2]-cz[0];
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
        glBegin(GL_QUADS);
        glNormal3fv(n);
        glVertex3f(cx[0],cy[0],cz[0]);
        glVertex3f(cx[1],cy[1],cz[1]);
        glVertex3f(cx[2],cy[2],cz[2]);
        glVertex3f(cx[3],cy[3],cz[3]);
        glEnd();
      }
      
      else {
        /* display triangle */
        for (l=0; l<3; l++) {
          k1 = nbneg == 3 ? neg[l] : pos[l];
          k2 = nbneg == 3 ? pos[0] : neg[0];
          p0 = &mesh->point[pt->v[k1]];
          p1 = &mesh->point[pt->v[k2]];
          cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
          cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
          cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
          cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
        }
        /* compute face normal */
        ax = cx[1]-cx[0]; ay = cy[1]-cy[0]; az = cz[1]-cz[0];
        bx = cx[2]-cx[0]; by = cy[2]-cy[0]; bz = cz[2]-cz[0];
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
        glBegin(GL_TRIANGLES);
        glNormal3fv(n);
          glVertex3f(cx[0],cy[0],cz[0]);
          glVertex3f(cx[1],cy[1],cz[1]);
          glVertex3f(cx[2],cy[2],cz[2]);
        glEnd();
      } 
      k = pt->nxt;
    }
    if ( transp ) {
      glDepthMask(GL_TRUE);
      glDisable(GL_BLEND);
    }
  }

  glEndList();
  return(dlist);
}


/* build list for capping map */
GLuint capTetraMap(pMesh mesh) {
  pScene     sc;
  pClip      clip;
  GLuint     dlist = 0;
  pTetra     pt;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     dd1[6],d,ax,ay,az,bx,by,bz;
  double     cx[4],cy[4],cz[4],cc;
  float      n[3],sol[4];
  int        m,k,k1,k2,l,pos[6],neg[6],nbpos,nbneg,nbnul;
  static int tn[4] = {0,0,1,1};
  static int tp[4] = {0,1,1,0};
  triangle   t1,t2;

  /* default */
  if ( !mesh->ntet || !mesh->nbb )  return(0);
  if ( ddebug ) printf("create capping map list / TETRA\n");
  sc  = cv.scene[currentScene()];
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
  clip = sc->clip;
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  glBegin(GL_TRIANGLES);
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;
 
    while ( k != 0 ) {
      pt = &mesh->tetra[k];
      if ( !pt->v[0] || !pt->clip ) {
        k = pt->nxt;
        continue;
      }
      nbpos = nbneg = nbnul = 0;
      for (l=0; l<4; l++) {
        p0 = &mesh->point[pt->v[l]];
        if ( p0->clip == 2 )      pos[nbpos++] = l;
        else if (p0->clip == 1 )  neg[nbneg++] = l;
        else                      nbnul++;
        dd1[l] = p0->c[0]*clip->eqn[0] + p0->c[1]*clip->eqn[1] \
               + p0->c[2]*clip->eqn[2] + clip->eqn[3];
      }

      if ( nbneg == 2 && nbpos == 2 ) {
        /* display quadrilateral */
        for (l=0; l<4; l++) {
          k1 = neg[tn[l]];
          k2 = pos[tp[l]];
          p0 = &mesh->point[pt->v[k1]];
          p1 = &mesh->point[pt->v[k2]];
          cc = 1.0f;
          if ( dd1[k2]-dd1[k1] != 0.0f )
            cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
          cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
          cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
          cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
          if ( mesh->typage == 2 ) {
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            sol[l] = ps0->bb+cc*(ps1->bb-ps0->bb);
          }
          else {
            ps0 = &mesh->sol[k];
            sol[l] = ps0->bb;
          }
        }
        
        /* compute face normal */
        ax = cx[1]-cx[0]; ay = cy[1]-cy[0]; az = cz[1]-cz[0];
        bx = cx[2]-cx[0]; by = cy[2]-cy[0]; bz = cz[2]-cz[0];
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

        /* store triangles */
        t1.a[0] = t2.a[0] = cx[0];
        t1.a[1] = t2.a[1] = cy[0];
        t1.a[2] = t2.a[2] = cz[0];

        t1.b[0] = cx[1];
        t1.b[1] = cy[1];
        t1.b[2] = cz[1];

        t1.c[0] = t2.b[0] = cx[2]; 
        t1.c[1] = t2.b[1] = cy[2]; 
        t1.c[2] = t2.b[2] = cz[2]; 

        t2.c[0] = cx[3]; 
        t2.c[1] = cy[3]; 
        t2.c[2] = cz[3]; 
        
        /* store normals */
        memcpy(t1.na,n,3*sizeof(float));
	memcpy(t1.nb,n,3*sizeof(float));
	memcpy(t1.nc,n,3*sizeof(float));
	memcpy(t2.na,n,3*sizeof(float));
	memcpy(t2.nb,n,3*sizeof(float));
	memcpy(t2.nc,n,3*sizeof(float));

        /* store solutions */
        t1.va = t2.va = sol[0];
	t1.vb = sol[1];
	t1.vc = t2.vb = sol[2];
	t2.vc = sol[3];
        /* color interpolation */
        cutTriangle(sc,t1);
        cutTriangle(sc,t2);
      }

      else {
        /* display triangle */
        for (l=0; l<3; l++) {
          k1 = nbneg == 3 ? neg[l] : pos[l];
          k2 = nbneg == 3 ? pos[0] : neg[0];
          p0 = &mesh->point[pt->v[k1]];
          p1 = &mesh->point[pt->v[k2]];
          cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
          cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
          cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
          cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
          if ( mesh->typage == 2 ) {
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            sol[l] = ps0->bb+cc*(ps1->bb-ps0->bb);
          }
          else {
            ps0 = &mesh->sol[k];
            sol[l] = ps0->bb;
          }
        }

        /* compute face normal */
        ax = cx[1]-cx[0]; ay = cy[1]-cy[0]; az = cz[1]-cz[0];
        bx = cx[2]-cx[0]; by = cy[2]-cy[0]; bz = cz[2]-cz[0];
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

        /* store triangle */
        t1.a[0] = cx[0];
        t1.a[1] = cy[0];
        t1.a[2] = cz[0];

        t1.b[0] = cx[1];
        t1.b[1] = cy[1];
        t1.b[2] = cz[1];

        t1.c[0] = cx[2];
        t1.c[1] = cy[2]; 
        t1.c[2] = cz[2]; 

        /* store normal */
        memcpy(t1.na,n,3*sizeof(float));
        memcpy(t1.nb,n,3*sizeof(float));
        memcpy(t1.nc,n,3*sizeof(float));
 
        /* store solutions */
        t1.va = sol[0];
        t1.vb = sol[1];
        t1.vc = sol[2];
        
        /* color interpolation */
        cutTriangle(sc,t1);
      }     
      k = pt->nxt;
    }
  }
  
  glEnd();
  glEndList();
  return(dlist);
}


/* build list for capping isolines */
GLuint capTetraIso(pMesh mesh) {
  pScene     sc;
  pClip      clip;
  GLuint     dlist = 0;
  pTetra     pt;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     dd1[6];
  double     rgb[3],cx[4],cy[4],cz[4],ccx,ccy,ccz,cc;
  float      iso,kc,sol[4];
  int        i,m,k,k1,k2,l,l1,nc,pos[6],neg[6],nbpos,nbneg,nbnul,ncol;
  static int tn[4] = {0,0,1,1};
  static int tp[4] = {0,1,1,0};
  static double hsv[3]   = { 0.0f, 1.0f, 0.20f };
  static int idirt[5] = {0,1,2,0,1};
  static int idirq[6] = {0,1,2,3,0,1};

  /* default */
  if ( !mesh->ntet || !mesh->nbb || mesh->typage == 1 )  return(0);
  if ( ddebug ) printf("create capping iso list / TETRA\n");
  sc = cv.scene[currentScene()];
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
  clip  = sc->clip;
  if ( !(clip->active & C_ON) )  return(0);
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  glBegin(GL_LINES);
  ncol = NBCOL;
  for (i=0; i<=ncol*(MAXISO-1); i++) {
    if ( i < ncol*(MAXISO-1) ) {
      l   = i / ncol;
      kc  = (i % ncol) / (float)ncol;
      iso = sc->iso.val[l]*(1.0-kc)+sc->iso.val[l+1]*kc;
      hsv[0] = sc->iso.col[l]*(1.0-kc)+sc->iso.col[l+1]*kc;
    }
    else {
      iso    = sc->iso.val[MAXISO-1];
      hsv[0] = sc->iso.col[MAXISO-1];
    }
    
    hsvrgb(hsv,rgb);
    glColor3dv(rgb);

    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTets];
      if ( !k || pm->flag )  continue;
 
      while ( k != 0 ) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] || !pt->clip ) {
          k = pt->nxt;
          continue;
        }
        nbpos = nbneg = nbnul = 0;
        for (l=0; l<4; l++) {
          p0 = &mesh->point[pt->v[l]];
          if ( p0->clip == 2 )      pos[nbpos++] = l;
          else if (p0->clip == 1 )  neg[nbneg++] = l;
          else                      nbnul++;
          dd1[l] = p0->c[0]*clip->eqn[0] + p0->c[1]*clip->eqn[1] \
                 + p0->c[2]*clip->eqn[2] + clip->eqn[3];
        }

        if ( nbneg == 2 && nbpos == 2 ) {
          /* analyze quadrilateral */
          for (l=0; l<4; l++) {
            k1 = neg[tn[l]];
            k2 = pos[tp[l]];
            p0 = &mesh->point[pt->v[k1]];
            p1 = &mesh->point[pt->v[k2]];
            cc = 1.0f;
            if ( dd1[k2]-dd1[k1] != 0.0f )
              cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
            cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            sol[l] = ps0->bb+cc*(ps1->bb-ps0->bb);
          }
          for (l=0; l<4; l++) {
            l1  = idirq[l+1];
            if ( (sol[l] > iso && sol[l1] <= iso) ||
                 (sol[l] < iso && sol[l1] >= iso) ) {
              cc = 0.0f;
 
              if ( fabs(sol[l1]-sol[l]) > 0.0f )
                cc = (iso-sol[l]) / (sol[l1]-sol[l]);
              if ( cc == 0.0f || cc == 1.0f )  continue;
              ccx = cx[l]+cc*(cx[l1]-cx[l]);
              ccy = cy[l]+cc*(cy[l1]-cy[l]);
              ccz = cz[l]+cc*(cz[l1]-cz[l]);
              glVertex3f(ccx,ccy,ccz);
              nc++;
            }
            else if ( sol[l] == iso && sol[l1] == iso ) {
              nc = 2;
              glVertex3f(cx[l],cy[l],cz[l]);
              glVertex3f(cx[l1],cy[l1],cz[l1]);
              break;
            }
          }
        }
        else {
          /* analyze triangle */
          for (l=0; l<3; l++) {
            k1 = nbneg == 3 ? neg[l] : pos[l];
            k2 = nbneg == 3 ? pos[0] : neg[0];
            p0 = &mesh->point[pt->v[k1]];
            p1 = &mesh->point[pt->v[k2]];
            cc = fabs(dd1[k1] / (dd1[k2]-dd1[k1]));
            cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            sol[l] = ps0->bb+cc*(ps1->bb-ps0->bb);
          }
          
          nc = 0;
          ccx = ccy = ccz = 0.0f;
          for (l=0; l<3; l++) {
            l1  = idirt[l+1];
            if ( (sol[l] > iso && sol[l1] <= iso) ||
                 (sol[l] < iso && sol[l1] >= iso) ) {
              cc = 0.0f;
 
              if ( fabs(sol[l1]-sol[l]) > 0.0f )
                cc = (iso-sol[l]) / (sol[l1]-sol[l]);
              if ( cc == 0.0f || cc == 1.0f )  continue;
              ccx = cx[l]+cc*(cx[l1]-cx[l]);
              ccy = cy[l]+cc*(cy[l1]-cy[l]);
              ccz = cz[l]+cc*(cz[l1]-cz[l]);
              glVertex3f(ccx,ccy,ccz);
              nc++;
            }
            else if ( sol[l] == iso && sol[l1] == iso ) {
              nc = 2;
              glVertex3f(cx[l],cy[l],cz[l]);
              glVertex3f(cx[l1],cy[l1],cz[l1]);
              break;
            }
          }
          if ( nc > 0 && nc != 2 )
            glVertex3f(ccx,ccy,ccz);
        }
        k = pt->nxt;
      }
    }
  }
  
  glEnd();
  glEndList();
  return(dlist);
}
