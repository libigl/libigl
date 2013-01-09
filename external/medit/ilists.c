#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define BASETR   0.2
#define NBCOL    9


/* build lists for iso-surfaces */
GLuint listTriaIso(pScene sc,pMesh mesh) {
  GLuint     dlist = 0;
  pTriangle  pt;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     rgb[3];
  float      iso,cx,cy,cz,cc,kc;
  int        m,k,i,l,l1,nc,ncol;
  static double hsv[3]  = { 0.0, 1.0, 0.9 };
  static int    idir[5] = {0,1,2,0,1};

  /* default */
  if ( !mesh->nt || !mesh->nbb || mesh->typage == 1 )  return(0);
  if ( ddebug ) printf("create iso-values map list / TRIA\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* create display list */
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
      k  = pm->depmat[LTria];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tria[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        /* analyze edges */
        nc = 0;
        cx = cy = cz = 0.0;
        for (l=0; l<3; l++) {
          l1  = idir[l+1];
          p0  = &mesh->point[pt->v[l]];
          p1  = &mesh->point[pt->v[l1]];
          ps0 = &mesh->sol[pt->v[l]];
          ps1 = &mesh->sol[pt->v[l1]];
          if ( (ps0->bb > iso && ps1->bb <= iso) ||
               (ps0->bb < iso && ps1->bb >= iso) ) {
            cc = 0.0;
            if ( fabs(ps1->bb-ps0->bb) > 0.0 )
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
            if ( cc == 0.0 || cc == 1.0 )  continue;
            cx = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            nc++;
            if ( mesh->dim == 2 )
              glVertex2f(cx,cy);
            else {
              cz = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
              glVertex3f(cx,cy,cz);
            }
          }
          else if ( ps0->bb == iso && ps1->bb == iso ) {
            nc = 2;
            if ( mesh->dim == 2 ) {
              glVertex2f(p0->c[0],p0->c[1]);
              glVertex2f(p1->c[0],p1->c[1]);
              break;
            }
            else {
              glVertex3f(p0->c[0],p0->c[1],p0->c[2]);
              glVertex3f(p1->c[0],p1->c[1],p1->c[2]);
              break;
            }
          }
        }
        if ( nc > 0 && nc != 2 ) {
          if ( mesh->dim ==2 )  glVertex2f(cx,cy);
          else                  glVertex3f(cx,cy,cz);
        }
        k = pt->nxt;
      }
    }
  }
  glEnd();
  glEndList();
  return(dlist);
}

/* build lists for iso-surfaces */
GLuint listQuadIso(pScene sc,pMesh mesh) {
  GLuint     dlist = 0;
  pQuad      pq;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     rgb[3];
  float      iso,cx,cy,cz,cc,kc;
  int        m,k,i,l,l1,ncol;
  static double hsv[3]   = { 0.0f, 1.0f, 0.8f };
  static int idir[6] = {0,1,2,3,0,1};

  /* default */
  if ( !mesh->nq || !mesh->nbb || mesh->typage == 1 )  return(0);
  if ( ddebug ) printf("create iso-values map list / QUAD\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
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
      k  = pm->depmat[LQuad];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pq = &mesh->quad[k];
        if ( !pq->v[0] ) {
          k = pq->nxt;
          continue;
        }

        /* analyze edges */
        for (l=0; l<4; l++) {
          p0  = &mesh->point[pq->v[l]];
          ps0 = &mesh->sol[pq->v[l]];
          l1  = idir[l+1];
          p1  = &mesh->point[pq->v[l1]];
          ps1 = &mesh->sol[pq->v[l1]];
          if ( (ps0->bb > iso && ps1->bb <= iso) ||
               (ps0->bb < iso && ps1->bb >= iso) ) {
            cc = 0.0f;
            if ( fabs(ps1->bb-ps0->bb) > 0.0f )
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
            if ( cc == 0.0f || cc == 1.0f )  continue;
            cx = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            if ( mesh->dim == 2 )
              glVertex2f(cx,cy);
            else {
              cz = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
              glVertex3f(cx,cy,cz);
            }
          }
        }
        k = pq->nxt;
      }
    }
  }
  glEnd();
  glEndList();
  return(dlist);
}


/* build lists for iso-surfaces */
GLuint listTetraIso(pScene sc,pMesh mesh) {
  FILE      *outv,*outf;
  GLuint     dlist = 0;
  pTetra     pt;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     delta,rgb[4],d,ax,ay,az,bx,by,bz;
  float      iso,n[3],cx[4],cy[4],cz[4],cc;
  int        m,k,k1,k2,i,l,pos[4],neg[4],nbpos,nbneg,nbnul,nv,nf;
  static double hsv[3]   = { 0.0f, 1.0f, 0.80f };
  static int tn[4] = {0,0,1,1};
  static int tp[4] = {0,1,1,0};

  /* default */
  if ( !mesh->ntet || !mesh->nbb || mesh->typage == 1 )  return(0);
  if ( ddebug ) printf("create iso-values map list / TETRA\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  delta = sc->iso.val[MAXISO-1] - sc->iso.val[0];

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);

  if (ddebug) {
	outv = fopen("vertex.mesh","w");
	fprintf(outv,"MeshVersionFormatted 1\n Dimension\n 3\n\nVertices\n \n");
	outf = fopen("faces.mesh2","w");
	fprintf(outv,"Triangles\n \n");
  }
  nv = nf = 0;

  glBegin(GL_TRIANGLES);
  for (i=MAXISO-1; i>=0; i--) {
    iso = sc->iso.val[i];

    /* base color */
    /*hsv[0] = 240.0f*(1.0f - (iso-sc->iso.val[0])/delta);*/
    hsv[0] = sc->iso.col[i];
    hsvrgb(hsv,rgb);
    rgb[0] = MEDIT_MIN(1.0,rgb[0]+BASETR);
    rgb[1] = MEDIT_MIN(1.0,rgb[1]+BASETR);
    rgb[2] = MEDIT_MIN(1.0,rgb[2]+BASETR);
    rgb[3] = BASETR + (float)(i-1)/(float)MAXISO*(1.0-BASETR);
    /*rgb[3] = 0.5; */  
    glColor4dv(rgb);

    if ( i == MAXISO-1 )  iso -= 0.001*fabs(iso)/delta;
    else if ( i == 0 )    iso += 0.001*fabs(iso)/delta;

    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTets];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        /* analyze vertices */
        nbpos = nbneg = nbnul = 0;
        for (l=0; l<4; l++) {
          p0  = &mesh->point[pt->v[l]];
          ps0 = &mesh->sol[pt->v[l]];
          /*if ( ps0->bb < sc->iso.val[0] )  ps0->bb = sc->iso.val[0];*/
          
          if ( ps0->bb > iso )      pos[nbpos++] = l;
          else if ( ps0->bb < iso ) neg[nbneg++] = l;
          else                      nbnul++;
        }
        if ( nbneg == 4 || nbpos == 4 ) {
          k = pt->nxt;
          continue;
        }

        if ( nbneg == 2 && nbpos == 2 ) {
          for (l=0; l<4; l++) {
            k1  = neg[tn[l]];
            k2  = pos[tp[l]];
            p0  = &mesh->point[pt->v[k1]];
            p1  = &mesh->point[pt->v[k2]];
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            cc = 0.0f;
            if ( fabs(ps1->bb-ps0->bb) > 0.0f )
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
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
          glNormal3fv(n);
          glVertex3f(cx[0],cy[0],cz[0]);
          glVertex3f(cx[1],cy[1],cz[1]);
          glVertex3f(cx[2],cy[2],cz[2]);
          
	      glNormal3fv(n);
          glVertex3f(cx[0],cy[0],cz[0]);
          glVertex3f(cx[2],cy[2],cz[2]);
          glVertex3f(cx[3],cy[3],cz[3]);

          if ( ddebug ) {
            fprintf(outv,"%f %f %f 0\n",cx[0],cy[0],cz[0]);
            fprintf(outv,"%f %f %f 0\n",cx[1],cy[1],cz[1]);
            fprintf(outv,"%f %f %f 0\n",cx[2],cy[2],cz[2]);
            fprintf(outv,"%f %f %f 0\n",cx[3],cy[3],cz[3]);

            fprintf(outf,"%d %d %d 0\n",nv+1,nv+2,nv+3);
            fprintf(outf,"%d %d %d 0\n",nv+1,nv+3,nv+4);
          }
          nv+= 4;
          nf+= 2;
        }
        else if ( !nbnul ) {
          for (l=0; l<3; l++) {
            k1 = nbneg == 3 ? neg[l] : pos[l];
            k2 = nbneg == 3 ? pos[0] : neg[0];
            p0 = &mesh->point[pt->v[k1]];
            p1 = &mesh->point[pt->v[k2]];
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            cc = 0.0f;
            if ( fabs(ps1->bb-ps0->bb) > 0.0f ) 
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
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
          glNormal3fv(n);
          glVertex3f(cx[0],cy[0],cz[0]);
          glVertex3f(cx[1],cy[1],cz[1]);
          glVertex3f(cx[2],cy[2],cz[2]);

          if ( ddebug ) {
            fprintf(outv,"%f %f %f 0\n",cx[0],cy[0],cz[0]);
            fprintf(outv,"%f %f %f 0\n",cx[1],cy[1],cz[1]);
            fprintf(outv,"%f %f %f 0\n",cx[2],cy[2],cz[2]);
            fprintf(outf,"%d %d %d 0\n",nv+1,nv+2,nv+3);
          }
          nv += 3;
          nf += 1;
        }
        k = pt->nxt;
      }
    }
  }
  glEnd();
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);

  if ( ddebug ) {
    fclose(outv);
    fclose(outf);
  }
  printf("  Vertices %d   Triangles %d\n",nv,nf);
  glEndList();
  return(dlist);
}


int tetraIsoPOVray(pScene sc,pMesh mesh) {
  FILE      *isofil;
  pTetra     pt;
  pPoint     p0,p1;
  pMaterial  pm;
  pSolution  ps0,ps1;
  double     delta;
  float      iso,cx[4],cy[4],cz[4],cc;
  int        m,k,k1,k2,i,l,pos[4],neg[4],nbpos,nbneg,nbnul;
  char      *ptr,data[128];
  static int tn[4] = {0,0,1,1};
  static int tp[4] = {0,1,1,0};

  /* default */
  if ( !mesh->ntet || !mesh->nbb || mesh->typage == 1 )  return(0);
  if ( ddebug ) printf("create isosurfaces POVray\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  delta = sc->iso.val[MAXISO-1] - sc->iso.val[0];

  strcpy(data,mesh->name);
  ptr = strstr(data,".mesh");
  if ( ptr )  ptr = '\0';
  strcat(data,".pov"); 
  if ( ddebug )  fprintf(stdout,"  Writing POVRay file %s\n",data);
  isofil = fopen(data,"w");
  if ( !isofil )  return(0);

  for (i=MAXISO-1; i>=0; i--) {
    iso = sc->iso.val[i];

    if ( i == MAXISO-1 )  iso -= 0.001*fabs(iso)/delta;
    else if ( i == 0 )    iso += 0.001*fabs(iso)/delta;

    fprintf(isofil,"\n#declare isosurf%d = mesh {\n",i);

    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      k  = pm->depmat[LTets];
      if ( !k || pm->flag )  continue;

      while ( k != 0 ) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] ) {
          k = pt->nxt;
          continue;
        }

        /* analyze vertices */
        nbpos = nbneg = nbnul = 0;
        for (l=0; l<4; l++) {
          p0  = &mesh->point[pt->v[l]];
          ps0 = &mesh->sol[pt->v[l]];
          
          if ( ps0->bb > iso )      pos[nbpos++] = l;
          else if ( ps0->bb < iso ) neg[nbneg++] = l;
          else                      nbnul++;
        }
        if ( nbneg == 4 || nbpos == 4 ) {
          k = pt->nxt;
          continue;
        }

        if ( nbneg == 2 && nbpos == 2 ) {
          for (l=0; l<4; l++) {
            k1  = neg[tn[l]];
            k2  = pos[tp[l]];
            p0  = &mesh->point[pt->v[k1]];
            p1  = &mesh->point[pt->v[k2]];
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            cc = 0.0f;
            if ( fabs(ps1->bb-ps0->bb) > 0.0f )
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
            cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
          }

	  fprintf(isofil,"triangle {\n");
	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[0]+mesh->xtra,cy[0]+mesh->ytra,cz[0]+mesh->ztra);
 	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[1]+mesh->xtra,cy[1]+mesh->ytra,cz[1]+mesh->ztra);
	  fprintf(isofil,"  <%f,%f,%f>\n" ,
	  cx[2]+mesh->xtra,cy[2]+mesh->ytra,cz[2]+mesh->ztra);
          fprintf(isofil,"}\n");

	  fprintf(isofil,"triangle {\n");
	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[0]+mesh->xtra,cy[0]+mesh->ytra,cz[0]+mesh->ztra);
 	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[2]+mesh->xtra,cy[2]+mesh->ytra,cz[2]+mesh->ztra);
	  fprintf(isofil,"  <%f,%f,%f>\n" ,
	  cx[3]+mesh->xtra,cy[3]+mesh->ytra,cz[3]+mesh->ztra);
          fprintf(isofil,"}\n");
        }
        else if ( !nbnul ) {
          for (l=0; l<3; l++) {
            k1 = nbneg == 3 ? neg[l] : pos[l];
            k2 = nbneg == 3 ? pos[0] : neg[0];
            p0 = &mesh->point[pt->v[k1]];
            p1 = &mesh->point[pt->v[k2]];
            ps0 = &mesh->sol[pt->v[k1]];
            ps1 = &mesh->sol[pt->v[k2]];
            cc = 0.0f;
            if ( fabs(ps1->bb-ps0->bb) > 0.0f ) 
              cc = (iso-ps0->bb) / (ps1->bb-ps0->bb);
            cx[l] = p0->c[0]+cc*(p1->c[0]-p0->c[0]);
            cy[l] = p0->c[1]+cc*(p1->c[1]-p0->c[1]);
            cz[l] = p0->c[2]+cc*(p1->c[2]-p0->c[2]);
          }
	  fprintf(isofil,"triangle {\n");
	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[0]+mesh->xtra,cy[0]+mesh->ytra,cz[0]+mesh->ztra);
 	  fprintf(isofil,"  <%f,%f,%f>,\n",
	  cx[1]+mesh->xtra,cy[1]+mesh->ytra,cz[1]+mesh->ztra);
	  fprintf(isofil,"  <%f,%f,%f>\n",
	  cx[2]+mesh->xtra,cy[2]+mesh->ytra,cz[2]+mesh->ztra);
          fprintf(isofil,"}\n");
        }
        k = pt->nxt;
      }
    }
    fprintf(isofil,"}\n");
  }

  fclose(isofil);
  return(1);
}
