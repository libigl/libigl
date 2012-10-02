#include "medit.h"
#include "extern.h"
#include "sproto.h"

typedef struct saddle {
  double r1,r2;
  float  pt[3];
  float  vp1[2],vp2[2];
  int    k;
} Saddle;


#define EPSD  1.e-14
#define EPS3  1.e-03
#define MAXPS 100
#define HSIZ  0.001

static int idir[5] = {0,1,2,0,1};


int closedBall(pMesh mesh,int depart,ubyte i) {
  pTriangle  pt;
  int        adj,iadr;
  ubyte      voy;

  pt   = &mesh->tria[depart];
  voy  = idir[i+1];
  iadr = 3*(depart-1)+1;
  adj  = mesh->adja[iadr+voy];

  /* search triangles in ball */
  while ( adj && adj != depart ) {
    voy  = mesh->voy[iadr+voy];
    iadr = 3*(adj-1)+1;
    voy  = idir[voy+2];
    adj  = mesh->adja[iadr+voy];
  }

  return( adj == depart );
}


GLuint listCritPoint(pScene sc,pMesh mesh) {
  pTriangle   pt;
  pPoint      p0,p1,p2;
  pSolution   s0,s1,s2;
  pMaterial   pm;
  Saddle      sad[MAXPS];
  GLuint      dlist;
  double      aire,ux,uy,vx,vy,dd,cb0[3],cb1[3],cb2[3],vv[3][2],bc[3];
  double      rgb[3],a0,a1,delta,rr1,rr2,aa,dmin;
  float       p[3];
  int        *adj,iadr,i,i1,i2,k,m,ncp,ps,ifilt;
  ubyte       typ,tag;
  static double hsv[3] = {0.0f, 1.0f, 0.80f};
  time_t      t;

  if ( !mesh->nbb || mesh->nfield != mesh->dim ) return(0);
  if ( mesh->nt && !hashTria(mesh) )  return(0);
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug )  printf("find critical points\n");

  /* build list */
  typ = 0;
  ncp = 0;
  ps  = 0;
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
  glPointSize(4.0);

  dmin   = sc->dmax * EPS;
  dmin  *= dmin;
  ifilt  = 0;
  hsv[0] = sc->iso.col[0];
  hsvrgb(hsv,rgb);

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

      p0 = &mesh->point[pt->v[0]];
      p1 = &mesh->point[pt->v[1]];
      p2 = &mesh->point[pt->v[2]];
      s0 = &mesh->sol[pt->v[0]];
      s1 = &mesh->sol[pt->v[1]];
      s2 = &mesh->sol[pt->v[2]];
      
      ux = p1->c[0] - p0->c[0];
      uy = p1->c[1] - p0->c[1];
      vx = p2->c[0] - p0->c[0];
      vy = p2->c[1] - p0->c[1];
      
      aire = ux*vy - uy*vx;
      if ( aire == 0.0 ) {
        k = pt->nxt;
        continue;
      }
      else if ( aire < 0.0 ) {
        p1 = &mesh->point[pt->v[2]];
        p2 = &mesh->point[pt->v[1]];
        s1 = &mesh->sol[pt->v[2]];
        s2 = &mesh->sol[pt->v[1]];
        aire = -aire;
      }

      /* coef des barycentriques */
      aire   = 1.0 / aire;
      cb0[0] =   p1->c[1] - p2->c[1];
      cb0[1] = -(p1->c[0] - p2->c[0]);
      cb0[2] =   p1->c[0]*p2->c[1] - p1->c[1]*p2->c[0];

      cb1[0] =   p2->c[1] - p0->c[1];
      cb1[1] = -(p2->c[0] - p0->c[0]);
      cb1[2] =   p2->c[0]*p0->c[1] - p2->c[1]*p0->c[0];

      cb2[0] =   p0->c[1] - p1->c[1];
      cb2[1] = -(p0->c[0] - p1->c[0]);
      cb2[2] =   p0->c[0]*p1->c[1] - p0->c[1]*p1->c[0];
      for (i=0; i<3; i++) {
        vv[i][0] = aire * (cb0[i]*s0->m[0] + cb1[i]*s1->m[0] + cb2[i]*s2->m[0]);
        vv[i][1] = aire * (cb0[i]*s0->m[1] + cb1[i]*s1->m[1] + cb2[i]*s2->m[1]);
      }

      dd = vv[0][0]*vv[1][1] - vv[0][1]*vv[1][0];
      if ( fabs(dd) < EPSD ) {
        k = pt->nxt;
        continue;
      }

      dd   = 1.0 / dd;
      p[0] = dd * (vv[1][0]*vv[2][1] - vv[2][0]*vv[1][1]);
      p[1] = dd * (vv[0][1]*vv[2][0] - vv[0][0]*vv[2][1]);
      p[2] = 0.0;
  
      if ( p[0] < mesh->xmin-mesh->xtra || p[0] > mesh->xmax-mesh->xtra
        || p[1] < mesh->ymin-mesh->ytra || p[1] > mesh->ymax-mesh->ytra ) {
        k = pt->nxt;
        continue;
      }
      else if ( !inTria(mesh,k,p,bc) ) {
        k = pt->nxt;
        continue;
      }

      /* filtering boundary points */
      tag  = 0;
      for (i=0; i<3; i++)  tag |= (bc[i] < EPS3) << i;
      if ( tag ) {
        iadr = 3*(k-1)+1;
        adj  = &mesh->adja[iadr];
        ifilt ++;
	switch (tag) {
          case 1:
	    if ( !adj[0] ) {
	      k = pt->nxt;
	      continue;
	    }
	    break;
	  case 2:
	    if ( !adj[1] ) {
	      k = pt->nxt;
	      continue;
	    }
	    break;
	  case 4:
	    if ( !adj[2] ) {
	      k = pt->nxt;
	      continue;
	    }
	    break;
	  
	  case 3:
            if ( !closedBall(mesh,k,2) ) {
	      k = pt->nxt;
	      continue;
	    }
	    break;
	  case 5:
            if ( !closedBall(mesh,k,1) ) {
	      k = pt->nxt;
	      continue;
	    }	  
	    break;  
	  case 6:
            if ( !closedBall(mesh,k,0) ) {
	      k = pt->nxt;
	      continue;
	    }
	    break;
	}
      }

      /* eigenvalues of jacobian */
      a1 = -(vv[0][0] + vv[1][1]);
      a0 = vv[0][0] *vv[1][1] - vv[0][1]*vv[1][0];
      delta = a1*a1 - 4*a0;
      i1 = i2 = 0;
      if ( delta >= 0.0 ) {
        delta = sqrt(delta);
        rr1 = 0.5 * (-a1+delta);
        rr2 = 0.5 * (-a1-delta);
      }
      else {
        delta = sqrt(fabs(delta));
        rr1 = rr2 = -0.5 * a1;
        i1  = i2  =  0.5 * delta;
      }

      /* classification */
      if ( i1 && i2 ) {
        glColor3f(0.0,1.0,0.5);
        if ( rr1 == 0.0f && rr2 == 0.0f )
          output2(p[0],p[1],"Cp");   /* center */
        else if ( rr1 > 0.0f && rr2 > 0.0f )
          output2(p[0],p[1],"Rf");   /* repelling focus */
        else if ( rr1 < 0.0f && rr2 < 0.0f )
          output2(p[0],p[1],"Af");   /* attracting focus */
      }
      else if ( !i1 && !i2 ) {
        glColor3f(1.0,0.5,0.0);
        if ( rr1 > 0.0f && rr2 > 0.0f )
          output2(p[0],p[1],"Rn");   /* repelling node */
        else if ( rr1 < 0.0f && rr2 < 0.0f )
          output2(p[0],p[1],"An");   /* attracting node */
        else if ( rr1*rr2 < 0.0f ) {
          output2(p[0],p[1],"Sp");   /* Saddle point */
          if ( ddebug ) 
            printf("  saddle point %f %f\n",p[0]+mesh->xtra,p[1]+mesh->ytra);
          if ( ps < MAXPS-5 ) {
            ++ps;
            sad[ps].pt[0] = p[0];
            sad[ps].pt[1] = p[1];
            sad[ps].pt[2] = 0.0f;

            /* eigenvalues */
            sad[ps].r1 = rr1;
            sad[ps].r2 = rr2;
            
            /* eigenvectors */
            aa = vv[0][0]*vv[0][0];
            dd = vv[1][1]*vv[1][1];
            delta = sqrt(aa-2.0*vv[0][0]*vv[1][1]+dd+4.0*vv[0][1]*vv[1][0]);
            sad[ps].vp1[0] = -0.5*(-vv[0][0]+vv[1][1]-delta);
            sad[ps].vp1[1] = vv[0][1];

            sad[ps].vp2[0] = -0.5*(-vv[0][0]+vv[1][1]+delta);
            sad[ps].vp2[1] = vv[0][1];
            
            sad[ps].k = k;
          }
        }
      }

      /* point color */
      glBegin(GL_POINTS);
        glColor3dv(rgb);
        glVertex2fv(p);
      glEnd();
      pt->cpt--;
      ++ncp;

      k = pt->nxt;
    }
  }
  glPointSize(1.0);
  glEndList();

  if ( ncp )
    fprintf(stdout,"   %d critical points identified (%d filtered)\n",ncp,ifilt);

return(dlist);
  
  if ( ps ) {
    fprintf(stdout," Building streamline(s)");
    fflush(stdout);
    t = clock();

    if ( !sc->slist ) {
      sc->stream = createStream(sc,mesh);
      if ( !sc->stream )  return(dlist);
    }
    for (k=1; k<=ps; k++) {
      if ( ddebug )  printf("  eigenv1 %f %f\n",sad[k].vp1[0],sad[k].vp1[1]);      
      
      listSaddleStream(sc,mesh,sad[k].k,sad[k].pt,sad[k].vp1,sad[k].r1);
      sad[k].vp1[0] = -sad[k].vp1[0];
      sad[k].vp1[1] = -sad[k].vp1[1];
      listSaddleStream(sc,mesh,sad[k].k,sad[k].pt,sad[k].vp1,sad[k].r1);
      
      if ( ddebug )  printf("  eigenv2 %f %f\n",sad[k].vp2[0],sad[k].vp2[1]);
      listSaddleStream(sc,mesh,sad[k].k,sad[k].pt,sad[k].vp2,sad[k].r2);
      sad[k].vp2[0] = -sad[k].vp2[0];
      sad[k].vp2[1] = -sad[k].vp2[1];
      listSaddleStream(sc,mesh,sad[k].k,sad[k].pt,sad[k].vp2,sad[k].r2);
    }
    sc->isotyp |= S_STREAML;
    fprintf(stdout,": %d lines",ps);
    t = clock() - t;
    fprintf(stdout," %6.2f sec.\n",t/(float)CLOCKS_PER_SEC);
  }
  
  return(dlist);
}
