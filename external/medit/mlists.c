#include "medit.h"
#include "extern.h"
#include "sproto.h"


static int ct[4][3] = { {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3} };
static int ch[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, 
            {1,2,6,5}, {2,3,7,6}, {0,3,7,4} };

#ifdef IGL
#include <igl/jet.h>
#include <igl/rgb_to_hsv.h>
void IGLParams::rgb(double x, double * rgb)
{

  double hsv[3];
  igl::rgb_to_hsv(this->easter_red,hsv);
  //hsv[0] = 0;
  //hsv[1] = this->easter_s;
  //hsv[2] = this->easter_v;
  switch(color_map)
  {
    case COLOR_MAP_JET:
      igl::jet(x,rgb);
      return;
    case COLOR_MAP_EASTER:
    {
      hsv[0] = x*360;
      hsvrgb(hsv,rgb);
      return;
    }
    case COLOR_MAP_WINDING_THEN_EASTER:
    {
      if(x == 0)
      {
        rgb[0] = 1;
        rgb[1] = 0;
        rgb[2] = 0;
      }else
      {
        hsv[0] = x*360;
        hsvrgb(hsv,rgb);
      }
      return;
    }
    case COLOR_MAP_DEFAULT:
    default:
    {
      double def[3] = { 0.0, 1.0, 0.80 };
      def[0] = 240.0-x*240.0;
      hsvrgb(def,rgb);
      return;
    }
  }
}
#endif

#ifdef IGL
double filter(pScene sc, double ss)
{
  double s = ss;
  if(sc->igl_params->fade_flip)
  {
    s = 1.0-s;
  }
  const double & max_s = sc->igl_params->fade_max_s;
  const double & min_s = sc->igl_params->fade_min_s;
  const double & max_v = sc->igl_params->fade_max_v;
  const double & min_v = sc->igl_params->fade_min_v;
  if(s>max_s)
  {
    return max_v;
  }else if(s<=min_s)
  {
    return min_v;
  }else
  {
    double f = (s-min_s)/(max_s-min_s);
    f = -2*f*f*f+3*f*f;
    return f*(max_v-min_v)+min_v;
  }
}
#endif

/* recursively subdivide a triangle */
void cutTriangle(pScene sc,triangle t) {
  triangle      t1,t2;
  double        kc,x,dd,rgb[4],maxe;
  int           i,ia,ib,ic,iedge;
  // Alec: hard coded pallete
  static double hsv[3] = { 0.0, 1.0, 0.80 };
#ifdef IGL
  rgb[3] = sc->igl_params->alpha_holder;
#endif

  /* analyze triangle edges */
  if ( t.va < sc->iso.val[0] )  
    t.va = sc->iso.val[0];
  else if ( t.va > sc->iso.val[MAXISO-1] )
    t.va = sc->iso.val[MAXISO-1];
  if ( t.vb < sc->iso.val[0] )  
    t.vb = sc->iso.val[0];
  else if ( t.vb > sc->iso.val[MAXISO-1] )
    t.vb = sc->iso.val[MAXISO-1];
  if ( t.vc < sc->iso.val[0] )  
    t.vc = sc->iso.val[0];
  else if ( t.vc > sc->iso.val[MAXISO-1] )
    t.vc = sc->iso.val[MAXISO-1];

  for (ia=0; ia<MAXISO-1; ia++)
    if ( t.va < sc->iso.val[ia] )  break;

  for (ib=0; ib<MAXISO-1; ib++)
    if ( t.vb < sc->iso.val[ib] )  break;

  for (ic=0; ic<MAXISO-1; ic++)
    if ( t.vc < sc->iso.val[ic] )  break;

  /* search longest edge */
  maxe  = fabs(t.va-t.vb);
  iedge = 1;
  if ( maxe < fabs(t.vb-t.vc) ) {
    maxe  = fabs(t.vb-t.vc);
    iedge = 2;
  }
  if ( maxe < fabs(t.va-t.vc) ) {
    maxe = fabs(t.va-t.vc);
    iedge = 3;
  }

  /* split longest edge */
  if ( maxe > 0.0 ) {
    switch( iedge ) {
    case 1:  /* edge a-b */
      x  = ia < ib ? sc->iso.val[ia] : sc->iso.val[ib];
      dd = (x-t.va) / (t.vb-t.va);
      if ( dd > 0.001 && dd < 0.999 ) {
        memcpy(&t1,&t,sizeof(struct triangle));
        memcpy(&t2,&t,sizeof(struct triangle));
        for (i=0; i<3; i++) {
          t1.b[i]  = t2.a[i]  = t.a[i]  + dd*(t.b[i] -t.a[i]);
          t1.nb[i] = t2.na[i] = t.na[i] + dd*(t.nb[i]-t.na[i]);
        }
        t1.vb = t2.va = x;
        cutTriangle(sc,t1);
        cutTriangle(sc,t2);
        return;
      }
      break;

    case 2:  /* edge b-c */
      x  = ib < ic ? sc->iso.val[ib] : sc->iso.val[ic];
      dd = (x-t.vb) / (t.vc-t.vb);
      if ( dd > 0.001f && dd < 0.999f ) {
        memcpy(&t1,&t,sizeof(struct triangle));
        memcpy(&t2,&t,sizeof(struct triangle));
        for (i=0; i<3; i++) {
          t1.c[i]  = t2.b[i]  = t.b[i]  + dd*(t.c[i] -t.b[i]);
          t1.nc[i] = t2.nb[i] = t.nb[i] + dd*(t.nc[i]-t.nb[i]);
        }
        t1.vc = t2.vb = x;
        cutTriangle(sc,t1);
        cutTriangle(sc,t2);
        return;
      }
      break;
    case 3:  /* edge c-a */
      x  = ia < ic ? sc->iso.val[ia] : sc->iso.val[ic];
      dd = (x-t.va) / (t.vc-t.va);
      if ( dd > 0.001f && dd < 0.999f ) {
        memcpy(&t1,&t,sizeof(struct triangle));
        memcpy(&t2,&t,sizeof(struct triangle));
        for (i=0; i<3; i++) {
          t1.c[i]  = t2.a[i]  = t.a[i]  + dd*(t.c[i] -t.a[i]);
          t1.nc[i] = t2.na[i] = t.na[i] + dd*(t.nc[i]-t.na[i]);
        }
        t1.vc = t2.va = x;
        cutTriangle(sc,t1);
        cutTriangle(sc,t2);
        return;
      }
      break;
    }
  }

  /* draw triangle */
  if ( t.va < sc->iso.val[0] ) 
    t.va = sc->iso.val[0];  
  else if ( t.va > sc->iso.val[MAXISO-1] )
    t.va = sc->iso.val[MAXISO-1];
  kc = (t.va-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
#ifdef IGL
  {
    double s = 1.0-(sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc)/240.0;
    // Alpha mess
    //rgb[3] = filter(sc,s)*sc->igl_params->alpha_holder;
    sc->igl_params->rgb(s,rgb);
  }
#else
  hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
  hsvrgb(hsv,rgb);
#endif

  glColor4dv(rgb);
  glNormal3fv(t.na);
  glVertex3fv(t.a);

  if ( t.vb < sc->iso.val[0] ) 
    t.vb = sc->iso.val[0];
  else if ( t.vb > sc->iso.val[MAXISO-1] )
    t.vb = sc->iso.val[MAXISO-1];
  kc = (t.vb-sc->iso.val[ib-1]) / (sc->iso.val[ib] - sc->iso.val[ib-1]);
#ifdef IGL
  {
    double s  = 1.0-(sc->iso.col[ib-1]*(1.0-kc)+sc->iso.col[ib]*kc)/240.0;
    // Alpha mess
    //rgb[3] = filter(sc,s)*sc->igl_params->alpha_holder;
    sc->igl_params->rgb( s,rgb);
  }
#else
  hsv[0] = sc->iso.col[ib-1]*(1.0-kc)+sc->iso.col[ib]*kc;
  hsvrgb(hsv,rgb);
#endif

  glColor4dv(rgb);
  glNormal3fv(t.nb);
  glVertex3fv(t.b);

  if ( t.vc < sc->iso.val[0] ) 
    t.vc = sc->iso.val[0];
  else if ( t.vc > sc->iso.val[MAXISO-1] )
    t.vc = sc->iso.val[MAXISO-1];
  kc = (t.vc-sc->iso.val[ic-1]) / (sc->iso.val[ic] - sc->iso.val[ic-1]);
#ifdef IGL
  {
    double s = 1.0-(sc->iso.col[ic-1]*(1.0-kc)+sc->iso.col[ic]*kc)/240.0;
    // Alpha mess
    // rgb[3] = filter(sc,s)*sc->igl_params->alpha_holder;
    sc->igl_params->rgb( s,rgb);
  }
#else
  hsv[0] = sc->iso.col[ic-1]*(1.0-kc)+sc->iso.col[ic]*kc;
  hsvrgb(hsv,rgb);
#endif

  glColor4dv(rgb);
  glNormal3fv(t.nc);
  glVertex3fv(t.c);
}


/* metric map: use linear interpolation on values
   rather than color interpolation ! */
GLuint listTriaMap(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pTriangle  pt;
  pPoint     p0,p1,p2;
  pSolution  ps0,ps1,ps2;
  GLint      dlist;
  double     ax,ay,az,bx,by,bz,dd;
  float      cx,cy,cz,n[3];
  int        k,m,is0,is1,is2;
  triangle   t;

  /* default */
  if ( !mesh->nt )  return(0);
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("create display list map / TRIA\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
#ifdef IGL
  bool transp = sc->material->dif[3] < 0.999;
  int old_depth_func =0;
  glGetIntegerv(GL_DEPTH_FUNC,&old_depth_func);
  if ( transp )
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_ALWAYS);
    sc->igl_params->alpha_holder = sc->material->dif[3];
  }else
  {
    sc->igl_params->alpha_holder = 1.0;
  }
#endif

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTria];
    if ( !k || pm->flag )  continue;

    if ( sc->type & S_FLAT ) {
      glBegin(GL_TRIANGLES);
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
        if ( dd > 0.0f ) {
          dd = 1.0f / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }
        memcpy(t.na,n,3*sizeof(float));
        memcpy(t.nb,n,3*sizeof(float));
        memcpy(t.nc,n,3*sizeof(float));

        if ( sc->shrink < 1.0 ) {
          cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.0;
          cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.0;
          cz = (p0->c[2] + p1->c[2] + p2->c[2]) / 3.0;
          t.a[0] = sc->shrink*(p0->c[0]-cx)+cx;
          t.a[1] = sc->shrink*(p0->c[1]-cy)+cy;
          t.a[2] = sc->shrink*(p0->c[2]-cz)+cz;
          t.b[0] = sc->shrink*(p1->c[0]-cx)+cx;
          t.b[1] = sc->shrink*(p1->c[1]-cy)+cy;
          t.b[2] = sc->shrink*(p1->c[2]-cz)+cz;
          t.c[0] = sc->shrink*(p2->c[0]-cx)+cx;
          t.c[1] = sc->shrink*(p2->c[1]-cy)+cy;
          t.c[2] = sc->shrink*(p2->c[2]-cz)+cz;
        }
        else {
          t.a[0] = p0->c[0];
          t.a[1] = p0->c[1];
          t.a[2] = p0->c[2];
          t.b[0] = p1->c[0];
          t.b[1] = p1->c[1];
          t.b[2] = p1->c[2];
          t.c[0] = p2->c[0];
          t.c[1] = p2->c[1];
          t.c[2] = p2->c[2];
        }
        if ( mesh->typage == 2 ) {
          ps0  = &mesh->sol[pt->v[0]];
          ps1  = &mesh->sol[pt->v[1]];
          ps2  = &mesh->sol[pt->v[2]];
          t.va = ps0->bb;
          t.vb = ps1->bb;
          t.vc = ps2->bb;
        }
        else {
          ps0 = &mesh->sol[k];
          t.va = t.vb = t.vc = ps0->bb;
        }
#ifdef IGL
        if(pt->ref == 1)
        {
          // Self-intersection
          float red[4] = {1.0, 0.0, 0.0, 1.0};
          red[3] = sc->material->dif[3];
          glColor4fv(red);
          glNormal3fv(n);
          glVertex3fv(t.a);
          glColor4fv(red);
          glNormal3fv(n);
          glVertex3fv(t.b);
          glColor4fv(red);
          glNormal3fv(n);
          glVertex3fv(t.c);
        }else{
#endif
        cutTriangle(sc,t);
#ifdef IGL
        }
#endif
        k = pt->nxt;
      }
      glEnd();
    }
    else {
      glBegin(GL_TRIANGLES);
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
        if ( dd > 0.0f ) {
      dd = 1.0f / sqrt(dd);
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
          t.a[0] = sc->shrink*(p0->c[0]-cx)+cx;
          t.a[1] = sc->shrink*(p0->c[1]-cy)+cy;
          t.a[2] = sc->shrink*(p0->c[2]-cz)+cz;
          t.b[0] = sc->shrink*(p1->c[0]-cx)+cx;
          t.b[1] = sc->shrink*(p1->c[1]-cy)+cy;
          t.b[2] = sc->shrink*(p1->c[2]-cz)+cz;
          t.c[0] = sc->shrink*(p2->c[0]-cx)+cx;
          t.c[1] = sc->shrink*(p2->c[1]-cy)+cy;
          t.c[2] = sc->shrink*(p2->c[2]-cz)+cz;
        }
        else {
          t.a[0] = p0->c[0];
          t.a[1] = p0->c[1];
          t.a[2] = p0->c[2];
          t.b[0] = p1->c[0];
          t.b[1] = p1->c[1];
          t.b[2] = p1->c[2];
          t.c[0] = p2->c[0];
          t.c[1] = p2->c[1];
          t.c[2] = p2->c[2];
        }
        if ( !is0 )
          memcpy(t.na,n,3*sizeof(float));
        else {
          t.na[0] = mesh->extra->n[3*(is0-1)+1];
          t.na[1] = mesh->extra->n[3*(is0-1)+2];
          t.na[2] = mesh->extra->n[3*(is0-1)+3];
        }
        if ( !is1 )
          memcpy(t.nb,n,3*sizeof(float));
        else {
          t.nb[0] = mesh->extra->n[3*(is1-1)+1];
          t.nb[1] = mesh->extra->n[3*(is1-1)+2];
          t.nb[2] = mesh->extra->n[3*(is1-1)+3];
        }
        if ( !is2 )
          memcpy(t.nc,n,3*sizeof(float));
        else {
          t.nc[0] = mesh->extra->n[3*(is2-1)+1];
          t.nc[1] = mesh->extra->n[3*(is2-1)+2];
          t.nc[2] = mesh->extra->n[3*(is2-1)+3];
        }
        if ( mesh->typage == 2 ) {
      ps0  = &mesh->sol[pt->v[0]];
      ps1  = &mesh->sol[pt->v[1]];
      ps2  = &mesh->sol[pt->v[2]];
      t.va = ps0->bb;
      t.vb = ps1->bb;
      t.vc = ps2->bb;
        }
        else {
      ps0 = &mesh->sol[k];
      t.va = t.vb = t.vc = ps0->bb;
        }
        cutTriangle(sc,t);
        k = pt->nxt;
      }
      glEnd();
    }
  }
#ifdef IGL
  if(transp)
  {
    glDepthFunc(old_depth_func);
    glDisable(GL_BLEND);
  }
#endif

  glEndList();
  return(dlist);
}


/* build list of quadrilaterals */
GLuint listQuadMap(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pQuad      pq;
  pPoint     p0,p1,p2,p3;
  pSolution  ps0,ps1,ps2,ps3;
  GLint      dlist = 0;
  double     ax,ay,az,bx,by,bz,dd;
  float      cx,cy,cz,n[3];
  int        k,m,is0,is1,is2,is3;
  triangle   t1,t2;

  /* default */
  if ( !mesh->nq ) return(0);
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("create display list map / QUADS\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LQuad];
    if ( !k || pm->flag )  continue;

    glBegin(GL_TRIANGLES);
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
      if ( dd > 0.0f ) {
        dd = 1.0f / sqrt(dd);
        n[0] *= dd;
        n[1] *= dd;
        n[2] *= dd;
      }

      if ( sc->shrink < 1.0 ) {
        cx = 0.25 * (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]);
        cy = 0.25 * (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]);
        cz = 0.25 * (p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]);
        t1.a[0] = t2.a[0] = sc->shrink*(p0->c[0]-cx)+cx;
        t1.a[1] = t2.a[1] = sc->shrink*(p0->c[1]-cy)+cy;
        t1.a[2] = t2.a[2] = sc->shrink*(p0->c[2]-cz)+cz;

        t1.b[0] = sc->shrink*(p1->c[0]-cx)+cx;
        t1.b[1] = sc->shrink*(p1->c[1]-cy)+cy;
        t1.b[2] = sc->shrink*(p1->c[2]-cz)+cz;

        t1.c[0] = t2.b[0] = sc->shrink*(p2->c[0]-cx)+cx; 
        t1.c[1] = t2.b[1] = sc->shrink*(p2->c[1]-cy)+cy; 
        t1.c[2] = t2.b[2] = sc->shrink*(p2->c[2]-cz)+cz; 

        t2.c[0] = sc->shrink*(p3->c[0]-cx)+cx; 
        t2.c[1] = sc->shrink*(p3->c[1]-cy)+cy; 
        t2.c[2] = sc->shrink*(p3->c[2]-cz)+cz; 
      }
      else {
        t1.a[0] = t2.a[0] = p0->c[0];
        t1.a[1] = t2.a[1] = p0->c[1];
        t1.a[2] = t2.a[2] = p0->c[2];

        t1.b[0] = p1->c[0];
        t1.b[1] = p1->c[1];
        t1.b[2] = p1->c[2];

        t1.c[0] = t2.b[0] = p2->c[0]; 
        t1.c[1] = t2.b[1] = p2->c[1]; 
        t1.c[2] = t2.b[2] = p2->c[2]; 

        t2.c[0] = p3->c[0]; 
        t2.c[1] = p3->c[1]; 
        t2.c[2] = p3->c[2]; 
      }
      if ( sc->type & S_FLAT ) {
        memcpy(t1.na,n,3*sizeof(float));
        memcpy(t1.nb,n,3*sizeof(float));
        memcpy(t1.nc,n,3*sizeof(float));
        memcpy(t2.na,n,3*sizeof(float));
        memcpy(t2.nb,n,3*sizeof(float));
        memcpy(t2.nc,n,3*sizeof(float));
      }
      else {
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
        
        if ( !is0 )
          memcpy(t1.na,n,3*sizeof(float));
        else {
          t1.na[0] = t2.na[0] = mesh->extra->n[3*(is0-1)+1];
          t1.na[1] = t2.na[1] = mesh->extra->n[3*(is0-1)+2];
          t1.na[2] = t2.na[2] = mesh->extra->n[3*(is0-1)+3];
        }
        if ( !is1 )
          memcpy(t1.nb,n,3*sizeof(float));
        else {
          t1.nb[0] = mesh->extra->n[3*(is1-1)+1];
          t1.nb[1] = mesh->extra->n[3*(is1-1)+2];
          t1.nb[2] = mesh->extra->n[3*(is1-1)+3];
        }
        if ( !is2 )
          memcpy(t1.nc,n,3*sizeof(float));
        else {
          t1.nc[0] = t2.nb[0] = mesh->extra->n[3*(is2-1)+1];
          t1.nc[1] = t2.nb[1] = mesh->extra->n[3*(is2-1)+2];
          t1.nc[2] = t2.nb[2] = mesh->extra->n[3*(is2-1)+3];
        }
        if ( !is3 )
          memcpy(t1.nc,n,3*sizeof(float));
        else {
          t2.nc[0] = mesh->extra->n[3*(is3-1)+1];
          t2.nc[1] = mesh->extra->n[3*(is3-1)+2];
          t2.nc[2] = mesh->extra->n[3*(is3-1)+3];
        }
      }

      if ( mesh->typage == 2 ) {
        /* solutions at vertices */
        ps0 = &mesh->sol[pq->v[0]];
        ps1 = &mesh->sol[pq->v[1]];
        ps2 = &mesh->sol[pq->v[2]];
        ps3 = &mesh->sol[pq->v[3]];
        t1.va = t2.va = ps0->bb;
        t1.vb = ps1->bb;
        t1.vc = t2.vb = ps2->bb;
        t2.vc = ps3->bb;
      }
      else {
        /* solution at element */
        ps0 = &mesh->sol[k];
        t1.va = t1.vb = t1.vc = ps0->bb;
        t2.va = t2.vb = t2.vc = ps0->bb;
      }
      /* color interpolation */
      cutTriangle(sc,t1);
      cutTriangle(sc,t2);
      k = pq->nxt;
    }
    glEnd();
  }

  glEndList();
  return(dlist);
}


/* build list of tetrahedra */
GLuint listTetraMap(pScene sc,pMesh mesh,ubyte clip) {
  pMaterial  pm;
  pTetra     pt;
  pPoint     p0,p1,p2;
  pSolution  ps0,ps1,ps2;
  GLint      dlist = 0;
  float      cx,cy,cz,ax,ay,az,bx,by,bz,d,n[3];
  int        k,l,m;
  triangle   t;

  /* default */
  if ( !mesh->ntet )  return(0);
  if ( ddebug ) printf("create display list map / TETRA\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  // By Leo: get number of triangles to render tet colors correctly
  int boundary_faces = mesh->nt;

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
#ifdef IGL
  bool transp = sc->igl_params->tet_color[3] < 0.999;
  int old_depth_func =0;
  glGetIntegerv(GL_DEPTH_FUNC,&old_depth_func);
  if ( transp )
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_ALWAYS);
    sc->igl_params->alpha_holder = sc->igl_params->tet_color[3];
  }else
  {
    sc->igl_params->alpha_holder = 1.0;
  }
#endif

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LTets];
    if ( !k || pm->flag )  continue;

    glBegin(GL_TRIANGLES);
    while ( k != 0 ) {
      pt = &mesh->tetra[k];
      if ( !pt->v[0] || (clip && !pt->clip) ) {
        k = pt->nxt;
        continue;
      }
      /* build 4 faces */
      cx = cy = cz = 0.0f;
      for (l=0; l<4; l++) {
    p0  = &mesh->point[pt->v[l]];
    cx += p0->c[0];
        cy += p0->c[1];
        cz += p0->c[2];
      }
      cx /= 4.;
      cy /= 4.;
      cz /= 4.;

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

        /* store triangle */
        t.a[0] = sc->shrink*(p0->c[0]-cx)+cx;
        t.a[1] = sc->shrink*(p0->c[1]-cy)+cy;
        t.a[2] = sc->shrink*(p0->c[2]-cz)+cz;

        t.b[0] = sc->shrink*(p1->c[0]-cx)+cx;
        t.b[1] = sc->shrink*(p1->c[1]-cy)+cy;
        t.b[2] = sc->shrink*(p1->c[2]-cz)+cz;

        t.c[0] = sc->shrink*(p2->c[0]-cx)+cx;
        t.c[1] = sc->shrink*(p2->c[1]-cy)+cy; 
        t.c[2] = sc->shrink*(p2->c[2]-cz)+cz; 
        
        /* store normals */
        memcpy(t.na,n,3*sizeof(float));
        memcpy(t.nb,n,3*sizeof(float));
        memcpy(t.nc,n,3*sizeof(float));

        if ( mesh->typage == 2 ) {
          /* solutions at vertices */
      ps0 = &mesh->sol[pt->v[ct[l][0]]];
      ps1 = &mesh->sol[pt->v[ct[l][1]]];
      ps2 = &mesh->sol[pt->v[ct[l][2]]];
          t.va = ps0->bb;
      t.vb = ps1->bb;
      t.vc = ps2->bb;
        }
        else {
          /* solution at element */  
          ps0 = &mesh->sol[k+boundary_faces];
      t.va = t.vb = t.vc = ps0->bb;
        }
        /* color interpolation */
        cutTriangle(sc,t);
      }
      k = pt->nxt;
    }
    glEnd();
  }
#ifdef IGL
  if(transp)
  {
    glDepthFunc(old_depth_func);
    glDisable(GL_BLEND);
  }
#endif

  glEndList();
  return(dlist);
}


/* build list of hexahedra */
GLuint listHexaMap(pScene sc,pMesh mesh,ubyte clip) {
  pMaterial  pm;
  pHexa      ph;
  pPoint     p0,p1,p2,p3;
  pSolution  ps0,ps1,ps2,ps3;
  GLint      dlist = 0;
  double     ax,ay,az,bx,by,bz,d;
  float      n[3],cx,cy,cz;
  int        k,l,m;
  triangle   t1,t2;

  if ( !mesh->nhex )  return(0);
  if ( ddebug ) printf("create display list map / HEXA\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  /* build list */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    k  = pm->depmat[LHexa];
    if ( !k || pm->flag )  continue;

    glBegin(GL_TRIANGLES);
    while ( k != 0 ) {
      ph = &mesh->hexa[k];
      if ( !ph->v[0] || (clip && !ph->clip) ) {
        k = ph->nxt;
        continue;
      }
      cx = cy = cz = 0.0f;
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

        /* store triangles */
        t1.a[0] = t2.a[0] = sc->shrink*(p0->c[0]-cx)+cx;
        t1.a[1] = t2.a[1] = sc->shrink*(p0->c[1]-cy)+cy;
        t1.a[2] = t2.a[2] = sc->shrink*(p0->c[2]-cz)+cz;

        t1.b[0] = sc->shrink*(p1->c[0]-cx)+cx;
        t1.b[1] = sc->shrink*(p1->c[1]-cy)+cy;
        t1.b[2] = sc->shrink*(p1->c[2]-cz)+cz;

        t1.c[0] = t2.b[0] = sc->shrink*(p2->c[0]-cx)+cx; 
        t1.c[1] = t2.b[1] = sc->shrink*(p2->c[1]-cy)+cy; 
        t1.c[2] = t2.b[2] = sc->shrink*(p2->c[2]-cz)+cz; 

        t2.c[0] = sc->shrink*(p3->c[0]-cx)+cx; 
        t2.c[1] = sc->shrink*(p3->c[1]-cy)+cy; 
        t2.c[2] = sc->shrink*(p3->c[2]-cz)+cz; 

        /* store normals */
    memcpy(t1.na,n,3*sizeof(float));
    memcpy(t1.nb,n,3*sizeof(float));
    memcpy(t1.nc,n,3*sizeof(float));
    memcpy(t2.na,n,3*sizeof(float));
    memcpy(t2.nb,n,3*sizeof(float));
    memcpy(t2.nc,n,3*sizeof(float));

        if ( mesh->typage == 2 ) {
          /* solutions at vertices */
      ps0 = &mesh->sol[ph->v[ch[l][0]]];
      ps1 = &mesh->sol[ph->v[ch[l][1]]];
      ps2 = &mesh->sol[ph->v[ch[l][2]]];
      ps3 = &mesh->sol[ph->v[ch[l][3]]];
      t1.va = t2.va = ps0->bb;
      t1.vb = ps1->bb;
      t1.vc = t2.vb = ps2->bb;
      t2.vc = ps3->bb;
        }
        else {
          /* solution at element */
      ps0 = &mesh->sol[k];
      t1.va = t1.vb = t1.vc = ps0->bb;
      t2.va = t2.vb = t2.vc = ps0->bb;
        }
        /* color interpolation */
        cutTriangle(sc,t1);
        cutTriangle(sc,t2);
     }
     k = ph->nxt;
   }
   glEnd();
  }
  
  glEndList();
  return(dlist);
}


GLuint alt2dList(pScene sc,pMesh mesh,int geomtype,float shrink,float altcoef) {
  pTriangle  pt,pt1;
  pMaterial  pm;
  pQuad      pq;
  pPoint     p0,p1,p2,p3;
  pSolution  ps0,ps1,ps2,ps3;
  GLuint     dlist;
  double     ax,ay,az,bx,by,bz,dd,kc,rgb[4];
  float      cx,cy,cz,n[3];
  int       *adj,k,m,ia,iadr;
  ubyte     *voy;
  triangle   t,t1,t2;
  static double hsv[3] = { 0.0, 1.0, 0.80 };
  static float  nn[3] = {1.0, 0.0, 0.0 };

  /* default */
  if ( ddebug ) printf("create 2d elevation map list\n");

  if ( geomtype == LTria && !mesh->nt )  return(0);
  if ( geomtype == LQuad && !mesh->nq )  return(0);
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);

  mesh->zmin = altcoef*mesh->bbmin;
  mesh->zmax = altcoef*mesh->bbmax;
  if ( mesh->bbmin*mesh->bbmax < 0.0 ) {
    mesh->ztra = mesh->zmin;
  }
  else {
    mesh->ztra = 0.95 * mesh->zmin;
  }

  switch (geomtype) {
  case LTria:
    if ( ddebug ) printf("create triangle list %d\n",mesh->nt);
    
    if ( mesh->typage == 1 ) {
      if ( mesh->nt && !hashTria(mesh) )    return(0);
    }

    glBegin(GL_TRIANGLES);

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
      
        if ( mesh->typage == 1 )
          ps0 = ps1 = ps2 = &mesh->sol[k];
        else {
          ps0  = &mesh->sol[pt->v[0]];
          ps1  = &mesh->sol[pt->v[1]];
          ps2  = &mesh->sol[pt->v[2]];
        }
        cx = (p0->c[0] + p1->c[0] + p2->c[0]) / 3.0;
        cy = (p0->c[1] + p1->c[1] + p2->c[1]) / 3.0;
        cz = (ps0->bb + ps1->bb + ps2->bb) / 3.0;

        t.a[0] = shrink*(p0->c[0]-cx) + cx;
        t.a[1] = shrink*(p0->c[1]-cy) + cy;
        t.a[2] = shrink*(altcoef*ps0->bb-cz) + cz - 0.25*mesh->ztra;

        t.b[0] = shrink*(p1->c[0]-cx) + cx;
        t.b[1] = shrink*(p1->c[1]-cy) + cy;
        t.b[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
        
        t.c[0] = shrink*(p2->c[0]-cx) + cx;
        t.c[1] = shrink*(p2->c[1]-cy) + cy;
        t.c[2] = shrink*(altcoef*ps2->bb-cz) + cz - 0.25*mesh->ztra;
        
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

        memcpy(t.na,n,3*sizeof(float));
        memcpy(t.nb,n,3*sizeof(float));
        memcpy(t.nc,n,3*sizeof(float));

        t.va = ps0->bb;
        t.vb = ps1->bb;
        t.vc = ps2->bb;

        if ( mesh->typage == 2 ) 
          cutTriangle(sc,t);

        else {
          if ( t.va < sc->iso.val[0] ) 
            t.va = sc->iso.val[0];  
          else if ( t.va > sc->iso.val[MAXISO-1] )
            t.va = sc->iso.val[MAXISO-1];
          for (ia=0; ia<MAXISO-1; ia++)
            if ( t.va < sc->iso.val[ia] )  break;
          kc = (t.va-sc->iso.val[ia-1]) / (sc->iso.val[ia] - sc->iso.val[ia-1]);
          hsv[0] = sc->iso.col[ia-1]*(1.0-kc)+sc->iso.col[ia]*kc;
          hsvrgb(hsv,rgb);

          glColor4dv(rgb);
          glNormal3fv(t.na);
          glVertex3fv(t.a);
          glVertex3fv(t.b);
          glVertex3fv(t.c);

          /* add quads to sides (thanks to F. Lagoutiere) */
          iadr = 3*(k-1)+1;
          adj  = &mesh->adja[iadr];
          voy  = &mesh->voy[iadr];

          if ( adj[0] && adj[0] < k ) {
            pt1 = &mesh->tria[ adj[0] ]; 
            p3  = &mesh->point[ pt1->v[voy[0]] ];
            ps1 = &mesh->sol[ adj[0] ];
            
            cx = (p1->c[0] + p2->c[0] + p3->c[0]) / 3.0;
            cy = (p1->c[1] + p2->c[1] + p3->c[1]) / 3.0;
            cz = ps1->bb;

            memcpy(t1.a,t.b,3*sizeof(float));
            memcpy(t1.b,t.c,3*sizeof(float));
            memcpy(t1.c,t.b,3*sizeof(float));
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps0->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);

            memcpy(t1.a,t.c,3*sizeof(float));
            memcpy(t1.b,t.c,3*sizeof(float));
            memcpy(t1.c,t.b,3*sizeof(float));
            t1.b[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps1->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);
          }
          if ( adj[1] && adj[1] < k ) {
            pt1 = &mesh->tria[ adj[1] ]; 
            p3  = &mesh->point[ pt1->v[voy[1]] ];
            ps1 = &mesh->sol[ adj[1] ];
            
            cx = (p0->c[0] + p2->c[0] + p3->c[0]) / 3.0;
            cy = (p0->c[1] + p2->c[1] + p3->c[1]) / 3.0;
            cz = ps1->bb;

            memcpy(t1.a,t.a,3*sizeof(float));
            memcpy(t1.b,t.c,3*sizeof(float));
            memcpy(t1.c,t.a,3*sizeof(float));
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps0->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);

            memcpy(t1.a,t.c,3*sizeof(float));
            memcpy(t1.b,t.c,3*sizeof(float));
            memcpy(t1.c,t.a,3*sizeof(float));
            t1.b[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps1->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);
          }
          if ( adj[2] && adj[2] < k ) {
            pt1 = &mesh->tria[ adj[2] ]; 
            p3  = &mesh->point[ pt1->v[voy[2]] ];
            ps1 = &mesh->sol[ adj[2] ];
            
            cx = (p0->c[0] + p1->c[0] + p3->c[0]) / 3.0;
            cy = (p0->c[1] + p1->c[1] + p3->c[1]) / 3.0;
            cz = ps1->bb;

            memcpy(t1.a,t.a,3*sizeof(float));
            memcpy(t1.b,t.b,3*sizeof(float));
            memcpy(t1.c,t.a,3*sizeof(float));
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps0->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);

            memcpy(t1.a,t.b,3*sizeof(float));
            memcpy(t1.b,t.b,3*sizeof(float));
            memcpy(t1.c,t.a,3*sizeof(float));
            t1.b[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.c[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
            t1.va = ps0->bb;
            t1.vb = ps1->bb;
            t1.vc = ps1->bb;
            memcpy(t1.na,nn,3*sizeof(float));
            memcpy(t1.nb,nn,3*sizeof(float));
            memcpy(t1.nc,nn,3*sizeof(float));
            cutTriangle(sc,t1);
          }
        }

        k = pt->nxt;
      }
    }
    glEnd();
    break;
    
  case LQuad:
    if ( ddebug ) printf("create quadrilateral list %d\n",mesh->nq);

    glBegin(GL_TRIANGLES);
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
        
        p0 = &mesh->point[pq->v[0]];
        p1 = &mesh->point[pq->v[1]];
        p2 = &mesh->point[pq->v[2]];
        p3 = &mesh->point[pq->v[3]];
      
        if ( mesh->typage == 1 )
          ps0 = ps1 = ps2 = ps3 = &mesh->sol[k];
        else {
      ps0 = &mesh->sol[pq->v[0]];
      ps1 = &mesh->sol[pq->v[1]];
      ps2 = &mesh->sol[pq->v[2]];
      ps3 = &mesh->sol[pq->v[3]];
    }
        cx = 0.25 * (p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]);
        cy = 0.25 * (p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]);
        cz = 0.25 * (ps0->bb + ps1->bb + ps2->bb + ps3->bb);

        t.a[0] = t2.a[0] = shrink*(p0->c[0]-cx) + cx;
        t.a[1] = t2.a[1] = shrink*(p0->c[1]-cy) + cy;
        t.a[2] = t2.a[2] = shrink*(altcoef*ps0->bb-cz) + cz - 0.25*mesh->ztra;

        t.b[0] = shrink*(p1->c[0]-cx) + cx;
        t.b[1] = shrink*(p1->c[1]-cy) + cy;
        t.b[2] = shrink*(altcoef*ps1->bb-cz) + cz - 0.25*mesh->ztra;
        
        t.c[0] = t2.b[0] = shrink*(p2->c[0]-cx) + cx;
        t.c[1] = t2.b[1] = shrink*(p2->c[1]-cy) + cy;
        t.c[2] = t2.b[2] = shrink*(altcoef*ps2->bb-cz) + cz - 0.25*mesh->ztra;

        t2.c[0] = shrink*(p3->c[0]-cx) + cx;
        t2.c[1] = shrink*(p3->c[1]-cy) + cy;
        t2.c[2] = shrink*(altcoef*ps3->bb-cz) + cz - 0.25*mesh->ztra;

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
        memcpy(t.na,n,3*sizeof(float));
        memcpy(t.nb,n,3*sizeof(float));
        memcpy(t.nc,n,3*sizeof(float));
        memcpy(t2.na,n,3*sizeof(float));
        memcpy(t2.nb,n,3*sizeof(float));
        memcpy(t2.nc,n,3*sizeof(float));

        t.va = t2.va = ps0->bb;
        t.vb = ps1->bb;
        t.vc = t2.vb = ps2->bb;
        t2.vc = ps3->bb;

        cutTriangle(sc,t);
        cutTriangle(sc,t2);

        k = pq->nxt;
      }
    }
    glEnd();
    break;
  }
  glEndList();

  return(dlist);
}


/* setup color table */
void setupPalette(pScene sc,pMesh mesh) {
  double     delta;
  int        i;

  if ( ddebug ) printf("create palette %f %f\n",mesh->bbmin,mesh->bbmax);

  if ( !sc->iso.palette ) {
    delta = mesh->bbmax - mesh->bbmin;
    for (i=0; i<MAXISO; i++) {
      sc->iso.col[i] = 240.0 *(1.0 - (float)i/(MAXISO-1));
      sc->iso.val[i] = mesh->bbmin + i * delta/(MAXISO-1);
    }
    sc->iso.palette = 1;
  }
  else {
    for (i=0; i<MAXISO; i++)
      sc->iso.col[i] = 240.0 *(1.0 - (float)i/(MAXISO-1));
  }
}


/* build color palette */
GLuint drawPalette(pScene sc) {
  double     rgb[3];
  float      xpos,ypos,inc,top,bottom,left,right;
  int        i;
  static double hsv[3] = {1.0, 1.0, 0.80};

  if ( sc->iso.palette < 3 ) {
    if ( sc->iso.palette <= 2 ) {
      top     = sc->par.ys - 20;
      bottom  = top - 10;
      left    = sc->par.xs / 10;
      right   = sc->par.xs - left;
    }
    else if ( sc->iso.palette == 2 ) {
      top     = 40;
      bottom  = top - 10;
      left    = sc->par.xs / 10;
      right   = sc->par.xs - left;
    }

    inc  = (sc->par.xs-2*left) / (MAXISO-1);
    xpos = left;
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    for (i=1; i<MAXISO; i++) {
#ifdef IGL
    sc->igl_params->rgb( 1.0-(sc->iso.col[i-1])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i-1];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(xpos,bottom);
      glVertex2f(xpos,top);

#ifdef IGL
    sc->igl_params->rgb( 1.0-(sc->iso.col[i])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(xpos+inc,top);
      glVertex2f(xpos+inc,bottom);
      xpos += inc;
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(1.0);
    glColor3fv(sc->par.line);
    glRectf(left,bottom,xpos,top);

    /* graduate  */
    xpos = left / 2;
    bottom -= 10;
    top    += 5;
    output2(xpos,bottom,"%.4E",sc->iso.val[0]);
    xpos = left;
    for (i=1; i<MAXISO-1; i++) {
      xpos += inc;
      if ( (i % 2) == 0 ) {
        output2(xpos-25.0,bottom,"%.4E",sc->iso.val[i]);
        glBegin(GL_LINES);
          glVertex2f(xpos,bottom+15);
          glVertex2f(xpos,bottom+10);
        glEnd();
      }
      else {
        output2(xpos-25.0,top,"%.4E",sc->iso.val[i]);
        glBegin(GL_LINES);
          glVertex2f(xpos,top-10);
          glVertex2f(xpos,top-5);
        glEnd();
      }
    }
    output2(right-left/2,bottom,"%.4E",sc->iso.val[MAXISO-1]);
  }

  else if ( sc->iso.palette == 3 ) {
    bottom  = sc->par.ys / 10;
    top     = sc->par.ys - bottom;
    left    = 10;
    right   = left + 15;

    inc  = (sc->par.ys-2*bottom) / (MAXISO-1);
    ypos = bottom;
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    for (i=1; i<MAXISO; i++) {
#ifdef IGL
    sc->igl_params->rgb( 1.0-(sc->iso.col[i-1])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i-1];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(left,ypos);
      glVertex2f(right,ypos);

#ifdef IGL
      sc->igl_params->rgb( 1.0-(sc->iso.col[i])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(right,ypos+inc);
      glVertex2f(left,ypos+inc);
      ypos += inc;
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(1.0);
    glColor3fv(sc->par.line);
    glRectf(left,bottom,right,top);

    /* graduate  */
    ypos    = bottom;
    right  += 5;
    output2(right,bottom,"%.4E",sc->iso.val[0]);
    for (i=1; i<MAXISO-1; i++) {
      ypos += inc;
      output2(right,ypos,"%.4E",sc->iso.val[i]);
      glBegin(GL_LINES);
        glVertex2f(right-5,ypos);
        glVertex2f(right-13,ypos);
      glEnd();
    }
    output2(right,top,"%.4E",sc->iso.val[MAXISO-1]);
  }

  else {
    bottom  = sc->par.ys / 10;
    top     = sc->par.ys - bottom;
    right   = sc->par.xs - 10;
    left    = right - 15;
    inc  = (sc->par.ys-2*bottom) / (MAXISO-1);

    ypos = bottom;
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    for (i=1; i<MAXISO; i++) {
#ifdef IGL
      sc->igl_params->rgb( 1.0-(sc->iso.col[i-1])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i-1];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(left,ypos);
      glVertex2f(right,ypos);

#ifdef IGL
      sc->igl_params->rgb( 1.0-(sc->iso.col[i])/240.0,rgb);
#else
      hsv[0] = sc->iso.col[i];
      hsvrgb(hsv,rgb);
#endif
      glColor3dv(rgb);
      glVertex2f(right,ypos+inc);
      glVertex2f(left,ypos+inc);
      ypos += inc;
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(1.0);
    glColor3fv(sc->par.line);
    glRectf(left,bottom,right,top);

    /* graduate  */
    ypos  = bottom;
    left -= 65;
    output2(left,bottom,"%.4E",sc->iso.val[0]);
    for (i=1; i<MAXISO-1; i++) {
      ypos += inc;
      output2(left,ypos,"%.4E",sc->iso.val[i]);
      glBegin(GL_LINES);
        glVertex2f(left+65,ypos);
        glVertex2f(left+75,ypos);
      glEnd();
    }
    output2(left,top,"%.4E",sc->iso.val[MAXISO-1]);
  }

  return(1);
}

