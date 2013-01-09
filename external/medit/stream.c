#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define MAX_PTS   100000
#define MAX_CPT   5000
#define MAX_LST   1024

#define COS175  -0.996194698
#define COS178  -0.99939083
#define HSIZ     0.03
#define EPST    -1.e-14
#define EPSR     1.e+14


extern int       reftype,refitem;
extern GLfloat   altcoef;
int    nbar=0;
enum   {Euler=1,RK4=2};


/* find tetra containg p, starting nsdep */
int locateTetra(pMesh mesh,int nsdep,int base,float *p,double *cb) {
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   epsra,vol1,vol2,vol3,vol4,dd; 
  int     *adj,iadr,it,nsfin;

  it    = 0;
  nsfin = nsdep;
  /*  printf("locateTetra: searching for %f %f %f\n",p[0],p[1],p[2]);*/
  do {
    if ( !nsfin )  return(0);
    pt = &mesh->tetra[nsfin];
    if ( !pt->v[0] )  return(0);
    if ( pt->mark == base )  return(0);
    pt->mark = base;
    iadr = 4*(nsfin-1)+1;
    adj  = &mesh->adja[iadr];
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];

    /* barycentric */
    bx  = p1->c[0] - p0->c[0];
    by  = p1->c[1] - p0->c[1];
    bz  = p1->c[2] - p0->c[2];
    cx  = p2->c[0] - p0->c[0];
    cy  = p2->c[1] - p0->c[1];
    cz  = p2->c[2] - p0->c[2];
    dx  = p3->c[0] - p0->c[0];
    dy  = p3->c[1] - p0->c[1];
    dz  = p3->c[2] - p0->c[2];

    /* test volume */
    vx  = cy*dz - cz*dy;
    vy  = cz*dx - cx*dz;
    vz  = cx*dy - cy*dx;

    epsra = EPST*(bx*vx + by*vy + bz*vz);
    apx = p[0] - p0->c[0];
    apy = p[1] - p0->c[1];
    apz = p[2] - p0->c[2];

    /* p in 2 */
    vol2  = apx*vx + apy*vy + apz*vz;
    if ( epsra > vol2 ) {
      nsfin = adj[1];
      continue;
    }

    /* p in 3 */
    vx  = by*apz - bz*apy;
    vy  = bz*apx - bx*apz;
    vz  = bx*apy - by*apx;
    vol3 = dx*vx + dy*vy + dz*vz;
    if ( epsra > vol3 ) {
      nsfin = adj[2];
      continue;
    }
    
    /* p in 4 */
    vol4 = -cx*vx - cy*vy - cz*vz;
    if ( epsra > vol4 ) {
      nsfin = adj[3];
      continue;
    }
    
    /* p in 1 */
    vol1 = -epsra * EPSR - vol2 - vol3 - vol4;
    if ( epsra > vol1 ) {
      nsfin = adj[0];
      continue;
    }

    dd = vol1+vol2+vol3+vol4;
    if ( dd != 0.0 )  dd = 1.0 / dd;
    cb[0] = vol1 * dd;
    cb[1] = vol2 * dd;
    cb[2] = vol3 * dd;
    cb[3] = vol4 * dd; 
    
    pt->cpt++;
    return(nsfin);
  }
  while ( ++it <= mesh->ntet );

  return(0);
}


/* return 1 if point in tetra, adjacent #, if not */
int inSubTetra(pPoint pt[4],float *p,double *cb) {
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   epsra,vol1,vol2,vol3,vol4,dd; 

  /* barycentric */
  bx  = pt[1]->c[0] - pt[0]->c[0];
  by  = pt[1]->c[1] - pt[0]->c[1];
  bz  = pt[1]->c[2] - pt[0]->c[2];
  cx  = pt[2]->c[0] - pt[0]->c[0];
  cy  = pt[2]->c[1] - pt[0]->c[1];
  cz  = pt[2]->c[2] - pt[0]->c[2];
  dx  = pt[3]->c[0] - pt[0]->c[0];
  dy  = pt[3]->c[1] - pt[0]->c[1];
  dz  = pt[3]->c[2] - pt[0]->c[2];

  /* test volume */
  vx  = cy*dz - cz*dy;
  vy  = cz*dx - cx*dz;
  vz  = cx*dy - cy*dx;

  epsra = EPST*(bx*vx + by*vy + bz*vz);
  apx = p[0] - pt[0]->c[0];
  apy = p[1] - pt[0]->c[1];
  apz = p[2] - pt[0]->c[2];

  /* p in 2 */
  vol2 = apx*vx + apy*vy + apz*vz;
  if ( epsra > vol2 )  return(-1);

  /* p in 3 */
  vx  = by*apz - bz*apy;
  vy  = bz*apx - bx*apz;
  vz  = bx*apy - by*apx;
  vol3 = dx*vx + dy*vy + dz*vz;
  if ( epsra > vol3 )  return(-2);
  
  /* p in 4 */
  vol4 = -cx*vx - cy*vy - cz*vz;
  if ( epsra > vol4 )  return(-3);

  /* p in 1 */
  vol1 = -epsra * EPSR - vol2 - vol3 - vol4;
  if ( epsra > vol1 )  return(0);

  dd = vol1+vol2+vol3+vol4;
  if ( dd != 0.0f )  dd = 1.0 / dd;
  cb[0] = vol1 * dd;
  cb[1] = vol2 * dd;
  cb[2] = vol3 * dd;
  cb[3] = vol4 * dd; 
  
  return(1);
}

int locateHexa(pMesh mesh,int nsdep,int base,float *p,double *cb,pPoint pt[4]) {
  pHexa    ph;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   epsra,vol1,vol2,vol3,vol4,dd; 
  int     *adj,iadr,it,nsfin,in;

  it    = 0;
  nsfin = nsdep;
  /*printf("locateHexa: searching for %f %f %f\n",p[0],p[1],p[2]);*/
  do {
    if ( !nsfin )  return(0);
    ph = &mesh->hexa[nsfin];
/*printf("\nnsfin %d  base %d  mark %d\n",nsfin,base,ph->mark);*/
    if ( !ph->v[0] || ph->mark == base )  return(0);
    ph->mark = base;
    iadr = 6*(nsfin-1)+1;
    adj  = &mesh->adja[iadr];
/*printf("adj %d %d %d %d %d %d\n",adj[0],adj[1],adj[2],adj[3],adj[4],adj[5]);*/

    /* tetra1: 0,2,3,7 : 3 external faces */
/*printf("tet1: %d %d %d %d\n",ph->v[0],ph->v[2],ph->v[3],ph->v[7]);*/
    pt[0] = &mesh->point[ph->v[0]];
    pt[1] = &mesh->point[ph->v[2]];
    pt[2] = &mesh->point[ph->v[3]];
    pt[3] = &mesh->point[ph->v[7]];
    in = inSubTetra(pt,p,cb);
printf("tet1 : on sort en %d\n",in);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == 0 ) {
      nsfin = adj[4];
      continue;
    }
    else if ( in == -1 ) {
      nsfin = adj[5];
      continue;
    }
    else if ( in == -3 ) {
      nsfin = adj[0];
      continue;
    }

    /* tetra2: 1,4,5,6 : 3 external faces */
    pt[0] = &mesh->point[ph->v[1]];
    pt[1] = &mesh->point[ph->v[4]];
    pt[2] = &mesh->point[ph->v[5]];
    pt[3] = &mesh->point[ph->v[6]];
    in = inSubTetra(pt,p,cb);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == 0 ) {
      nsfin = adj[1];
      continue;
    }
    else if ( in == -1 ) {
      nsfin = adj[3];
      continue;
    }
    else if ( in == -3 ) {
      nsfin = adj[2];
      continue;
    }

    /* tetra3: 0,4,6,7 : 2 external faces */
    pt[0] = &mesh->point[ph->v[0]];
    pt[1] = &mesh->point[ph->v[4]];
    pt[2] = &mesh->point[ph->v[6]];
    pt[3] = &mesh->point[ph->v[7]];
    in = inSubTetra(pt,p,cb);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == 0 ) {
      nsfin = adj[1];
      continue;
    }
    else if ( in == -2 ) {
      nsfin = adj[5];
      continue;
    }

    /* tetra4: 0,1,2,6 : 2 external faces */
    pt[0] = &mesh->point[ph->v[0]];
    pt[1] = &mesh->point[ph->v[1]];
    pt[2] = &mesh->point[ph->v[2]];
    pt[3] = &mesh->point[ph->v[6]];
    in = inSubTetra(pt,p,cb);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == 0 ) {
      nsfin = adj[3];
      continue;
    }
    else if ( in == -3 ) {
      nsfin = adj[0];
      continue;
    }
    
    /* tetra5: 0,6,2,7 : 1 external face */
    pt[0] = &mesh->point[ph->v[0]];
    pt[1] = &mesh->point[ph->v[6]];
    pt[2] = &mesh->point[ph->v[2]];
    pt[3] = &mesh->point[ph->v[7]];
    in = inSubTetra(pt,p,cb);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == 0 ) {
      nsfin = adj[4];
      continue;
    }

    /* tetra6: 0,4,1,6 : 1 external face */
    pt[0] = &mesh->point[ph->v[0]];
    pt[1] = &mesh->point[ph->v[4]];
    pt[2] = &mesh->point[ph->v[1]];
    pt[3] = &mesh->point[ph->v[6]];
    in = inSubTetra(pt,p,cb);
    if ( in > 0 ) {
      ph->cpt++;
      return(nsfin);
    }
    else if ( in == -3 ) {
      nsfin = adj[2];
      continue;
    }
puts("PROBLEME");
exit(1);
  }
  while ( ++it <= mesh->nhex );

  return(0);
}

int locateTria(pMesh mesh,int nsdep,int base,float *p,double *cb) {
  pTriangle pt;
  pPoint    p0,p1,p2;
  double    ax,ay,bx,by,cx,cy;
  double    epsra,aire1,aire2,aire3,dd; 
  int      *adj,iadr,it,isign,nsfin;

  it    = 0;
  nsfin = nsdep;
  /*printf("locateTria: searching for %f %f\n",p[0],p[1]);*/
  do {
    pt = &mesh->tria[nsfin];
    if ( !pt->v[0] )  return(0);
    if ( pt->mark == base )  return(0);
    pt->mark = base;
    iadr = 3*(nsfin-1)+1;
    adj  = &mesh->adja[iadr];

    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    isign= dd > 0 ? 1 : -1;
    epsra = isign > 0 ? EPST*dd : -(EPST*dd);
    /* barycentric */
    bx = p1->c[0] - p[0];
    by = p1->c[1] - p[1];
    cx = p2->c[0] - p[0];
    cy = p2->c[1] - p[1];
    /* p in 1 */
    aire1 = isign*(bx*cy - by*cx);
    if ( epsra > aire1 ) {
      nsfin = adj[0];
      continue;
    }

    ax = p0->c[0] - p[0];
    ay = p0->c[1] - p[1];
    aire2 = isign*(cx*ay - cy*ax);
    if ( epsra > aire2 ) {
      nsfin = adj[1];
      continue;
    }

    aire3 = -epsra*EPSR - aire1 - aire2;
    if ( epsra > aire3 ) {
      nsfin = adj[2];
      continue;
    }

    dd = aire1+aire2+aire3;
    if ( dd != 0.0f )  dd = 1.0 / dd;
    cb[0] = aire1 * dd;
    cb[1] = aire2 * dd;
    cb[2] = aire3 * dd;
    
    pt->cpt++;
    return(nsfin);
  }
  while ( ++it <= mesh->nt );

  return(0);
}

/* point in tetra */
int inTetra(pMesh mesh,int nsdep,float *p,double *cb) {
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
  double   epsra,vol1,vol2,vol3,vol4,dd; 

  pt = &mesh->tetra[nsdep];
  if ( !pt->v[0] )  return(0);

  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];
  p3 = &mesh->point[pt->v[3]];

  /* barycentric */
  bx  = p1->c[0] - p0->c[0];
  by  = p1->c[1] - p0->c[1];
  bz  = p1->c[2] - p0->c[2];
  cx  = p2->c[0] - p0->c[0];
  cy  = p2->c[1] - p0->c[1];
  cz  = p2->c[2] - p0->c[2];
  dx  = p3->c[0] - p0->c[0];
  dy  = p3->c[1] - p0->c[1];
  dz  = p3->c[2] - p0->c[2];

  /* test volume */
  vx  = cy*dz - cz*dy;
  vy  = cz*dx - cx*dz;
  vz  = cx*dy - cy*dx;

  epsra = EPST*(bx*vx + by*vy + bz*vz);
  apx = p[0] - p0->c[0];
  apy = p[1] - p0->c[1];
  apz = p[2] - p0->c[2];

  /* p in 2 */
  vol2  = apx*vx + apy*vy + apz*vz;
  if ( epsra > vol2 )  return(0);

  /* p in 3 */
  vx  = by*apz - bz*apy;
  vy  = bz*apx - bx*apz;
  vz  = bx*apy - by*apx;
  vol3 = dx*vx + dy*vy + dz*vz;
  if ( epsra > vol3 )  return(0);
    
  /* p in 4 */
  vol4 = -cx*vx - cy*vy - cz*vz;
  if ( epsra > vol4 )  return(0);
  
  /* p in 1 */
  vol1 = -epsra * EPSR - vol2 - vol3 - vol4;
  if ( epsra > vol1 )  return(0);

  dd = vol1+vol2+vol3+vol4;
  if ( dd != 0.0f )  dd = 1.0 / dd;
  cb[0] = vol1 * dd;
  cb[1] = vol2 * dd;
  cb[2] = vol3 * dd;
  cb[3] = vol4 * dd; 

  pt->cpt++;
  return(1);
}

int inHexa(pMesh mesh,int nsdep,float *p,double *cb,pPoint pt[4]) {
  return(0);
}

int inTria(pMesh mesh,int nsdep,float *p,double *cb) {
  pTriangle  pt;
  pPoint     p0,p1,p2;
  double     ax,ay,bx,by,cx,cy;
  double     epsra,dd,aire1,aire2,aire3;
  int        isign;

  pt = &mesh->tria[nsdep];
  if ( !pt->v[0] )  return(0);

  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];

  ax = p1->c[0] - p0->c[0];
  ay = p1->c[1] - p0->c[1];
  bx = p2->c[0] - p0->c[0];
  by = p2->c[1] - p0->c[1];
  dd = ax*by - ay*bx;
  isign= dd > 0 ? 1 : -1;
  epsra = isign > 0 ? EPST*dd : -(EPST*dd);

  /* barycentric */
  bx = p[0] - p1->c[0];
  by = p[1] - p1->c[1];
  cx = p[0] - p2->c[0];
  cy = p[1] - p2->c[1];
  aire1 = isign*(bx*cy - by*cx);
  if ( epsra > aire1 )  return(0);

  ax = p[0] - p0->c[0];
  ay = p[1] - p0->c[1];
  aire2 = isign*(cx*ay - cy*ax);
  if ( epsra > aire2 )  return(0);
  
  aire3 = -epsra*EPSR - aire1 - aire2;
  if ( epsra > aire3 )  return(0);

  dd = aire1+aire2+aire3;
  if ( dd != 0.0f )  dd = 1.0 / dd;
  cb[0] = aire1 * dd;
  cb[1] = aire2 * dd;
  cb[2] = aire3 * dd;

  pt->cpt++;
  return(1);
}


/* return size of tetra */
double sizeTetra(pMesh mesh,int k) {
  pTetra   pt;
  pPoint   p[4];
  double   ax,ay,az,dd;
  double   hmin;
  int      i;
  static int  idire[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};

  pt = &mesh->tetra[k];
  for (i=0; i<4; i++)
    p[i] = &mesh->point[pt->v[i]];
  hmin = FLT_MAX;
  for (i=1; i<6; i++) {
    ax = p[idire[i][0]]->c[0] - p[idire[i][1]]->c[0];
    ay = p[idire[i][0]]->c[1] - p[idire[i][1]]->c[1];
    az = p[idire[i][0]]->c[2] - p[idire[i][1]]->c[2];
    dd = ax*ax + ay*ay + az*az;
    hmin = MEDIT_MIN(dd,hmin);
  }
  return(sqrt(hmin));
}

double sizeHexa(pMesh mesh,int k) {
  pHexa    ph;
  pPoint   p[8];
  double   ax,ay,az,dd;
  double   hmin;
  int      i;
  static int idire[12][2] = {{0,1}, {1,2}, {2,3}, {0,3}, {4,5}, {5,6},
                             {6,7}, {4,7}, {0,4}, {1,5}, {2,6}, {3,7}};

  ph = &mesh->hexa[k];
  for (i=0; i<8; i++)
    p[i] = &mesh->point[ph->v[i]];
  hmin = FLT_MAX;
  for (i=1; i<12; i++) {
    ax = p[idire[i][0]]->c[0] - p[idire[i][1]]->c[0];
    ay = p[idire[i][0]]->c[1] - p[idire[i][1]]->c[1];
    az = p[idire[i][0]]->c[2] - p[idire[i][1]]->c[2];
    dd = ax*ax + ay*ay + az*az;
    hmin = MEDIT_MIN(dd,hmin);
  }
  return(sqrt(hmin));
}

double sizeTria(pMesh mesh,int k) {
  pTriangle pt;
  pPoint    p0,p1;
  double    ax,ay,dd;
  double    hmin;
  int       i;
  static int idir[5] = {0,1,2,0,1};

  pt = &mesh->tria[k];
  hmin = FLT_MAX;
  for (i=0; i<3; i++) {
    p0 = &mesh->point[pt->v[i]];
    p1 = &mesh->point[pt->v[idir[i+1]]];
    ax = p0->c[0] - p1->c[0];
    ay = p0->c[1] - p1->c[1];
    dd = ax*ax + ay*ay;
    hmin = MEDIT_MIN(dd,hmin);
  }
  return(sqrt(hmin));
}

double sizeQuad(pMesh mesh,int k) {
  pQuad     pq;
  pPoint    p0,p1;
  double    ax,ay,dd;
  double    hmin;
  int       i;
  static int idir[7] = {0,1,2,3,0,1,2};

  pq = &mesh->quad[k];
  hmin = FLT_MAX;
  for (i=0; i<4; i++) {
    p0 = &mesh->point[pq->v[i]];
    p1 = &mesh->point[pq->v[idir[i+1]]];
    ax = p0->c[0] - p1->c[0];
    ay = p0->c[1] - p1->c[1];
    dd = ax*ax + ay*ay;
    hmin = MEDIT_MIN(dd,hmin);
  }
  return(sqrt(hmin));
}

/* vector interpolation */
double field3DInterp(pMesh mesh,int iel,double *cb,double *v) {
  pTetra     pt;
  pSolution  ps0,ps1,ps2,ps3;
  double     dd;

  pt  = &mesh->tetra[iel];
  ps0 = &mesh->sol[pt->v[0]];
  ps1 = &mesh->sol[pt->v[1]];
  ps2 = &mesh->sol[pt->v[2]];
  ps3 = &mesh->sol[pt->v[3]];

  v[0] = cb[0]*ps0->m[0] + cb[1]*ps1->m[0] + \
         cb[2]*ps2->m[0] + cb[3]*ps3->m[0];
  v[1] = cb[0]*ps0->m[1] + cb[1]*ps1->m[1] + \
         cb[2]*ps2->m[1] + cb[3]*ps3->m[1];
  v[2] = cb[0]*ps0->m[2] + cb[1]*ps1->m[2] + \
         cb[2]*ps2->m[2] + cb[3]*ps3->m[2];
  dd = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if ( dd > 0.0f ) {
    v[0] /= dd;
    v[1] /= dd;
    v[2] /= dd;
  }
  return(dd);
}

/* vector interpolation */
double vector3DInterp(pMesh mesh,pPoint pt[4],double *cb,double *v) {
  pSolution  ps0,ps1,ps2,ps3;
  double     dd;

  ps0 = &mesh->sol[pt[0]-&mesh->point[0]];
  ps1 = &mesh->sol[pt[1]-&mesh->point[0]];
  ps2 = &mesh->sol[pt[2]-&mesh->point[0]];
  ps3 = &mesh->sol[pt[3]-&mesh->point[0]];

  v[0] = cb[0]*ps0->m[0] + cb[1]*ps1->m[0] + \
         cb[2]*ps2->m[0] + cb[3]*ps3->m[0];
  v[1] = cb[0]*ps0->m[1] + cb[1]*ps1->m[1] + \
         cb[2]*ps2->m[1] + cb[3]*ps3->m[1];
  v[2] = cb[0]*ps0->m[2] + cb[1]*ps1->m[2] + \
         cb[2]*ps2->m[2] + cb[3]*ps3->m[2];
  dd = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if ( dd > 0.0f ) {
    v[0] /= dd;
    v[1] /= dd;
    v[2] /= dd;
  }
  return(dd);
}


double field2DInterp(pMesh mesh,int iel,double *cb,double *v) {
  pTriangle  pt;
  pSolution  ps0,ps1,ps2;
  double     dd;

  pt  = &mesh->tria[iel];
  ps0 = &mesh->sol[pt->v[0]];
  ps1 = &mesh->sol[pt->v[1]];
  ps2 = &mesh->sol[pt->v[2]];

  v[0] = cb[0]*ps0->m[0] + cb[1]*ps1->m[0] + cb[2]*ps2->m[0];
  v[1] = cb[0]*ps0->m[1] + cb[1]*ps1->m[1] + cb[2]*ps2->m[1];
  v[2] = 0.0f;
  dd = sqrt(v[0]*v[0] + v[1]*v[1]);
  if ( dd > 0.0f ) {
    v[0] /= dd;
    v[1] /= dd;
  }
  return(dd);
}


/* add point to display list, if needed */
int filterPoint(pScene sc,Stream *st,float *p,ubyte color) {
  double  norm,rgb[3],kc,ux,uy,uz,vx,vy,vz,dd;
  int     i;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* store point */
  memcpy(st->stpt[++st->stnp],p,3*sizeof(float));
  nbar++;

  /* point color */
  norm = st->norm;
  if ( !color ) {
    if ( norm < sc->iso.val[0] ) 
      norm = sc->iso.val[0];  
    else if ( norm > sc->iso.val[MAXISO-1] )
      norm = sc->iso.val[MAXISO-1];
    for (i=0; i<MAXISO-1; i++)
      if ( norm < sc->iso.val[i] )  break;
    kc = (norm-sc->iso.val[i-1]) / (sc->iso.val[i] - sc->iso.val[i-1]);
    st->stcol[st->stnp] = sc->iso.col[i-1]*(1.0-kc)+sc->iso.col[i]*kc;
    st->stiso[st->stnp] = i;
  }
  else {
    st->stcol[st->stnp] = sc->iso.col[MAXISO-1];
    st->stiso[st->stnp] = MAXISO-1;
  }

  if ( sc->mode & S_ALTITUDE )
    st->stpt[st->stnp][2] = altcoef*norm;

  if ( !color && (st->stiso[0] != st->stiso[st->stnp]) ) {
    hsv[0] = st->stcol[st->stnp];
    hsvrgb(hsv,rgb);
    glColor3dv(rgb);
    glVertex3fv(st->stpt[st->stnp]);

    memcpy(st->stpt[0],st->stpt[st->stnp],3*sizeof(float));
    st->stcol[0] = st->stcol[st->stnp];
    st->stiso[0] = st->stiso[st->stnp];
    st->stnp = 0;
    return(1);
  }

  if ( st->stnp < 2 )  return(0);

  /* filtering point */
  ux = st->stpt[0][0] - st->stpt[1][0];
  uy = st->stpt[0][1] - st->stpt[1][1];
  uz = st->stpt[0][2] - st->stpt[1][2];
  dd = ux*ux + uy*uy + uz*uz;
  
  if ( dd > 0.0 ) {  
    dd = 1.0f / sqrt(dd);
    ux *= dd;
    uy *= dd;
    uz *= dd;
  }
  else {
    memcpy(st->stpt[0],st->stpt[1],2*3*sizeof(float));
    st->stnp--;
    return(0);
  }

  vx = st->stpt[2][0] - st->stpt[1][0];
  vy = st->stpt[2][1] - st->stpt[1][1];
  vz = st->stpt[2][2] - st->stpt[1][2];
  dd = vx*vx + vy*vy + vz*vz;
  if ( dd > 0.0 ) {
    dd = 1.0f / sqrt(dd);
    vx *= dd;
    vy *= dd;
    vz *= dd;
  }
  else {
    memcpy(st->stpt[1],st->stpt[2],3*sizeof(float));
    st->stnp--;
    return(0);
  }

  dd = ux*vx + uy*vy + uz*vz;
  if ( dd > COS178 ) {
    hsv[0] = st->stcol[st->stnp];
    hsvrgb(hsv,rgb);
    glColor3dv(rgb);
    glVertex3fv(st->stpt[st->stnp]);

    memcpy(st->stpt[0],st->stpt[st->stnp],3*sizeof(float));
    st->stcol[0] = st->stcol[st->stnp];
    st->stiso[0] = st->stiso[st->stnp];
    st->stnp = 0;
    return(1);
  }
  else {
    memcpy(st->stpt[st->stnp-1],st->stpt[st->stnp],3*sizeof(float));
    st->stcol[st->stnp-1] = st->stcol[st->stnp];
    st->stiso[st->stnp-1] = st->stiso[st->stnp];
    st->stnp--;
    return(0);
  }
}


/* add vertex to display list */
void addPoint(pScene sc,Stream *st,float *p,ubyte color) {
  double  norm,rgb[3],kc;
  int     i;
  static double hsv[3] = { 0.0f, 1.0f, 0.80f };

  /* point color */
  norm = st->norm;
  i = MAXISO-1;
  if ( !color ) {
    norm = st->norm;
    if ( norm < sc->iso.val[0] )
      norm = sc->iso.val[0];  
    else if ( norm > sc->iso.val[MAXISO-1] )
      norm = sc->iso.val[MAXISO-1];
    for (i=0; i<MAXISO-1; i++)
      if ( norm < sc->iso.val[i] )  break;
    kc = (norm-sc->iso.val[i-1]) / (sc->iso.val[i] - sc->iso.val[i-1]);
    hsv[0] = sc->iso.col[i-1]*(1.0-kc)+sc->iso.col[i]*kc;
  }

  if ( sc->mode & S_ALTITUDE )
    p[2] = altcoef*norm;

  hsvrgb(hsv,rgb);
  glColor3dv(rgb);
  glVertex3fv(p);
  st->stnp = 0;
  memcpy(st->stpt[st->stnp],p,3*sizeof(float));
  st->stcol[st->stnp] = hsv[0];
  st->stiso[st->stnp] = i;
  st->stnp++;
}


int nxtPoint3D(pMesh mesh,int nsdep,float *p,float step,double *v) {
  pTetra     pt;
  double     norm,h6,cb[4],v1[3],v2[3],v3[3];
  float      xp1[3],xp2[3],xp3[3];
  int        k;

  /* 4th order Runge-Kutta */
  xp1[0] = p[0] + 0.5*step*v[0];
  xp1[1] = p[1] + 0.5*step*v[1];
  xp1[2] = p[2] + 0.5*step*v[2];

  k = locateTetra(mesh,nsdep,++mesh->mark,xp1,cb);
  if ( !k )  return(0);
  norm = field3DInterp(mesh,k,cb,v1);
  pt = &mesh->tetra[k];
  pt->cpt--;

  xp2[0] = p[0] + 0.5*step*v1[0];
  xp2[1] = p[1] + 0.5*step*v1[1];
  xp2[2] = p[2] + 0.5*step*v1[2];

  k = locateTetra(mesh,k,++mesh->mark,xp2,cb);
  if ( !k )  return(0);
  norm = field3DInterp(mesh,k,cb,v2);
  pt = &mesh->tetra[k];
  pt->cpt--;

  xp3[0] = p[0] + step*v2[0];
  xp3[1] = p[1] + step*v2[1];
  xp3[2] = p[2] + step*v2[2];

  k = locateTetra(mesh,k,++mesh->mark,xp3,cb);
  if ( !k )  return(0);
  norm = field3DInterp(mesh,k,cb,v3);
  pt = &mesh->tetra[k];
  pt->cpt--;

  h6    = step / 6.0;
  p[0] += h6 * (v[0] + 2*(v1[0] + v2[0]) + v3[0]);
  p[1] += h6 * (v[1] + 2*(v1[1] + v2[1]) + v3[1]);
  p[2] += h6 * (v[2] + 2*(v1[2] + v2[2]) + v3[2]);

  return(1);
}

int nxtPoint2D(pMesh mesh,int nsdep,float *p,float step,double *v) {
  pTriangle  pt;
  double     norm,h6,cb[3],v1[3],v2[3],v3[3];
  float      xp1[3],xp2[3],xp3[3];
  int        k;

  /* 4th order Runge-Kutta */
  xp1[0] = p[0] + 0.5*step*v[0];
  xp1[1] = p[1] + 0.5*step*v[1];

  k = locateTria(mesh,nsdep,++mesh->mark,xp1,cb);
  if ( !k )  return(0);
  norm = field2DInterp(mesh,k,cb,v1);
  pt = &mesh->tria[k];
  pt->cpt--;
  
  xp2[0] = p[0] + 0.5*step*v1[0];
  xp2[1] = p[1] + 0.5*step*v1[1];

  k = locateTria(mesh,k,++mesh->mark,xp2,cb);
  if ( !k )  return(0);
  norm = field2DInterp(mesh,k,cb,v2);
  pt = &mesh->tria[k];
  pt->cpt--;

  xp3[0] = p[0] + step*v2[0];
  xp3[1] = p[1] + step*v2[1];

  k = locateTria(mesh,k,++mesh->mark,xp3,cb);
  if ( !k )  return(0);
  norm = field2DInterp(mesh,k,cb,v3);
  pt = &mesh->tria[k];
  pt->cpt--;

  h6 = step / 6.0;
  p[0] += h6 * (v[0] + 2*(v1[0] + v2[0]) + v3[0]);
  p[1] += h6 * (v[1] + 2*(v1[1] + v2[1]) + v3[1]);

  return(1);
}


/* read streamlines origins */
int parseStream(pScene sc,pMesh mesh) {
  FILE    *in;
  pStream  st;
  float    x,y,z;
  int      i,k,nbp,ret;
  char    *ptr,data[128],key[256],tmp[128];

  /* input file */
  strcpy(tmp,mesh->name);
  ptr = (char*)strstr(tmp,".mesh");
  if ( ptr )  *ptr = '\0';

  sprintf(data,"%s.iso",tmp);
  in = fopen(data,"r");
  if ( !in ) {
    sscanf(data,"DEFAULT.iso");
    in = fopen(data,"r");
    if ( !in )  return(0);
  }
  if ( !quiet )  fprintf(stdout,"  Reading %s\n",data);
  sc->stream = createStream(sc,mesh); 
  st = sc->stream;

  while ( !feof(in) ) {
    fscanf(in,"%s",key);
    for (i=0; i<strlen(key); i++) key[i] = tolower(key[i]);

    if ( !strcmp(key,"nblines") ) {
      fscanf(in,"%d",&nbp);
      st->nbstl = nbp;
      if ( mesh->dim == 3 )
        for (k=1; k<=3*st->nbstl; k+=3) {
          ret = fscanf(in,"%f %f %f\n",&x,&y,&z);
printf("x %f %f %f\n",x,y,z);
          st->listp[k]   = x - mesh->xtra;
          st->listp[k+1] = y - mesh->ytra;
          st->listp[k+2] = z - mesh->ztra;
printf("x %f %f %f\n",st->listp[k],st->listp[k+1],st->listp[k+2]);
        }
      else
        for (k=1; k<=2*st->nbstl; k+=2) {
          ret = fscanf(in,"%f %f\n",&x,&y);
          st->listp[k]   = x - mesh->xtra;
          st->listp[k+1] = y - mesh->ytra;
        }
    }
    else if ( !strcmp(key,"euler") ) {
      st->typtrack = Euler;
    }
    else if ( !strcmp(key,"box") ) {
      ret = fscanf(in,"%f %f",&x,&y);
      if ( ret != 2 )  break; 
      st->xmin = x - mesh->xtra;
      st->xmax = y - mesh->xtra;
      ret = fscanf(in,"%f %f",&x,&y);
      if ( ret != 2 )  break;
      st->ymin = x - mesh->ytra;
      st->ymax = y - mesh->ytra; 
      if ( mesh->dim == 3 ) {
        ret = fscanf(in,"%f %f",&x,&y);
        if ( ret != 2 )  break;
        st->zmin = x - mesh->ztra;
        st->zmax = y - mesh->ztra;
      }
    }
    else if ( key[0] == '#' ) {
      fgets(key,255,in);
    }
  }
  fclose(in);

  if ( !st->nbstl ) {
    fprintf(stderr,"   ## No data found.\n");
    return(0);
  }
k = 1;
printf("fin proc %f %f %f\n",sc->stream->listp[k],sc->stream->listp[k+1],sc->stream->listp[k+2]);

  return(1);
}


/* build lists for streamlines */
int listTetraStream(pScene sc,pMesh mesh,float *pp,int squiet) {
  pTetra     pt;
  pStream    st;
  double     dd,cb[4],cbdep[4],v[4],vdep[4],sizedep,normdep;
  float      step,p[3],ox,oy,oz,ldt;
  int        i,k,exh,depart,nsdep,nsfin,nsold,nbp,maxpts;
  clock_t    ct;
  FILE      *out;

  /* default */
  if ( !mesh->ntet )  return(0);
  else if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("\n create streamlines list / TETRA\n");

  if ( !squiet && !ddebug ) {
    fprintf(stdout," Building streamline(s)");
    fflush(stdout);
  }
  ct = clock();

  /* build display list */
  st = sc->stream;
  if ( st->nbstl > MAX_LST-1 )  return(0);
  sc->slist[st->nbstl] = glGenLists(1);
  glNewList(sc->slist[st->nbstl],GL_COMPILE);
  if ( glGetError() )  return(0);

  k = st->nbstl*3 + 1;
  st->listp[k+0] = pp[0];
  st->listp[k+1] = pp[1];
  st->listp[k+2] = pp[2];
  st->nbstl++;
printf("\n%d: pp = %f %f %f\n",st->nbstl,st->listp[k+0],st->listp[k+1],st->listp[k+2]);

  maxpts = MEDIT_MAX(MAX_PTS,5*mesh->ntet);
  glLineWidth(2.0);

  /* compute streamline */
  nbp   = 0;
  exh   = 0;
  nsdep = mesh->ntet / 2;
  step  = 0.0;
  nbar  = 0;
  if ( ddebug ) 
    printf("   start point %d: %f %f %f\n",
           3*k/3,st->listp[k],st->listp[k+1],st->listp[k+2]);

  for (i=1; i<mesh->ntet; i++) {
    pt = &mesh->tetra[i];
    pt->cpt = 0;
  }

  /* find enclosing tet */
  memcpy(p,pp,3*sizeof(float));
  depart = locateTetra(mesh,nsdep,++mesh->mark,p,cb);
printf("depart = %d\n",depart);
  if ( !depart ) {
    for (depart=1; depart<=mesh->ntet; depart++) {
      pt = &mesh->tetra[depart];
      if ( pt->mark != mesh->mark && inTetra(mesh,depart,p,cb) )
        break;
    }
    if ( depart > mesh->ntet ) {
      glEndList();
      return(0);
    }
  }

  st->norm = field3DInterp(mesh,depart,cb,v);
  memcpy(cbdep,cb,4*sizeof(double));
  memcpy(vdep,v,4*sizeof(double));
  st->size = sizeTetra(mesh,depart);
  sizedep  = st->size;
  normdep  = st->norm;
  ldt      = 0.0;

  if ( st->size == 0.0 )
    step = EPS*sc->dmax;
  else
    step = HSIZ * MEDIT_MIN(st->size,st->norm);

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  if ( sc->par.maxtime < FLT_MAX ) {
    sc->par.cumtim = 0.0;
    step = MEDIT_MIN(0.05*sc->par.dt,step);
    out  = fopen("particules.dat","a+");
    assert(out);
    fprintf(out,"\n%8.2f  %f %f %f\n",
            sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra,p[2]+mesh->ztra);
  }

  do {
    ox = p[0];
    oy = p[1];
    oz = p[2];

    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint3D(mesh,nsdep,p,step,v) ) {
      p[0] += step*v[0];
      p[1] += step*v[1];
      p[2] += step*v[2];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin || p[2]<st->zmin ||
         p[0]>st->xmax || p[1]>st->ymax || p[2]>st->zmax )
      break;
    else if ( sc->par.maxtime < FLT_MAX ) {
      ox -= p[0];
      oy -= p[1];
      oz -= p[2];
      dd  = sqrt(ox*ox + oy*oy + oz*oz) / st->norm;
      ldt += dd;
      sc->par.cumtim     += dd;
      if ( sc->par.cumtim  > sc->par.maxtime )  break;
      
      if ( ldt > sc->par.dt  ) {
        fprintf(out,"%8.2f  %f %f %f\n",
	        sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra,p[2]+mesh->ztra);
    	ldt = fabs(sc->par.dt - ldt);
      }
    }

    /* find tet containing p */
    nsfin = locateTetra(mesh,nsdep,++mesh->mark,p,cb);
    if ( !nsfin )  break;
    nsdep = nsfin;
    pt = &mesh->tetra[nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeTetra(mesh,nsdep);
      nsold = nsdep;
    }

    /* vector field interpolation */
    st->norm = field3DInterp(mesh,nsdep,cb,v);
    step     = HSIZ*MEDIT_MIN(st->size,st->norm);
    if ( sc->par.maxtime < FLT_MAX )
      step = MEDIT_MIN(0.05*sc->par.dt,step);
    if ( step == 0.0 )  break;

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  addPoint(sc,st,p,0);
  glEnd();

  if ( nbp >= maxpts || sc->par.maxtime < FLT_MAX ) {
    glLineWidth(1.0);
    glEndList();
    if ( !squiet && !ddebug ) {
      fprintf(stdout,": %d (%d, %.2f) / %d lines",
              nbar,nbp,(float)nbp/nbar,k/3);
      ct = difftime(clock(),ct);
      fprintf(stdout," %6.2f sec.\n",ct);
    }
    if ( sc->par.maxtime < FLT_MAX ) {
      fprintf(out,"%8.2f  %f %f %f\n",
              sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra,p[2]+mesh->ztra);
      fclose(out);
    }
    return(1);
  }

  /* reverse orientation */
  memcpy(p,pp,3*sizeof(float));
  memcpy(cb,cbdep,4*sizeof(double));
  memcpy(v,vdep,4*sizeof(double));
  st->norm = normdep;
  st->size = sizedep;
  if ( st->size == 0.0 )
    step = EPS * sc->dmax;
  else
    step = HSIZ * MEDIT_MIN(st->size,st->norm);

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  do {
    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint3D(mesh,nsdep,p,-step,v) ) {
      p[0] -= step*v[0];
      p[1] -= step*v[1];
      p[2] -= step*v[2];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin || p[2]<st->zmin ||
         p[0]>st->xmax || p[1]>st->ymax || p[2]>st->zmax )
      break;

    /* find tet containing p */
    nsfin = locateTetra(mesh,nsdep,++mesh->mark,p,cb);
    if ( !nsfin )  break;
    nsdep = nsfin;
    pt = &mesh->tetra[nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeTetra(mesh,nsdep);
      nsold = nsdep;
    }

    /* vector field interpolation */
    st->norm = field3DInterp(mesh,nsdep,cb,v);
    step = HSIZ * MEDIT_MIN(st->size,st->norm);
    if ( step == 0.0 )  break;

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  addPoint(sc,st,p,0);
  glEnd();
  glLineWidth(1.0);
  glEndList();

  if ( !nbp ) {
    st->nbstl--;
    if ( !squiet && !ddebug )  fprintf(stdout,"..empty\n");
    return(0);
  }

  if ( !squiet && !ddebug ) {  
    if ( nbar )
      fprintf(stdout,": %d (%d, %.2f) / %d lines",nbar,nbp,(float)nbp/nbar,k/3);
    ct = difftime(clock(),ct);
    fprintf(stdout," %6.2f sec.\n",ct);
  }

  return(1);
}


int listHexaStream(pScene sc,pMesh mesh,float *pp,int squiet) {
  pHexa      ph;
  pPoint     pt[4];
  pStream    st;
  double     cbdep[4],cb[4],v[6],vdep[6],sizedep,normdep;
  float      step,p[3];
  int        i,k,exh,depart,nsdep,nsfin,nsold,nbp,maxpts;
  mytime     tt;

  /* default */
  if ( !mesh->nhex )  return(0);
  else if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( ddebug ) printf("create streamlines list / HEXA\n");

  if ( !squiet && !ddebug ) {
    fprintf(stdout," Building streamline(s)");
    fflush(stdout);
    chrono(ON,&tt);
  }

  /* build display list */
  st = sc->stream;
  if ( st->nbstl > MAX_LST-1 )  return(0);
  sc->slist[st->nbstl] = glGenLists(1);
  glNewList(sc->slist[st->nbstl],GL_COMPILE);
  if ( glGetError() )  return(0);

  st->nbstl++;
  k = st->nbstl*3;
  st->listp[k]   = pp[0];
  st->listp[k+1] = pp[1];
  st->listp[k+2] = pp[2];

  maxpts = MEDIT_MAX(MAX_PTS,5*mesh->nhex);
  glLineWidth(2.0);

  /* compute streamline */
  nbp   = 0;
  exh   = 0;
  nsdep = mesh->nhex / 2;
  step  = 0.0f;
  nbar  = 0;
  if ( ddebug ) 
    printf("   start point %d: %f %f %f\n",
           3*k/3,st->listp[k],st->listp[k+1],st->listp[k+2]);

  for (i=1; i<mesh->nhex; i++) {
    ph = &mesh->hexa[i];
    ph->cpt  = 0;
  }

  /* find enclosing tet */
  memcpy(p,pp,3*sizeof(float));
  depart = locateHexa(mesh,nsdep,++mesh->mark,p,cb,pt);
printf("DEPART %d\n",depart);
ph = &mesh->hexa[depart];
printf("sommets %d %d %d %d %d %d %d %d\n",
ph->v[0],ph->v[1],ph->v[2],ph->v[3],ph->v[4],ph->v[5],ph->v[6],ph->v[7]);
  
  if ( !depart ) {
    for (depart=1; depart<=mesh->nhex; depart++) {
      ph = &mesh->hexa[depart];
      if ( ph->mark != mesh->mark && inHexa(mesh,depart,p,cb,pt) )
        break;
    }
    if ( depart > mesh->nhex )  return(0);
  }

  st->norm = vector3DInterp(mesh,pt,cb,v);
  memcpy(cbdep,cb,4*sizeof(double));
  memcpy(vdep,v,4*sizeof(double));
  st->size = sizeHexa(mesh,depart);
  sizedep  = st->size;
  normdep  = st->norm;
  if ( st->size == 0.0f )
    step = EPS*sc->dmax;
  else
    step = HSIZ * st->size;

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  do {
    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint3D(mesh,nsdep,p,step,v)) {
      p[0] += step*v[0];
      p[1] += step*v[1];
      p[2] += step*v[2];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin || p[2]<st->zmin ||
         p[0]>st->xmax || p[1]>st->ymax || p[2]>st->zmax )
      break;

    /* find tet containing p */
    nsfin = locateHexa(mesh,nsdep,++mesh->mark,p,cb,pt);
    if ( !nsfin )  break;
    nsdep = nsfin;
    ph = &mesh->hexa[nsdep];
    if ( ph->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeHexa(mesh,nsdep);
      step  = HSIZ*st->size;
      nsold = nsdep;
    }

    /* vector field interpolation */
    st->norm = vector3DInterp(mesh,pt,cb,v);
    if ( st->norm < EPS*step )  break;
    step = MEDIT_MIN(step,st->norm);
    if ( step == 0.0f )  break; /*step = 1.0e-06*sc->dmax;*/

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  glEnd();

  if ( nbp >= maxpts ) {
    glLineWidth(1.0);
    glEndList();
    if ( !squiet && !ddebug ) {
      fprintf(stdout,": %d (%d, %.2f) / %d lines",
              nbar,nbp,(float)nbp/nbar,k/3);
      chrono(OFF,&tt);
      fprintf(stdout," %6.2f sec.\n",gttime(tt));
    }
    return(1);
  }

  /* reverse orientation */
  memcpy(p,pp,3*sizeof(float));
  memcpy(cb,cbdep,4*sizeof(double));
  memcpy(v,vdep,4*sizeof(double));
  st->norm = normdep;
  st->size = sizedep;
  if ( st->size == 0.0f )
    step = EPS * sc->dmax;
  else
    step = HSIZ*st->size;

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  do {
    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint3D(mesh,nsdep,p,-step,v) ) {
      p[0] -= step*v[0];
      p[1] -= step*v[1];
      p[2] -= step*v[2];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin || p[2]<st->zmin ||
         p[0]>st->xmax || p[1]>st->ymax || p[2]>st->zmax )
      break;

    /* find tet containing p */
    nsfin = locateHexa(mesh,nsdep,++mesh->mark,p,cb,pt);
    if ( !nsfin )  break;
    nsdep = nsfin;
    ph = &mesh->hexa[nsdep];
    if ( ph->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeHexa(mesh,nsdep);
      step  = HSIZ * st->size;
      nsold = nsdep;
    }

    /* vector field interpolation */
    st->norm = vector3DInterp(mesh,pt,cb,v);
    if ( st->norm < EPS*step )   break;
    step = MEDIT_MIN(step,st->norm);
    if ( step == 0.0f )  break; /*step = 1.e-06 * sc->dmax;*/

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  glEnd();
  glLineWidth(1.0);
  glEndList();

  if ( nbp >= maxpts ) {
    glLineWidth(1.0);
    glEndList();
    if ( !squiet && !ddebug ) {
      fprintf(stdout,": %d (%d, %.2f) / %d lines",
              nbar,nbp,(float)nbp/nbar,k/3);
      chrono(OFF,&tt);
      fprintf(stdout," %6.2f sec.\n",gttime(tt));
    }
    return(1);
  }

  if ( !nbp ) {
    st->nbstl--;
    if ( !squiet && !ddebug )  fprintf(stdout,"..empty\n");
    return(0);
  }

  if ( !squiet && !ddebug ) {  
    if ( nbar )
      fprintf(stdout,": %d (%d, %.2f) / %d lines",nbar,nbp,(float)nbp/nbar,k/3);
    chrono(OFF,&tt);
    fprintf(stdout," %6.2f sec.\n",gttime(tt));
  }

  return(1);
}


int listTriaStream(pScene sc,pMesh mesh,float *pp) {
  pTriangle  pt;
  pStream    st;
  double     dd,cb[3],cbdep[3],v[3],vdep[3],sizedep,normdep;
  float      step,p[3],ox,oy,ldt;
  int        i,k,exh,depart,nsdep,nsfin,nsold,nbp,maxpts;
  clock_t    ct;
  FILE      *out;

  /* default */
  if ( !mesh->nt )  return(0);
  if ( ddebug ) printf("create streamlines / TRIA\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  fprintf(stdout," Building streamline(s)");
  fflush(stdout);
  ct = clock();

  /* build display list */
  st = sc->stream;
  if ( st->nbstl > MAX_LST-1 )  return(0);
  sc->slist[st->nbstl] = glGenLists(1);
  glNewList(sc->slist[st->nbstl],GL_COMPILE);
  if ( glGetError() )  return(0);

  st->nbstl++;
  k = st->nbstl*3;
  st->listp[k]   = pp[0];
  st->listp[k+1] = pp[1];
  st->listp[k+2] = pp[2];

  maxpts = MEDIT_MAX(MAX_PTS,5*mesh->nt);
  glLineWidth(2.0);

  /* compute streamlines */
  nbp   = 0;
  exh   = 0;
  nsdep = mesh->nt / 2;
  step  = 0.0;
  nbar  = 0;
  if ( ddebug )
    printf("   start point %d: %f %f\n",3*k/3,st->listp[k],st->listp[k+1]);

  for (i=1; i<=mesh->nt; i++) {
    pt = &mesh->tria[i];
    pt->cpt = 0;
  }

  /* find enclosing triangle */
  memcpy(p,pp,3*sizeof(float));
  depart = locateTria(mesh,nsdep,++mesh->mark,p,cb);
  if ( !depart ) {
    if ( ddebug )  printf("exhaustif search\n");
    for (depart=1; depart<=mesh->nt; depart++) {
      pt = &mesh->tria[depart];
      if ( pt->mark != mesh->mark && inTria(mesh,depart,p,cb) )
        break;
    }
    if ( depart > mesh->nt ) {
      st->nbstl--;
      glEndList();
      glLineWidth(1.0);
      fflush(stdout);
      return(0);
    }
  }

  st->norm = field2DInterp(mesh,depart,cb,v);
  memcpy(cbdep,cb,3*sizeof(double));
  memcpy(vdep,v,3*sizeof(double));
  st->size   = sizeTria(mesh,depart);
  sizedep    = st->size;
  normdep    = st->norm;
  ldt        = 0.0;

  if ( st->size == 0.0 )
    step = EPS * sc->dmax;
  else 
    step = HSIZ * MEDIT_MIN(st->size,st->norm);

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  if ( sc->par.maxtime < FLT_MAX ) {
    sc->par.cumtim = 0.0;
    step = MEDIT_MIN(0.05*sc->par.dt,step);
    out  = fopen("particules.dat","a+");
    assert(out);
    fprintf(out,"\n%8.2f  %f %f\n",
            sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra);
  }

  do {
    ox = p[0];
    oy = p[1];

    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint2D(mesh,nsdep,p,step,v) ) {
      p[0] += step*v[0];
      p[1] += step*v[1];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin ||
         p[0]>st->xmax || p[1]>st->ymax )
      break;
    else if ( sc->par.maxtime < FLT_MAX ) {
      ox -= p[0];
      oy -= p[1];
      dd  = sqrt(ox*ox + oy*oy) / st->norm;
      ldt += dd;
      sc->par.cumtim     += dd;
      if ( sc->par.cumtim >= sc->par.maxtime )  break;
      if ( ldt > sc->par.dt ) {
        fprintf(out,"%8.2f  %f %f\n",
                sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra);
        ldt = fabs(sc->par.dt - ldt);
      }
    }

    /* find tet containing p */
    nsfin = locateTria(mesh,nsdep,++mesh->mark,p,cb);
    if ( !nsfin )   break;
    nsdep = nsfin;
    pt = &mesh->tria[nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeTria(mesh,nsdep);
      nsold    = nsdep;
    }

    /* vector field interpolation */
    st->norm = field2DInterp(mesh,nsdep,cb,v);
    step     = HSIZ * MEDIT_MIN(st->size,st->norm);
    if ( sc->par.maxtime < FLT_MAX )
      step = MEDIT_MIN(0.05*sc->par.dt,step);
    if ( step == 0.0 )   break;

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  addPoint(sc,st,p,0);
  glEnd();

  if (  nbp >= maxpts || sc->par.maxtime < FLT_MAX ) {
    glLineWidth(1.0);
    glEndList();
    fprintf(stdout,": %d (%d, %.2f) / %d lines",nbar,nbp,(float)nbp/nbar,k/3);
    ct = difftime(clock(),ct);
    fprintf(stdout," %6.2f sec.\n",ct);
    if ( sc->par.maxtime < FLT_MAX ) {
      fprintf(out,"%8.2f  %f %f\n",
              sc->par.cumtim,p[0]+mesh->xtra,p[1]+mesh->ytra);
      fclose(out);
    }
    return(1);
  }

  /* reverse orientation */
  memcpy(p,pp,3*sizeof(float));
  memcpy(cb,cbdep,3*sizeof(double));
  memcpy(v,vdep,3*sizeof(double));
  st->norm = normdep;
  st->size = sizedep;
  if ( st->size == 0.0 )
    step = EPS * sc->dmax;
  else
    step = HSIZ* st->size;

  /* build display list incrementally */
  nsdep = nsold = depart;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,p,0);
  nbp++;

  do {
    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint2D(mesh,nsdep,p,-step,v) ) {
      p[0] -= step*v[0];
      p[1] -= step*v[1];
    }
    if ( p[0]<st->xmin || p[1]<st->ymin || p[0]>st->xmax || p[1]>st->ymax )
      break;

    /* find tet containing p */
    nsfin = locateTria(mesh,nsdep,++mesh->mark,p,cb);
    if ( !nsfin )  break;
    nsdep = nsfin;
    pt = &mesh->tria[nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeTria(mesh,nsdep);
      nsold    = nsdep;
    }

    /* vector field interpolation */
    st->norm = field2DInterp(mesh,nsdep,cb,v);
    step = HSIZ * MEDIT_MIN(st->size,st->norm);
    if ( step == 0.0 )  break;

    nbp += filterPoint(sc,st,p,0);
  }
  while ( nbp < maxpts );
  addPoint(sc,st,p,0);
  glEnd();
  glLineWidth(1.0);
  glEndList();

  if ( !nbp ) {
    st->nbstl--;
    fprintf(stdout,".. empty.\n");
    return(0);
  }

  /*if ( ddebug ) */
  if ( nbar )
    fprintf(stdout,": %d (%d, %.2f) / %d lines",nbar,nbp,(float)nbp/nbar,k/3);
  ct = difftime(clock(),ct);
  fprintf(stdout," %6.2f sec.\n",ct);

  return(1);
}


int listSaddleStream(pScene sc,pMesh mesh,int depart,
                     float *pp,float *vv,double lambda) {
  pTriangle  pt;
  pStream    st;
  double     cb[3],v[3];
  float      sens,step,p[3];
  int        i,k,nsdep,nsfin,nsold,nbp,maxpts;

  /* default */
  if ( !mesh->nt )  return(0);
  if ( ddebug ) printf("create streamlines for saddle point\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  /* build display list */
  st = sc->stream;
  if ( st->nbstl > MAX_LST-1 )  return(0);
  sc->slist[st->nbstl] = glGenLists(1);
  glNewList(sc->slist[st->nbstl],GL_COMPILE);
  if ( glGetError() )  return(0);
  maxpts = MEDIT_MAX(MAX_PTS,5*mesh->nt);
  glLineWidth(2.0);

  st->nbstl++;
  k = st->nbstl*3;
  st->listp[k]   = p[0] = pp[0];
  st->listp[k+1] = p[1] = pp[1];
  st->listp[k+2] = p[2] = 0.0f;
  
  for (i=1; i<=mesh->nt; i++) {
    pt = &mesh->tria[i];
    pt->cpt = 0;
  }
 
  /* compute streamlines */
  nsold = nsdep = depart;
  nbp   = nbar   = 0;
  glBegin(GL_LINE_STRIP);
  addPoint(sc,st,pp,1);

  st->size = sizeTria(mesh,depart);
  st->norm = sqrt(vv[0]*vv[0] + vv[1]*vv[1]);
  if ( st->size == 0.0f )
    step = EPS * sc->dmax;
  else
    step = HSIZ* st->size;
  if ( st->norm > 0.0f ) {
    v[0] = vv[0] / st->norm;
    v[1] = vv[1] / st->norm;
  }
  sens = lambda < 0.0f ? -1. : 1.;

  /* build display list incrementally */
  do {
    /* move to next point */
    if ( st->typtrack == Euler || !nxtPoint2D(mesh,nsdep,p,step,v) ) {
      p[0] += sens*step*v[0];
      p[1] += sens*step*v[1];
    }

    if ( p[0]<st->xmin || p[1]<st->ymin ||
         p[0]>st->xmax || p[1]>st->ymax )  break;

    /* find tet containing p */
    nsfin = locateTria(mesh,nsdep,++mesh->mark,p,cb);
    if ( !nsfin )  break;
    nsdep = nsfin;
    pt = &mesh->tria[nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( nsdep != nsold ) {
      st->size = sizeTria(mesh,nsdep);
      step  = HSIZ * st->size;
      nsold = nsdep;
    }

    /* vector field interpolation */
    st->norm = field2DInterp(mesh,nsdep,cb,v);
    if ( st->norm < EPS*step )  break;
    step = MEDIT_MIN(step,st->norm);
    if ( step == 0.0f )  break;

    nbp += filterPoint(sc,st,p,1);
  }
  while ( nbp < maxpts );
  glEnd();
  glLineWidth(1.0);
  glEndList();

  if ( !nbp ) {
    glDeleteLists(sc->slist[st->nbstl--],1);
    return(0);
  }

  return(1);
}


pStream createStream(pScene sc,pMesh mesh) {
  pStream    st;

  /* hash simplices */
  if ( mesh->dim == 2 ) {
    if ( mesh->nt && !hashTria(mesh) )    return(0);
  }
  else {
    if ( mesh->ntet && !hashTetra(mesh) ) return(0);
    if ( mesh->nhex && !hashHexa(mesh) )  return(0);
  }

  st = (pStream)calloc(1,sizeof(struct sstream));
  if ( !st )  return(0);

  /*st->typtrack = Euler;*/
  st->typtrack = RK4;
  st->stnp     = 0;
  st->nbstl    = 0;

  /* bounding box */
  st->xmin = mesh->xmin - mesh->xtra;
  st->ymin = mesh->ymin - mesh->ytra;
  st->zmin = mesh->zmin - mesh->ztra;
  st->xmax = mesh->xmax - mesh->xtra;
  st->ymax = mesh->ymax - mesh->ytra;
  st->zmax = mesh->zmax - mesh->ztra;

  /* init list */
  st->listp = (float*)malloc((MAX_LST*3+1)*sizeof(float));
  assert(st->listp);

  sc->slist = (GLuint*)calloc(MAX_LST,sizeof(GLuint));
  if ( !sc->slist )  return(0);

  return(st);
}

/* create from point picking */
int streamRefPoint(pScene sc,pMesh mesh) {
  pPoint    ppt;
  float     s[3];

  ppt = &mesh->point[refitem];
  if ( ppt->flag )  return(0);
  s[0] = ppt->c[0];
  s[1] = ppt->c[1];
  s[2] = ppt->c[2];
  if ( mesh->dim == 2 )
    listTriaStream(sc,mesh,s);
  else {
    if ( mesh->ntet )      listTetraStream(sc,mesh,s,0);
    else if ( mesh->nhex ) listHexaStream(sc,mesh,s,0);
  }
  ppt->flag = 1;

  return(1);
}


/* read starting point in file.iso */
int streamIsoPoint(pScene sc,pMesh mesh) {
  pStream   st;
  int       k,nbp,nbstl;
  time_t    t;

  if ( !parseStream(sc,mesh) )  return(0);
  t = clock();
  fprintf(stdout," Building streamline(s)");
  fflush(stdout);

  st    = sc->stream;
  nbstl = st->nbstl;
  nbp   = 0;
  st->nbstl = 0;
  if ( mesh->dim == 3 ) {
    if ( !mesh->ntet )  return(0);
    nbp = 0;
    for (k=1; k<=3*nbstl; k+=3) {
printf("\n ici: %f %f %f\n",st->listp[k],st->listp[k+1],st->listp[k+2]);
      nbp += listTetraStream(sc,mesh,&st->listp[k],1);
    }
  }

  if ( !nbp )    return(0);
  if ( ddebug )  printf("stream start: %d points  ",nbp);
  fprintf(stdout,": %d lines",nbp);
  t = clock() - t;
  fprintf(stdout," %6.2f sec.\n",t/(float)CLOCKS_PER_SEC);

  return(1);
}


int streamRefTria(pScene sc,pMesh mesh) {
  pMaterial   pm;
  pTriangle   pt;
  pPoint      ppt;
  float       s[3];
  int         i,k,nmat,nbp,base;
  time_t      t;

  /* build list */
  base = ++mesh->mark;
  nbp  = 0;
  pt   = &mesh->tria[refitem];

  nmat = matRef(sc,pt->ref);
  /*nmat = !pt->ref ? DEFAULT_MAT : 1+(pt->ref-1)%(sc->par.nbmat-1);*/
  pm = &sc->material[nmat];
  k  = pm->depmat[LTria];
  if ( !k || pm->flag )  return(0);

  t = clock();
  fprintf(stdout," Building streamline(s)");
  fflush(stdout);

  for (i=1; i<=mesh->np; i++) {
    ppt = &mesh->point[i];
    ppt->mark = base;
  }
  ++base;
  while ( k != 0 ) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) {
      k = pt->nxt;
      continue;
    }
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->mark != base ) {
        ppt->mark = base;
        s[0] = ppt->c[0];
        s[1] = ppt->c[1];
        s[2] = ppt->c[2];
        if ( ++nbp > MAX_LST-1 ) break;
        listTetraStream(sc,mesh,s,1);
        ppt->flag = 1;
      }
    }
    k = pt->nxt;
  }
  if ( !nbp )  return(0);
  if ( ddebug )  printf("stream start: %d points  ",nbp);
  fprintf(stdout,": %d lines",nbp);
  t = clock() - t;
  fprintf(stdout," %6.2f sec.\n",t/(float)CLOCKS_PER_SEC);

  return(1);
}

int streamRefQuad(pScene sc,pMesh mesh) {
  pMaterial   pm;
  pQuad       pq;
  pPoint      ppt;
  float       s[3];
  int         i,k,nmat,nbp,base;
  time_t      t;

  /* build list */
  base = ++mesh->mark;
  nbp  = 0;
  pq   = &mesh->quad[refitem];

  nmat = matRef(sc,pq->ref);
  /*nmat = !pt->ref ? DEFAULT_MAT : 1+(pt->ref-1)%(sc->par.nbmat-1);*/
  pm = &sc->material[nmat];
  k  = pm->depmat[LQuad];
  if ( !k || pm->flag )  return(0);

  t = clock();
  fprintf(stdout," Building streamline(s)");
  fflush(stdout);

  for (i=1; i<=mesh->np; i++) {
    ppt = &mesh->point[i];
    ppt->mark = base;
  }
  ++base;
  while ( k != 0 ) {
    pq = &mesh->quad[k];
    if ( !pq->v[0] ) {
      k = pq->nxt;
      continue;
    }
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pq->v[i]];
      if ( ppt->mark != base ) {
        ppt->mark = base;
        s[0] = ppt->c[0];
        s[1] = ppt->c[1];
        s[2] = ppt->c[2];
        if ( ++nbp > MAX_LST-1 ) break;
        listHexaStream(sc,mesh,s,1);
        ppt->flag = 1;
      }
    }
    k = pq->nxt;
  }
  if ( !nbp )  return(0);
  if ( ddebug )  printf("stream start: %d points  ",nbp);
  fprintf(stdout,": %d lines",nbp);
  t = clock() - t;
  fprintf(stdout," %6.2f sec.\n",t/(float)CLOCKS_PER_SEC);

  return(1);
}
