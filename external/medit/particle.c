#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define MAX_PRT    10
#define MAX_LST    1024
#define MAX_CPT    5000
#define HSIZ       0.03

extern int     reftype,refitem;


typedef struct spart {
  double    cb[4];
  float     pos[MAX_PRT+1][3],col[MAX_PRT+1][3],size,norm,ct,step;
  int       nsdep,cur;
  ubyte     flag;
} Particle;
typedef Particle *pParticle;

Particle  *tp;


void colorParticle(pScene sc,pParticle pp) {
  double        rgb[3],norm,kc;
  int           i;
  static double hsv[3] = { 0.0, 1.0, 0.80 };

  norm = pp->norm;
  if ( norm < sc->iso.val[0] )
    norm = sc->iso.val[0];  
  else if ( norm > sc->iso.val[MAXISO-1] )
    norm = sc->iso.val[MAXISO-1];
  for (i=0; i<MAXISO-1; i++)
    if ( norm < sc->iso.val[i] )  break;
  kc = (norm-sc->iso.val[i-1]) / (sc->iso.val[i] - sc->iso.val[i-1]);
  hsv[0] = sc->iso.col[i-1]*(1.0-kc)+sc->iso.col[i]*kc;

  hsvrgb(hsv,rgb);
  pp->col[pp->cur][0] = rgb[0];
  pp->col[pp->cur][1] = rgb[1];
  pp->col[pp->cur][2] = rgb[2];
}


void drawParticle(pScene sc,pParticle pp) {
  float         radius;
  int           i;

  radius = 0.005*sc->dmax;

  for (i=1; i<=sc->par.nbpart; i++) {
    glPushMatrix();
    glTranslatef(pp->pos[i][0],pp->pos[i][1],pp->pos[i][2]);
    glColor3fv(pp->col[i]);
    glutSolidSphere(radius,10,10);
    glPopMatrix();
  }
}


void computeTetraParticle(pScene sc,pMesh mesh,int k) {
  pTetra       pt;
  pStream      st;
  pParticle    pp;
  double       dd,cb[4],v[4];
  float        ux,uy,uz,pos[3],ldt;
  int          cur,nsfin,nsold,nbp;

  st   = sc->stream;
  pp   = &tp[k];
  if ( pp->ct > sc->par.pertime )  return;

  ldt   = 0.0;
  nbp   = 0;
  nsold = pp->nsdep;

  pos[0] = pp->pos[pp->cur][0];
  pos[1] = pp->pos[pp->cur][1];
  pos[2] = pp->pos[pp->cur][2];

  do {
    ux = pos[0];
    uy = pos[1];
    uz = pos[2];
    if ( st->typtrack == 1 || !nxtPoint3D(mesh,pp->nsdep,pos,pp->step,v) ) {
      pos[0] += pp->step*v[0];
      pos[1] += pp->step*v[1];
      pos[2] += pp->step*v[2];
    }
    if ( pos[0]<st->xmin || pos[0]>st->xmax ||
         pos[1]<st->ymin || pos[1]>st->ymax ||
	     pos[2]<st->zmin || pos[2]>st->zmax ) {
      pp->flag = 0;
      break;
    }
    ux -= pos[0];
    uy -= pos[1];
    uz -= pos[2];

    dd      = sqrt(ux*ux + uy*uy + uz*uz) / pp->norm;
    ldt    += dd;
    pp->ct += dd;

    if ( pp->ct >= sc->par.pertime ) {
      pp->flag = 0;
      sc->par.cumpertime = sc->par.pertime+1.e-06;
      break; /*return;*/
    }
    else if ( ldt >= sc->par.dt  ) {
      pp->cur = 1 + (pp->cur % sc->par.nbpart);
      pp->pos[pp->cur][0] = pos[0];
      pp->pos[pp->cur][1] = pos[1];
      pp->pos[pp->cur][2] = pos[2];
      colorParticle(sc,pp);
      break; /*return;*/
    }
 
    /* find tet containing p */
    nsfin = locateTetra(mesh,pp->nsdep,++mesh->mark,pos,cb);
    if ( !nsfin ) {
      pp->flag = 0;
      break; /*return;*/
    }
    pp->nsdep = nsfin;
    pt = &mesh->tetra[pp->nsdep];
    if ( pt->cpt > MAX_CPT )  break;

    /* adjust local stepsize */
    if ( pp->nsdep != nsold ) {
      pp->size = sizeTetra(mesh,pp->nsdep);
      nsold    = pp->nsdep;
    }

    /* vector field interpolation */
    pp->norm = field3DInterp(mesh,pp->nsdep,pp->cb,v);
    pp->step = HSIZ*MEDIT_MIN(pp->size,pp->norm);
    if ( sc->par.maxtime < FLT_MAX )
      pp->step = MEDIT_MIN(0.05*sc->par.dt,pp->step);
    if ( pp->step == 0.0 ) {
      pp->flag = 0;
      return;
    }
    nbp++;
  }
  while ( ldt <= sc->par.dt );

  cur = (pp->cur % MAX_PRT) + 1;
  pp->pos[cur][0] = pp->pos[pp->cur][0];
  pp->pos[cur][1] = pp->pos[pp->cur][1];
  pp->pos[cur][2] = pp->pos[pp->cur][2];
  pp->cur = cur;
  colorParticle(sc,pp);
}


int displayParticle(pScene sc,pMesh mesh) {
  pParticle   pp;
  pStream     st;
  int         k;

  st = sc->stream;
  if ( sc->par.advtim )
    for (k=1; k<=st->nbstl; k++)
      computeTetraParticle(sc,mesh,k);

  /* redraw particles */
  for (k=1; k<=st->nbstl; k++) {
    pp = &tp[k];
    if ( pp->flag )
      drawParticle(sc,pp);
  }

  return(1);
}


int createParticle(pScene sc,pMesh mesh) {
  pParticle   pp;
  pStream     st;
  pMaterial   pm;
  pTetra      pt1;
  pTriangle   pt;
  pPoint      ppt;
  double      v[4],cx,cy,cz;
  int         i,j,k,l,nmat,nbp,base;

  if ( ddebug )  printf("create particles\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);

  if ( sc->stream ) { 
    st = sc->stream;
    if ( st->listp )  free(st->listp);
    free(sc->stream);
  }
  sc->stream     = createStream(sc,mesh);
  sc->par.cumtim = 0.0;
  sc->par.advtim = 0;
  assert(sc->stream);
  st = sc->stream;

  fprintf(stdout," Creating particles :");  fflush(stdout);
  st->nbstl = 0;

  base = ++mesh->mark;
  pt   = &mesh->tria[refitem];
  nmat = matRef(sc,pt->ref);
  pm   = &sc->material[nmat];
  k    = pm->depmat[LTria];
  if ( !k || pm->flag )  return(0);

  /* point positions */
  for (i=1; i<=mesh->np; i++) {
    ppt = &mesh->point[i];
    ppt->mark = base;
  }
  if ( sc->par.nbpart >= MAX_PRT )
    sc->par.nbpart = MAX_PRT;

  ++base;
  l   = 1;
  nbp = 0;
  while ( k != 0 && st->nbstl < MAX_LST-1 ) {
    pt = &mesh->tria[k];
    if ( pt->v[0] ) {
      cx = cy = cz = 0.0;
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt->v[i]];
        cx += 0.33 * ppt->c[0];
        cy += 0.33 * ppt->c[1];
        cz += 0.33 * ppt->c[2];
        ppt->flag = 1;
      }  
      st->listp[l++] = cx;
      st->listp[l++] = cy;
      st->listp[l++] = cz;
      ++st->nbstl;
      if ( st->nbstl > MAX_LST-1 ) break;
      k = pt->nxt;
    }
  }
  fprintf(stdout,"%d\n",st->nbstl);
  if ( !st->nbstl )  return(0);

  /* init positions */
#ifdef IGL
  tp = static_cast<pParticle>(calloc((st->nbstl+1),sizeof(Particle)));
#else
  tp = calloc((st->nbstl+1),sizeof(Particle));
#endif
  assert(tp);
  for (k=1; k<=st->nbstl; k++)
    tp[k].nsdep = mesh->ntet / 2;

  for (k=1; k<=mesh->ntet; k++)
    mesh->tetra[k].cpt = 0;

  l = 1;
  for (k=1; k<=st->nbstl; k++) {
    pp = &tp[k];
    pp->pos[1][0] = st->listp[l++];
    pp->pos[1][1] = st->listp[l++];
    pp->pos[1][2] = st->listp[l++];

    tp[k].nsdep = locateTetra(mesh,pp->nsdep,++mesh->mark,pp->pos[1],pp->cb);

    if ( !pp->nsdep ) {
      for (j=1; j<=mesh->ntet; j++) {
        pt1 = &mesh->tetra[j];
        if ( pt1->mark != mesh->mark && inTetra(mesh,j,pp->pos[1],pp->cb) )
          break;
      }
      if ( j > mesh->ntet )  return(0);
      else
        pp->nsdep = j;
    }
    pp->norm  = field3DInterp(mesh,tp[k].nsdep,tp[k].cb,v);
    pp->size  = sizeTetra(mesh,tp[k].nsdep);
    if ( pp->size == 0.0 )
      pp->step = EPS*sc->dmax;
    else
      pp->step = HSIZ * MEDIT_MIN(pp->size,pp->norm);
    pp->step = MEDIT_MIN(0.05*sc->par.dt,pp->step);
    pp->flag =  1;
    pp->cur  =  1;
    colorParticle(sc,pp);

    for (i=2; i<=sc->par.nbpart; i++) {
      pp->pos[i][0] = pp->pos[1][0];
      pp->pos[i][1] = pp->pos[1][1];
      pp->pos[i][2] = pp->pos[1][2];
    }
  }

  return(1);
}


int advectParticle(pScene sc,pMesh mesh) {
  pParticle   pp;
  pStream     st;
  pTetra      pt1;
  pPoint      ppt;
  double      v[4];
  int         i,j,k,l,base;

  if ( ddebug )  printf("advect particles\n");
  if ( egal(sc->iso.val[0],sc->iso.val[MAXISO-1]) )  return(0);
  if ( !sc->stream )  return(0);
  st = sc->stream;

  if ( mesh->ntet && !hashTetra(mesh) ) return(0);

  st->xmin = mesh->xmin - mesh->xtra;
  st->ymin = mesh->ymin - mesh->ytra;
  st->zmin = mesh->zmin - mesh->ztra;
  st->xmax = mesh->xmax - mesh->xtra;
  st->ymax = mesh->ymax - mesh->ytra;
  st->zmax = mesh->zmax - mesh->ztra;

  sc->par.cumpertime = 0.0;
  sc->par.advtim     = 0;

  /* init list */
  sc->slist = (GLuint*)calloc(MAX_LST,sizeof(GLuint));
  if ( !sc->slist )  return(0);

  /* point positions */
  base = ++mesh->mark;
  for (i=1; i<=mesh->np; i++) {
    ppt = &mesh->point[i];
    ppt->mark = base;
  }
  for (k=1; k<=mesh->ntet; k++)
    mesh->tetra[k].cpt = 0;

  l = 1;
  for (k=1; k<=st->nbstl; k++) {
    pp = &tp[k];
    pp->ct = 0.0;
    st->listp[l++] = pp->pos[pp->cur][0];
    st->listp[l++] = pp->pos[pp->cur][1];
    st->listp[l++] = pp->pos[pp->cur][2];
    pp->cur = 1;
  }

  ++base;
  l = 1;
  for (k=1; k<=st->nbstl; k++) {
    pp = &tp[k];
    pp->pos[1][0] = st->listp[l++];
    pp->pos[1][1] = st->listp[l++];
    pp->pos[1][2] = st->listp[l++];

    tp[k].nsdep = locateTetra(mesh,pp->nsdep,++mesh->mark,pp->pos[1],pp->cb);
    if ( !pp->nsdep ) {
      for (j=1; j<=mesh->ntet; j++) {
        pt1 = &mesh->tetra[j];
        if ( pt1->mark != mesh->mark && inTetra(mesh,j,pp->pos[1],pp->cb) )
          break;
      }
      if ( j > mesh->ntet )  continue;
      else
        pp->nsdep = j;
    }
    pp->norm  = field3DInterp(mesh,tp[k].nsdep,tp[k].cb,v);
    pp->size  = sizeTetra(mesh,tp[k].nsdep);
    if ( pp->size == 0.0 )
      pp->step = EPS*sc->dmax;
    else
      pp->step = HSIZ * MEDIT_MIN(pp->size,pp->norm);
    pp->step = MEDIT_MIN(0.05*sc->par.dt,pp->step);
    pp->flag =  1;
    pp->cur  =  1;
    colorParticle(sc,pp);

    for (i=2; i<=sc->par.nbpart; i++) {
      pp->pos[i][0] = pp->pos[1][0];
      pp->pos[i][1] = pp->pos[1][1];
      pp->pos[i][2] = pp->pos[1][2];
    }
  }

  return(1);
}

