#include "medit.h"
#include "extern.h"
#include "sproto.h"
#include <time.h>

#define NC     32
#define NC2    (NC*NC)
#define NC3    (NC2*NC)

#define KA     31
#define KB     57
#define KC     79

static int ch[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, 
                        {1,2,6,5}, {2,3,7,6}, {0,3,7,4} };
static int idir[5]  = {0,1,2,0,1};
static int idirt[7] = {0,1,2,3,0,1,2};


/* very sioux! (09/2002) */
int hashTetra(pMesh mesh) {
  pTetra    pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int      *hcode,*link,inival,hsize;
  char     *hvoy;
  ubyte     i,ii,i1,i2,i3;
  unsigned  int  key;

  /* avoid building */
  if ( mesh->adja )  return(1);
  if ( 4*sizeof(char) != sizeof(int) )  exit(1);

  /* default */
  if ( ddebug) {
    fprintf(stdout,"  Setting topology.");
    fflush(stdout);
  }

  /* memory alloc */
  hcode = (int*)M_calloc(  MEDIT_MAX(100,mesh->ntet+1),sizeof(int),"hash.tetra");
  link  = (int*)M_calloc(4*MEDIT_MAX(100,mesh->ntet+1),sizeof(int),"hash.tetra");
  hsize = MEDIT_MAX(100,mesh->ntet);
  assert(hcode);
  assert(link);

  hvoy = (char*)hcode;

  /* init */
  inival = 2<<30;
  for (k=0; k<=mesh->ntet; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->ntet; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<4; i++) {
      i1 = idirt[i+1];
      i2 = idirt[i+2];
      i3 = idirt[i+3];
      mins = MEDIT_MIN(pt->v[i1],pt->v[i2]);
      mins = MEDIT_MIN(mins,pt->v[i3]);
      maxs = MEDIT_MAX(pt->v[i1],pt->v[i2]);
      maxs = MEDIT_MAX(maxs,pt->v[i3]);

      /* compute key */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = KA*mins + KB*maxs + KC*sum;
      key = key % hsize + 1;

      /* insert */
      iadr = 4*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  if ( ddebug ) {
    fprintf(stdout,".");
    fflush(stdout);
  }

  /* set adjacency */
  for (l=4*mesh->ntet; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 4 + 1;
    i = (l-1) % 4;
    i1 = idirt[i+1];
    i2 = idirt[i+2];
    i3 = idirt[i+3];
    pt = &mesh->tetra[k];

    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];
    mins = MEDIT_MIN(pt->v[i1],pt->v[i2]);
    mins = MEDIT_MIN(mins,pt->v[i3]);
    maxs = MEDIT_MAX(pt->v[i1],pt->v[i2]);
    maxs = MEDIT_MAX(maxs,pt->v[i3]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = idirt[ii+1];
      i2 = idirt[ii+2];
      i3 = idirt[ii+3];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = MEDIT_MIN(pt1->v[i1],pt1->v[i2]);
        mins1 = MEDIT_MIN(mins1,pt1->v[i3]);
        if ( mins1 == mins ) {
          maxs1 = MEDIT_MAX(pt1->v[i1],pt1->v[i2]);
          maxs1 = MEDIT_MAX(maxs1,pt1->v[i3]);
          if ( maxs1 == maxs ) {
            /* adjacent found */
            if ( pp != 0 )  link[pp] = link[ll];
            link[l] = kk;
            hvoy[l] = ii;
            link[ll]= k;
            hvoy[ll]= i;
            break;
          }
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  mesh->adja = (int*)link;
  mesh->voy  = (ubyte*)hcode;

  if ( ddebug )
    fprintf(stdout,"..\n");

  return(1);
}


/* very sioux! (09/2002) */
int hashHexa(pMesh mesh) {
  pHexa     ph,ph1;
  int       k,kk,iadr,pp,l,ll,v;
  int       imin,mins,mins1,opps,opps1;
  int      *hcode,*link,inival,hsize;
  char     *hvoy;
  ubyte     i,i1,ii;
  unsigned  int  key;

  /* avoid building again! */
  if ( mesh->adja )  return(1);
  if ( 4*sizeof(char) != sizeof(int) )  exit(1);

  /* default */
  if ( ddebug ) {
    fprintf(stdout,"  Setting topology.");
    fflush(stdout);
  }

  /* memory alloc */
/* bug fixe: 17/04/2007
  hcode = (int*)M_calloc(max(11,mesh->nhex+1),sizeof(int),"hash.hexa");
  link  = (int*)M_calloc(6*max(11,mesh->nhex+1),sizeof(int),"hash.hexa");
  hsize = MEDIT_MAX(10,mesh->nhex);
  if ( !hcode || !link ) {
    myerror.coderr = 1000;
    return(0);
  }
  hvoy = (char*)hcode;
*/
  hcode = (int*)M_calloc(MEDIT_MAX(10,6*mesh->nhex/4+1),sizeof(int),"hash.hexa");
  assert(hcode);
  link  = (int*)M_calloc(MEDIT_MAX(10,6*mesh->nhex+1),sizeof(int),"hash.hexa");
  assert(link);
  hsize = MEDIT_MAX(2,mesh->nhex);
  hvoy  = (char*)hcode;

  /* init */
  inival = 2 << 30;
  for (k=0; k<=6*mesh->nhex/4; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nhex; k++) {
    ph = &mesh->hexa[k];
    if ( !ph->v[0] )  continue;
    for (i=0; i<6; i++) {
      mins = ph->v[ch[i][0]];
      imin = 0;
      for (v=1; v<4; v++)
        if ( ph->v[ch[i][v]] < mins ) {
          mins = ph->v[ch[i][v]];
          imin = v;
        }
      i1   = (imin+2) % 4;
      opps = ph->v[ch[i][i1]];

      /* compute key */
      key = KA*mins + KB*opps;
      key = key % hsize + 1;

      /* insert */
      iadr = 6*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  if ( ddebug ) {
    fprintf(stdout,".");
    fflush(stdout);
  }

  /* set adjacency */
  for (l=6*mesh->nhex; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 6 + 1;
    i = (l-1) % 6;
    ph   = &mesh->hexa[k];
    mins = ph->v[ch[i][0]];
    imin = 0;
    for (v=1; v<4; v++)
      if ( ph->v[ch[i][v]] < mins ) {
        mins = ph->v[ch[i][v]];
        imin = v;
      }
    i1   = (imin+2) % 4;
    opps = ph->v[ch[i][i1]];

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 6 +1;
      ii = (ll-1) % 6;
      ph1   = &mesh->hexa[kk];
      mins1 = ph1->v[ch[ii][0]];
      imin  = 0;
      for (v=1; v<4; v++)
        if ( ph1->v[ch[ii][v]] < mins1 ) {
          mins1 = ph1->v[ch[ii][v]];
          imin  = v;
        }
      i1    = (imin+2) % 4;
      opps1 = ph1->v[ch[ii][i1]];

      /* adjacent found */
      if ( mins1 == mins && opps1 == opps ) {
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = kk;
        hvoy[l] = ii;
        link[ll]= k;
        hvoy[ll]= i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  mesh->adja = (int*)link;
  mesh->voy  = (ubyte*)hcode;

  if ( ddebug )
    fprintf(stdout,"..\n");

  return(1);
}


/* very sioux! (09/2002) */
int hashTria(pMesh mesh) {
  pTriangle pt,pt1;
  int       k,kk,l,ll,mins,maxs,mins1,maxs1,hsize;
  int      *hcode,*link,inival,iadr,pp;
  char     *hvoy;
  ubyte     i,i1,i2,ii;
  unsigned int key;

  /* avoid building again! */
  if ( mesh->adja )  return(1);
  if ( 4*sizeof(char) != sizeof(int) )  exit(1);

  /* default */
  if ( ddebug) {
    fprintf(stdout,"  Setting topology.");
    fflush(stdout);
  }

  /* memory alloc */
  hcode = (int*)M_calloc(MEDIT_MAX(1,3*mesh->nt/4)+1,sizeof(int),"hash.tria");
  link  = (int*)M_calloc(3*mesh->nt+1,sizeof(int),"hash.tria");
  hsize = MEDIT_MAX(2,3*mesh->nt/4-1);
  assert(hcode);
  assert(link);
  hvoy = (char*)hcode;

  /* init */
  inival = 2 << 30;
  for (k=0; k<=3*mesh->nt/4; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = idir[i+1];
      i2 = idir[i+2];
      mins = MEDIT_MIN(pt->v[i1],pt->v[i2]);
      maxs = MEDIT_MAX(pt->v[i1],pt->v[i2]);

      /* compute key */
      key = KA*mins + KB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }
  if ( ddebug ) {
    fprintf(stdout,".");
    fflush(stdout);
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = idir[i+1];
    i2 = idir[i+2];
    pt = &mesh->tria[k];

    mins = MEDIT_MIN(pt->v[i1],pt->v[i2]);
    maxs = MEDIT_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = idir[ii+1];
      i2 = idir[ii+2];
      pt1   = &mesh->tria[kk];
      mins1 = MEDIT_MIN(pt1->v[i1],pt1->v[i2]);
      maxs1 = MEDIT_MAX(pt1->v[i1],pt1->v[i2]);
      
      /* adjacent found */
      if ( mins1 == mins && maxs1 == maxs ) {
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = kk;
        hvoy[l] = ii;
        link[ll]= k;
        hvoy[ll]= i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  mesh->adja = (int*)link;
  mesh->voy  = (ubyte*)hcode;

  if ( ddebug )
    fprintf(stdout,".\n");

  return(1);
}
