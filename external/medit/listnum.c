#include "medit.h"
#include "extern.h"
#include "sproto.h"


void updatePoints(pScene sc,pMesh mesh,int refmat) {
  pPoint    ppt;
  pTriangle pt;
  pQuad     pq;
  pTetra    ptet;
  pHexa     ph;
  pMaterial pm;
  int       i,k,m;

  /* mark all points */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tag = M_UNUSED;
  }

  /* unmark points */
  for (m=0; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    if ( pm->flag )  continue;

    /*triangles */
    k = pm->depmat[LTria];
    while ( k > 0 ) {
      pt = &mesh->tria[k];
      if ( pt->v[0] ) {
        for (i=0; i<3; i++) {
          ppt = &mesh->point[pt->v[i]];
          ppt->tag &= ~M_UNUSED;
        }
      }
      k = pt->nxt;
    }
    
    /* quads */
    k  = pm->depmat[LQuad];
    while ( k > 0 ) {
      pq = &mesh->quad[k];
      if ( pq->v[0] ) {
        for (i=0; i<4; i++) {
          ppt = &mesh->point[pq->v[i]];
          ppt->tag &= ~M_UNUSED;
        }
      }
      k = pq->nxt;
    }
    
    /* tetras */
    k = pm->depmat[LTets];
    while ( k > 0 ) {
      ptet = &mesh->tetra[k];
      if ( ptet->v[0] ) {
        for (i=0; i<4; i++) {
          ppt = &mesh->point[ptet->v[i]];
          ppt->tag &= ~M_UNUSED;
        }
      }
      k = ptet->nxt;
    }
    
    /* hexas */
    k = pm->depmat[LHexa];
    while ( k > 0 ) {
       ph = &mesh->hexa[k];
       if ( ph->v[0] ) {
         for ( i=0; i<8; i++) {
           ppt = &mesh->point[ph->v[i]];
           ppt->tag &= ~M_UNUSED;
         }
       }
       k = ph->nxt;
    }
  }
}

void listNum(pScene sc,pMesh mesh) {
  pMaterial  pm;
  pPoint     ppt;
  pTriangle  pt;
  pQuad      pq;
  double     cx,cy,cz;
  int        k,i,m;

  glDisable(GL_LIGHTING);
  /* point numbers */
  if ( sc->item & S_NUMP ) {

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( ppt->tag == M_UNUSED ) continue;
      if ( ppt->ref ) {
        m  = matRef(sc,ppt->ref);
        pm = &sc->material[m];
        glColor3f(pm->dif[0],pm->dif[1],pm->dif[2]);
      }
      else
        glColor3fv(sc->par.line);
      output3(ppt->c[0],ppt->c[1],ppt->c[2],"%d",k);
    }
  }

  /* facets numbers */
  if ( sc->item & S_NUMF ) {
    glColor4fv(sc->par.line);
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      if ( pm->flag )  continue;
    
      /*triangles */
      k = pm->depmat[LTria];
      while ( k > 0 ) {
        pt = &mesh->tria[k];
        if ( pt->v[0] ) {
          cx = cy = cz = 0.0f;
          for (i=0; i<3; i++) {
            ppt = &mesh->point[pt->v[i]];
            cx += 0.333 * ppt->c[0];
            cy += 0.333 * ppt->c[1];
            cz += 0.333 * ppt->c[2];
          }
          output3(cx,cy,cz,"%d",k);
        }
        k = pt->nxt;
      }
      
      /* quads */
      k = pm->depmat[LQuad];
      while ( k > 0 ) {
        pq = &mesh->quad[k];
        if ( pq->v[0] ) {
          cx = cy = cz = 0.0f;
          for (i=0; i<4; i++) {
            ppt = &mesh->point[pq->v[i]];
            cx += 0.25 * ppt->c[0];
            cy += 0.25 * ppt->c[1];
            cz += 0.25 * ppt->c[2];
          }
          output3(cx,cy,cz,"%d",k);
        }
        k = pq->nxt;
      }
    }
  }
}
