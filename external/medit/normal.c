#include "medit.h"
#include "extern.h"
#include "sproto.h"


GLuint drawNormals(pMesh mesh,pScene sc) {
  pTriangle   pt;
  pQuad       pq;
  pEdge       pr;
  pPoint      ppt;
  GLuint      dlist = 0;
  double      nn[3];
  float      *n,p[3],scal;
  int         i,k,kk,ki,iadr;

  /* default */
  if ( !mesh->nvn )  return(0);
  if ( ddebug ) printf("display normals, tangents\n");

  /* build display list */
  dlist = glGenLists(1);
  glNewList(dlist,GL_COMPILE);
  if ( glGetError() )  return(0);
  if ( sc->type & S_OPPOS )
    scal = -0.025;
  else
    scal =  0.025;
  scal *= sc->dmax;

  glColor3f(0.8,0.2,0.6);
  glLineWidth(2.0);
  glBegin(GL_LINES);

  /* normals at vertices */
  if ( mesh->extra->iv && mesh->extra->nv ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      ki  = mesh->extra->nv[k];
      if ( ki > 0 ) {
	iadr  = 3*(ki-1)+1;
        n     = &mesh->extra->n[ iadr ];
        p[0]  = ppt->c[0];
        p[1]  = ppt->c[1];
        p[2]  = ppt->c[2];
        nn[0] = n[0];
        nn[1] = n[1];
        nn[2] = n[2];
        drawVector3D(p,nn,scal);
      }
    }
  }

  /* normal at triangle vertices */
  if ( mesh->extra->it && mesh->extra->nt ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !pt->v[0] )  continue;
      for ( i=0; i<3; i++) {
        kk = 3*(k-1)+i+1;
        ki = mesh->extra->nt[kk];
	if ( ki > 0 ) {
	  ppt   = &mesh->point[pt->v[i]];
	  iadr  = 3*(ki-1)+1;
          n     = &mesh->extra->n[ iadr ];
          p[0]  = ppt->c[0];
          p[1]  = ppt->c[1];
          p[2]  = ppt->c[2];
          nn[0] = n[0];
          nn[1] = n[1];
          nn[2] = n[2];
          drawVector3D(p,nn,scal);
	}
      }
    }
  }

  /* normal at quad vertices */
  if ( mesh->extra->iq && mesh->extra->nq ) {
    for (k=1; k<=mesh->nq; k++) {
      pq = &mesh->quad[k];
      if ( !pq->v[0] )  continue;
      for ( i=0; i<4; i++) {
   	kk = 4*(k-1)+i+1;
        ki = mesh->extra->nq[kk];
	if ( ki > 0 ) {
	  ppt   = &mesh->point[pq->v[i]];
          iadr  = 3*(ki-1)+1;
	  n     = &mesh->extra->n[ iadr ];
          p[0]  = ppt->c[0];
          p[1]  = ppt->c[1];
          p[2]  = ppt->c[2];
          nn[0] = n[0];
          nn[1] = n[1];
          nn[2] = n[2];
          drawVector3D(p,nn,scal);
	}
      }
    }
  }


  /* tangents at vertices */
  if ( mesh->extra->jv && mesh->extra->tv ) {
    glColor3f(0.2,0.6,0.8);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      ki  = mesh->extra->tv[k];
      if ( ki > 0 ) {
	iadr  = 3*(ki-1)+1;
	n     = &mesh->extra->t[ iadr ];
        p[0]  = ppt->c[0];
        p[1]  = ppt->c[1];
        p[2]  = ppt->c[2];
        nn[0] = n[0];
        nn[1] = n[1];
        nn[2] = n[2];
        drawVector3D(p,nn,scal);
      }
    }
  }

  /* tangent at edge vertices */
  if ( mesh->extra->je && mesh->extra->te ) {
    for (k=1; k<=mesh->na; k++) {
      pr = &mesh->edge[k];
      if ( !pr->v[0] )  continue;
      for ( i=0; i<2; i++) {
        kk = 2*(k-1)+i+1;
        ki = mesh->extra->te[kk];
        if ( ki > 0 ) {
	       ppt   = &mesh->point[pr->v[i]];
          iadr  = 3*(ki-1)+1;
          n     = &mesh->extra->t[ iadr ];
          p[0]  = ppt->c[0];
          p[1]  = ppt->c[1];
          p[2]  = ppt->c[2];
          nn[0] = n[0];
          nn[1] = n[1];
          nn[2] = n[2];
          drawVector3D(p,nn,scal);
        }
      }
    }
  }

  glEnd();
  glLineWidth(1.0);
  glEndList();
  return(dlist);
}
