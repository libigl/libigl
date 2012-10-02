#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "medit.h"
#include "extern.h"
#include "sproto.h"


/* compute circumradius and center */
int cenrad(pMesh mesh,int iel,double *c,double *rad) {
  pTetra      pt;
  pPoint      p1,p2,p3,p4;
  double      dd,ux,uy,uz,n1[3],n2[3],n3[3],c1,c2,c3,pl1,pl2,pl3;

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return(0);

  p1 = &mesh->point[pt->v[0]];
  p2 = &mesh->point[pt->v[1]];
  p3 = &mesh->point[pt->v[2]];
  p4 = &mesh->point[pt->v[3]];
  
  ux = p4->c[0] - p1->c[0];
  uy = p4->c[1] - p1->c[1];
  uz = p4->c[2] - p1->c[2];
  dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
  n1[0] = ux*dd;
  n1[1] = uy*dd;
  n1[2] = uz*dd;
  /* plan: vecteur directeur passant par milieu(1,4) */
  pl1 = n1[0]*(p4->c[0]+p1->c[0]) \
      + n1[1]*(p4->c[1]+p1->c[1]) + n1[2]*(p4->c[2]+p1->c[2]);

  ux = p4->c[0] - p2->c[0];
  uy = p4->c[1] - p2->c[1];
  uz = p4->c[2] - p2->c[2];
  dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
  n2[0] = ux*dd;
  n2[1] = uy*dd;
  n2[2] = uz*dd;
  pl2 = n2[0]*(p4->c[0]+p2->c[0]) \
      + n2[1]*(p4->c[1]+p2->c[1]) + n2[2]*(p4->c[2]+p2->c[2]);

  ux = p4->c[0] - p3->c[0];
  uy = p4->c[1] - p3->c[1];
  uz = p4->c[2] - p3->c[2];
  dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
  n3[0] = ux*dd;
  n3[1] = uy*dd;
  n3[2] = uz*dd;
  pl3 = n3[0]*(p4->c[0]+p3->c[0]) \
      + n3[1]*(p4->c[1]+p3->c[1]) + n3[2]*(p4->c[2]+p3->c[2]);

  /* center = intersection of 3 planes */
  ux = n2[1]*n3[2] - n2[2]*n3[1];
  uy = n1[2]*n3[1] - n1[1]*n3[2];
  uz = n1[1]*n2[2] - n1[2]*n2[1];

  dd = n1[0]*ux + n2[0]*uy + n3[0]*uz;
  dd = 0.5 / dd;

  c1 = ux*pl1 + uy*pl2 + uz*pl3;
  c2 = pl1 * (n2[2]*n3[0] - n2[0]*n3[2]) \
     + pl2 * (n1[0]*n3[2] - n3[0]*n1[2]) \
     + pl3 * (n2[0]*n1[2] - n2[2]*n1[0]);
  c3 = pl1 * (n2[0]*n3[1] - n2[1]*n3[0]) \
     + pl2 * (n3[0]*n1[1] - n3[1]*n1[0]) \
     + pl3 * (n1[0]*n2[1] - n2[0]*n1[1]);

  c[0] = dd * c1;
  c[1] = dd * c2;
  c[2] = dd * c3;

  /* radius (squared) */
  *rad = (c[0] - p4->c[0]) * (c[0] - p4->c[0]) \
       + (c[1] - p4->c[1]) * (c[1] - p4->c[1]) \
       + (c[2] - p4->c[2]) * (c[2] - p4->c[2]);

  return(1);
}
