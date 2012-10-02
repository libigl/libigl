#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define MAX_MORPH   100

int     imstep,imreverse;

int modeMorphing() {
  pScene   scene;
  pMesh    mesh1,mesh2;
  pPoint   ppt1,ppt2;
  int      k;
  clock_t  ct;

  /* default */
  cv.nbs    = 1;
  imstep    = 1;
  imreverse = 1;
  if ( ddebug ) printf("morphing: create window\n");
  fprintf(stdout,"\n Building scene\n");
  ct = clock();

  /* create grafix */
  scene = cv.scene[0];
  mesh1 = cv.mesh[0];
  iniopt(scene,mesh1);
  parsop(scene,mesh1);
  meshRef(scene,mesh1);
  matSort(scene);

  mesh2 = cv.mesh[1];
  parsop(scene,mesh2);
  meshRef(scene,mesh2);
  for (k=1; k<=mesh2->np; k++) {
    ppt1 = &mesh1->point[k]; 
    ppt2 = &mesh2->point[k];
    ppt2->c[0] -= ppt1->c[0]; 
    ppt2->c[1] -= ppt1->c[1];
    ppt2->c[2] -= ppt1->c[2];
  }

  if ( !createScene(scene,0) ) {
    fprintf(stderr,"  ## Unable to create scene\n");
    return(0);
  }
  ct = difftime(clock(),ct);
  fprintf(stdout,"  Scene seconds:      %.2f\n",
          (double)ct/(double)CLOCKS_PER_SEC);
  return(1);
}


int morphMesh(pScene sc,pMesh mesh1) {
  pMesh       mesh2;
  pPoint      ppt1,ppt2;
  int         k;
  static  float  dt = 1.0 / MAX_MORPH;

  imstep++;
  if ( imstep == 0 ) {
    imstep = 2;
    dt = -dt;
    if ( !imreverse ) {
      glutIdleFunc(NULL);
      morphing = 0;
      return(0);
    }
  }
  else if ( imstep == MAX_MORPH+1 ) {
    imstep = -MAX_MORPH+1;
    dt = -dt;
    if ( !imreverse ) {
      glutIdleFunc(NULL);
      morphing = 0;
      return(0);
    }
  }

  mesh2 = cv.mesh[1];

  for (k=1; k<=mesh1->np; k++) {
    ppt1 = &mesh1->point[k];
    ppt2 = &mesh2->point[k];
    ppt1->c[0] += dt*ppt2->c[0]; 
    ppt1->c[1] += dt*ppt2->c[1]; 
    ppt1->c[2] += dt*ppt2->c[2]; 
  }

  doLists(sc,mesh1);
  return(1);
}
