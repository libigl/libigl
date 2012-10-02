#include "medit.h"
#include "extern.h"
#include "sproto.h"


/* display cube */
void updateCube(pCube cube,pMesh mesh) {
  pScene      sc;
  pTransform  cubetr,view;
  GLfloat     inv[16],axis[4],trans[4];
  int         idw;

  /* default */
  if ( ddebug ) printf("updateCube\n");

  /* retrieve context */
  idw    = currentScene();
  sc     = cv.scene[idw];
  view   = sc->view;
  cubetr = cube->cubetr;

  /* compute cube transform */
  if ( cube->active & C_EDIT ) {
    invertMatrix(view->matrix,inv);
    inv[3] = inv[7] = inv[11] = 0.0f;  inv[15] = 1.0;
 
    /* rotation cumulative */
    if ( cubetr->angle != 0.0 ) {
      transformVector(trans,cubetr->axis,inv);
      glPushMatrix();
      glLoadIdentity();
      glRotatef(cubetr->angle,trans[0],trans[1],trans[2]);
      glMultMatrixf(cubetr->rot);
      glGetFloatv(GL_MODELVIEW_MATRIX,cubetr->rot);
      glPopMatrix();
    }

    /* translation cumulative */
    axis[0] = cubetr->panx;
    axis[1] = cubetr->pany;
    axis[2] = 0.0;
    axis[3] = 1.0;
    transformVector(trans,axis,inv);

    cubetr->tra[12] = trans[0];
    cubetr->tra[13] = trans[1];
    cubetr->tra[14] = trans[2];

   /* final transformation */
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixf(cubetr->tra);
    glMultMatrixf(cubetr->rot);
    glGetFloatv(GL_MODELVIEW_MATRIX,cubetr->matrix);
    glPopMatrix();
  }
  if ( !cubetr->manim ) {
    cubetr->angle = 0.0;
    cube->active ^= C_UPDATE;
  }
}


void dumpCube(pScene sc,pMesh mesh,pCube cube) {
  float  *tr,u[4];
  double  v[4];
  int     i;
  FILE   *out;

  tr = cube->cubetr->matrix;

  out = fopen("tr.data","w");
  for (i=0; i<4; i++)
    fprintf(out,"%f %f %f %f\n",tr[i],tr[4+i],tr[8+i],tr[12+i]);
  
  u[0] = cube->cmi[0] - mesh->xtra;
  u[1] = cube->cmi[1] - mesh->ytra;
  u[2] = cube->cmi[2] - mesh->ztra;
  u[3] = 1.0;
  /*printf("avant %f %f %f %f\n",u[0],u[1],u[2],u[3]);*/
  transformPoint2(v,u,tr);
  fprintf(out,"\n%f %f %f   %f\n",
          v[0]+mesh->xtra,v[1]+mesh->ytra,v[2]+mesh->ztra,v[3]);

  u[0] = cube->cma[0] - mesh->xtra;
  u[1] = cube->cma[1] - mesh->ytra;
  u[2] = cube->cma[2] - mesh->ztra;
  /*printf("avant %f %f %f %f\n",u[0],u[1],u[2],u[3]);*/
  transformPoint2(v,u,tr);
  fprintf(out,"%f %f %f   %f\n",
          v[0]+mesh->xtra,v[1]+mesh->ytra,v[2]+mesh->ztra,v[3]);

  fprintf(out,"\n%f %f %f\n",
          cube->cubetr->tra[12],cube->cubetr->tra[13],cube->cubetr->tra[14]);

  fprintf(out,"%f   %f %f %f\n",
    cube->cubetr->angle,cube->cubetr->axis[0],cube->cubetr->axis[1],
    cube->cubetr->axis[2]);
  fclose(out);
}


void resetCube(pScene sc,pCube cube,pMesh mesh) {

  resetTransform(cube->cubetr);

  cube->active |= C_REDO;
  cube->cmi[0]  = mesh->xmin;
  cube->cmi[1]  = mesh->ymin;
  cube->cmi[2]  = mesh->zmin;
  cube->cma[0]  = mesh->xmax;
  cube->cma[1]  = mesh->ymax;
  cube->cma[2]  = mesh->zmax;
}


pCube createCube(pScene sc,pMesh mesh) {
  pCube   cube;

  cube = (pCube)M_calloc(1,sizeof(struct cube),"cube");
  assert(cube);
  
  cube->cubetr = (pTransform)M_calloc(1,sizeof(struct transform),"cube");
  if ( !cube->cubetr )  return(0);

  resetCube(sc,cube,mesh);
  return(cube);	
}
