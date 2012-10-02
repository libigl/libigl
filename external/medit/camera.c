#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "extern.h"
#include "sproto.h"


double Azimuth(pCamera c) {
  double  dd,azim,cosazim,sinazim;

  dd      = sqrt((double)c->speed[0]*c->speed[0]+(double)c->speed[2]*c->speed[2]);
  cosazim = c->speed[2] / dd;
  sinazim = c->speed[0] / dd;
  azim = atan2(sinazim,cosazim) * RTOD;
  return(azim);
}

double Elevation(pCamera c) {
  double  elev;

  elev  = sqrt((double)c->speed[0]*c->speed[0]+(double)c->speed[2]*c->speed[2]);
  elev  = atan2(c->speed[1],elev) * RTOD;
  return(elev);
}

/* compute new sun position */
void updateSun(pScene sc,pCamera c) {
  double    dd;
  GLfloat   axe[3],sunf[4];
  GLdouble  speed[3],sunp[4],matrix[16];

  axe[0] = c->speed[2];
  axe[1] = 0.0f;
  axe[2] = -c->speed[0];
  dd = sqrt(axe[0]*axe[0] + axe[2]*axe[2]);
  if ( dd != 0.0f ) {
    axe[0] /= dd;
    axe[2] /= dd;
  }
  
  speed[0] = c->speed[0];
  speed[1] = c->speed[1];
  speed[2] = c->speed[2];
  
  glPushMatrix();
  glLoadIdentity();
    glRotatef(-30.0f,axe[0],axe[1],axe[2]);
    glGetDoublev(GL_MODELVIEW_MATRIX,matrix);
  glPopMatrix();
  transformPointd(sunp,speed,matrix);
  sunf[0] = -sc->dmax*sunp[0];
  sunf[1] = -sc->dmax*sunp[1];
  sunf[2] = -sc->dmax*sunp[2];
  sunf[3] = 0.0;

  glLightfv(GL_LIGHT0,GL_POSITION,sunf);

  if ( ddebug ) {
    printf("    speed  %g %g %g\n",c->speed[0],c->speed[1],c->speed[2]);
    printf("    axe    %g %g %g\n",axe[0],axe[1],axe[2]);
    printf("    sunpos %g %g %g\n",sunp[0],sunp[1],sunp[2]);
  }
}


void updateCamera(pScene sc,pCamera c,double azim,double elev) {
  double     d,lazim,lelev;

  /* compute speed vector */
  if ( elev > 89.0f )       elev =  89.0;
  else if ( elev < -89.0f ) elev = -89.0;
  lazim = azim * DTOR;
  lelev = elev * DTOR;
  c->speed[0] = sin(lazim)*cos(lelev);
  c->speed[1] = sin(lelev);
  c->speed[2] = cos(lazim)*cos(lelev);

  d = (double)c->speed[0]*sc->par.sunpos[0] + c->speed[1]*sc->par.sunpos[1] +
              c->speed[2]*sc->par.sunpos[2];
  d = d / sqrt((double)sc->par.sunpos[0]*sc->par.sunpos[0]+
               sc->par.sunpos[1]*sc->par.sunpos[1]+
               sc->par.sunpos[2]*sc->par.sunpos[2]);
  d = acos(d);
  if ( fabs(d) > 0.10*sc->persp->fovy*DTOR )
    updateSun(sc,c);
}


pCamera initCamera(pScene sc,int up) {
  pCamera  c;
  pMesh    mesh;
  double   dd;
  double   look[3];

  if ( ddebug ) printf("    initCamera dmax %g\n",sc->dmax);
  if ( sc->camera ) c = sc->camera;
  else {
    c = (pCamera)M_calloc(1,sizeof(struct camera),"camera");
    if ( !c ) {
      printf("  ## unable to allocate memory / camera\n");
      exit(1);
    }
  }

  /* adjust coeffs */
  mesh = cv.mesh[sc->idmesh];
  c->eye[0] = c->eye[1] = 0.0;
  c->eye[2] = sc->dmax;

  c->vecup  = up;

  c->speed[0] = 0.0;
  c->speed[1] = 0.0;
  c->speed[2] = c->eye[2];
  dd = -1.0 / sqrt(c->speed[2]*c->speed[2]);
  c->speed[2] *= dd;
  c->spmod  = 0.01*sc->dmax;
  c->altinc = 0.01*(mesh->ymax-mesh->ymin);

  /* set sun position */
  updateSun(sc,c);

  if ( ddebug ) {
    look[0] = c->eye[0] + sc->dmax*c->speed[0];
    look[1] = c->eye[1] + sc->dmax*c->speed[1];
    look[2] = c->eye[2] + sc->dmax*c->speed[2];
    printf("    eye   %g %g %g\n",c->eye[0],c->eye[1],c->eye[2]);
    printf("    speed %g %g %g\n",c->speed[0],c->speed[1],c->speed[2]);
    printf("    look  %g %g %g\n",look[0],look[1],look[2]);
  }

  return(c);
}

#ifdef __cplusplus
}
#endif
