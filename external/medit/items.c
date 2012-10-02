#include "medit.h"
#include "extern.h"
#include "sproto.h"


void drawRulers(pScene sc) {
  pMesh  mesh = cv.mesh[sc->idmesh];
  
  if ( ddebug ) printf("draw rulers\n");
  
  glPushMatrix();
  glTranslatef(-mesh->xtra,-mesh->ytra,-mesh->ztra);
  glColor3fv(sc->par.line);
  glEnable(GL_LINE_STIPPLE);
  glLineStipple(1,0x0444);
  glLineWidth(2.);
  glBegin(GL_LINES);
    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmin);
    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmin);

    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmax);
    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmax);

    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmin,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmin);
    glVertex3f(mesh->xmin,mesh->ymax,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymin,mesh->zmax);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmin);
    glVertex3f(mesh->xmax,mesh->ymax,mesh->zmax);
  glEnd();
  glLineWidth(1.);
  glDisable(GL_LINE_STIPPLE);
  glPopMatrix();

}

void drawAxis(pScene sc,int dim) {
  pMesh  mesh;

  /* default */
  if ( ddebug ) printf("draw axis\n");
  mesh = cv.mesh[sc->idmesh];

  glPushMatrix();
  glTranslatef(1.01*(mesh->xmin-mesh->xtra),
               1.01*(mesh->ymin-mesh->ytra),
               1.01*(mesh->zmin-mesh->ztra));
  glScalef(0.6*sc->dmin,0.6*sc->dmin,0.6*sc->dmin);
  glLineWidth(max(2,sc->par.linewidth));
  glColor3f(1.0,0.,0.);
  if ( mesh->dim == 2 ) {
    glBegin(GL_LINE_STRIP);
    glVertex2f(0.0, 0.0);
    glVertex2f(1.0, 0.0);
    glVertex2f(0.95, 0.01);
    glVertex2f(0.95, -0.01);
    glVertex2f(1.0, 0.0);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex2f(0.0, 0.0);
    glVertex2f(0.0, 1.0);
    glVertex2f(-0.01, 0.95);
    glVertex2f(0.01, 0.95);
    glVertex2f(0.0, 1.0);
    glEnd();
  
    /*glColor3f(0.0f,1.0f,0.0f);*/
    glColor3f(1.-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
    glRasterPos3f(1.02, 0.0, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'x');
    glRasterPos3f(0.0, 1.02, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'y');
  }
  else {
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glVertex3f(0.95, 0.01, 0.0);
    glVertex3f(0.95, -0.01, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    glVertex3f(0.95, 0.0, 0.01);
    glVertex3f(0.95, 0.0, -0.01);
    glVertex3f(1.0, 0.0, 0.0);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.95, 0.01);
    glVertex3f(0.0, 0.95, -0.01);
    glVertex3f(0.0, 1.0, 0.0);
    glVertex3f(0.01, 0.95, 0.0);
    glVertex3f(-0.01, 0.95, 0.0);
    glVertex3f(0.0, 1.0, 0.0);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glVertex3f(0.01, 0.0, 0.95);
    glVertex3f(-0.01, 0.0, 0.95);
    glVertex3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.01, 0.95);
    glVertex3f(0.0, -0.01, 0.95);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();
    
    /*glColor3f(0.0f,1.0f,0.0f);*/
    glColor3f(1.-sc->par.back[0],1.0-sc->par.back[1],1.0-sc->par.back[2]);
    glRasterPos3f(1.02, 0.0, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'x');
    glRasterPos3f(0.0, 1.02, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'y');
    glRasterPos3f(0.0, 0.0, 1.02);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'z');
  }
  
  glLineWidth(1.0);
  glPopMatrix();
}


void drawBox(pScene sc,pMesh mesh,int mode) {
  pMaterial  pm;
  float      cx,cy,cz;
  int        i,k,m;

  glDisable(GL_LIGHTING);
  glPushMatrix();
  glScalef(1.01 * fabs(mesh->xmax-mesh->xmin),
           1.01 * fabs(mesh->ymax-mesh->ymin),
           1.01 * fabs(mesh->zmax-mesh->zmin));
  glColor3f(1.0,0.0,0.5);
  glutWireCube(1.0);
  glPopMatrix();

  /* one box per sub-domain */
  if ( mode ) {
    for (m=0; m<sc->par.nbmat; m++) {
      pm = &sc->material[m];
      for (i=0; i<MAX_LIST; i++) {
        k  = pm->depmat[i];
        if ( !k || pm->flag )  continue;
        cx = 0.5 * (pm->ext[3]+pm->ext[0]);
        cy = 0.5 * (pm->ext[4]+pm->ext[1]);
        cz = 0.5 * (pm->ext[5]+pm->ext[2]);
        glPushMatrix();
        glColor3fv(pm->dif);
        glTranslatef(cx,cy,cz);
        glScalef(pm->ext[3]-pm->ext[0],pm->ext[4]-pm->ext[1],
                 pm->ext[5]-pm->ext[2]);
        glutWireCube(1.0);
        glPopMatrix();
      }
    }
  }
}


void drawCube(pScene sc,pMesh mesh) {
  pTransform  cubetr;
  pCube       cube;
  float       x1,y1,z1,x2,y2,z2,xd,yd,zd;

  cube   = sc->cube;
  cubetr = cube->cubetr;

  if ( cube->active & C_UPDATE )  updateCube(cube,mesh);

  glDisable(GL_LIGHTING);
  glPushMatrix();
  glMultMatrixf(cubetr->matrix);

  glLineWidth(3.0);
  if ( cube->active & C_EDIT )        glColor3f(1.0,0.0,1.0);
  else if ( cube->active & C_FREEZE ) glColor3f(0.0,0.6,0.9);
  else                                glColor3f(0.0,1.0,0.0);
  x1 = cube->cmi[0] - mesh->xtra;
  y1 = cube->cmi[1] - mesh->ytra;
  z1 = cube->cmi[2] - mesh->ztra;
  x2 = cube->cma[0] - mesh->xtra;
  y2 = cube->cma[1] - mesh->ytra;
  z2 = cube->cma[2] - mesh->ztra;
  xd = cube->cma[0] - cube->cmi[0];
  yd = cube->cma[1] - cube->cmi[1];
  zd = cube->cma[2] - cube->cmi[2];
  
  glBegin(GL_QUADS);
    glVertex3f(x1,y1,z1);
    glVertex3f(x1+xd,y1,z1);
    glVertex3f(x1+xd,y1+yd,z1);
    glVertex3f(x1,y1+yd,z1);

    glVertex3f(x1,y1,z2);
    glVertex3f(x1+xd,y1,z2);
    glVertex3f(x1+xd,y1+yd,z2);
    glVertex3f(x1,y1+yd,z2);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(x1,y1,z1);
    glVertex3f(x1,y1,z2);

    glVertex3f(x1+xd,y1,z1);
    glVertex3f(x1+xd,y1,z2);

    glVertex3f(x1+xd,y1+yd,z1);
    glVertex3f(x1+xd,y1+yd,z2);

    glVertex3f(x1,y1+yd,z1);
    glVertex3f(x1,y1+yd,z2);
  glEnd();

  glLineWidth(1.0);

  glPopMatrix();
}


void drawGrid(pScene sc,pMesh mesh) {
  int k;

  /* default */
  if ( ddebug ) printf("draw grid + graduation\n");

  if ( !sc->grid ) {
    sc->grid = glGenLists(1);
    glNewList(sc->grid,GL_COMPILE);
    glBegin(GL_LINES);
    for (k=0; k<5; k++) {
      glVertex3f(k*0.25,0.,0.);  glVertex3f(k*0.25,1.,0.);
      glVertex3f(0.,k*0.25,0.);  glVertex3f(1.,k*0.25,0.);
      glVertex3f(0.,k*0.25,0.);  glVertex3f(0.,k*0.25,1.);
      glVertex3f(0.,0.,k*0.25);  glVertex3f(0.,1.,k*0.25);
      glVertex3f(k*0.25,0.,0.);  glVertex3f(k*0.25,0.,1.);
      glVertex3f(0.,0.,k*0.25);  glVertex3f(1.,0.,k*0.25);
     }
    glEnd();
    glEndList();
  }

  /* call display list */
  glPushMatrix();
  glTranslatef(-0.3*fabs(sc->dmax),0.,-4.7*sc->dmax);
  glRotatef(-60.,1.,0.,0.);
  glRotatef(-120.,0.,0.,1.);
  glScalef(2.5*sc->dmax,2.5*sc->dmax,2.5*sc->dmax);
  glDisable(GL_LIGHTING);

  glColor3f(0.4,0.4,0.4);
  glCallList(sc->grid);

  glColor3fv(sc->par.line);
  output3(0.0,0.0,0.0,"%.2f",mesh->xmin);
  output3(1.1,0.0,0.0,"%.2f",mesh->xmax);
  output3(0.0,1.01,0.0,"%.2f",mesh->ymax);
  output3(0.0,0.0,1.01,"%.2f",mesh->zmax);
  glEnable(GL_LIGHTING);
  glPopMatrix();
}


void drawBase(pScene sc,pMesh mesh) {
  int  k;

  /* default */
  if ( ddebug ) printf("draw base\n");

  if ( !sc->grid ) {
    sc->grid = glGenLists(1);
    glNewList(sc->grid,GL_COMPILE);
    if ( glGetError() )  return;
    glColor3f(0.5,0.5,0.5);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    for (k=0; k<21; k+=5) {
      glVertex3f(k*0.05,0.,0.);  glVertex3f(k*0.05,1.,0.);
      glVertex3f(0.,k*0.05,0.);  glVertex3f(1.,k*0.05,0.);
    }
    glEnd();
    glColor3f(0.6,0.6,0.6);
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (k=0; k<21; k++) {
      if ( k%5 == 0 ) continue;
      glVertex3f(k*0.05,0.,0.);  glVertex3f(k*0.05,1.,0.);
      glVertex3f(0.,k*0.05,0.);  glVertex3f(1.,k*0.05,0.);
    }
    glEnd();
    glEndList();
  }

  glPushMatrix();
  glTranslatef(-1.5*sc->dmax,-1.5*sc->dmax,-0.5*(mesh->zmax-mesh->zmin));
  glScalef(3*sc->dmax,3*sc->dmax,3*sc->dmax);
  glDisable(GL_LIGHTING);
    glCallList(sc->grid);
  glPopMatrix();
}


/* draw HUD system for flight */
void drawHUD(pScene sc) {
  pCamera  c;
  pMesh    mesh;
  GLfloat  xm,ym,x,y,dx,dy,alt;
  double   azim,elev;
  int      i,j;

  if ( ddebug )  fprintf(stdout,"drawHUD\n");
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(-0.5,639.5,-0.5,479.5);

  c  = sc->camera;
  mesh = cv.mesh[sc->idmesh];
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  glColor3f(0.0,0.5,0.0);
  glLineWidth(1.0f);
  xm = sc->par.xs / 2.0f;
  ym = sc->par.ys / 2.0f;
  glRecti(200,160,440,320);

  /* altitude */
  glColor3f(1.0,0.5,0.0);
  output2(230,330,"speed: %6.1f",1000.0f*c->spmod/sc->dmax);

  switch(c->vecup) {
    case X_AXIS: alt = (c->eye[1]+mesh->xtra)/(mesh->xmax-mesh->xmin); break;
    case Y_AXIS: alt = (c->eye[1]+mesh->ytra)/(mesh->ymax-mesh->ymin); break;
    case Z_AXIS: alt = (c->eye[1]+mesh->ztra)/(mesh->zmax-mesh->zmin); break;
    default: alt = 0.0;  break;
  }
  glColor3f(1.0,0.5,0.0);
  output2(350,330,"alt: %9.1f",1000.0f*alt);

  /* horiz rulers */
  output2(310,139.0f,"azim");
  glColor3f(0.0,1.0,0.0);
  output2(197,150,"-180");
  output2(257,150," -90");
  output2(317,150,"0");
  output2(377,150," 90");
  output2(437,150,"180");
  x = 200.0f;
  glBegin(GL_LINES);
  for (i=1; i<8; i++) {
    x += 240.0 / 8.0;
    glVertex2f(x,158.0);
    glVertex2f(x,162.0);
  }
  glEnd();
  
  /* vert rulers */
  glColor3f(0.0,1.0,0.0);
  output2(185,160,"-90");
  output2(185,200,"-45");
  output2(185,240,"0");
  output2(185,280,"45");
  output2(185,320,"90");
  y = 160.0f;
  glBegin(GL_LINES);
  for (i=1; i<8; i++) {
    y += 160.0 / 8.0;
    glVertex2f(198,y);
    glVertex2f(202,y);
  }
  glEnd();

  /* azimuth */
  azim = Azimuth(c);
  if ( azim > 0.0f )      azim =  180.0f - azim;
  else if ( azim < 0.0f ) azim = -180.0f - azim;
  x = 2.0/3.0*azim + 320.0f;
  glColor3f(1.0,0.0,0.0);
  glLineWidth(1.0);
  output2(x,143.0,"%d",azim>0 ? (int)(azim+0.5) : (int)(azim-0.5));
  glBegin(GL_LINES);
    glVertex2f(x,166.0);
    glVertex2f(x,318.0);
  glEnd();
  y  = 160.0f;
  dy = 160.0 / 8.0;
  glBegin(GL_LINES);
  for (i=0; i<8; i++) {
    glVertex2f(x-4,y);
    glVertex2f(x+4,y);
    for (j=0; j<5; j++) {
      glVertex2f(x-2,y+j*dy/5.0);
      glVertex2f(x+2,y+j*dy/5.0);
    }
    y += dy;
  }
  glEnd();
  
  /* elevation */
  elev = Elevation(c);
  if ( elev > 90.0f )       y = 320.0f;
  else if ( elev < -90.0f ) y = 160.0f;
  else y = 8.0/9.0 * elev + 240.0;
  glColor3f(1.0,0.0,0.0);
  output2(175.0,y,"%5.1f",elev);
  glBegin(GL_LINES);
    glVertex2f(206.0,y);
    glVertex2f(438.0,y);
  glEnd();
  x  = 200.0f;
  dx = 240.0f / 8.0f;
    glBegin(GL_LINES);
  for (i=1; i<=8; i++) {
    glVertex2f(x,y-4);
    glVertex2f(x,y+4);
    for (j=0; j<5; j++) {
      glVertex2f(x+j*dx/5.0,y-2);
      glVertex2f(x+j*dx/5.0,y+2);
    }
    x += dx;
  }
  glEnd();

  /* horizon */
  glLineWidth(2.0f);
  glColor3f(0.0,0.0,1.0);
  glBegin(GL_LINES);
    glVertex2f(200.0,240.0);
    glVertex2f(440.0,240.0);
  glEnd();
  glLineWidth(1.0f);

  /* HUD */
  glColor3f(1.0,0.0,0.0);
  glBegin(GL_LINES);
    glVertex2f(310,230);
    glVertex2f(330,230);
    glVertex2f(310,250);
    glVertex2f(330,250);
    
    glVertex2f(310,230);
    glVertex2f(310,234);
    glVertex2f(330,230);
    glVertex2f(330,234);
    
    glVertex2f(310,246);
    glVertex2f(310,250);
    glVertex2f(330,246);
    glVertex2f(330,250);
  glEnd();
  /*glRecti(318,238,322,242);*/
  glColor3f(0.0,1.0,0.0);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
    glVertex2f(320.,235.);
    glVertex2f(320.,245.);
    glVertex2f(315.,240.);
    glVertex2f(325.,240.);
  glEnd();
  glLineWidth(1.0f);
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glEnable(GL_DEPTH_TEST);
}
