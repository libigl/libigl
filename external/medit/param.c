#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "extern.h"
#include "sproto.h"

#ifndef  ON
#define  ON     1
#define  OFF    0
#endif

extern void ortho2D(pScene ,ubyte );

/* globals */
typedef struct sparval {
  int   arg;
} Parval;
typedef Parval * pParval;


static void parMotion(int x,int y) {
  pScene    sc;

  sc = cv.scene[currentScene()];
  glEnable(GL_COLOR_LOGIC_OP);
  glLogicOp(GL_XOR);
  
  glColor3ub(255,255,0);
  setFont("helvetica",10);
  drwstr(10,sc->par.ys-120,"Vector length");
  glColor3ub(0,255,128);
  drwstr(150,sc->par.ys-120,"%g",10.1);
  glFlush();
  glDisable(GL_COLOR_LOGIC_OP);
}

static void parMouse(int button,int state,int x,int y) {
  pScene   sc;
  if ( button != GLUT_LEFT_BUTTON )  return;
  sc = cv.scene[currentScene()];

  if ( state == GLUT_DOWN ) {
    glColor3ub(0,255,128);
    glDrawBuffer(GL_FRONT);
    ortho2D(sc,ON);
    glutMotionFunc(parMotion);
  }
  else {
    glDrawBuffer(GL_BACK);
    ortho2D(sc,OFF);
    glutMotionFunc(parMotion);
  }
}



void parEdit(pScene sc) {
  glutMouseFunc(parMouse);
}


#ifdef __cplusplus
}
#endif
