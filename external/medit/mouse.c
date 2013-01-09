#include "medit.h"
#include "extern.h"
#include "sproto.h"
	
#ifndef  ON
#define  ON     1
#define  OFF    0
#endif


GLuint    lasttime;
GLboolean tracking = GL_FALSE,ctracking = GL_FALSE;
GLboolean picking = GL_FALSE;
int       cbutton  = 0;
int       startx,starty,curx,cury;
int       rxi,ryi,rx1,rx2,ry1,ry2; 

#define MinWH   0.1
#define MaxWH   0.9

#ifndef GLUT_BUTTON_3
#define GLUT_BUTTON_3   2
#define GLUT_BUTTON_4   3
#endif

// Patch for machines without a middle mouse button. Interprets
// CTRL+ALT+LEFT_CLICK and ALT+RIGHT_CLICK as a MIDDLE_CLICK
//
// Input:
//   button  current button value
// Output:
//   button  new button value
//   
int alt_ctrl_left_to_middle(const int button)
{
#ifdef __APPLE__
  int keyact = glutGetModifiers();
  if (keyact & GLUT_ACTIVE_ALT && 
    ((keyact & GLUT_ACTIVE_CTRL && button == GLUT_LEFT_BUTTON) || 
    button == GLUT_RIGHT_BUTTON))
  {
    return GLUT_MIDDLE_BUTTON;
  }
#endif
  return button;
}

/* project x,y, onto a hemi-sphere */
static void point2Vect(int x,int y,int w,int h,float *v) {
  double   d,a,areax,areay;

  areax = (w-startx) / w;
  areay = (h-starty) / h;
  if ( areax > MinWH && areax < MaxWH && areay > MinWH && areay < MaxWH ) {
    v[0] = (2.0 * x - w) / w;
    v[1] = (h - 2.0 * y) / h;
    v[2] = 1.0f;
  }
  else {
    v[0] =  2.0f*(x-startx) / w;
    v[1] = -2.0f*(y-starty) / h;
    v[2] = 1.0f;
  }
  d = v[0]*v[0]+v[1]*v[1];
  if ( d == 0.0f )  return;
  d = sqrt(d);

  v[2] = cos(M_PI_2 * ((d < 1.0) ? d : 1.0));
  d    = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  a    = 1.0f / sqrt(d);
  v[0] *= a;  
  v[1] *= a;  
  v[2] *= a;
}


void ortho2D(pScene sc,ubyte mode) {
  if ( mode == ON ) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,sc->par.xs,0.,sc->par.ys);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
  }
  else if ( mode == OFF ) {
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.);
  }
}

static void drawRubberBand(int xa,int ya,int xb,int yb) {
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  glLineWidth(2.0);
  glColor3f(0.5,0.5,0.5);
  glRecti(xa,ya,xb,yb);
}

static void rubberMotion(int x,int y) {
  pScene      sc;
  pPersp      p;

  sc = cv.scene[currentScene()];
  p  = sc->persp;

  glEnable(GL_COLOR_LOGIC_OP);
  glLogicOp(GL_XOR);

  /* draw frame */
  drawRubberBand(rxi,ryi,rx1,ry1);
  rx2 = x;
  ry2 = sc->par.ys-y;
  drawRubberBand(rxi,ryi,rx2,ry2);
  glFlush();

  glLogicOp(GL_COPY);
  glDisable(GL_COLOR_LOGIC_OP);

  /* keep old coords */
  rx1 = rx2;
  ry1 = ry2;
}


void zoomMotion(int x,int y) {
  pScene      sc;
  pPersp      p;
  int         dy;

  sc = cv.scene[currentScene()];
  p  = sc->persp;

  dy = starty - y;
  if ( dy > 0 )
    if ( p->fovy < 1.0e-02 )  
      return;
    else 
      p->fovy = MEDIT_MAX(0.95*p->fovy,1e-05);
  else if ( p->fovy > 160.0 )  return;
  else
    p->fovy = MEDIT_MIN(1.1*p->fovy,179.0);

  farclip(1);
  starty = y;
  glutPostRedisplay();
}


void mouse(int button,int state,int x,int y) {
  button = alt_ctrl_left_to_middle(button);
  pScene      sc;
  pTransform  tr;
  pPersp      p;
  int         keyact,idw = currentScene();
  static int  olds = -1;
#ifdef IGL
  // Tweakbar has precedence over everything else
  if(TwEventMouseButtonGLUT(button,state,x,y) && state == GLUT_DOWN)
  {
    return;
  }
#endif

  /* default */
  if ( ddebug ) printf("control mouse %d\n",state);
  sc = cv.scene[idw];
  p  = sc->persp;
  if ( sc->cube->active & C_EDIT )
    tr = sc->cube->cubetr;
  else if ( sc->clip->active & C_EDIT ) 
    tr = sc->clip->cliptr;
  else
    tr = sc->view;
  tr->mstate  = state;
  tr->mbutton = button;

  /* check if ctrl-shift-alt pressed */
  keyact = glutGetModifiers();

  if ( state == GLUT_DOWN ) {
  
  picking=GL_FALSE;
    tracking = GL_TRUE;
    lasttime = glutGet(GLUT_ELAPSED_TIME);

    if ( button == GLUT_LEFT_BUTTON ) {
      if ( keyact & GLUT_ACTIVE_SHIFT ) {
        /* entity designation */
        picking = GL_TRUE;
	if ( sc->picklist ) glDeleteLists(sc->picklist,1);
	sc->picklist = pickingScene(sc,x,y,0);
	return;
      }

      else if ( keyact & GLUT_ACTIVE_ALT ) {
	    /* zoom */
	    starty = y;
	    glutMotionFunc(zoomMotion);
	    return;
      }

      else if ( keyact & GLUT_ACTIVE_CTRL ) {
        /* rubberband selection */
        glutSetCursor(GLUT_CURSOR_CROSSHAIR);
        p->rubix  = p->rubfx = x;
        p->rubiy  = p->rubfy = sc->par.ys-y;
        rxi = rx1 = x;
        ryi = ry1 = sc->par.ys-y;
        p->rubber = 1;
        glDrawBuffer(GL_BACK_LEFT);
        ortho2D(sc,ON);
        glutMotionFunc(rubberMotion);
        return;
      }
    }
    
    else if ( button == GLUT_MIDDLE_BUTTON && keyact & GLUT_ACTIVE_SHIFT ) {
      picking = GL_TRUE;
      if ( sc->picklist ) glDeleteLists(sc->picklist,1);
      sc->picklist = pickingScene(sc,x,y,LPoint);
      return;
    }

    /* transformation */
    startx = x;
    starty = y;
    point2Vect(x,y,sc->par.xs,sc->par.ys,tr->pos);
    glutSetCursor(GLUT_CURSOR_INFO);
  }

  else if ( state == GLUT_UP ) {

    if ( button == GLUT_LEFT_BUTTON ) {
      
      if ( keyact & GLUT_ACTIVE_CTRL ) {
        /* rubberband selection */
        p->rubfx  = x;
        p->rubfy  = sc->par.ys-y;
        p->rubber = 2;
        glDrawBuffer(GL_BACK_LEFT);
        ortho2D(sc,OFF);
        glutMotionFunc(motion);
        return;
      }
      
      else if ( keyact & GLUT_ACTIVE_ALT ) {
        glutMotionFunc(motion);
	    return;
      }

      else if ( picking == GL_TRUE ) {
        picking = GL_FALSE;
        reshapeScene(sc->par.xs,sc->par.ys);
        glutPostRedisplay();
      }
    }
    
    glutMotionFunc(motion);
    if ( sc->clip->active & C_EDIT )
      sc->clip->active |= C_REDO;

    /* transformation */
    glutSetCursor(GLUT_CURSOR_INHERIT);
    tracking = GL_FALSE;
    if ( glutGet(GLUT_ELAPSED_TIME) >= lasttime ) { 
      if ( tr->manim == GL_TRUE )  glutIdleFunc(glutIdle);
      else  tr->angle = 0.0;
      /*if ( abs(startx-x) + abs(starty-y) > 0 )*/
        glutPostRedisplay();
    }
    else if ( tr->manim == GL_TRUE && olds == idw )  
      glutIdleFunc(NULL);
  }

  olds = idw;
}

#ifdef IGL
void passive_motion(int x,int y)
{
  // Tweakbar has precedence over everything else
  if(TwEventMouseMotionGLUT(x,y))
  {
    glutPostRedisplay();
  }
}
#endif

void motion(int x,int y) {
  pScene      sc;
  pTransform  tr;
  pPersp      p;
  GLuint      gtime;
  double      deltax,deltay;
  float       coeff,pos[3],dx,dy,dz;
  int         idw = currentScene();

  /* default */
  if ( picking )  return;
  if ( ddebug ) fprintf(stdout,"motion\n");

#ifdef IGL
  // Tweakbar has precedence over everything else
  if(TwEventMouseMotionGLUT(x,y))
  {
    if(!tracking)
    {
      glutPostRedisplay();
    }
  }
#endif

  if ( tracking == GL_FALSE )  return;
  sc = cv.scene[idw];
  p  = sc->persp;
  if ( p->rubber == 1 )  return;

  /* what is transformed ? */
  if ( sc->cube->active & C_EDIT ) 
    tr = sc->cube->cubetr;
  else if ( sc->clip->active & C_EDIT ) 
    tr = sc->clip->cliptr;
  else
    tr = sc->view;

  if ( tr->mstate != GLUT_DOWN )  return;
  if ( picking )  tr->angle = 0.0f;

  gtime = glutGet(GLUT_ELAPSED_TIME);
  if ( (animate || sc->type & S_FOLLOW) && gtime < lasttime+40 )  return;

  if ( tr->mbutton == GLUT_LEFT_BUTTON ) {
    /* calculate axis of rotation: cross product */
    point2Vect(x,y,sc->par.xs,sc->par.ys,pos);
    tr->axis[0] = tr->pos[1]*pos[2] - tr->pos[2]*pos[1];
    tr->axis[1] = tr->pos[2]*pos[0] - tr->pos[0]*pos[2];
    tr->axis[2] = tr->pos[0]*pos[1] - tr->pos[1]*pos[0];
 
    /* calculate angle to rotate by */
    if ( animate && saveimg )
      tr->angle = 2.0f;
    else {
      dx = pos[0] - tr->pos[0];
      dy = pos[1] - tr->pos[1];
      dz = pos[2] - tr->pos[2];
      tr->angle = 180.0*sqrt(dx*dx+dy*dy+dz*dz);
    }

    /* reset for next time */
    tr->pos[0] = pos[0];
    tr->pos[1] = pos[1];
    tr->pos[2] = pos[2];
    lasttime   = gtime;

    if ( sc->cube->active & C_ON && sc->cube->active & C_EDIT )
      sc->cube->active |= C_UPDATE;
    else if ( sc->clip->active & C_ON && 
	         (sc->clip->active & C_EDIT || sc->clip->active & C_FREEZE) )
      sc->clip->active |= C_UPDATE;

    glutPostRedisplay();
  }

  else if ( tr->mbutton == GLUT_MIDDLE_BUTTON ) {
    coeff  = tr->manim == GL_TRUE ? 0.2 : 2.0;
    deltax = coeff * (x-startx) / (float)sc->par.xs;
    deltay = coeff * (starty-y) / (float)sc->par.ys;

    if ( deltax != 0.0 )
      tr->panx += -deltax * p->depth * tan(p->fovy/360.*M_PI);
    if ( deltay != 0.0 )
      tr->pany += -deltay * p->depth * tan(p->fovy/360.*M_PI);
    tr->angle = 0.0;
    startx = x;
    starty = y;

    lasttime = gtime;
    if ( sc->cube->active & C_ON && sc->cube->active & C_EDIT )
      sc->cube->active |= C_UPDATE;
    else if ( sc->clip->active & C_ON && 
        (sc->clip->active & C_EDIT || sc->clip->active & C_FREEZE) )
      sc->clip->active |= C_UPDATE;

    glutPostRedisplay();
  }
}

void mouseCamera(int button,int state,int x,int y) {
  /* default */
  if ( ddebug ) printf("control mouse camera %d button %d\n",state,button);
  button = alt_ctrl_left_to_middle(button);

  cbutton = button;
  if ( state == GLUT_DOWN ) {
    ctracking = GL_TRUE;
    startx = x;
    starty = y;
    curx   = x;
    cury   = y;
  }
  else {
    startx = x;
    starty = y;
    ctracking = GL_FALSE;
  }
}

void motionCamera(int x,int y) {
  pScene   sc;
  pCamera  c;
  double   dazim,delev,azim,elev;
  float    cfelev,cfazim;

  /* keep current pos */
  curx = x;
  cury = y;

  if ( animate ) return;
  sc = cv.scene[currentScene()];
  c  = sc->camera;
  azim  = Azimuth(c);
  elev  = Elevation(c);
   switch (cbutton) {
    case GLUT_LEFT_BUTTON:
      cfelev = 50.0;
      cfazim = 50.0;
      delev  = cfelev * (y-starty)/(float)sc->par.ys;
      dazim  = cfazim * (x-startx)/(float)sc->par.xs;
      startx = x;
      starty = y;
      elev  += delev;
      azim  -= dazim;
      break;
    case GLUT_MIDDLE_BUTTON:
      break;
    case GLUT_BUTTON_3:
      puts("button3");
      break;
    case GLUT_BUTTON_4:
      puts("button4");
      break;
  }
  updateCamera(sc,c,azim,elev);
  reshapeScene(sc->par.xs,sc->par.ys);
  glutPostRedisplay();
}


void animateCamera() {
  pScene   sc;
  pCamera  c;
  double   dazim,delev,azim,elev;
  float    cfelev,cfazim;

  if ( !animate || !ctracking )  return;
  sc = cv.scene[currentScene()];
  c  = sc->camera;
  azim  = Azimuth(c);
  elev  = Elevation(c);

  switch (cbutton) {
    case GLUT_LEFT_BUTTON:
      cfelev = 3.0;
      cfazim = 3.0;
      delev  = 2.0*(cury-starty)/(float)sc->par.ys;
      dazim  = 2.0*(curx-startx)/(float)sc->par.xs;
      if ( delev >= 0.0 ) delev *= delev;
      else                delev = -delev*delev;
      if ( dazim >= 0.0 ) dazim *= dazim;
      else                dazim  = -dazim*dazim;
      elev += cfelev * delev;
      azim -= cfazim * dazim;
      break;
    
    case GLUT_MIDDLE_BUTTON:
      break;
    case GLUT_BUTTON_3:
      puts("button3");
      break;
    case GLUT_BUTTON_4:
      puts("button4");
      break;
  }
  updateCamera(sc,c,azim,elev);
  reshapeScene(sc->par.xs,sc->par.ys);
  glutPostRedisplay();
}
