#include "medit.h"
#include "extern.h"
#include "sproto.h"

extern int refmat,reftype,refitem;
int       *pilmat,ipilmat;
int        refpick = -1;

GLfloat  dazim = 1.0;
GLfloat  delev = 1.0;

extern void mouseCamera(int button,int state,int x,int y);
extern void motionCamera(int x,int y);


static void usage() {
  fprintf(stdout,"\n");
  fprintf(stdout,"-- Medit: OnLine Help --\n");

  fprintf(stdout,"** Rendering options (toggle):\n");
  fprintf(stdout,"f - facets        l - lines          g - const. entities\n");
  fprintf(stdout,"c - obj. color    e - material       b - back color\n");
  fprintf(stdout,"A - axis          B - box            G - grid\n");
  fprintf(stdout,"C - capping       r(R) - hide (show) refs.\n");
  fprintf(stdout,"n - smooth/flat shading\n");

  fprintf(stdout,"\n** Controls options:\n");
  fprintf(stdout,"i - reset view    a - animate        I - interactive\n");
  fprintf(stdout,"h - online help   q - quit           X - close window\n");

  fprintf(stdout,"\n** View controls (ALT + key):\n");
  fprintf(stdout,"c - copy   d - duplicate p - paste  l - link   u - unlink\n");
  fprintf(stdout,"J - toggle flight        y - change axis\n");
  
  fprintf(stdout,"\n** Misc. features\n");
  fprintf(stdout,"L - load prefs    W - save prefs\n");
  fprintf(stdout,"H - hardcopy PPM\n");
  fprintf(stdout,"F - face nums     P - point nums     # - select entity\n");
  fprintf(stdout,"N - normals       O - oppos. normals ! - select plane\n");
  fprintf(stdout,"m - data          o - iso-lines      w - tensor/vector\n");
  fprintf(stdout,"v - streamlines   k - toggle elev.   K - elev. coeff.\n");
  fprintf(stdout,"j - toggle deco   p - palette        S - play animation\n");
  
  fprintf(stdout,"+/- scale object  z/Z scale view\n");
  fprintf(stdout,"F1,F2,F3 - clipping: On/Off, Edit, Freeze\n");
}

/* special keys CAMERA mode */
void specCamera(pScene sc,int key) {
  pCamera     c;
  double      dd,azim,elev;
  GLfloat     axe[3];
  int         keyact;
  
  c  = sc->camera;
  keyact = glutGetModifiers();

  axe[0] = c->speed[2];
  axe[1] = 0.0;
  axe[2] = -c->speed[0];
  dd = sqrt(axe[0]*axe[0] + axe[2]*axe[2]);
  if ( dd != 0.0f ) {
    axe[0] /= dd;
    axe[2] /= dd;
  }

  switch (key) {
  case GLUT_KEY_LEFT:
    if ( keyact & GLUT_ACTIVE_SHIFT ) {
      c->eye[0] += axe[0]*c->altinc;
      c->eye[2] += axe[2]*c->altinc;
      reshapeScene(sc->par.xs,sc->par.ys);
      glutPostRedisplay();
      return;
    }
    else {
      azim = Azimuth(c);
      elev = Elevation(c);
      azim += dazim;
    }
    break;
  case GLUT_KEY_RIGHT:
    if ( keyact & GLUT_ACTIVE_SHIFT ) {
      c->eye[0] -= axe[0]*c->altinc;
      c->eye[2] -= axe[2]*c->altinc;
      reshapeScene(sc->par.xs,sc->par.ys);
      glutPostRedisplay();
      return;
    }
    else {
      azim = Azimuth(c);
      elev = Elevation(c);
      azim -= dazim;
    }
    break;
  case GLUT_KEY_UP:
    if ( keyact & GLUT_ACTIVE_SHIFT ) {
      c->eye[1] += c->altinc;
      reshapeScene(sc->par.xs,sc->par.ys);
      glutPostRedisplay();
      return;
    }
    else {
      azim = Azimuth(c);
      elev = Elevation(c);
      elev -= delev;
    }
    break;
  case GLUT_KEY_DOWN:
    if ( keyact & GLUT_ACTIVE_SHIFT ) {
      c->eye[1] -= c->altinc;
      reshapeScene(sc->par.xs,sc->par.ys);
      glutPostRedisplay();
      return;
    }
    else {
      azim = Azimuth(c);
      elev = Elevation(c);
      elev += delev;
    }
    break;
  default:
    return;
  }

  updateCamera(sc,c,azim,elev);

  /* refresh scene */
  reshapeScene(sc->par.xs,sc->par.ys);
  glutPostRedisplay();
}


/* change center of scene */
static void changeCenter(pScene sc,pMesh mesh) {
  pPoint     p0;
  pTriangle  pt;
  pQuad      pq;
  pTetra     ptt;
  pHexa      ph;
  float      cx,cy,cz;
  int        i;

  if ( !refitem )  return;
  cx = cy = cz = 0.0;
  switch(reftype) {
  case LPoint:
    p0 = &mesh->point[refitem];
    cx = p0->c[0];
    cy = p0->c[1];
    cz = p0->c[2];
    break;
  case LTria:
    pt = &mesh->tria[refitem];
    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      cx += 0.33 * p0->c[0];
      cy += 0.33 * p0->c[1];
      cz += 0.33 * p0->c[2];
    }
    break;
  case LQuad:
    pq = &mesh->quad[refitem];
    for (i=0; i<4; i++) {
      p0 = &mesh->point[pq->v[i]];
      cx += 0.25 * p0->c[0];
      cy += 0.25 * p0->c[1];
      cz += 0.25 * p0->c[2];
    }
    break;
  case LTets:
    ptt = &mesh->tetra[refitem];
    for (i=0; i<4; i++) {
      p0 = &mesh->point[ptt->v[i]];
      cx += 0.25 * p0->c[0];
      cy += 0.25 * p0->c[1];
      cz += 0.25 * p0->c[2];
    }
    break;
  case LHexa:
    ph = &mesh->hexa[refitem];
    for (i=0; i<8; i++) {
      p0 = &mesh->point[ph->v[i]];
      cx += 0.125 * p0->c[0];
      cy += 0.125 * p0->c[1];
      cz += 0.125 * p0->c[2];
    }
    break;
  }

  /* reset translation and move to center */
  sc->view->panx = sc->view->pany = 0.0;
  sc->cx = -cx;
  sc->cy = -cy;
  sc->cz = -cz;
}


/* special keys PERSPECTIVE mode */
void special(int key,int x,int y) {
  pTransform  view;
  pScene      sc;
  pMesh       mesh;
  pClip       clip;
  pCube       cube;
  float       pancoeff = 0.1f;
  int         keyact,idw = currentScene();
  ubyte       post = TRUE;
 
  /* default */
  if ( ddebug ) printf("special key  %d\n",key);
  sc   = cv.scene[idw];

  /* special mode camera */
  if ( sc->persp->pmode == CAMERA ) {
    specCamera(sc,key);
    return;
  }
  
  view = sc->view;
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;
  cube = sc->cube;
  keyact = glutGetModifiers();
  if ( keyact & GLUT_ACTIVE_SHIFT )
    pancoeff = 0.01;

  switch(key) {
  case GLUT_KEY_F1:  /* toggle clipping plane */
    if ( mesh->dim == 3 ) {
      if ( cube->active & C_ON)
    keyCube('C',0,0);
      else
    keyClip('C',0,0);
    }
    post = FALSE;
    break;
  case GLUT_KEY_F2:  /* edit clipping plane */
    if ( mesh->dim == 3 ) {
      if ( cube->active & C_ON)
        keyCube('E',0,0);     
      else if ( clip->active & C_ON ) 
    keyClip('E',0,0);
    }
    post = FALSE;
    break;
  case GLUT_KEY_F3:  /* freeze clipping plane */
    if ( cube->active & C_ON )
      keyCube('F',0,0);
    else if ( clip->active & C_ON ) 
      keyClip('F',0,0);
    post = FALSE;
    break;
  case GLUT_KEY_F4:  /* toggle Vclip */
    if ( mesh->dim == 3 )  keyClip('Z',0,0);
    post = FALSE;
    break;
  case GLUT_KEY_F5:  /* Toggle Shrink */
    if ( mesh->ntet+mesh->nhex > 0 )
      keyFeature('V',0,0);
    else
      keyFeature('S',0,0);
    break;
  case GLUT_KEY_F6:  /* Increase Shrink */
    keyFeature('I',0,0);
    break;
  case GLUT_KEY_F7:  /* Decrease Shrink */
    keyFeature('i',0,0);
    break;

  case GLUT_KEY_LEFT:  /* translate eyes or object */
    if ( clip->active & C_EDIT ) {
      clip->cliptr->panx -= 0.02 * sc->dmax;
      clip->cliptr->angle = 0.0;
      clip->active |= C_UPDATE + C_REDO;
    }
    else if ( keyact & GLUT_ACTIVE_CTRL ) {
      sc->par.eyesep *= 0.9;
      printf("eyesep %f\n",sc->par.eyesep);
    }
    else
      view->panx -= pancoeff * sc->dmax;
    break;
  case GLUT_KEY_RIGHT:
    if ( clip->active & C_EDIT ) {
      clip->cliptr->panx += 0.02 * sc->dmax;
      clip->cliptr->angle = 0.0;
      clip->active |= C_UPDATE + C_REDO;
    }
    else if ( keyact & GLUT_ACTIVE_CTRL ) {
      sc->par.eyesep *= 1.1;
      printf("eyesep %f\n",sc->par.eyesep);
    }
    else
      view->panx += pancoeff * sc->dmax;
    break;
  case GLUT_KEY_UP:
    if ( clip->active & C_EDIT ) {
      clip->cliptr->pany += 0.02 * sc->dmax;
      clip->cliptr->angle = 0.0;
      clip->active |= C_UPDATE + C_REDO;
    }
    else
      view->pany += pancoeff * sc->dmax;
    break;
  case GLUT_KEY_DOWN:
    if ( clip->active & C_EDIT ) {
      clip->cliptr->pany -= 0.02 * sc->dmax;
      clip->cliptr->angle = 0.0;
      clip->active |= C_UPDATE + C_REDO;
    }
    else
      view->pany -= pancoeff * sc->dmax;
    break;
  default:
    return;
  }
  if ( post )  glutPostRedisplay();
}


void keyScene(unsigned char key,int x,int y) {
  pMaterial   pm;
  pTetra      ptt;
  pHexa       ph;
  pTriangle   pt;
  pQuad       pq;
  pScene      sc,sc1;
  pMesh       mesh;
  pClip       clip;
  pCube       cube;
  pPersp      p;
  pCamera     cam;
  double      dd;
  float       a,b,c,d;
  int         k,keyact,numit,idw = currentScene();
  ubyte       post = FALSE,dolist = FALSE;


  if ( idw < 0 ) exit(0);

  /* ESC = end medit */
  if ( key == 'q' || key == 27 ) 
#ifdef IGL
  {
    deleteScene(cv.scene[idw]);
    exit(0);
  }
#else
    exit(0);
#endif
  else if ( key == 'h' || key == '?' )
    usage();

  /* default */
  sc   = cv.scene[idw];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;
  cube = sc->cube;
  p    = sc->persp;

#ifdef IGL
  // Tweakbar has precedence over everything else
  if(TwEventKeyboardGLUT(key,x,y))
  {
    dolist = TRUE;
    post = TRUE;
    goto keySceneFinish;
  }

  // Hack so that '=/+' acts like zoom in/out
  key = (key == '+' || key == '=' ? 'z' : key);
  key = (key == '-' || key == '_' ? 'Z' : key);
#endif

  keyact = glutGetModifiers();
  if ( key == ' ' ) {
    if ( option == MORPHING )
      morphMesh(sc,mesh);
    else if ( sc->isotyp & S_PARTICLE ) {
      glutIdleFunc(0);
    } 
    else {      
      cam = sc->camera;
      cam->eye[0] += cam->spmod*cam->speed[0];
      cam->eye[1] += cam->spmod*cam->speed[1];
      cam->eye[2] += cam->spmod*cam->speed[2];
      reshapeScene(sc->par.xs,sc->par.ys);
    }
    post = TRUE;
  }

  else if ( islower(key) ) {
    switch(key) {
    case 'a':  /* toggle animate */
      keyAnim('A',0,0);
      break;
    case 'b':  /* backcolor */
      keyColor('b',0,0);
      break;
    case 'c':
      if ( keyact & GLUT_ACTIVE_ALT )
    keyView('C',0,0);
      else {
    sc->mode ^= S_COLOR;
    post = TRUE;
      }
      break;
    case 'd':
      if ( keyact & GLUT_ACTIVE_ALT ) keyView('D',0,0);
      break;
    case 'e': 
      keyColor('e',0,0);
      break;
    case 'f':
      sc->mode ^= S_FILL;
      post=TRUE;
      break;
    case 'g': 
      keyItem('g',0,0);
      break;
    case 'h':
      usage();
      break;
    case 'i':
      keyView('R',0,0);
      if ( clip->active & C_ON )  resetClip(sc,clip,mesh);
      if ( cube->active & C_ON )  resetCube(sc,cube,mesh);
      sc1 = sc;
      while ( sc1->slave > -1 ) {
        sc1 = cv.scene[sc1->slave];
        glutSetWindow(sc1->idwin);
        keyScene('i',0,0);
      }
      glutSetWindow(sc->idwin);
      break;
    case 'j':
      sc->type ^= S_DECO;
      post = TRUE;
      break;
    case 'k':
      keyMetric('k',0,0);
      break;
    case 'l':
      if ( keyact & GLUT_ACTIVE_ALT ) 
        keyView('L',0,0);
      else
        sc->mode ^= S_BDRY;  
      post = TRUE;  
      break;
    case 'm':  /* toggle metric */
      if ( mesh->nbb )  keyMetric('m',0,0);
      keyMode('n',0,0);
      break;
    case 'o': /* iso-lines */
      if ( mesh->nbb ) keyMetric('l',0,0);
      break;
    case 'p':
      if ( keyact & GLUT_ACTIVE_ALT ) 
        keyView('P',0,0);
      else
        keyMetric('p',0,0);
      break;
    case 'r':
      if ( refmat<0 || ipilmat == sc->par.nbmat ) return;
      pilmat[++ipilmat] = refmat;
      pm = &sc->material[refmat];
      pm->flag = 1;
      updatePoints(sc,mesh,refmat);
      if ( sc->picklist ) glDeleteLists(sc->picklist,1);
      sc->picklist = 0;
      refmat = -1;
      dolist = TRUE;
      post   = TRUE;
      break;
    case 's':
      if ( !refitem ) return;
      switch(reftype) {
      case LTria:
        pt = &mesh->tria[refitem];
        pt->v[0] = 0;
        break;
      case LQuad:
        pq = &mesh->quad[refitem];
        pq->v[0] = 0;
        break;
      case LTets:
        ptt = &mesh->tetra[refitem];
        ptt->v[0] = 0;
        break;
      case LHexa:
        ph = &mesh->hexa[refitem];
        ph->v[0] = 0;
        break;
      }
      sc->picklist = 0;
      dolist = TRUE;
      post   = TRUE;
      break;
    case 't':  /* toggle texture */
      /*keyColor('t',0,0);*/
      break;
    case 'u':
      if ( keyact & GLUT_ACTIVE_ALT ) 
    keyView('U',0,0);
      break;
    case 'v':
      keyMetric('v',0,0);
      break;
    case 'w':
      keyMetric('w',0,0);
      break;
    case 'x':
      if ( cube->active & C_EDIT ) {
        cube->cma[0] += 0.1*sc->dmax;
        post = TRUE;
        break;
      }
    case 'y':
      if ( cube->active & C_EDIT ) {
        cube->cma[1] += 0.1*sc->dmax;
        post = TRUE;
        break;
      }
      if ( p->pmode & CAMERA ) {
        cam = sc->camera;
        cam->vecup = (cam->vecup+1) % 3;
        switch(cam->vecup) {
          case X_AXIS: cam->altinc = (mesh->xmax-mesh->xmin); break;
          case Y_AXIS: cam->altinc = (mesh->ymax-mesh->ymin); break;
          case Z_AXIS: cam->altinc = (mesh->zmax-mesh->zmin); break;
        }
        cam->altinc *= 0.005f;
        sc->type |= S_RESET;
        reshapeScene(sc->par.xs,sc->par.ys);
      }
      /* tilt cut plane */
      else if ( clip->active & C_EDIT )
        tiltClip(sc,clip);
      else return;
      post = TRUE;
      break;

    case 'z':  /* zoom in */
      if ( cube->active & C_EDIT ) {
        cube->cma[2] += 0.1*sc->dmax;
        post = TRUE;
        break;
      }
      else {
    if ( p->rubber == 2 )
      setPersp(sc,p,1);
        else 
          p->fovy = MEDIT_MAX(0.9*p->fovy,1e-05);
        farclip(1);
      }
      post = TRUE;

      /* get linked view */
      sc1 = sc;
      while ( sc1->slave > -1 ) {
        sc1 = cv.scene[sc1->slave];
        memcpy(sc1->view,sc->view,sizeof(struct transform));
        memcpy(sc1->persp,sc->persp,sizeof(struct sperspective));
        glutSetWindow(sc1->idwin);
        reshapeScene(sc1->par.xs,sc1->par.ys);
      }
      glutSetWindow(sc->idwin);
      break;
    }
  }

  else if ( isupper(key) ) {
    switch(key) {
    case 'A':  /* toggle axis */
      keyItem('A',0,0);
      break;
    case 'B':  /* toggle box */
      keyItem('B',0,0);
      break;
    case 'C':  /* toggle capping */
      keyClip('K',0,0);
      break;
    case 'D':  /* selection cube */
      if (mesh->dim == 3 )  keyItem('D',0,0);
        break;
    case 'E':
      keyColor('E',0,0);
      glutSetWindow(sc->idwin);
      break;
    case 'F': /* toggle face nums */
      keyItem('F',0,0);
      break;
    case 'G':  /* toggle grid */
      keyItem('G',0,0);
      break;
    case 'H': /* hardcopy PPM */
      keyFile('H',0,0);
      break;
    case 'I':  /* follows mouse */
      sc->type ^= S_FOLLOW;
      sc1 = sc;
      while ( sc1->slave > -1 ) {
        sc1 = cv.scene[sc1->slave];
        sc1->type = sc->type;
      }
      break;
    case 'J':
      if ( mesh->dim == 2 && !(sc->mode & S_ALTITUDE) )  break;
      if ( p->pmode == PERSPECTIVE ) {
        p->pmode = CAMERA;  
        glutMouseFunc(mouseCamera);
        glutMotionFunc(motionCamera);
      }
      else {
        p->pmode = PERSPECTIVE;
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
      }
      sc->type = S_RESET;
      reshapeScene(sc->par.xs,sc->par.ys);
      post = TRUE;
      break;
    case 'K':
      keyMetric('K',0,0);
      break;
    case 'L':
      keyFile('L',0,0);
      break;
    case 'M':
      morphing = 1-morphing;
      post = TRUE;
      break;
    case 'N':  /* draw normals */
    case 'O':
      if ( mesh->nvn )
        keyItem(key,0,0);
      break;
    case 'P': /* toggle point nums */
      keyItem('P',0,0);
      break;
    case 'Q':
      /*keyMetric('q',0,0);*/
      break;
    case 'R':
      if ( ipilmat < 1 ) return;
      refmat = pilmat[ipilmat--];
      pm = &sc->material[refmat];
      pm->flag = 0;
      updatePoints(sc,mesh,refmat);
      dolist = TRUE;
      post   = TRUE;
      break;
    case 'S':  /* softcopy */
      /*keyFile('S',0,0);*/
      keyAnim('S',0,0);
      break;
    case 'V':  /* change center of scene */
      if ( !refitem )  return;
      changeCenter(sc,mesh);
      reshapeScene(sc->par.xs,sc->par.ys);
      post = TRUE;
      /* linked view */
      sc1 = sc;
      while ( sc1->slave > -1 ) {
    sc1 = cv.scene[sc1->slave];
        sc1->cx = sc->cx;
        sc1->cy = sc->cy;
        sc1->cz = sc->cz;
        memcpy(sc1->view,sc->view,sizeof(struct transform));
        memcpy(sc1->persp,sc->persp,sizeof(struct sperspective));
    glutSetWindow(sc1->idwin);
    reshapeScene(sc1->par.xs,sc1->par.ys);
      }
      glutSetWindow(sc->idwin);
      break;
    case 'W':
      keyFile('W',0,0);
      break;
    case 'X':  /* close window */
      if ( cube->active & C_EDIT ) {
        cube->cma[0] -= 0.1*sc->dmax;
        post = TRUE;
        break;
      }
      if ( idw != cv.nbs ) {
    deleteScene(sc);
    for (k=idw+1; k<cv.nbs; k++)
      cv.scene[k-1] = cv.scene[k];
    cv.scene[cv.nbs-1] = 0;
      }
      glutHideWindow();
      if ( --cv.nbs == 0 ) exit(0);
      break;

    case 'Y':  /* close window */
      if ( cube->active & C_EDIT ) {
        cube->cma[1] -= 0.1*sc->dmax;
        post = TRUE;
        break;
      }
    case 'Z':  /* zoom out */
      if ( cube->active & C_EDIT ) {
        cube->cma[2] -= 0.1*sc->dmax;
        post = TRUE;
        break;
      }
      else {
        if ( p->rubber == 2 )
          setPersp(sc,p,0);
        else 
          p->fovy = MEDIT_MIN(1.1*p->fovy,179.0);
        farclip(1);
      }
      post = TRUE;

      /* get linked view */
      sc1 = sc;
      while ( sc1->slave > -1 ) {
        sc1 = cv.scene[sc1->slave];
        memcpy(sc1->view,sc->view,sizeof(struct transform));
        memcpy(sc1->persp,sc->persp,sizeof(struct sperspective));
        glutSetWindow(sc1->idwin);
        reshapeScene(sc1->par.xs,sc1->par.ys);
      }
      glutSetWindow(sc->idwin);
      break;
    }
  }
  
  else {
    switch (key) {
    case '-':
      if (keyact & GLUT_ACTIVE_ALT ) {
        keyAnim('p',0,0);
        break;
      }
      if ( cube->active & C_EDIT ) {
        cube->cma[0] *= 0.95;
        cube->cma[1] *= 0.95;
        cube->cma[2] *= 0.95;
        post = TRUE;
        break;
      }  
      if ( sc->persp->pmode == CAMERA ) {
        cam = sc->camera;
        cam->spmod -= 1.e-04 * sc->dmax;
        post = TRUE;
        break;
      }
      if ( p->depth < -20.0*sc->dmax )  break;
      p->depth -= 0.1*sc->dmax;
      farclip(1);
      post = TRUE;

      /*
      if ( p->rubber == 2 )
        setPersp(sc,p,0);
      else 
        p->fovy = MEDIT_MIN(1.1*p->fovy,179.0);
      farclip(1);
      post = TRUE;
      */
      /* get linked view */
      sc1 = sc;
      while ( sc1->slave > -1 ) {
        sc1 = cv.scene[sc1->slave];
        memcpy(sc1->view,sc->view,sizeof(struct transform));
        memcpy(sc1->persp,sc->persp,sizeof(struct sperspective));
        glutSetWindow(sc1->idwin);
        reshapeScene(sc1->par.xs,sc1->par.ys);
      }
      glutSetWindow(sc->idwin);
      break;

    case '+':
      if (keyact & GLUT_ACTIVE_ALT ) {
        keyAnim('n',0,0);
        break;
      }
      if ( cube->active & C_EDIT ) {
    cube->cma[0] *= 1.05;
    cube->cma[1] *= 1.05;
    cube->cma[2] *= 1.05;
        post = TRUE;
        break;
      }  
      if ( sc->persp->pmode == CAMERA ) {
        cam = sc->camera;
        cam->spmod += 1.e-04 * sc->dmax;
        post = TRUE;
        break;
      }
      if ( p->depth > 0.0 )  break;
      p->depth += 0.1*sc->dmax;
      farclip(1);
      post = TRUE;

/*
      if ( p->rubber == 2 )
    setPersp(sc,p,1);
      else 
      p->fovy = MEDIT_MAX(0.9*p->fovy,1e-05);
      farclip(1);
      post = TRUE;
*/      
      /* update linked view */
      sc1 = sc;
      while ( sc1->slave > -1 ) {
    sc1 = cv.scene[sc1->slave];
        memcpy(sc1->view,sc->view,sizeof(struct transform));
        memcpy(sc1->camera,sc->camera,sizeof(struct camera));
        memcpy(sc1->persp,sc->persp,sizeof(struct sperspective));
    glutSetWindow(sc1->idwin);
    reshapeScene(sc1->par.xs,sc1->par.ys);
      }
      glutSetWindow(sc->idwin);
      break;

    case '#':  /* select entity */
      fprintf(stdout,"ENTITY NUMBER: "); fflush(stdout);
      fflush(stdin);  fscanf(stdin,"%d",&numit);
      if ( sc->picklist )  glDeleteLists(sc->picklist,1);
      if ( numit > 0 )
    sc->picklist = pickItem(mesh,sc,numit);
      post = TRUE;
      break;
    case '!':  /* clip plane */
    return;
      if ( !(clip->active & C_ON) )  return;
      dd = clip->eqn[3]-clip->eqn[0]*mesh->xtra \
         - clip->eqn[1]*mesh->ytra-clip->eqn[2]*mesh->ztra;
      fprintf(stdout,"\nCurrent plane: %gx %+gy %+gz %+g = 0\n",
              clip->eqn[0],clip->eqn[1],clip->eqn[2],dd);
      fprintf(stdout,"Plane coeffs : "); fflush(stdout);
      fflush(stdin);  fscanf(stdin,"%f %f %f %f",&a,&b,&c,&d);
      resetClip(sc,clip,mesh);
      clip->eqn[0] = a;
      clip->eqn[1] = b;
      clip->eqn[2] = c;
      clip->eqn[3] = d;

      fprintf(stdout,"New plane eq.: ");
      if ( clip->eqn[0] )
        fprintf(stdout,"%+gx",clip->eqn[0]);
      if ( clip->eqn[1] )
        fprintf(stdout," %+gy",clip->eqn[1]);
      if ( clip->eqn[2] )  
        fprintf(stdout," %+gz",clip->eqn[2]);
      if ( clip->eqn[3] ) 
        fprintf(stdout," %+g",clip->eqn[3]);
      fprintf(stdout," = 0\n");
      clip->eqn[3] += (a*mesh->xtra+b*mesh->ytra+c*mesh->ztra);
      post   = TRUE;
      dolist = TRUE;
      break;
    
    case '@': /* add trajectoire point */
      if ( p->pmode == CAMERA )
        pathAdd(sc,x,y);
      break;
    
    case '%':
      fprintf(stdout,"reference (%d): ",refpick); fflush(stdout);
      fflush(stdin);
      fscanf(stdin,"%d",&refpick);
      break;
    
    case '&':
      puts("ADJUST");
      parEdit(sc);
      break;
    }
  }
#ifdef IGL
keySceneFinish:
#endif
  
  if ( dolist ) {
    doLists(sc,mesh);
    doMapLists(sc,mesh,1);
    doIsoLists(sc,mesh,1);
  }
  if ( post )   glutPostRedisplay();
}

