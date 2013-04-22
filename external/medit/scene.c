#include "medit.h"
#include "extern.h"
#include "sproto.h"

#ifdef IGL
#  include "IGLParams.h"
#  include <igl/canonical_quaternions.h>
#  include <igl/quat_to_mat.h>
#  include <igl/mat_to_quat.h>
#  include <igl/snap_to_canonical_view_quat.h>
#  include <igl/C_STR.h>
#  include <igl/png/render_to_png_async.h>
#endif

extern GLboolean  hasStereo;
extern int       *pilmat,ipilmat,refmat,reftype,refitem;
extern short      schw,schh;
extern ubyte      quiet,fullscreen,tiling,stereoMode;


/* return current active scene */
int currentScene() {
  int  k,idw;

  idw = glutGetWindow();
  for (k=0; k<MAX_SCENE; k++) {
    if ( cv.scene[k] && idw == cv.scene[k]->idwin )
      return(k);
  }
  return(0);
}

/* check for OpenGL error */
void checkErrors(void) {
  GLenum error;

  while ( (error = glGetError()) != GL_NO_ERROR ) {
    fprintf(stderr,"  ## ERROR: %d: %s\n",
            (int)error,(char*)gluErrorString(error));
    exit(1);
  }
}

void farclip(GLboolean b) {
  pScene  sc;
  pPersp  p;
  pCamera c;
  float   look[3],ratio,units;
  int     idw = currentScene();
  static  GLfloat up[3] = { 0.0, 1.0, 0.0};

  /* default */
  sc    = cv.scene[idw];
  p     = sc->persp;
  c     = sc->camera;
  ratio = (GLfloat)sc->par.xs / sc->par.ys;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  switch (p->pmode) {
  case ORTHO:
    glOrtho(-1.,1.,-1.,0.1,0.01,0.01);
    break;

  case PERSPECTIVE:
    if ( b )
      gluPerspective(p->fovy,ratio,sc->dmax,4.0*sc->dmax);
    else
      gluPerspective(p->fovy,ratio,0.01,10000.0*sc->dmax);

    units = 1.e-02;
    glPolygonOffset(1.0, units);
    break;

  case CAMERA:
    gluPerspective(p->fovy,ratio,0.001*sc->dmax,4.0*sc->dmax);
    look[0] = c->eye[0] + 0.001*sc->dmax*c->speed[0];
    look[1] = c->eye[1] + 0.001*sc->dmax*c->speed[1];
    look[2] = c->eye[2] + 0.001*sc->dmax*c->speed[2];
    gluLookAt(c->eye[0],c->eye[1],c->eye[2], 
              look[0],look[1],look[2],
              up[0],up[1],up[2]);
    break;
  }

  /* zoom transformation */
  if ( p->rubber == 2 ) {
    glPushMatrix();
    glLoadIdentity();
    glRotatef(-p->gamma,1.,0.,0.);
    glRotatef(p->alpha,0.,1.,0.);
    glMultMatrixf(p->matrix);
    glGetFloatv(GL_PROJECTION_MATRIX,p->matrix);
    glPopMatrix();
    p->rubber = 0;
  }

  /* apply transformation */
  glMultMatrixf(p->matrix);

  glMatrixMode(GL_MODELVIEW);
}


void reshapeScene(int width,int height) {
  pScene   sc;

  if ( ddebug ) printf("reshape scene\n");
  sc = cv.scene[currentScene()];
  sc->par.xs = width;
  sc->par.ys = height;

  glViewport(0,0,width,height);
  farclip(GL_TRUE);
#ifdef IGL
  // Tell AntTweakBar about new size
  TwWindowSize(width, height);
  // Keep AntTweakBar on right side of screen and height == opengl height
  // get the current position of a bar
  int size[2];
  sc->rebar.TwGetParam(NULL, "size", TW_PARAM_INT32, 2, size);
  int pos[2];
  // Place bar on right side of opengl rect (padded by 10 pixels)
  pos[0] = MEDIT_MAX(10,(int)width - size[0] - 10);
  // place bar at top (padded by 10 pixels)
  pos[1] = 10;
  // Set height to new height of window (padded by 10 pixels on bottom)
  size[1] = height-pos[1]-10;
  sc->rebar.TwSetParam( NULL, "position", TW_PARAM_INT32, 2, pos);
  sc->rebar.TwSetParam( NULL, "size", TW_PARAM_INT32, 2,size);
#endif
}


static void drawList(pScene sc,int clip,int map) {
  pMesh   mesh = cv.mesh[sc->idmesh];
  ubyte   elev = sc->mode & S_ALTITUDE;

  if ( ddebug ) printf("drawList %p %d %d\n",sc,clip,map);
  if ( mesh->dim == 2 && !elev ) glDisable(GL_DEPTH_TEST);

  glLineWidth(1.0);
  if ( clip ) {
    if ( map ) {
      if ( sc->cmlist[LTets] ) glCallList(sc->cmlist[LTets]);
      if ( sc->cmlist[LHexa] ) glCallList(sc->cmlist[LHexa]);
    }
    else {
      if ( sc->clist[LTets] ) glCallList(sc->clist[LTets]);
      if ( sc->clist[LHexa] ) glCallList(sc->clist[LHexa]);
    }
  }
  else if ( map ) {
    if ( mesh->nt+mesh->nq ) {
//#ifdef IGL
//      // Hack for winding figures to show usual mesh
//      if ( sc->dlist[LTria] ) glCallList(sc->dlist[LTria]);
//#else
      if ( sc->mlist[LTria] ) glCallList(sc->mlist[LTria]);
//#endif
      if ( sc->mlist[LQuad] ) glCallList(sc->mlist[LQuad]);
    }
    else {
      if ( sc->mlist[LTets] ) glCallList(sc->mlist[LTets]);
      if ( sc->mlist[LHexa] ) glCallList(sc->mlist[LHexa]);
    }
  }
  else {
    if ( mesh->nt+mesh->nq ) {
      if ( sc->dlist[LTria] ) glCallList(sc->dlist[LTria]);
      if ( sc->dlist[LQuad] ) glCallList(sc->dlist[LQuad]);
    }
    else {
      if ( sc->dlist[LTets] ) glCallList(sc->dlist[LTets]);
      if ( sc->dlist[LHexa] ) glCallList(sc->dlist[LHexa]);
    }
  }
  if ( mesh->dim == 2 && !elev ) glEnable(GL_DEPTH_TEST);
}

#ifdef ppc
void bogusQuad(pScene sc) {
  /* bogus polygon (nvidia) */
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLineWidth(1.0);
    glColor3fv(sc->par.back);
    glBegin(GL_QUADS);
      glVertex3f(0., 0.,-sc->persp->depth);
      glVertex3f(0., 0.,-sc->persp->depth);
      glVertex3f(0., 0.,-sc->persp->depth);
      glVertex3f(0., 0.,-sc->persp->depth);
    glEnd();
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}
#endif

static void displayScene(pScene sc,int mode,int clip) {
  int     map;
  
  map = mode & S_MAP;


  switch(mode) {
  case FILL:  /* solid fill */
    if ( ddebug ) printf("solid fill\n");
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glDisable(GL_POLYGON_OFFSET_FILL);
      drawList(sc,clip,0);
    glDisable(GL_LIGHTING);
    break;

  case WIRE:  /* basic wireframe */
  case WIRE+S_MATERIAL:
    if ( ddebug ) printf("wireframe\n");
#ifdef ppc
    bogusQuad(sc);
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glColor4fv(sc->par.line);
    glDisable(GL_POLYGON_OFFSET_FILL);
      drawList(sc,clip,0);
    break;

  case DEPTH:  /* depth wireframe */
  case DEPTH + S_MATERIAL:
    if ( ddebug ) printf("depth wireframe\n");
    glEnable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      drawList(sc,clip,0);
    glDisable(GL_LIGHTING);
    break;

  case HIDDEN: /* hidden lines removal */
  case HIDDEN + S_MATERIAL:
    if ( ddebug ) printf("hidden lines\n");
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) glDisable(GL_POLYGON_OFFSET_FILL);
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glColor3fv(sc->par.back);
      drawList(sc,clip,0);
    glDisable(GL_POLYGON_OFFSET_FILL);
#ifdef ppc
    bogusQuad(sc);
#endif
#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) break;
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glColor4fv(sc->par.line);
      drawList(sc,clip,0);
    break;

  case SHADED: /* shaded polygons */
  case SHADED+S_MATERIAL:
    if ( ddebug ) printf("shaded polygons\n");
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) glDisable(GL_POLYGON_OFFSET_FILL);
#endif
      drawList(sc,clip,0);
    glDisable(GL_LIGHTING);
    glDisable(GL_POLYGON_OFFSET_FILL);
#ifdef ppc
    bogusQuad(sc);
#endif

#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) break;
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glColor4fv(sc->par.line);
      drawList(sc,clip,0);
    break;

  case SIZEMAP: /* display metric map */
  case SIZEMAP+S_MATERIAL:
    if ( ddebug ) printf("display sizemap\n");
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      drawList(sc,clip,map);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_POLYGON_OFFSET_FILL);
#ifdef ppc
    bogusQuad(sc);
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glColor4fv(sc->par.line);
    glLineWidth(1.0);
      drawList(sc,clip,0);
    if ( sc->mode & S_ALTITUDE ) {
      glColor4fv(sc->par.line);
      if ( sc->mlist[LTets] ) glCallList(sc->mlist[LTets]);
      if ( sc->mlist[LHexa] ) glCallList(sc->mlist[LHexa]);
    }
    break;

  default: /* other modes */
    if ( ddebug ) printf("rendering mode %d\n",sc->mode);
    /* interior */
    if ( sc->mode & S_FILL ) {
      if ( sc->mode & S_COLOR )  glEnable(GL_LIGHTING);
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glEnable(GL_POLYGON_OFFSET_FILL);
#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) glDisable(GL_POLYGON_OFFSET_FILL);
#endif
#ifdef IGL
      if ( sc->mode & S_MAP && (~(sc->clip->active & C_ON) || clip) ) {
        glEnable(GL_COLOR_MATERIAL);
      }
#else
    if ( sc->mode & S_MAP ) {
      glEnable(GL_COLOR_MATERIAL);
#endif
#ifdef IGL
    if(!sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) glDisable(GL_LIGHTING);
#endif
	drawList(sc,clip,map);
        glDisable(GL_COLOR_MATERIAL);
      }
      else {
 	glColor4fv(sc->par.back);
	drawList(sc,clip,0);
      }

    /* boundary */
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_POLYGON_OFFSET_FILL);
#ifdef ppc
    bogusQuad(sc);
#endif
#ifdef IGL
    if(
        //(mode==31) // Seems to be data+cap
        //&&
        !sc->igl_params->lines_on_cap && (clip && sc->clip->active & C_CAP)) break;
#endif
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    if ( !(sc->mode & S_BDRY) )  break;
    if ( sc->mode & S_COLOR && !(sc->mode & S_FILL) ) {
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_LIGHTING);
    }
    if ( sc->mode & S_MAP) {
      if ( sc->mode & S_FILL ) {
        glColor4fv(sc->par.line);
        drawList(sc,clip,0);
      }
      else
        drawList(sc,clip,map);
    }
    else if ( sc->mode & S_ALTITUDE ) {
      glColor4fv(sc->par.line);
      drawList(sc,clip,map);
      if ( sc->mlist[LTets] ) glCallList(sc->mlist[LTets]);
      if ( sc->mlist[LHexa] ) glCallList(sc->mlist[LHexa]);
    }
    else {
      glColor4fv(sc->par.line);
      drawList(sc,clip,0);
    }
  }
}

static void displayData(pScene sc,pMesh mesh) {
  int  kk;
  
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);

  /* iso-lines */
  if ( sc->isotyp & S_ISOLINE ) {
    glLineWidth(1.0);
    if ( sc->ilist[LTria] )  glCallList(sc->ilist[LTria]);
    if ( sc->ilist[LQuad] )  glCallList(sc->ilist[LQuad]);
    glDisable(GL_CLIP_PLANE0);
    if ( sc->ilist[LTets] && sc->clip->active & C_ON )
      glCallList(sc->ilist[LTets]);
  }

  /* vector field */
  if ( sc->isotyp & S_VECTOR ) {
    if ( mesh->dim == 2 ) {
      if ( sc->vlist[LTria] )  glCallList(sc->vlist[LTria]);
      if ( sc->vlist[LQuad] )  glCallList(sc->vlist[LQuad]);
    }
    else {
      if ( sc->clip->active & C_ON ) {
        glDisable(GL_CLIP_PLANE0);
        if ( sc->vlist[LTets] )  glCallList(sc->vlist[LTets]);
        if ( sc->vlist[LHexa] )  glCallList(sc->vlist[LHexa]);
      }
      else if (mesh->ntet+mesh->nhex == 0 )
        if ( sc->vlist[LTria] )  glCallList(sc->vlist[LTria]);
    }
  }

  /* streamlines */
  if ( sc->isotyp & S_CRITP && sc->cplist )
    glCallList(sc->cplist);
  if ( sc->isotyp & S_STREAML ) {
    for (kk=0; kk<sc->stream->nbstl; kk++)
      glCallList(sc->slist[kk]);
  }
  else if ( sc->isotyp & S_PARTICLE ) {
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    displayParticle(sc,mesh);
    glDisable(GL_COLOR_MATERIAL);
  }

  /* iso-surfaces */
  if ( sc->isotyp & S_ISOSURF ) {
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    if ( sc->ilist[LTets] )  glCallList(sc->ilist[LTets]);
    if ( sc->ilist[LHexa] )  glCallList(sc->ilist[LHexa]);
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  }
  
}


void setupView(pScene sc) {
  pScene       slave;
  pMesh        mesh;
  pTransform   view;
  pPersp       p;
  pCamera      c;
  int          clvol;

  /* default */
  if ( ddebug )  fprintf(stdout,"setupView\n");
  mesh = cv.mesh[sc->idmesh];
  view = sc->view;
  p    = sc->persp;
  c    = sc->camera;

  /* init transformation matrix */
  if ( sc->type & S_RESET ) {
    glPushMatrix();
    glLoadIdentity();
    if ( p->pmode != CAMERA ) {
      if ( mesh->dim == 3 || sc->mode & S_ALTITUDE ) {
        glRotatef(-60.,1.,0.,0.);
        glRotatef(-120.,0.,0.,1.);
      }
    }
    else {
      if ( c->vecup == X_AXIS )
        glRotatef(90.0,0.0,0.0,1.0);
      else if ( c->vecup == Z_AXIS )
        glRotatef(-90.0,1.0,0.0,0.0);
    }
    glGetFloatv(GL_MODELVIEW_MATRIX,view->matrix);
    glPopMatrix();
    sc->type ^= S_RESET;
  }

  /* keep old transformation */
  memcpy(view->oldmat,view->matrix,16*sizeof(float));

  /* compute new transformation */
  glPushMatrix();
  glLoadIdentity();
  if ( p->pmode != CAMERA ) {
    glTranslatef(view->panx,view->pany,0.0);
    if(ddebug)
    {
      printf("%g %g %g\n", -view->panx,-view->pany,0.);
    }
    if ( mesh->dim == 3 || sc->mode & S_ALTITUDE )
    {
      glRotatef(view->angle,view->axis[0],view->axis[1],view->axis[2]);
      if(ddebug)
      {
        printf("%g, %g %g %g\n",
            view->angle,view->axis[0],view->axis[1],view->axis[2]);
      }
    }
    glTranslatef(-view->opanx,-view->opany,0.);
    if(ddebug)
    {
      printf("%g %g %g\n",
          -view->opanx,-view->opany,0.);
    }
    glMultMatrixf(view->matrix);
    glGetFloatv(GL_MODELVIEW_MATRIX,view->matrix);
    if(ddebug)
    {
      int i,j;
      for(i = 0;i<4;i++)
      {
        for(j = 0;j<4;j++)
        {
          printf("%g ",view->matrix[i+4*j]);
        }
        printf("\n");
      }
    }
  }
  else if ( animate ) {
    c->eye[0] += c->spmod*c->speed[0];
    c->eye[1] += c->spmod*c->speed[1];
    c->eye[2] += c->spmod*c->speed[2];
    animateCamera();
    reshapeScene(sc->par.xs,sc->par.ys);
  }
  glPopMatrix();

  /* keep old translation */
  view->opanx = view->panx;
  view->opany = view->pany;

  /* copy views */
  if ( !animate && sc->slave > -1 ) {
    slave = cv.scene[sc->slave];
    memcpy(slave->view,sc->view,sizeof(struct transform));
    memcpy(slave->camera,sc->camera,sizeof(struct camera));
    slave->view->angle = 0.0f;
    clvol = slave->clip->active & C_VOL;
    memcpy(slave->clip,sc->clip,sizeof(struct clip));
    if ( clvol )  slave->clip->active |= C_VOL;
  }
}


void drawModel(pScene sc) {
  pMesh        mesh;
  pTransform   view;
  pClip        clip;
  ubyte        sstatic;

  /* default */
  mesh = cv.mesh[sc->idmesh];
  view = sc->view;
  clip = sc->clip;
  if ( ddebug ) printf("\n-- redraw scene %d, mesh %d\n",sc->idwin,sc->idmesh);

#ifdef IGL
  const bool hot_dog_view = sc->igl_params->hot_dog_view;
  const double width = sc->igl_params->width(mesh);
  const double hot_dog_ratio = sc->igl_params->hot_dog_ratio;
  int eff_slices = 1;
  if(hot_dog_view)
  {
    eff_slices = sc->igl_params->num_hot_dog_slices+1;
  }
  for(int h = 0; h < eff_slices ;h+=2)
  {
#endif

  glDisable(GL_LIGHTING);

  /* draw clipping plane */
  if ( clip->active & C_ON ) {
#ifdef IGL
    if(hot_dog_view)
    {
      drawClip(sc,clip,mesh,0);
      // first clipping plane
      GLdouble eqn1[4];
      for(int c = 0;c<3;c++) eqn1[c] = clip->eqn[c];
      eqn1[3] = clip->eqn[3] + width*h;
      
      glClipPlane(GL_CLIP_PLANE0,eqn1);
      glEnable(GL_CLIP_PLANE0);

      // second clipping plane
      if(h>0)
      {
        GLdouble eqn2[4];
        for(int c = 0;c<3;c++) eqn2[c] = -clip->eqn[c];
        eqn2[3] = -clip->eqn[3] - width*(h) + (1.-hot_dog_ratio)*width*2;
        //eqn2[3] = -clip->eqn[3] - width*(h-1);

        glClipPlane(GL_CLIP_PLANE1,eqn2);
        glEnable(GL_CLIP_PLANE1);
      }else
      {
        glDisable(GL_CLIP_PLANE1);
      }

    }else
    {
#endif
    drawClip(sc,clip,mesh,0);
    glClipPlane(GL_CLIP_PLANE0,clip->eqn);
    glEnable(GL_CLIP_PLANE0);
#ifdef IGL
    }
#endif
  }
  else
  {
#ifdef IGL
    glDisable(GL_CLIP_PLANE1);
#endif
    glDisable(GL_CLIP_PLANE0);
  }

  /* draw object if static scene */
  sstatic = view->mstate > 0 && clip->cliptr->mstate > 0;
  if ( sstatic || sc->type & S_FOLLOW ) {
    displayScene(sc,sc->mode,0);
    if ( sc->item & S_NUMP || sc->item & S_NUMF )  listNum(sc,mesh);
    /* draw normals */
    if ( sc->type & S_NORMAL ) {
      if ( !sc->nlist ) sc->nlist = drawNormals(mesh,sc);
      glCallList(sc->nlist);
    }
    /* draw data */
    if ( sstatic )  displayData(sc,mesh);
  }
  else if ( !(sc->item & S_BOX) )
    drawBox(sc,mesh,0);

  /* draw ridges, corners, etc. */
  if ( (sc->item & S_GEOM) && sc->glist ) {
    glDisable(GL_LIGHTING);
    if ( !mesh->ne )
      glPointSize(1);
    else
      glPointSize(5);
    glDisable(GL_COLOR_MATERIAL);
    glCallList(sc->glist);
  }
#ifdef IGL
  }
  glDisable(GL_CLIP_PLANE1);
#endif

  glDisable(GL_CLIP_PLANE0);
#if IGL
  if ( sc->item & S_BOX )
#else
  if ( clip->active & C_EDIT || sc->item & S_BOX )
#endif
    drawBox(sc,mesh,0);
  if ( sc->item & S_AXIS )
    drawAxis(sc,mesh->dim);
  if ( (mesh->dim == 3 || sc->mode & S_ALTITUDE) && sc->item & S_GRID )
    drawBase(sc,mesh);
  if ( sc->cube->active & C_ON )
    drawCube(sc,mesh);

  //glColorMask(false, false, false, true);//This ensures that only alpha will be effected
  //glClearColor(0, 0, 0, 1);//alphaValue - Value to which you need to clear
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //glColorMask(true, true, true, true);//This ensures that only alpha will be effected

  sstatic |= tiling;
  if ( sstatic && clip->active & C_ON && clip->active & C_VOL )
    displayScene(sc,sc->mode,1);

  if ( sc->picklist && !(sc->isotyp & S_PARTICLE) ) {
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glDisable(GL_POLYGON_OFFSET_FILL);
      glCallList(sc->picklist);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }

  /* show path, if any */
  if ( sc->type & S_PATH && sc->path.tlist )
    glCallList(sc->path.tlist);
}


void redrawMorphing(pScene sc) {
  pMesh   mesh;
  
  if ( morphing ) {
    mesh = cv.mesh[sc->idmesh];
    if ( !morphMesh(sc,mesh) )  return;
  }
}


/* animate */
void glutIdle(void) {
  static float timePassed = 0.0;
  static float timeRedraw = 0.02;  /*1.0 / 50.0*/
  float  clk,timeElaps;

  clk = clock();
  timeElaps  = (clk - timePassed) / (float)CLOCKS_PER_SEC;
  if ( timeElaps >= timeRedraw ) {
    timePassed = clk;
    glutPostRedisplay();
  }
}


void streamIdle(void) {
  pScene         sc;
  pMesh          mesh;
  float          elp;
  clock_t        tim;
  static clock_t timbase= 0;
  static float   maxtim = 0.;

  sc  = cv.scene[currentScene()];
  sc->par.cumtim += sc->par.dt;
  sc->par.advtim  = 1;
  if ( sc->par.cumtim > sc->par.maxtime ) {
    sc->par.cumtim     = 0.0;
    sc->par.cumpertime = 0.0;
    saveimg = 0;
    sc->isotyp &= ~S_PARTICLE;
    glutIdleFunc(0);
    printf("\nfin");
  }
  else if ( sc->par.cumpertime >= sc->par.pertime ) {
    mesh = cv.mesh[sc->idmesh];
    sc->par.cumpertime = 0.0;
    if ( !animParticle(sc,mesh) ) {
      sc->par.cumtim     = 0.0;
      sc->par.cumpertime = 0.0;
      saveimg = 0;
      glutIdleFunc(0);
      printf("\nfin");
    }
    else
      glutPostRedisplay();
  }
  else
    glutPostRedisplay();

  return;

  if ( timbase < 1 ) {
    timbase = clock();
    maxtim  = 0.0;
    sc->par.advtim = 0;
    glutPostRedisplay();
  }
  else {
    tim = clock();
    elp = (tim - timbase) / (float)CLOCKS_PER_SEC;
    if ( elp > sc->par.dt ) {
      timbase = tim;
      maxtim += elp;  
      sc->par.advtim  = 1;
      if ( maxtim > sc->par.maxtime )
        glutIdleFunc(0);
    }
    sc->par.cumtim += sc->par.dt;
    glutPostRedisplay();
  }
}


/* OpenGL callbacks */
void redrawScene() {
  pScene       sc,slave;
  pTransform   view;
  pPersp       p;
  pCamera      c;
  double       ndfl,ratio,top,bottom,left,right,nnear,ffar;

  sc   = cv.scene[currentScene()];
  view = sc->view;
  p    = sc->persp;
  c    = sc->camera;

  if ( stereoMode == MONO || !hasStereo ) {
    glDrawBuffer(GL_BACK_LEFT);
    glClearColor(sc->par.back[0],sc->par.back[1],
                 sc->par.back[2],0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.,0.,-p->depth, 0.,0.,0., 0.0,1.0,0.0);

    setupView(sc);
    glMultMatrixf(view->matrix);
    glTranslatef(sc->cx,sc->cy,sc->cz);
    drawModel(sc);
    if ( sc->type & S_DECO )  redrawStatusBar(sc);
  }

  else {
    nnear   = -p->depth - 0.5 * sc->dmax;
    if ( nnear < 0.1 )  nnear = 0.1;
    ffar    = -p->depth + 0.5 * sc->dmax;
    ratio   = sc->par.xs / (double)sc->par.ys;
    top     = nnear * tan(DTOR * 0.5 * p->fovy);
    ndfl    = nnear / p->depth;
    if ( sc->par.eyesep < 0.0 )
      sc->par.eyesep = fabs(p->depth / 20.0);

    /* left view */
    glDrawBuffer(GL_BACK_LEFT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    left  = -ratio * top + 0.5 * sc->par.eyesep * ndfl;
    right =  ratio * top + 0.5 * sc->par.eyesep * ndfl;
    bottom= -top;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(left,right,top,bottom,nnear,ffar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(-sc->par.eyesep,0.,-p->depth, 
              sc->par.eyesep/3.0,0.,0., 
              0.0,1.0,0.0);

    setupView(sc);
    glMultMatrixf(view->matrix);
    glTranslatef(sc->cx,sc->cy,sc->cz);
    drawModel(sc);
    if ( sc->type & S_DECO )  redrawStatusBar(sc);

    /* right view */
    glDrawBuffer(GL_BACK_RIGHT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    left  = -ratio * top - 0.5 * sc->par.eyesep * ndfl;
    right =  ratio * top - 0.5 * sc->par.eyesep * ndfl;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(left,right,top,bottom,nnear,ffar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(sc->par.eyesep,0.,-p->depth, 
              sc->par.eyesep/3.0,0.,0., 
              0.0,1.0,0.0);

    setupView(sc);
    glMultMatrixf(view->matrix);
    glTranslatef(sc->cx,sc->cy,sc->cz);
    drawModel(sc);
    if ( sc->type & S_DECO )  redrawStatusBar(sc);
  }
#ifdef IGL
  if(sc->igl_params->render_on_next)
  {
    printf("FUCK SHIT\n");
    char path[1024];
    pMesh mesh;
    mesh = cv.mesh[sc->idmesh];
    strcpy(path,mesh->name);
    static int numdep = 0;
#define PNG "png"
    numdep = filnum(path,numdep,PNG);
    if ( numdep == -1)
    {
      glClearColor(1,0,0,1);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }else
    {
      sprintf(path,"%s.%.3d." PNG,path,numdep);
      igl::render_to_png_async(std::string(path),sc->par.xs,sc->par.ys,true,false);
    }
    sc->igl_params->render_on_next = false;
  }
  TwDraw();
#endif

  /* refresh screen */
  if ( saveimg && animate )
    glFlush();
  else
    glutSwapBuffers();

  if ( ddebug ) checkErrors();

  if ( saveimg && !(sc->type & S_SCISSOR) )  keyFile('H',0,0);

  /* redraw linked scene */
  if ( !animate && sc->slave > -1 ) {
    slave = cv.scene[sc->slave];
    glutSetWindow(slave->idwin);
    redrawScene();
  }
}


/* OpenGL callbacks */
void redrawSchnauzer() {
  pScene  sc = cv.scene[currentScene()];
  pMesh   mesh;
  char   *ptr,data[256];

  mesh = cv.mesh[sc->idmesh];
  strcpy(data,mesh->name);
  ptr = (char*)strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  ptr = (char*)strstr(data,".gis");
  if ( ptr ) *ptr = '\0';

  strcat(data,".ppm");
  redrawScene();
  imgHard(sc,data,'H');

  exit(0);
}


void deleteScene(pScene sc) {
  // Alec: This isn't actually called on exit
  /* default */
  if ( ddebug) printf("deleteScene\n");

  // try to save rebar
#ifdef IGL
  sc->rebar.save("rebar.rbr");
#endif

  M_free(sc->view);
  M_free(sc->clip);
  M_free(sc->persp);
  M_free(sc->camera);
  M_free(sc->material);
  M_free(sc->matsort);
  M_free(sc);
}

#ifdef IGL
void initIGL(pScene sc,pMesh mesh)
{
  sc->igl_params = new IGLParams();
}

// No-op setter, does nothing
void TW_CALL no_op(const void *value, void *clientData){}
// No-op getter, does nothing
void TW_CALL no_op(void *value, void *clientData){}

// Snap XY plane
void TW_CALL view_xy_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::XY_PLANE_QUAT_F, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->view->matrix[i*4+j] = mat[i*4+j];
}
// Snap XZ plane
void TW_CALL view_xz_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::XZ_PLANE_QUAT_F, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->view->matrix[i*4+j] = mat[i*4+j];
}
// Snap YZ plane
void TW_CALL view_yz_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::YZ_PLANE_QUAT_F, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->view->matrix[i*4+j] = mat[i*4+j];
}
void TW_CALL snap_rotation_to_canonical_view_quat(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float q[4];
  igl::mat4_to_quat(sc->view->matrix,q);
  igl::snap_to_canonical_view_quat<float>(q,1.0,q);
  float mat[16];
  igl::quat_to_mat(q, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->view->matrix[i*4+j] = mat[i*4+j];
}
void TW_CALL set_rotation(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  const float * q = (const float*)(value);
  float mat[16];
  igl::quat_to_mat(q, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->view->matrix[i*4+j] = mat[i*4+j];
}
void TW_CALL get_rotation(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float * q = (float *)(value);
  igl::mat4_to_quat(sc->view->matrix,q);
}

// Update clip after it changed
void update_clip(pScene sc)
{
  if ( sc->clip->active & C_ON )
    sc->clip->active |= C_UPDATE;
}
// Snap XY plane
void TW_CALL clip_xy_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::CANONICAL_VIEW_QUAT_F[5], mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->clip->cliptr->rot[i*4+j] = mat[i*4+j];
  update_clip(sc);
}
// Snap XZ plane
void TW_CALL clip_xz_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::XZ_PLANE_QUAT_F, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->clip->cliptr->rot[i*4+j] = mat[i*4+j];
  update_clip(sc);
}
// Snap YZ plane
void TW_CALL clip_yz_plane(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float mat[16];
  igl::quat_to_mat( igl::YZ_PLANE_QUAT_F, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->clip->cliptr->rot[i*4+j] = mat[i*4+j];
  update_clip(sc);
}
void TW_CALL snap_rotation_to_canonical_clip_quat(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float q[4];
  igl::mat4_to_quat(sc->clip->cliptr->rot,q);
  igl::snap_to_canonical_view_quat<float>(q,1.0,q);
  float mat[16];
  igl::quat_to_mat(q, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->clip->cliptr->rot[i*4+j] = mat[i*4+j];
  update_clip(sc);
}

void TW_CALL set_clip_rotation(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  const float * q = (const float*)(value);
  float mat[16];
  igl::quat_to_mat(q, mat);
  // Copy 3x3 linear block
  for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++)
      sc->clip->cliptr->rot[i*4+j] = mat[i*4+j];
  update_clip(sc);
}
void TW_CALL get_clip_rotation(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float * q = (float *)(value);
  igl::mat4_to_quat(sc->clip->cliptr->rot,q);
}
void TW_CALL set_clip_tra12(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  const float * tra12 = (const float*)(value);
  sc->clip->cliptr->tra[12] = tra12[0];
  sc->clip->cliptr->tra[12] = tra12[1];
  sc->clip->cliptr->tra[13] = tra12[2];
  update_clip(sc);
}
void TW_CALL get_clip_tra12(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  float * tra12 = (float *)(value);
  tra12[0] = sc->clip->cliptr->tra[12];
  tra12[1] = sc->clip->cliptr->tra[12];
  tra12[2] = sc->clip->cliptr->tra[13];
}

void TW_CALL invert_clip(void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  invertClip(sc,sc->clip);
  sc->clip->active |= C_REDO;
}

void TW_CALL set_hot_dog_view(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  sc->igl_params->hot_dog_view = *(const bool*)(value);
  printf("set_hot_dog_view\n");
  sc->clip->active |= C_REDO;
}
void TW_CALL get_hot_dog_view(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  *(bool *)(value) = sc->igl_params->hot_dog_view;
  sc->clip->active |= C_REDO;
}

void TW_CALL set_num_hot_dog_slices(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  sc->igl_params->num_hot_dog_slices = *(const int*)(value);
  sc->clip->active |= C_REDO;
}
void TW_CALL get_num_hot_dog_slices(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  *(int *)(value) = sc->igl_params->num_hot_dog_slices;
  sc->clip->active |= C_REDO;
}

void TW_CALL set_hot_dog_ratio(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  sc->igl_params->hot_dog_ratio = *(const double*)(value);
  sc->clip->active |= C_REDO;
}
void TW_CALL get_hot_dog_ratio(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  *(double *)(value) = sc->igl_params->hot_dog_ratio;
  sc->clip->active |= C_REDO;
}

void TW_CALL set_C_HIDE(const void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  bool v = *(const bool*)(value);
  if(v != (0!=(sc->clip->active & C_HIDE)))
  {
    sc->clip->active ^= C_HIDE;
  }
}
void TW_CALL get_C_HIDE(void *value, void *clientData)
{
  pScene sc = static_cast<const pScene>(clientData);
  *(bool *)(value) = sc->clip->active & C_HIDE;
}

void initAntTweakBar(pScene sc,pMesh mesh)
{
  printf("initAntTweakBar\n");
  TwInit(TW_OPENGL, NULL);
  sc->rebar.TwNewBar("bar");
  TwDefine("bar label='Release' size='200 550' text=light alpha='200' color='68 68 68'");
  TwGLUTModifiersFunc(glutGetModifiers);

  sc->rebar.TwAddVarRW(
    "camera_eye",
    TW_TYPE_DIR3F,
    sc->camera->eye,
    "group='View' ");
  sc->rebar.TwAddVarRW(
    "back",
    TW_TYPE_COLOR4F,
    sc->par.back,
    "group='View' "
    "label='background color' "
    "open "
    "colormode=hls ");
  sc->rebar.TwAddVarRW(
    "line",
    TW_TYPE_COLOR4F,
    sc->par.line,
    "group='View' "
    "label='line color' "
    "open "
    "colormode=hls ");

  sc->rebar.TwAddVarRW(
    "lines_on_cap",
    TW_TYPE_BOOLCPP,
    &(sc->igl_params->lines_on_cap),
    "group='View' "
    "help='Display lines on cap when showing lines' ");

  sc->rebar.TwAddVarCB(
    "hot_dog_view",
    TW_TYPE_BOOLCPP,
    set_hot_dog_view,
    get_hot_dog_view,
    sc,
    "group='View' "
    "key=H "
    "help='Display many evenly spaced cuts' ");

  sc->rebar.TwAddVarCB(
    "num_hot_dog_slices",
    TW_TYPE_INT32,
    set_num_hot_dog_slices,
    get_num_hot_dog_slices,
    sc,
    C_STR("group='View' "
    "min=1 max="<<MAX_HOT_DOG_SLICES<<" step=1"
    "help='number of hot_dog_slices' "));

  sc->rebar.TwAddVarCB(
    "hot_dog_ratio",
    TW_TYPE_DOUBLE,
    set_hot_dog_ratio,
    get_hot_dog_ratio,
    sc,
    "group='View' "
    "min=0 max=1 step=0.1"
    "help='ratio of in to out hot dog slice sizes' ");

  sc->rebar.TwAddButton(
    "view_xy",
    view_xy_plane,
    sc,
    "group='View' "
    "key='z' "
    "help='Set view rotation to view xy plane' ");
  sc->rebar.TwAddButton(
    "view_xz",
    view_xz_plane,
    sc,
    "group='View' "
    "key='y' "
    "help='Set view rotation to view xz plane' ");
  sc->rebar.TwAddButton(
    "view_yz",
    view_yz_plane,
    sc,
    "group='View' "
    "key='x' "
    "help='Set view rotation to view yz plane' ");
  sc->rebar.TwAddButton(
    "snap",
    snap_rotation_to_canonical_view_quat,
    sc,
    "group='View' "
    "key='Z' "
    "help='Snap view rotation to nearest canonical view' ");
  sc->rebar.TwAddVarCB(
    "rotation",
    TW_TYPE_QUAT4F,
    set_rotation,
    get_rotation,
    sc,
    "group='View' open ");
  sc->rebar.TwAddVarRW(
    "fov_y",
    TW_TYPE_FLOAT,
    &(sc->persp->fovy),
    "group='View' readonly=true");

  sc->rebar.TwAddVarRW(
    "open_color",
    TW_TYPE_COLOR3F,
    &(sc->igl_params->open_color),
    "group=View "
    "help='Color of open boundaries.' colormode=hls");

  sc->rebar.TwAddVarRW(
    "line_width",
    TW_TYPE_FLOAT,
    &(sc->par.linewidth),
    "group=View "
    "help='Width of lines.' ");

  sc->rebar.TwAddVarRW(
    "dot_size",
    TW_TYPE_FLOAT,
    &(sc->igl_params->dot_size),
    "group=View "
    "help='size of dots.' ");

  sc->rebar.TwAddVarRW(
    "nme_color",
    TW_TYPE_COLOR3F,
    &(sc->igl_params->nme_color),
    "group=View "
    "help='Color of nonmanifold edges.' colormode=hls");

  sc->rebar.TwAddVarRW(
    "tet_color",
    TW_TYPE_COLOR4F,
    &(sc->igl_params->tet_color),
    "group=View "
    "help='Color of tets.' colormode=hls");

  sc->rebar.TwAddVarRW(
    "fade_flip",
    TW_TYPE_BOOLCPP,
    &(sc->igl_params->fade_flip),
    "group=Fade ");
  sc->rebar.TwAddVarRW(
    "fade_max_s",
    TW_TYPE_DOUBLE,
    &(sc->igl_params->fade_max_s),
    "step=0.01 max=1 min=0 group=Fade ");
  sc->rebar.TwAddVarRW(
    "fade_min_s",
    TW_TYPE_DOUBLE,
    &(sc->igl_params->fade_min_s),
    "step=0.01 max=1 min=0 group=Fade ");
  sc->rebar.TwAddVarRW(
    "fade_max_v",
    TW_TYPE_DOUBLE,
    &(sc->igl_params->fade_max_v),
    "step=0.01 max=1 min=0 group=Fade ");
  sc->rebar.TwAddVarRW(
    "fade_min_v",
    TW_TYPE_DOUBLE,
    &(sc->igl_params->fade_min_v),
    "step=0.01 max=1 min=0 group=Fade ");

  TwEnumVal ColorMapTypeEV[NUM_COLOR_MAP] = 
  {
    {COLOR_MAP_DEFAULT,"DEFAULT"},
    {COLOR_MAP_JET,"JET"},
    {COLOR_MAP_EASTER,"EASTER"},
    {COLOR_MAP_WINDING_THEN_EASTER,"WN_EASTER"},
  };
  TwType ColorMapTypeTW = 
    igl::ReTwDefineEnum(
      "ColorMapType", 
      ColorMapTypeEV, 
      NUM_COLOR_MAP);
  sc->rebar.TwAddVarRW(
    "color_map",
    ColorMapTypeTW,
    &(sc->igl_params->color_map),
    "group=View "
    "keydecr='<' keyincr='>' "
    "help='Colormap for data vis.' ");
  sc->rebar.TwAddVarRW(
    "easter_red",
    TW_TYPE_COLOR3F,
    &(sc->igl_params->easter_red),
    "group=View "
    "help='Colormap easter hsv for red.' colormode=hls");
  //sc->rebar.TwAddVarRW(
  //  "easter_s",
  //  TW_TYPE_DOUBLE,
  //  &(sc->igl_params->easter_s),
  //  "group=View "
  //  "help='Colormap easter hsv for red.' colormode=hls");
  //sc->rebar.TwAddVarRW(
  //  "easter_v",
  //  TW_TYPE_DOUBLE,
  //  &(sc->igl_params->easter_v),
  //  "group=View "
  //  "help='Colormap easter hsv for red.' ");

  sc->rebar.TwAddButton(
    "clip_xy",
    clip_xy_plane,
    sc,
    "group='Clip' "
    "key='ALT+x' "
    "help='Set clip rotation to clip xy plane' ");
  sc->rebar.TwAddButton(
    "clip_xz",
    clip_xz_plane,
    sc,
    "group='Clip' "
    "key='ALT+z' "
    "help='Set clip rotation to clip xz plane' ");
  sc->rebar.TwAddButton(
    "clip_yz",
    clip_yz_plane,
    sc,
    "group='Clip' "
    "key='ALT+y' "
    "help='Set clip rotation to clip yz plane' ");
  sc->rebar.TwAddButton(
    "invert_clip",
    invert_clip,
    sc,
    "group='Clip' "
    "key='I' "
    "help='invert clip direction' ");
  sc->rebar.TwAddVarCB(
    "clip_rotation",
    TW_TYPE_QUAT4F,
    set_clip_rotation,
    get_clip_rotation,
    sc,
    "group='Clip' ");
  sc->rebar.TwAddButton(
    "clip_snap",
    snap_rotation_to_canonical_clip_quat,
    sc,
    "group='Clip' "
    "key='ALT+Z' "
    "help='Snap clip rotation to nearest canonical view' ");
  sc->rebar.TwAddVarCB(
    "clip_tra12",
    TW_TYPE_DIR3F,
    set_clip_tra12,
    get_clip_tra12,
    sc,
    "help='Translation component of clip (dont touch)' "
    "group='Clip' readonly=true");
  sc->rebar.TwAddVarCB(
    "clip_C_HIDE",
    TW_TYPE_BOOLCPP,
    set_C_HIDE,
    get_C_HIDE,
    sc,
    "help='Hide clipping plane ' "
    "key='P' "
    "group='Clip' ");


  sc->rebar.TwAddVarRW(
    "mat_amb",
    TW_TYPE_COLOR4F,
    sc->material->amb,
    "group='Material' "
    "label='Ambient' "
    "colormode=hls ");

  //sc->rebar.TwAddVarRW(
  //  "mat_emi",
  //  TW_TYPE_COLOR4F,
  //  sc->material->emi,
  //  "group='Material (double-tap e)' "
  //  "label='Emission' "
  //  "colormode=hls ");

  sc->rebar.TwAddVarRW(
    "mat_dif",
    TW_TYPE_COLOR4F,
    sc->material->dif,
    "group='Material (double-tap e)' "
    "label='Diffuse' "
    "colormode=hls ");

  sc->rebar.TwAddVarRW(
    "mat_spe",
    TW_TYPE_COLOR4F,
    sc->material->spe,
    "group='Material (double-tap e)' "
    "label='Specular' "
    "colormode=hls ");

  sc->rebar.TwAddVarRW(
    "mat_shininess",
    TW_TYPE_FLOAT,
    &(sc->material->shininess),
    "group='Material (double-tap e)' "
    "label='shininess' ");

  sc->rebar.TwAddVarRW(
    "scene_type",
    TW_TYPE_UINT8,
    &(sc->type),
    "group='Scene' "
    "readonly=true "
    "label='type' ");

  sc->rebar.TwAddVarRW(
    "scene_mode",
    TW_TYPE_UINT8,
    &(sc->mode),
    "group='Scene' "
    "readonly=true "
    "label='mode' ");

  sc->rebar.TwAddVarRW(
    "clip",
    TW_TYPE_UINT8,
    &(sc->clip->active),
    "group='Scene' "
    "readonly=true "
    "label='clip->active' ");

  sc->rebar.TwAddVarRW(
    "render_on_C_UPDATE",
    TW_TYPE_BOOLCPP,
    &(sc->igl_params->render_on_C_UPDATE),
    "group='Render' "
    "help='Write a png screenshot each time the clip is moved' ",
    false
    );

  // Try to load previous bar
  if(!sc->rebar.load("rebar.rbr"))
  {
    fprintf(stderr,"Error: could not reload rebar.rbr\n");
  }else
  {
    // Redo clipping to get proper clip
    ubyte old_ca = sc->clip->active;
    sc->clip->active |= C_REDO;
    sc->clip->active |= C_UPDATE;
    sc->clip->active |= C_EDIT;
    updateClip(sc->clip,mesh);
    sc->clip->active = old_ca;

    // Recompile lists to get proper colors
    doLists(sc,mesh);
    doMapLists(sc,mesh,true);
  }
}
#endif


void initGrafix(pScene sc,pMesh mesh) {
  GLfloat  lightamb[4] = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat  lightdif[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat  lightpos[4] = { 0.0, 0.0, 1.0, 0.0 };

  if ( ddebug )  printf("initGrafix\n");
  glEnable(GL_DEPTH_TEST);
  //glDepthFunc(GL_LEQUAL);
  glDepthFunc(GL_LESS);
  glPolygonOffset(1.0, 1.0 / (float)0x10000);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

  glShadeModel(GL_SMOOTH);
  glDisable(GL_NORMALIZE);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_DITHER);
  glDisable(GL_CULL_FACE);
  if ( mesh->typ == 2 ) {
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
  }
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);

  /* lighting */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
#if ( !defined(GL_VERSION_1_1) )
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL,GL_SEPARATE_SPECULAR_COLOR);
#endif

  glLightfv(GL_LIGHT0,GL_DIFFUSE,lightdif);
  glLightfv(GL_LIGHT0,GL_AMBIENT,lightamb);
  glEnable(GL_LIGHTING);
  if ( stereoMode != MONO ) {
    lightpos[2] = -1.0;
    if ( sc->par.sunpos )
      sc->par.sunpos[2] = -fabs(sc->par.sunpos[2]);
  }
  if ( sc->par.sunp )
    glLightfv(GL_LIGHT0,GL_POSITION,sc->par.sunpos);
  else
    glLightfv(GL_LIGHT0,GL_POSITION,lightpos);
  glEnable(GL_LIGHT0);
}


/* new scene */
int createScene(pScene sc,int idmesh) {
  pMesh     mesh;
  char      data[128];

  /* default */
  mesh = cv.mesh[idmesh];
  if ( !quiet ) fprintf(stdout,"   Computing 3D scene\n");

  /* set default mode */
  sc->idmesh = idmesh;
  sc->par.xi = sc->par.yi = 10;
  if ( option == SCHNAUZER ) {
    sc->par.xs = schw;
    sc->par.ys = schh;
  }
  else {
    if ( sc->par.xs == 0 )  sc->par.xs = 600;
    if ( sc->par.ys == 0 )  sc->par.ys = 600;
  }
  if ( !sc->mode )  sc->mode = HIDDEN;

  sc->item   = 0;
  sc->shrink = 1.0;
  sc->slave  = sc->master = -1;
  sc->picked = 0;
  if ( mesh->nvn == 0 )  sc->type ^= S_FLAT;
  if ( mesh->ne == 0 )   sc->item |= S_GEOM;

  /* compute scene depth */
  sc->dmax = sc->dmin= mesh->xmax - mesh->xmin;
  sc->dmax = MEDIT_MAX(sc->dmax,mesh->ymax - mesh->ymin);
  sc->dmin = MEDIT_MIN(sc->dmin,mesh->ymax - mesh->ymin);
  if ( mesh->dim == 3 ) {
    sc->dmax = MEDIT_MAX(sc->dmax,mesh->zmax - mesh->zmin);
    sc->dmin = MEDIT_MIN(sc->dmin,mesh->zmax - mesh->zmin);
  }
  sc->dmax = fabs(sc->dmax);
  sc->dmin = fabs(sc->dmin);
  if ( !sc->par.sunp ) {
    sc->par.sunpos[0] *= 2.0*sc->dmax;
    sc->par.sunpos[1] *= 2.0*sc->dmax;
    sc->par.sunpos[2] *= 2.0*sc->dmax;
  }
  sc->par.sunpos[3] = 1.0;

  /* create window */
  glutInitWindowSize(sc->par.xs,sc->par.ys);
  glutInitWindowPosition(sc->par.pxs,sc->par.pys);
  sc->idwin = glutCreateWindow("");
  assert(sc->idwin != 0);
  if ( fullscreen ) {
    glutFullScreen();
    sc->par.xs = glutGet(GLUT_SCREEN_WIDTH);
    sc->par.ys = glutGet(GLUT_SCREEN_HEIGHT);
  }

  /* set window name */
  sprintf(data,"Medit - [%s] #%d",mesh->name,sc->idwin);
  glutSetWindowTitle(data);
  glutSetIconTitle(data);
 
 /* required! to change background color */
  glClearColor(sc->par.back[0],sc->par.back[1],
               sc->par.back[2],0);

  /* init perspective */
  sc->persp  = initPersp(0,sc->dmax);
  sc->camera = (pCamera)initCamera(sc,Y_AXIS);
  if ( mesh->typ == 2 ) {
    sc->persp->pmode = CAMERA;
    sc->persp->depth *= 0.5;
  }

  /* create default view */
  sc->view = (pTransform)createTransform();
  if ( !sc->view )  return(0);
  sc->type |= S_RESET + S_DECO;
  sc->clip  = (pClip)createClip(sc,mesh);
  if ( !sc->clip )  return(0);
  sc->cube  = (pCube)createCube(sc,mesh);
  if ( !sc->cube )  return(0);

  /* create menus */
  if ( !createMenus(sc,mesh) )  return(0);

  /* assign callbacks */
  if ( sc->type & S_SCISSOR)
    glutDisplayFunc(scissorScene);
  else if ( option == SCHNAUZER )
    glutDisplayFunc(redrawSchnauzer);
  else if ( sc->persp->pmode == CAMERA ) {
    glutMouseFunc(mouseCamera);
    glutMotionFunc(motionCamera);
    glutDisplayFunc(redrawScene);
  }
  else {
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutDisplayFunc(redrawScene);
  }
#ifdef IGL
  glutPassiveMotionFunc(passive_motion); // same as MouseMotion
#endif
  glutReshapeFunc(reshapeScene);
  glutKeyboardFunc(keyScene);
  glutSpecialFunc(special);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

#ifndef IGL
  /* create display lists by geom type */
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  doLists(sc,mesh);
  sc->glist = geomList(sc,mesh);
#endif
  sc->type |= S_FOLLOW;

  /* local stack */
  if ( !pilmat ) {
    pilmat = (int*)calloc(sc->par.nbmat+2,sizeof(int));
    if ( !pilmat ) return(0);
  }

  /* color list */
  setupPalette(sc,mesh);
  sc->stream = NULL;

  initGrafix(sc,mesh);
#ifdef IGL
  initIGL(sc,mesh);
  initAntTweakBar(sc,mesh);
#endif

#ifdef IGL
  /* create display lists by geom type */
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  doLists(sc,mesh);
#endif

  return(1);
}
