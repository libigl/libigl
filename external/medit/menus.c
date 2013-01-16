#include "medit.h"
#include "extern.h"
#include "sproto.h"

GLfloat    altcoef=0.0;
static int nbimg;
extern int refitem,refmat,reftype,imstep,imreverse;
extern int *pilmat,ipilmat;

#define MAX_LST   128


/* rebuild display lists */
void doLists(pScene sc,pMesh mesh) {
  int     k;

  /*default */
  if ( ddebug ) printf("build display lists\n");

  if ( !morphing )  glutSetCursor(GLUT_CURSOR_WAIT);
  for (k=0; k<MAX_LIST; k++) {
    if ( sc->dlist[k] )  glDeleteLists(sc->dlist[k],1);
    if ( sc->clist[k] )  glDeleteLists(sc->clist[k],1);
    sc->dlist[k] = sc->clist[k] = (GLuint)0;
  }

  /* poly lists */
  sc->dlist[LTria] = listTria(sc,mesh);
  sc->dlist[LQuad] = listQuad(sc,mesh);
  if ( !mesh->nt && !mesh->nq ) {
    sc->dlist[LTets] = listTetra(sc,mesh,0);
    sc->dlist[LHexa] = listHexa(sc,mesh,0);
  }
  if ( sc->clip->active & C_VOL )  sc->clip->active |= C_REDO;
  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors();

#ifdef IGL
  /* create display lists by geom type */
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  sc->glist = geomList(sc,mesh);
#endif
}


/* build metric list */
void doMapLists(pScene sc,pMesh mesh,int reset) {
  int   k;

  /*default */
  if ( !mesh->sol )  return;
  if ( ddebug ) printf("build map lists\n");
  glutSetCursor(GLUT_CURSOR_WAIT);

  /* delete old lists */
  if ( reset ) {
    for (k=0; k<MAX_LIST; k++) {
      if ( sc->mlist[k] )  glDeleteLists(sc->mlist[k],1);
      sc->mlist[k] = (GLuint)0;
    }
  }

  /* create map list */
  if ( sc->mode & S_ALTITUDE ) {
    sc->mlist[LTria] = alt2dList(sc,mesh,LTria,sc->shrink,altcoef);
    sc->mlist[LQuad] = alt2dList(sc,mesh,LQuad,sc->shrink,altcoef);
  }
  else if ( sc->mode & S_MAP ) {
    if ( mesh->nt && !sc->mlist[LTria] )
      sc->mlist[LTria] = listTriaMap(sc,mesh);
    if ( mesh->nq && !sc->mlist[LQuad] )
      sc->mlist[LQuad] = listQuadMap(sc,mesh);
    if ( !mesh->nt && mesh->ntet && !sc->mlist[LTets] )
      sc->mlist[LTets] = listTetraMap(sc,mesh,0);
    if ( !mesh->nq && mesh->nhex && !sc->mlist[LHexa] )
      sc->mlist[LHexa] = listHexaMap(sc,mesh,0);
  }

  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors();
}


/* rebuild iso-values lists */
void doIsoLists(pScene sc,pMesh mesh,int reset) {
  pPoint  ppt;
  int     k,kk,ret;
  
  /*default */
  if ( ddebug ) printf("build iso lists\n");
  glutSetCursor(GLUT_CURSOR_WAIT);

  /* delete old lists */
  if ( reset ) {
    for (k=0; k<MAX_LIST; k++) {
      if ( sc->ilist[k] )  glDeleteLists(sc->ilist[k],1);
      if ( sc->vlist[k] )  glDeleteLists(sc->vlist[k],1);
      if ( sc->cplist )    glDeleteLists(sc->cplist,1);
      sc->ilist[k] = sc->vlist[k] = (GLuint)0;
      sc->cplist   = (GLuint)0;
    }
    if ( sc->stream && sc->stream->nbstl ) {
      for (kk=0; kk<sc->stream->nbstl; kk++) {
        if ( sc->slist[kk] )  glDeleteLists(sc->slist[kk],1);
        sc->slist[kk] = (GLuint)0;
      }
      sc->slist = (GLuint)0;
      if ( reset < 2 )   sc->stream->nbstl = 0;
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        ppt->flag = 0;
      }
    }
  }

  /* iso-lines */
  if ( sc->isotyp & S_ISOLINE ) {
    if ( mesh->nt && !sc->ilist[LTria] )
      sc->ilist[LTria] = listTriaIso(sc,mesh);
    if ( mesh->nq && !sc->ilist[LQuad] )
      sc->ilist[LQuad] = listQuadIso(sc,mesh);
  }
  /* iso-surfaces */
  if ( sc->isotyp & S_ISOSURF ) {
    if ( mesh->ntet && !sc->ilist[LTets] )
      sc->ilist[LTets] = listTetraIso(sc,mesh);
  }

  /* vector */
  if ( sc->isotyp & S_VECTOR ) {
    if ( mesh->dim == 2 ) {
      if ( mesh->nt && !sc->vlist[LTria] )
        sc->vlist[LTria] = listTria2dVector(mesh);
      if ( mesh->nq && !sc->vlist[LQuad] )
        sc->vlist[LQuad] = listQuad2dVector(mesh);
    }
    else {
      if ( mesh->ntet+mesh->nhex == 0 ) {
        sc->vlist[LTria] = listTria3dVector(mesh);
      }
      else {
        if ( mesh->ntet && !sc->vlist[LTets] )
          sc->vlist[LTets] = listClipTetraVector(mesh);
        else if ( mesh->nhex && !sc->vlist[LHexa] )
          sc->vlist[LHexa] = listClipHexaVector(mesh);
      }
    }
  }

  /* streamlines */
  if ( sc->isotyp & S_STREAML ) {
    if ( !mesh->adja || !sc->slist ) {
      sc->stream = createStream(sc,mesh);
      if ( !sc->stream )
        sc->isotyp &= ~S_STREAML;
    }
    if ( sc->stream ) {
      if ( reftype == LPoint )
        ret = streamRefPoint(sc,mesh);
      else if ( mesh->dim == 3 ) {
        if ( reftype == LTria )
          ret = streamRefTria(sc,mesh);
        else if ( reftype == LQuad )
          ret = streamRefQuad(sc,mesh);
        if ( sc->picklist )  glDeleteLists(sc->picklist,1);
        sc->picklist = 0;
      }
    }
  }

  /* critical points */
  if ( sc->isotyp & S_CRITP ) {
    if ( mesh->dim == 2 )
      if ( mesh->nt && !sc->cplist )
        sc->cplist = listCritPoint(sc,mesh);
  }
  
  glutSetCursor(GLUT_CURSOR_INHERIT);
  checkErrors();
}


void resetLists(pScene sc,pMesh mesh) {
  int    kk;

  for (kk=0; kk<MAX_LIST; kk++) {
    if ( sc->dlist[kk] )  glDeleteLists(sc->dlist[kk],1);
    if ( sc->mlist[kk] )  glDeleteLists(sc->mlist[kk],1);
    if ( sc->ilist[kk] )  glDeleteLists(sc->ilist[kk],1);
    if ( sc->clist[kk] )  glDeleteLists(sc->clist[kk],1);
    if ( sc->cmlist[kk] ) glDeleteLists(sc->cmlist[kk],1);
    if ( sc->vlist[kk] )  glDeleteLists(sc->vlist[kk],1);
    if ( sc->cplist )     glDeleteLists(sc->cplist,1);
    sc->dlist[kk] = sc->clist[kk] = sc->mlist[kk] = (GLuint)0;
    sc->ilist[kk] = sc->vlist[kk] = sc->cplist    = (GLuint)0;
  }
  if ( sc->glist )  glDeleteLists(sc->glist,1);
  if ( sc->nlist )  glDeleteLists(sc->nlist,1);
  sc->glist = sc->nlist = (GLuint)0;
}


void keyFile(unsigned char key,int x,int y) {
  pScene       sc;
  pMesh        mesh;
  pTransform   view;
  char        *ptr,data[128];
  ubyte        post  = FALSE,clipon;
  static int   nfree = 0;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  view = sc->view;

  switch(key) {
  case 'L':  /* load prefs */
    if ( parsop(sc,mesh) ) {
      refitem = reftype = 0;
      setupPalette(sc,mesh);
      initGrafix(sc,mesh);
      doLists(sc,mesh);
      if ( sc->mode & S_MAP )  doMapLists(sc,mesh,1);
      if ( sc->isotyp )        doIsoLists(sc,mesh,1);
      glClearColor(sc->par.back[0],sc->par.back[1],sc->par.back[2],sc->par.back[3]);
      post = GL_TRUE;
    }
    break;
  case 'W':  /* save prefs */
    saveMeditFile(mesh->name,sc);
    break;

  case 'R':  /* reload mesh */
    refitem = reftype = 0;
    if ( !meshUpdate(sc,mesh) ) exit(1);
    meshRef(sc,mesh);
    matSort(sc);
    setupPalette(sc,mesh);
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,1);
    if ( sc->isotyp )       doIsoLists(sc,mesh,1);
    glClearColor(sc->par.back[0],sc->par.back[1],sc->par.back[2],sc->par.back[3]);
    post = GL_TRUE;
    break;
  
  case 'S':  /* save mesh */
    strcpy(data,mesh->name);
    ptr = (char*)strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';
    sprintf(data,"%s.d.mesh",data);
    clipon = sc->clip->active & C_ON;
    if ( clipon )
      clipVertices(mesh,sc,sc->clip);
    saveMesh(sc,mesh,data,clipon);
    break;
    
  case 'B':
  case 'C':
  case 'G':
  case 'H':  /* hardcopy */
  case 'T':
    strcpy(data,mesh->name);
    ptr = (char*)strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';
    ptr = (char*)strstr(data,".gis");
    if ( ptr ) *ptr = '\0';

    if ( option != SEQUENCE ) {
      if ( key == 'H' )
        nfree = filnum(data,nfree,"ppm");
      else
        nfree = filnum(data,nfree,"ps");
      if ( nfree == -1 )  break;
      sprintf(data,"%s.%.3d",data,nfree);
    }
    else
      sprintf(data,"%s.ppm",data);
    if (!saveimg )  glutSetCursor(GLUT_CURSOR_WAIT);
    if ( key == 'B' || key == 'C' || key == 'G' )
      imgTiling(sc,data,key);
    else
      imgHard(sc,data,key);
    if ( animate && !(sc->persp->pmode & CAMERA) && ++nbimg > 179 ) {
      view->angle = 0.0;
      nbimg       = 0;
      keyAnim('A',0,0);
    }
    if (!saveimg ) glutSetCursor(GLUT_CURSOR_INHERIT);
    break;
  case 's':  /* softcopy */
    sftcpy(sc,mesh);
    post = FALSE;
    break;
  }
  if ( post ) glutPostRedisplay();
}

void menuFile(int item) {
  keyFile((unsigned char)item,0,0);
}


void keyItem(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  pCube   cube;
  ubyte   post = TRUE;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  cube = sc->cube;

  switch(key) {
  case 'A':  /* draw axis */
    sc->item ^= S_AXIS;
    break;
  case 'B':  /* draw bounding box */
    sc->item ^= S_BOX;
    break;
  case 'D':  /*draw selection cube */
    if ( cube->active & C_ON ) {
      cube->active &= ~(C_ON+C_EDIT);
      dumpCube(sc,mesh,cube);
    }
    else {
      cube->active |= C_ON;
    }
    break;
  case 'G':  /* draw grid */
    sc->item ^= S_GRID;
    break;
  case 'j':  /* toggle info */
    sc->type ^= S_DECO;
    break;
  case 'P':  /* toggle point nums */
    sc->item ^= S_NUMP;
    break;
  case 'F':  /* toggle face nums */
    sc->item ^= S_NUMF;
    break;
  case 'g':  /* const. items */
    if ( !sc->glist )  post = FALSE;
    sc->item ^= S_GEOM;
#ifdef IGL
    if(sc->item | S_GEOM)
    {
      doLists(sc,mesh);
    }
#endif
    break;
  case 'N':
    if ( mesh->nvn )
      sc->type ^= S_NORMAL;
    else
      post = FALSE;
    break;
  case 'O':
    if ( mesh->nvn ) {
      sc->type ^= S_OPPOS;
      if ( sc->nlist ) {
	glDeleteLists(sc->nlist,1);
	sc->nlist = 0;
      }
      else 
	post = FALSE;
    }
    else
      post = FALSE;
    break;
  }
  if ( post )  glutPostRedisplay();
}

void menuItem(int item) {
  keyItem((unsigned char)item,0,0);
}


void keyAnim(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  pClip   clip;
  ubyte   post = TRUE;
  char   *ptr,base[256];
  static int depart = -1;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;
  if ( depart == -1 )
    depart = animdep;

  switch(key) {
  case 'A':  /* switch animation */
    if ( animate ) {
      sc->view->manim     = GL_FALSE;
      clip->cliptr->manim = GL_FALSE;
      sc->view->angle = 0.0f;
      glutIdleFunc(NULL);
    }
    else {
      sc->view->manim     = GL_TRUE;
      clip->cliptr->manim = GL_TRUE;
      if ( sc->persp->pmode == CAMERA )
        glutIdleFunc(glutIdle);
    }
    animate = 1-animate;
    post = FALSE;
    break;
  case 'I':  /* save image */
    saveimg = 1 - saveimg;
    nbimg   = 0;
    post    = GL_FALSE;
    break;

  case 'M':  /* morphing */
    if ( morphing )
      glutIdleFunc(NULL);
    else
      glutIdleFunc(glutIdle);
    morphing = 1-morphing;
    post = FALSE;
    break;
  case 'R':  /* next morph */
    imreverse = 1-imreverse;
    break;

  case 'S':  /* start animation */
    if ( ddebug ) fprintf(stdout,"debut sequence %d a %d\n",animdep,animfin);
    glutSetWindow(sc->idwin);
    if ( option == SEQUENCE )
      playAnim(sc,mesh,animdep,animfin);
    if ( !saveimg )  post = FALSE;
    break;

  case 'f': /* first mesh */
    /* get basename */
    ptr = (char *)strrchr(mesh->name,'.');
    if ( ptr )  *ptr = '\0';
    strcpy(base,mesh->name);
    resetLists(sc,mesh);
    if ( !loadNextMesh(mesh,animdep,0) )  break;
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,0);
    if ( sc->isotyp )       doIsoLists(sc,mesh,0);
    strcpy(mesh->name,base);
    break;

  case 'l':
    /* get basename */
    ptr = (char *)strrchr(mesh->name,'.');
    if ( ptr )  *ptr = '\0';
    strcpy(base,mesh->name);
    resetLists(sc,mesh);
    if ( !loadNextMesh(mesh,animfin,0) )  break;
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,0);
    if ( sc->isotyp )       doIsoLists(sc,mesh,0);
    strcpy(mesh->name,base);
    break;

  case 'n':
    /* get basename */
    if ( ++depart > animfin )  depart = animdep;
    ptr = (char *)strrchr(mesh->name,'.');
    if ( ptr )  *ptr = '\0';
    strcpy(base,mesh->name);
    resetLists(sc,mesh);
    if ( !loadNextMesh(mesh,depart,0) )  break;
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,0);
    if ( sc->isotyp )       doIsoLists(sc,mesh,0);
    strcpy(mesh->name,base);
    break;

  case 'p':
    /* get basename */
    if ( --depart < animdep )  depart = animfin;
    ptr = (char *)strrchr(mesh->name,'.');
    if ( ptr )  *ptr = '\0';
    strcpy(base,mesh->name);
    resetLists(sc,mesh);
    if ( !loadNextMesh(mesh,depart,0) )  break;
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,0);
    if ( sc->isotyp )       doIsoLists(sc,mesh,0);
    strcpy(mesh->name,base);
    break;
  }
  
  if ( post == TRUE )   glutPostRedisplay();
}

void menuAnim(int item) {
  keyAnim((unsigned char)item,0,0);
}

void keyTrajet(unsigned char key,int x,int y) {
  pScene   sc;
  pMesh    mesh;

  sc = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  switch(key) {
  case 'C': /* add control point */
    pathAdd(sc,x,y);
    break;
  case 'S': /* show all points */
    sc->type ^= S_PATH;
    /*if ( sc->path.tlist )  glDeleteLists(sc->path.tlist,1);*/
    if ( sc->type & S_PATH && !sc->path.tlist )
      sc->path.tlist = pathList(sc);
    break;
  case 'F': /* follow path */
    break;
  case 'L': /* load path */
    pathLoad(mesh->name,sc);
    break;
  case 'W': /* save path */
    pathSave(mesh->name,sc);
    break;
  }
  glutPostRedisplay();
}

void menuTrajet(int item) {
  keyTrajet((unsigned char)item,0,0);
}

void keyMode(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  ubyte   post = TRUE,dolist = FALSE,material;

  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  material = sc->mode & S_MATERIAL;

  switch(key) {
  case 'D':  /* depth lines */
    if ( sc->mode == DEPTH ) break;
    sc->mode = DEPTH;
    break;
  case 'S': /* toggle smooth shading */
    if ( sc->mode == FILL )
      break;
    sc->mode = FILL;
    if ( material ) sc->mode |= S_MATERIAL;
    break;
  case 'P': /* toggle smooth shaded polys */
    if ( sc->mode == SHADED )
      break;
    sc->mode = SHADED;
    if ( material ) sc->mode |= S_MATERIAL;
    break;
  case 'H': /* hidden lines */
    if ( sc->mode == HIDDEN ) break;
    sc->mode = HIDDEN;
    if ( material ) sc->mode |= S_MATERIAL;
    break;
  case 'W':  /* wireframe */
    if ( sc->mode == WIRE ) break;
    sc->mode = WIRE;
    if ( material ) sc->mode |= S_MATERIAL;
    break;
  case 'n': /* toggle normals */
    if ( mesh->nvn == 0 ) {
      post = FALSE;
      break;
    }
    sc->type ^= S_FLAT;
    dolist = TRUE;
    break;
  default:
    post = FALSE;
    break;
  }
  if ( dolist == TRUE ) doLists(sc,mesh);
  if ( post == TRUE )   glutPostRedisplay();
}

void menuMode(int item) {
  keyMode((unsigned char)item,0,0);
}


void menuScene(int item) {
  keyScene((unsigned char)item,0,0);
}

void keyView(unsigned char key,int x,int y) {
  pScene  sc,sc1;
  pMesh   mesh;
  float   dmax;
  ubyte   post = FALSE;

  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  /*glutSetMenu(vmenu);*/

  switch(key) {
  case 'R':
    sc->type |= S_RESET;
    dmax = mesh->xmax - mesh->xmin;
    dmax = MEDIT_MAX(dmax,mesh->ymax - mesh->ymin);
    dmax = MEDIT_MAX(dmax,mesh->zmax - mesh->zmin);
    sc->cx = sc->cy = sc->cz = 0.0f;
    sc->dmax = fabs(dmax);
    if ( sc->persp->pmode == PERSPECTIVE ) {
      resetTransform(sc->view);
      initPersp(sc->persp,sc->dmax);
    }
    else if ( sc->persp->pmode == CAMERA ) {
      initPersp(sc->persp,sc->dmax);
      initCamera(sc,sc->camera->vecup); 
      sc->persp->pmode = CAMERA;
    }
    reshapeScene(sc->par.xs,sc->par.ys);
    post = TRUE;
    break;
  case 'C':
    copyView(sc->view,sc->camera,sc->persp);
    copyClip(sc->clip);
    break;
  case 'L':  /* link view */
    if ( !linkView(sc) )  return;
    reshapeScene(sc->par.xs,sc->par.ys);
    post = TRUE;
    break;
  case 'P':
    if ( pasteView(sc->view,sc->camera,sc->persp) && 
         pasteClip(sc->clip) ) {
      reshapeScene(sc->par.xs,sc->par.ys);
      post = TRUE;
    }
    break;
  case 'U':
    unlinkView(sc);
    break;

  case 'D':  /* duplicate view */
    if ( cv.nbs == MAX_SCENE )  break;
    if ( !cv.scene[++cv.nbs] ) {
      cv.scene[cv.nbs] = (pScene)M_calloc(1,sizeof(Scene),"menus.scene");
      if ( !cv.scene[cv.nbs] )  break;
      sc1 = cv.scene[cv.nbs];
      sc1->material = (pMaterial)calloc(2+sc->par.nbmat,sizeof(Material));
      assert(sc1->material);
    }
    sc1 = cv.scene[cv.nbs];
    memcpy(sc1,sc,sizeof(Scene));
    memcpy(sc1->material,sc->material,sc->par.nbmat*sizeof(Material));
    memcpy(&sc1->par,&sc->par,sizeof(Param));
    if ( !createScene(sc1,sc1->idmesh) ) {
      fprintf(stdout,"  ## Unable to create\n");
      return;
    }
    copyView(sc->view,sc->camera,sc->persp);
    copyClip(sc->clip);
    pasteView(sc1->view,sc1->camera,sc1->persp);
    pasteClip(sc1->clip);
    break;
  }
  if ( post )  glutPostRedisplay();
}

void menuView(int item) {
  keyView((unsigned char)item,0,0);
}

void keyColor(unsigned char key,int x,int y) {
  pScene     sc;
  pMesh      mesh;
  pMaterial  pm;
  int        i,k;
  ubyte      post = TRUE,dolist = FALSE;

  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];

  switch(key) {
  case 'b':  /* reverse backcolor */
    sc->par.back[0] = 1.0f - sc->par.back[0];
    sc->par.back[1] = 1.0f - sc->par.back[1];
    sc->par.back[2] = 1.0f - sc->par.back[2];
    glClearColor(sc->par.back[0],sc->par.back[1],sc->par.back[2],sc->par.back[3]);
    if ( !sc->par.linc ) {
      sc->par.line[0] = 1.0f - sc->par.line[0];
      sc->par.line[1] = 1.0f - sc->par.line[1];
      sc->par.line[2] = 1.0f - sc->par.line[2];
    }
    break;
  case 'e':  /* toggle matcolors */
    sc->mode ^= S_MATERIAL;
    dolist = TRUE;
    break;
  case 'E':  /* edit matcolors */
    post   = FALSE; 
    matEdit(sc);
    break;
  case 'r':
      if ( refmat<0 || ipilmat == sc->par.nbmat )  break;
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
  case 'R':  /* reset materials */
    for (k=0; k<sc->par.nbmat; k++) {
      pm = &sc->material[k];
      for (i=LTria; i<=LHexa; i++)
        pm->depmat[i] = abs(pm->depmat[i]);
    }
    dolist = TRUE;
    break;
  default:
    post = FALSE;
  }
  if ( dolist ) doLists(sc,mesh);
  if ( post )   glutPostRedisplay();
}

void menuColor(int item) {
  keyColor((unsigned char)item,0,0);
}


void keyClip(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  pClip   clip;
  ubyte   post = TRUE;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  clip = sc->clip;

  switch(key) {
  case 'C': /* toggle clipping */
    if ( clip->active & C_ON )
      clip->active &= ~(C_ON+C_EDIT);
    else {
      clip->active |= C_ON;
      if ( mesh->ntet+mesh->nhex )  clip->active |= C_VOL;
    }
    break;
  case 'E': /* edit clip plane */
    if ( !(clip->active & C_ON) )  return;
    if ( clip->active & C_FREEZE ) clip->active ^= C_FREEZE;
    clip->active ^= C_EDIT;
    break;
  case 'F':  /* freeze clip */
    if ( !(clip->active & C_ON) )  return;
    if ( clip->active & C_EDIT )   clip->active ^= C_EDIT;
    clip->active ^= C_FREEZE;
    break;
  case 'H':  /* toggle draw plane */
    clip->active ^= C_HIDE;
    break;
  case 'I': /* inverse orientation */
    invertClip(sc,clip);
    clip->active |= C_REDO;
    break;
  case 'K':  /* toggle capping */
    clip->active ^= C_CAP;
    clip->active |= C_REDO;
    break;
  case 'R':  /* reset clip */
    if ( !(clip->active & C_ON) ) break;
    resetClip(sc,clip,mesh);
    break;
  case 'Z':  /* toggle volclip */
    clip->active ^= C_VOL;
    if ( clip->active & C_VOL )
      clip->active |= C_REDO;
    break;
  }
  if ( post ) glutPostRedisplay();
}

void menuClip(int item) {
  keyClip((unsigned char)item,0,0);
}


void keyCube(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  pCube   cube;
  ubyte   post = TRUE;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];
  cube = sc->cube;

  switch(key) {
  case 'C': /* toggle clipping */
    if ( cube->active & C_ON )
      cube->active &= ~(C_ON+C_EDIT);
    else
      cube->active |= C_ON;
    break;
  case 'E': /* edit clip plane */
    if ( !(cube->active & C_ON) )  return;
    if ( cube->active & C_FREEZE ) cube->active ^= C_FREEZE;
    cube->active ^= C_EDIT;
    break;
  case 'F':  /* freeze cube */
    if ( !(cube->active & C_ON) )  return;
    if ( cube->active & C_EDIT )   cube->active ^= C_EDIT;
    cube->active ^= C_FREEZE;
    break;
  case 'R':  /* reset cube */
    if ( !(cube->active & C_ON) ) break;
    resetCube(sc,cube,mesh);
    break;
  }
  if ( post ) glutPostRedisplay();
}

void keyFeature(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  ubyte   post=TRUE,dolist=TRUE;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];

  switch(key) {
  case 'S':  /* shrink mode */
  case 'V':  /* volumic shrink */
    if ( sc->shrink < 0.99 )
      sc->shrink = 1.0f;
    else
      sc->shrink = 0.95;
    dolist = TRUE;
    break;
  case 'I':  /* increase shrink value */
    if ( sc->shrink > 0.99 ) break;
    sc->shrink -= .05;
    if ( sc->shrink < 0.1 ) sc->shrink = 0.1;
    dolist = TRUE;
    break;
  case 'i':  /* decrease shrink value */
    if ( sc->shrink > 0.99 ) break;
    sc->shrink += .05;
    if ( sc->shrink > 0.95 ) sc->shrink = 0.95;
    dolist = TRUE;
    break;
  case 's':  /* scissor mode */
    if ( mesh->dim == 2 )  break;
    sc->type ^= S_SCISSOR;
    if ( sc->type & S_SCISSOR )
      glutDisplayFunc(scissorScene);
    else {
      glutDisplayFunc(redrawScene);
      reshapeScene(sc->par.xs,sc->par.ys);
    }
    break;
  }
  if ( dolist == TRUE ) {
    doLists(sc,mesh);
    if ( sc->mode & S_MAP )  doMapLists(sc,mesh,1);
    /*if ( sc->isotyp )        doIsoLists(sc,mesh,1);*/
  }
  if ( post == TRUE )   glutPostRedisplay();
}

void menuFeature(int item) {
  keyFeature((unsigned char)item,0,0);
}


void menuImage(int item) {
  imgtype = item;
}

void keyMetric(unsigned char key,int x,int y) {
  pScene  sc;
  pMesh   mesh;
  pPoint  ppt;
  float   maxd;
  int     k,kk;
  ubyte   post=TRUE;

  /* default */
  sc   = cv.scene[currentScene()];
  mesh = cv.mesh[sc->idmesh];

  switch(key) {
  case 'c': /* critical points */
    if ( !mesh->nbb )  return;
    sc->isotyp ^= S_CRITP;
    doIsoLists(sc,mesh,0);
    break;
  case 'f':  /* flush streamlines */
    if ( sc->stream->nbstl ) {
      for (kk=0; kk<sc->stream->nbstl; kk++) {
        if ( sc->slist[kk] )  glDeleteLists(sc->slist[kk],1);
        sc->slist[kk] = (GLuint)0;
      }
      sc->stream->nbstl = 0;
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        ppt->flag = 0;
      }
    }
    break;
  case 'l': /* iso-lines */
    if ( !mesh->nbb )  return;
    sc->isotyp ^= S_ISOLINE;
    doIsoLists(sc,mesh,0);
    break;
  case 's': /* iso-surfaces */
    if ( !mesh->nbb ) return;
    sc->isotyp ^= S_ISOSURF;
    doIsoLists(sc,mesh,0);
    break;
  case 'u': /* field lines animation */
    if ( !mesh->nbb || mesh->nfield != mesh->dim )  return;
    if ( refitem ) {
      sc->isotyp |= S_PARTICLE;
      createParticle(sc,mesh);
    }
    glutIdleFunc(streamIdle);
    break;
  case 'd': /* displacement */
    if ( !mesh->nbb || mesh->nfield != mesh->dim )  return;
    sc->mode ^= S_DISPL;
	meshCoord(mesh,(sc->mode & S_DISPL) ? 1 : 0);
	meshBox(mesh,1);
    doLists(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,1);
    if ( sc->isotyp )       doIsoLists(sc,mesh,1);
    break;

  case 'v': /* streamlines */
    if ( !mesh->nbb || mesh->nfield != mesh->dim )  return;
    if ( refitem ) {
      sc->isotyp |= S_STREAML;
      doIsoLists(sc,mesh,0);
    }
    else {
      if ( !streamIsoPoint(sc,mesh) ) {
        post = FALSE;
        sc->isotyp &= ~S_STREAML;
      }
    }
    break;
  case 'w': /*vector/tensor */
    if ( mesh->nfield != mesh->dim ) {
      if ( mesh->dim==2 && mesh->nfield == 3 ) {
        if ( sc->picklist ) {
          glDeleteLists(sc->picklist,1);
          sc->picklist = 0;
          break;
        }
        else
          sc->picklist = drawAllEllipse(sc,mesh);
          break;
      }
      return;
    }
    sc->isotyp ^= S_VECTOR;
    if ( mesh->dim == 3 )
      if ( mesh->ntet+mesh->nhex && !(sc->clip->active & C_ON) ) 
        return;
    doIsoLists(sc,mesh,0);
    break;

  case 'k': /* toggle elevation */
    if ( mesh->dim != 2 )  return;
    sc->mode ^= S_ALTITUDE;
    if ( altcoef == 0.0 ) {
      maxd = MEDIT_MAX(mesh->xmax-mesh->xmin,mesh->ymax-mesh->ymin);
      altcoef = 0.3*maxd / mesh->bbmax;
    }
    if ( !(sc->mode & S_ALTITUDE) ) 
      sc->type |= S_RESET;
    doMapLists(sc,mesh,1);
    break;
  case 'K': /* elevation coeff */
    fprintf(stdout,"elevation coeff (%.2f): ",altcoef);
    fflush(stdout);
    fflush(stdin); fscanf(stdin,"%f",&altcoef);
    if ( altcoef == 0.0 ) sc->mode |= ~S_ALTITUDE;
    sc->type |= S_RESET;
    doMapLists(sc,mesh,1);
    doIsoLists(sc,mesh,1);
    break;

  case 'm': /* display metric */
    if ( !mesh->nbb ) return;
    sc->mode ^= S_MAP;
    doMapLists(sc,mesh,1);
    if ( sc->mode & S_MAP ) {
      if ( sc->clip->active & C_ON ) sc->clip->active |= C_REDO;
      if ( !(sc->item & S_PALETTE) )
        sc->item ^= S_PALETTE;
    }
    else if ( sc->item & S_PALETTE )
      sc->item ^= S_PALETTE;
    break;
  case 'p': /* toggle palette */
    if ( !mesh->nbb ) return;
    sc->item ^= S_PALETTE;
    break;

  default:
    post = FALSE;
    break;
  }
  if ( post )  glutPostRedisplay();
}

void menuMetric(int item) {
  keyMetric((unsigned char)item,0,0);
}


int createMenus(pScene sc,pMesh mesh) {
  int   menu,amenu,fmenu,femenu,vmenu,mmenu,smenu;
  int   clmenu,cmenu,vwmenu,trmenu;

  /* default */
  if ( ddebug )  printf("create menus\n");
  smenu = 0;

  /* File management menu */
  fmenu = glutCreateMenu(menuFile);
  glutAddMenuEntry("[L] Load prefs",'L');
  glutAddMenuEntry("[W] Save prefs",'W');
  glutAddMenuEntry("    Update mesh",'R');
  glutAddMenuEntry("    Save mesh",'S');
  glutAddMenuEntry("[H] Hardcopy PPM",'H');
  glutAddMenuEntry("    Hardcopy EPS (Color)",'C');
  glutAddMenuEntry("    Hardcopy EPS (Grey)",'G');
  glutAddMenuEntry("    Hardcopy EPS (B/W)",'B');
  glutAddMenuEntry("    Softcopy EPS",'s');

  /* rendering mode selector */
  mmenu = glutCreateMenu(menuMode);
  glutAddMenuEntry("Wireframe",'W');
  glutAddMenuEntry("Depth  lines",'D');
  glutAddMenuEntry("Hidden lines",'H');
  glutAddMenuEntry("Shading",'S');
  glutAddMenuEntry("Shading+lines",'P');
  if ( mesh->nvn > 0 )
    glutAddMenuEntry("[n] Toggle Normals",'n');

  /* color & material menu */
  cmenu = glutCreateMenu(menuColor);
  glutAddMenuEntry("[b] Toggle backcolor",'b');
  glutAddMenuEntry("[e] Toggle matcolors",'e');
  glutAddMenuEntry("[E] Edit   matcolors",'E');
  glutAddMenuEntry("[r] Hide   material",'r');
  glutAddMenuEntry("[R] Reset  materials",'R');

  /* metric */
  if ( mesh->nbb > 0 ) {
    smenu = glutCreateMenu(menuMetric);
    glutAddMenuEntry("[m] Toggle metric",'m');
    glutAddMenuEntry("[p] Toggle palette",'p');
    if ( mesh->typage == 2 )
      glutAddMenuEntry("[o] Toggle iso-lines",'l');
    if ( mesh->ntet+mesh->nhex > 0 && mesh->nfield == 1 )
      glutAddMenuEntry("    Toggle iso-surfaces",'s');
    if ( mesh->nfield == mesh->dim ) {
      glutAddMenuEntry("[w]  Toggle vector/tensor",'w');
	  glutAddMenuEntry("     Toggle displacement",'d');
      glutAddMenuEntry("[v]  Toggle streamlines",'v');
      glutAddMenuEntry("     Flush streamlines",'f');
      glutAddMenuEntry("     Particle advection",'u');
      if ( mesh->dim == 2 )
        glutAddMenuEntry("     Critical points",'c');
    }
    if ( mesh->dim == 2 ) {
      glutAddMenuEntry("[k]  Toggle elevation",'k');
      glutAddMenuEntry("[K]  Elevation coeff",'K');
    }
  }

  /* Show misc. items */
  vmenu = glutCreateMenu(menuItem);
  glutAddMenuEntry("[A] Axis",'A');
  glutAddMenuEntry("[B] Bounding box",'B');
  glutAddMenuEntry("[G] Grid ",'G');
  if ( sc->glist )
    glutAddMenuEntry("[g] Geometric items",'g');
  glutAddMenuEntry("[j] Toggle Info",'j');
  glutAddMenuEntry("[P] Toggle Point num",'P');
  glutAddMenuEntry("[F] Toggle Face num",'F');
  if ( mesh->nvn ) {
    glutAddMenuEntry("[N] Toggle normals",'N');
    glutAddMenuEntry("[O] Revert normals",'O');
  }

  /* clipping menu */
  clmenu = glutCreateMenu(menuClip);
  glutAddMenuEntry("[F1] Toggle clip",'C');
  glutAddMenuEntry("[F2] Edit clip",'E');
  glutAddMenuEntry("[F3] Freeze clip",'F');
  glutAddMenuEntry("     Inverse orient",'I');
  if ( mesh->ntet+mesh->nhex > 0 )
    glutAddMenuEntry("     Toggle capping",'K');
  glutAddMenuEntry("     Toggle plane",'H');
  glutAddMenuEntry(" --  Reset clip",'R');
  if ( mesh->ntet+mesh->nhex > 0 ) {
    sc->clip->active |= C_VOL;
    glutAddMenuEntry("[F4] Toggle Vclip",'Z');
  }

  /* feature menu */
  femenu = glutCreateMenu(menuFeature);
  if ( mesh->ntet+mesh->nhex > 0 )
    glutAddMenuEntry("[F5] Toggle Vshrink",'V');
  else
    glutAddMenuEntry("[F5] Toggle shrink",'S');
  glutAddMenuEntry("[F6] Increase shrink",'I');
  glutAddMenuEntry("[F7] Decrease shrink",'i');
  if ( mesh->dim == 3 )
    glutAddMenuEntry("Toggle splitview",'s');

  /* view handler menu */
  vwmenu = glutCreateMenu(menuView);
  glutAddMenuEntry("[i] Reset",'R');
  glutAddMenuEntry("[Alt-c] Copy",'C');
  glutAddMenuEntry("[Alt-p] Paste",'P');
  glutAddMenuEntry("[Alt+l] Link",'L');
  glutAddMenuEntry("[Alt+u] Unlink",'U');
  /*glutAddMenuEntry("[Alt+d] Duplicate",'D');*/

  /* animation menu */
  amenu = glutCreateMenu(menuAnim);
  glutAddMenuEntry("[a] Toggle Anim",'A');
  glutAddMenuEntry("Toggle ImgSave",'I');
  if ( option == SEQUENCE || option == SEQUENCE + PARTICLE ) {
    glutAddMenuEntry("Play  sequence",'S');
    glutAddMenuEntry("First mesh",'f');
    glutAddMenuEntry("Last  mesh",'l');
    glutAddMenuEntry("Next mesh",'n');
    glutAddMenuEntry("Prev mesh",'p');
  }
  else if ( option == MORPHING ) {
    glutAddMenuEntry("Start/Stop morph.",'M');
    glutAddMenuEntry("Toggle AutoReverse",'R');
  }

  /* trajectoire menu */
  if ( mesh->dim == 3 || mesh->nbb ) {
    trmenu = glutCreateMenu(menuTrajet);
    glutAddMenuEntry("New Ctrl point",'C');
    glutAddMenuEntry("Toggle path",'S');
    glutAddMenuEntry("Follow path",'F');
    glutAddMenuEntry("Load path",'L');
    glutAddMenuEntry("Save path",'W');
  }
  else
    trmenu = 0;

  /* main scene menu */
  menu = glutCreateMenu(menuScene);
  glutAddSubMenu("File",fmenu);
  glutAddSubMenu("Render mode",mmenu);
  glutAddSubMenu("Colors, Materials",cmenu);
  if ( mesh->nbb )
    glutAddSubMenu("Data",smenu);
  glutAddSubMenu("Items",vmenu);
  if ( mesh->dim == 3 )
    glutAddSubMenu("Clipping",clmenu);
  glutAddSubMenu("Features",femenu);
  glutAddSubMenu("View",vwmenu);
  glutAddSubMenu("Animation",amenu);
  /*if ( trmenu )
    glutAddSubMenu("Trajectory",trmenu);
  */
  glutAddMenuEntry("",'\0');
  glutAddMenuEntry("Close window",'X');
  glutAddMenuEntry("Quit",'q');


  return(1);
}
