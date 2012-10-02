#include <stdarg.h>
#include "medit.h"
#include "extern.h"
#include "sproto.h"

#define Width   630
#define Height  250
#define MAX_MAT  32

extern int refmat,reftype;

/* globals */
pScene    main_scene;
GLint     matwin=0,m_subwin;
GLint     selection;

static float colors[MAX_MAT+1][3] = {
  {0.1, 0.4, 0.9},    /* blue */
  {1.0, 0.0, 0.0},    /* red */
  {0.0, 1.0, 0.0},    /* green */
  {1.0, 1.0, 0.0},    /* yellow */
  {0.0, 1.0, 1.0},    /* cyan */ 
  {1.0, 0.5, 0.0},    /* orange */
  {0.5, 0.0, 1.0},    /* violet */ 
  {0.0, 0.0, 0.4},    /* dark blue */
  {0.0, 0.4, 0.0},    /* dark green */
  {0.4, 0.0, 0.0},    /* dark red */
  {1.0, 1.0, 0.5},
  {1.0, 0.5, 1.0},  
  {1.0, 0.5, 0.5},  
  {1.0, 0.5, 0.0},    /* orange */
  {1.0, 0.0, 1.0},  
  {1.0, 0.0, 0.5},  
  {0.5, 1.0, 1.0},  
  {0.5, 1.0, 0.5},  
  {0.5, 1.0, 0.0},
  {0.5, 0.5, 1.0},  
  {0.5, 0.5, 0.5},  
  {0.5, 0.5, 0.0},
  {0.5, 0.0, 0.5},  
  {0.5, 0.0, 0.0},
  {0.0, 1.0, 0.5},  
  {0.0, 0.5, 1.0},  
  {0.0, 0.5, 0.5},  
  {0.0, 0.5, 0.0},
  {0.0, 0.0, 0.5},
  {0.4, 0.4, 0.0},    /* dark yellow */
  {0.0, 0.4, 0.4},    /* dark cyan */ 
  {0.3, 0.7, 0.9},    /* default blue */
  {0.3, 0.7, 0.9}     /* default blue */
};

cell ambient[4] = {
    { 1, 120, 60, 0.0, 1.0, 0.0, 0.01,
      "Specifies R coordinate of ambient vector.", "%.2f" },
    { 2, 180, 60, 0.0, 1.0, 0.0, 0.01,
      "Specifies G coordinate of ambient vector.", "%.2f" },
    { 3, 240, 60, 0.0, 1.0, 0.0, 0.01,
      "Specifies B coordinate of ambient vector.", "%.2f" },
    { 4, 300, 60, 0.0, 1.0, 0.0, 0.01,
      "Specifies A coordinate of ambient vector.", "%.2f" },
};
cell diffuse[4] = {
    { 5, 120, 90, 0.0, 1.0, 0.0, 0.01,
      "Specifies R coordinate of diffuse vector.", "%.2f" },
    { 6, 180, 90, 0.0, 1.0, 0.0, 0.01,
      "Specifies G coordinate of diffuse vector.", "%.2f" },
    { 7, 240, 90, 0.0, 1.0, 0.0, 0.01,
      "Specifies B coordinate of diffuse vector.", "%.2f" },
    { 8, 300, 90, 0.0, 1.0, 0.0, 0.01,
      "Specifies A coordinate of diffuse vector.", "%.2f" },
};
cell specular[4] = {
    { 9, 120, 120, 0.0, 1.0, 0.0, 0.01,
      "Specifies R coordinate of specular vector.", "%.2f" },
    { 10, 180, 120, 0.0, 1.0, 0.0, 0.01,
      "Specifies G coordinate of specular vector.", "%.2f" },
    { 11, 240, 120, 0.0, 1.0, 0.0, 0.01,
      "Specifies B coordinate of specular vector.", "%.2f" },
    { 12, 300, 120, 0.0, 1.0, 0.0, 0.01,
      "Specifies A coordinate of specular vector.", "%.2f" },
};
cell emission[4] = {
    { 13, 120, 150, 0.0, 1.0, 0.0, 0.01,
      "Specifies R coordinate of emission vector.", "%.2f" },
    { 14, 180, 150, 0.0, 1.0, 0.0, 0.01,
      "Specifies G coordinate of emission vector.", "%.2f" },
    { 15, 240, 150, 0.0, 1.0, 0.0, 0.01,
      "Specifies B coordinate of emission vector.", "%.2f" },
    { 16, 300, 150, 0.0, 1.0, 0.0, 0.01,
      "Specifies A coordinate of emission vector.", "%.2f" },
};
cell shininess[1] = {
    { 17, 120, 180, 3.0, 128.0, 0.0, 0.5,
      "Specifies value of shininess.", "%.2f" },
};


void matInit(pScene sc) {
  pMaterial  pm;
  int        m,mm;

  /* default */
  if ( !sc->material ) {
    sc->material = (pMaterial)M_calloc(2+sc->par.nbmat,sizeof(Material),"matinit");
	assert(sc->material);
    sc->matsort = (int*)M_calloc(2+sc->par.nbmat,sizeof(int),"matinit");
	assert(sc->matsort);
  }

  /* store color in table */
  for (m=0; m<=sc->par.nbmat; m++) {
    pm = &sc->material[m];
    /* diffuse : primary color */
    mm = m % MAX_MAT;
    memcpy(pm->dif,colors[mm],3*sizeof(float));
    pm->dif[3] = 1.0;
    /* ambient : grey level */
    pm->amb[0] = pm->amb[1] = pm->amb[2] = 0.2;  pm->amb[3] = 1.0;
    /* emission: null (pas un neon!) */
    pm->emi[0] = pm->emi[1] = pm->emi[2] = 0.0;  pm->emi[3] = 1.0;
    /* specular: soleil blanc */
    pm->spe[0] = pm->spe[1] = pm->spe[2] = 0.4;  pm->spe[3] = 1.0;
    /* shininess: etalement des reflections spec. */
    pm->shininess = 80.0;
    if ( m != DEFAULT_MAT )
      sprintf(pm->name,"%s%.2d","MAT",m);
    else
      strcpy(pm->name,"DEFAULT_MAT");
    pm->flag = 0;
    if ( !pm->ref )  pm->ref = m; 
    pm->sort = m;
  }
}

void matSort(pScene sc) {
  pMaterial  pm,pm1;
  int        m,mm,transp;

  /* sort materials */
  if ( !quiet ) fprintf(stdout,"   Sorting %d materials\n",sc->par.nbmat);
  transp = 0;

  for (m=0; m<sc->par.nbmat; m++) {
    pm  = &sc->material[m];
    pm->sort       = m;
    sc->matsort[m] = m;
    if ( pm->dif[3] < 0.999 )  transp++;
  }
  if ( !transp )  return;

  /* sorting */
  if ( ddebug ) fprintf(stdout,"   %d translucent\n",transp);
  mm = sc->par.nbmat - 1;
  for (m=0; m<sc->par.nbmat-transp-1; m++) {
    pm = &sc->material[m];
    if ( pm->dif[3] < 0.995 ) {
      pm1 = &sc->material[mm];
      sc->matsort[mm] = m;
      sc->matsort[m]  = mm;
      pm->sort  = mm;
      pm1->sort = m;
      mm--;
    }
  }
}

int matRef(pScene sc,int ref) {
  pMaterial  pm;
  int        m;

  if ( !ref )  return(ref);
  for (m=1; m<sc->par.nbmat; m++) {
    pm = &sc->material[m];
    if ( pm->ref == ref )  return(m);
  }
  if ( sc->par.nbmat < 2 ) return(0);
  m = 1+(ref-1) % (sc->par.nbmat-1);
  return(m);
}

/* reshape material window */
void matReshape(int width,int height) {
  glutSetWindow(matwin);
  glViewport(0,0,Width,Height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,Width,Height,0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glutReshapeWindow(Width,Height);
}

static void cellInit(int mat) {
  pMaterial pm;
  int       i;

  pm = &main_scene->material[mat];
  for (i=0; i<4; i++) {
    ambient[i].value  = pm->amb[i];
    diffuse[i].value  = pm->dif[i];
    specular[i].value = pm->spe[i];
    emission[i].value = pm->emi[i];
  }
  shininess[0].value = pm->shininess;
  glutPostRedisplay();
}

static int cellCopy(int mat) {
  pMaterial pm;
  int       i,dosort;

  pm = &main_scene->material[mat];
  dosort = pm->dif[3] != diffuse[3].value;

  for (i=0; i<4; i++) {
    pm->amb[i] = ambient[i].value;
    pm->dif[i] = diffuse[i].value;
    pm->spe[i] = specular[i].value;
    pm->emi[i] = emission[i].value;
  }
  pm->shininess = shininess[0].value;
  return(dosort);
}

static void cellDraw(cell* cel) {
  if ( selection == cel->id ) {
    glColor3ub(255,255,0);
    drwstr(10, 240,cel->info);
    glColor3ub(255,0,0);
  }
  else
    glColor3ub(0,255,128);
  drwstr(cel->x,cel->y,cel->format,cel->value);
}

static int cellHit(cell* cel,int x,int y) {
  if ( x > cel->x    && x < cel->x + 60 &&
       y > cel->y-30 && y < cel->y+10 )
    return cel->id;
  else
    return 0;
}

static void cellUpdate(cell* cel,int update) {
  if ( selection != cel->id )
    return;

  cel->value += update * cel->step;

  if (cel->value < cel->min)
    cel->value = cel->min;
  else if (cel->value > cel->max) 
    cel->value = cel->max;
}


void matsubReshape(int x,int y) {
  GLfloat    sunpos[4] = {5.,2.,10.0,0.};

  glutSetWindow(m_subwin);
  glViewport(0,0,200,200);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(50.0,1.0,0.01f,100.0f);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
  glLightfv(GL_LIGHT0,GL_POSITION,sunpos);  
  glEnable(GL_LIGHT0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0,0.0,-5.0);
  glutSetWindow(matwin);
}

void matsubDisplay() {
  GLfloat    amb[4],dif[4],emi[4],spe[4];
  int        i,transp;

  glutSetWindow(m_subwin);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glClearColor(0.0f,0.0f,0.0f,1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glBegin(GL_QUADS);
    glColor3f(0.8,0.,0.);
    glVertex3f(-2.0,-1.5,-2.0);
    glVertex3f( 0.0,-1.5,-2.0);
    glVertex3f( 0.0,-1.5,5.0);
    glVertex3f(-2.0,-1.5,5.0);
    
    glColor3f(0.,0.8,0.);
    glVertex3f( 0.0,-1.5,-2.0);
    glVertex3f( 2.0,-1.5,-2.0);
    glVertex3f( 2.0,-1.5,5.0);
    glVertex3f( 0.0,-1.5,5.0);
    
    glColor3f(0.,0.,0.8);
    glVertex3f( 0.0,-1.5,-5.0);
    glVertex3f( 2.0,-1.5,-5.0);
    glVertex3f( 2.0,-1.5,-2.0);
    glVertex3f( 0.0,-1.5,-2.0);

    glColor3f(0.4,0.4,0.4);
    glVertex3f(-2.0,-1.5,-5.0);
    glVertex3f( 0.0,-1.5,-5.0);
    glVertex3f( 0.0,-1.5,-2.0);
    glVertex3f(-2.0,-1.5,-2.0);
  glEnd();  

  glEnable(GL_LIGHTING);
  transp = ambient[3].value < 0.999 || 
           diffuse[3].value < 0.999 || specular[3].value < 0.999;
  if ( transp ) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
  }
  for (i=0; i<4; i++) {
    amb[i] = ambient[i].value;
    dif[i] = diffuse[i].value;
    spe[i] = specular[i].value;
    emi[i] = emission[i].value;
  }
  glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,dif);
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,amb);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,spe);
  glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,emi);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&shininess[0].value);

  glutSolidSphere(1.5,20,10);
  if ( transp ) {
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
  }
  glutSwapBuffers();
  glutSetWindow(matwin);
}

/* display mat properties */
void matDisplay() {
  pMaterial  pm;

  glutSetWindow(matwin);
  glClearColor(0.1,0.1,0.1,1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  pm = &main_scene->material[refmat];

  glColor3ub(255,255,0);
  setFont("helvetica", 18);
  drwstr(10,20,pm->name);

  glColor3ub(255,255,255);
  setFont("helvetica", 18);
  drwstr(30,ambient[0].y,"Ambient(");
  drwstr(170,ambient[0].y,",");
  drwstr(230,ambient[1].y,","); 
  drwstr(290,ambient[2].y,",");
  drwstr(350,ambient[3].y,");");

  cellDraw(&ambient[0]);
  cellDraw(&ambient[1]);
  cellDraw(&ambient[2]);
  cellDraw(&ambient[3]);

  glColor3ub(255,255,255);
  drwstr(40,diffuse[0].y,"Diffuse(");
  drwstr(170,diffuse[0].y,",");
  drwstr(230,diffuse[1].y,","); 
  drwstr(290,diffuse[2].y,",");
  drwstr(350,diffuse[3].y,");");

  cellDraw(&diffuse[0]);
  cellDraw(&diffuse[1]);
  cellDraw(&diffuse[2]);
  cellDraw(&diffuse[3]);

  glColor3ub(255,255,255);
  drwstr(25,specular[0].y,"Specular(");
  drwstr(170,specular[0].y,",");
  drwstr(230,specular[1].y,","); 
  drwstr(290,specular[2].y,",");
  drwstr(350,specular[3].y,");");

  cellDraw(&specular[0]);
  cellDraw(&specular[1]);
  cellDraw(&specular[2]);
  cellDraw(&specular[3]);

  glColor3ub(255,255,255);
  drwstr(25,emission[0].y,"Emission(");
  drwstr(170,emission[0].y,",");
  drwstr(230,emission[1].y,","); 
  drwstr(290,emission[2].y,",");
  drwstr(350,emission[3].y,");");

  cellDraw(&emission[0]);
  cellDraw(&emission[1]);
  cellDraw(&emission[2]);
  cellDraw(&emission[3]);

  glColor3ub(255,255,255);
  drwstr(20,shininess[0].y,"Shininess(");
  drwstr(350,shininess[0].y,");");
  cellDraw(&shininess[0]);

  if ( !selection ) {
    glColor3ub(255,255,0);
    setFont("helvetica",18);
    drwstr(10,240,"Click and move the mouse to modify values ('A' to apply).");
  }

  glutSwapBuffers();
}

static int old_y = 0;

void matMouse(int button,int state,int x,int y) {
  glutSetWindow(matwin);
  selection = 0;
  if ( button != GLUT_LEFT_BUTTON ) return; 
  if ( state == GLUT_DOWN ) {
    /* mouse should only hit _one_ of the cells, so adding up all
       the hits just propagates a single hit. */
    selection += cellHit(&ambient[0],x,y);
    selection += cellHit(&ambient[1],x,y);
    selection += cellHit(&ambient[2],x,y);
    selection += cellHit(&ambient[3],x,y);
    selection += cellHit(&diffuse[0],x,y);
    selection += cellHit(&diffuse[1],x,y);
    selection += cellHit(&diffuse[2],x,y);
    selection += cellHit(&diffuse[3],x,y);
    selection += cellHit(&specular[0],x,y);
    selection += cellHit(&specular[1],x,y);
    selection += cellHit(&specular[2],x,y);
    selection += cellHit(&specular[3],x,y);
    selection += cellHit(&emission[0],x,y);
    selection += cellHit(&emission[1],x,y);
    selection += cellHit(&emission[2],x,y);
    selection += cellHit(&emission[3],x,y);
    selection += cellHit(&shininess[0],x,y);
  }
  else {
    /* update material */
    glutSetWindow(m_subwin);
  }
  old_y = y;
  glutPostRedisplay();
}

static void matMotion(int x,int y) {
  cellUpdate(&ambient[0],old_y-y);
  cellUpdate(&ambient[1],old_y-y);
  cellUpdate(&ambient[2],old_y-y);
  cellUpdate(&ambient[3],old_y-y);

  cellUpdate(&diffuse[0],old_y-y);
  cellUpdate(&diffuse[1],old_y-y);
  cellUpdate(&diffuse[2],old_y-y);
  cellUpdate(&diffuse[3],old_y-y);

  cellUpdate(&specular[0],old_y-y);
  cellUpdate(&specular[1],old_y-y);
  cellUpdate(&specular[2],old_y-y);
  cellUpdate(&specular[3],old_y-y);

  cellUpdate(&emission[0],old_y-y);
  cellUpdate(&emission[1],old_y-y);
  cellUpdate(&emission[2],old_y-y);
  cellUpdate(&emission[3],old_y-y);

  cellUpdate(&shininess[0],old_y-y);
  old_y = y;
  glutPostRedisplay();

  glutSetWindow(m_subwin);
  glutPostRedisplay();
}

void matKeyboard(unsigned char key,int x,int y) {

  if ( isdigit(key) )  return;
  glutSetWindow(matwin);

  switch(key) {
  case 'a':
  case 'A':
    if ( cellCopy(refmat) )
      matSort(main_scene);
    matwin = 0;
    glutHideWindow();
    /*glutPostRedisplay();*/
    glutSetWindow(main_scene->idwin);
    doLists(main_scene,cv.mesh[main_scene->idmesh]);
    glutPostRedisplay();
    break;
  case 'q':
  case 27:
    matwin = 0;
    glutHideWindow();
    glutPostRedisplay();
    break;
  }
}

void matEdit(pScene sc) {
  /* init window */
  if ( matwin || refmat < 0 )  return;
  matwin = m_subwin = 0;

  /* create 2d window */
  glutInitWindowSize(Width,Height);
  matwin = glutCreateWindow("Material Properties");
  if ( !matwin ) return;
  main_scene = sc;

  /* set callbacks */
  glutSetIconTitle("MatProp");
  glutReshapeFunc(matReshape);
  glutDisplayFunc(matDisplay); 
  glutKeyboardFunc(matKeyboard);
  glutMouseFunc(matMouse);
  glutMotionFunc(matMotion);

  /* create 3d sub-window */
  glutSetWindow(matwin);
  m_subwin = glutCreateSubWindow(matwin,400,20,200,200);
  if ( ddebug)   printf("window %d\n",(int)matwin);
  glutDisplayFunc(matsubDisplay);
  glutReshapeFunc(matsubReshape);

  glutSetWindow(matwin);
  cellInit(refmat);
}


#ifdef __cplusplus
}
#endif
