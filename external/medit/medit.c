/* medit.c      mesh visualization tool
 *
 * Written by Pascal Frey, LJLL
 * Copyright (c) Inria, 1999-2007. All rights reserved. */

#include "medit.h"
#include "compil.date"
#ifdef ppc
#include <unistd.h>
#endif


/* global variables (see extern.h) */
GLboolean hasStereo = 1;
Canvas    cv;
mytime    ctim[TIMEMAX]; 
ubyte     ddebug,animate,saveimg,imgtype,infogl,fullscreen;
ubyte     quiet,option,morphing,stereoMode;
int       menu,amenu,fmenu,femenu,vmenu,mmenu,smenu;
int       clmenu,cmenu,vwmenu,txmenu,trmenu;
int       animdep,animfin;

/**********************/
/*  Rajout pour popen */
ubyte     dpopen,dpopensol,dpopenbin;

static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  break;
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Program killed\n");  break;
  }
  exit(1);
}


static void endcod() {
  chrono(OFF,&ctim[0]);
  fprintf(stdout,"\n Total running seconds:  %.2f\n",gttime(ctim[0]));
  fprintf(stdout," Thank you for using Medit.\n");
}


static void grInfo(void) {
  GLboolean  b;
  GLint      i,win;

  win = glutCreateWindow("Info");
  fprintf(stdout,"Graphic info:\n");
  fprintf(stdout," GL Vendor:\t%s\n",glGetString(GL_VENDOR));
  fprintf(stdout," GL Version:\t%s\n",glGetString(GL_VERSION));
  fprintf(stdout," GL Renderer:\t%s\n\n",glGetString(GL_RENDERER));
  glGetBooleanv(GL_RGBA_MODE,&b);
  if ( b )  fprintf(stdout,"  RGBA Mode\n");
  glGetBooleanv(GL_DOUBLEBUFFER,&b);
  if ( b )  fprintf(stdout,"  Double Buffer\n");
  glGetBooleanv(GL_STEREO,&b);
  if ( b )  fprintf(stdout,"  Stereo\n");
  glGetIntegerv(GL_AUX_BUFFERS,&i);
  if ( i )  fprintf(stdout,"  Auxilary Buffers\t%2d\n",(int)i);
  glGetIntegerv(GL_INDEX_BITS,&i);
  if ( i )  fprintf(stdout,"  Index Bits\t\t%2d\n",(int)i);
  glGetIntegerv(GL_RED_BITS,&i);
  fprintf(stdout,"  RGBA Bits\t\t%2d",(int)i);
  glGetIntegerv(GL_GREEN_BITS,&i);
  fprintf(stdout,"\t%2d",(int)i);
  glGetIntegerv(GL_BLUE_BITS,&i);
  fprintf(stdout,"\t%2d",(int)i);
  glGetIntegerv(GL_ALPHA_BITS,&i);
  fprintf(stdout,"\t%2d\n",(int)i);
  glGetIntegerv(GL_ACCUM_RED_BITS,&i);
  fprintf(stdout,"  Accum RGBA Bits\t%2d",(int)i);
  glGetIntegerv(GL_ACCUM_GREEN_BITS,&i);
  fprintf(stdout,"\t%2d",(int)i);
  glGetIntegerv(GL_ACCUM_BLUE_BITS,&i);
  fprintf(stdout,"\t%2d",(int)i);
  glGetIntegerv(GL_ACCUM_ALPHA_BITS,&i);
  fprintf(stdout,"\t%2d\n",(int)i);
  glGetIntegerv(GL_DEPTH_BITS,&i);
  fprintf(stdout,"  Depth Bits\t\t%2d\n",(int)i);
  glGetIntegerv(GL_STENCIL_BITS,&i);
  fprintf(stdout,"  Stencil Bits\t\t%2d\n",(int)i);
  
  exit(1);
}


int medit0() {
  pMesh    mesh;
  char     data[128],*name;
  int      k,l,ret;
  clock_t  ct;

  /* default */
  //  fprintf(stdout," \n medit0() \n");
  fprintf(stdout," Loading data file(s)\n");
  ct = clock();

  /* enter name */
  if ( !cv.nbm ) {
    fprintf(stdout,"  File name(s) missing. Please enter : ");
    fflush(stdout); fflush(stdin);
    fgets(data,120,stdin);
    if ( !strlen(data) ) {
      fprintf(stdout,"  ## No data\n");
      return(0);
    }

    /* parse file name(s) */
    name = strtok(data," \n");
    while( name ) {
      if ( !cv.mesh[cv.nbm] ) {
        cv.mesh[cv.nbm] = (pMesh)M_calloc(1,sizeof(Mesh),"medit0.mesh");
        if ( !cv.mesh[cv.nbm] )  return(0);
      }
      /*(cv.mesh[cv.nbm])->name = calloc(strlen(name)+1,sizeof(char));*/
      strcpy(cv.mesh[cv.nbm]->name,name);
      name = strtok(NULL," \n\0");
      if ( ++cv.nbm == MAX_MESH )  break;
    }
    if ( !cv.nbm ) return(0);
  }

  if ( !cv.nbm ) { 
    fprintf(stdout,"  Number of mesh missing:. Please enter : "); 
    fflush(stdout); fflush(stdin);
    fgets(data,120,stdin);
    cv.nbm = atoi(data);
  }

  /* read mesh(es) */
  k = 0;
  do {
    if ( !cv.mesh[k] ) {
#ifdef IGL
      cv.mesh[k] = static_cast<pMesh>(M_calloc(1,sizeof(Mesh),"medit0.mesh"));
#else
      cv.mesh[k] = M_calloc(1,sizeof(Mesh),"medit0.mesh");
#endif
      if ( !cv.mesh[k] )  return(0);
    }
    mesh = cv.mesh[k];
    mesh->typ = 0;
    ret  = loadMesh(mesh);
    if ( ret < 0 ) {
      mesh->typ = 1;
      ret = inmsh2(mesh);
      if ( !ret ) {
        mesh->typ = 2;
        ret = loadGIS(mesh);
      }
    }
    if ( ret <= 0 ) {
      for (l=k+1; l<cv.nbm; l++)
	    cv.mesh[l-1] = cv.mesh[l];
      cv.nbm--;
      k--;
      continue;
    }

    /* compute mesh box */  
    if ( (mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq) )  
    {
      fprintf(stderr,"Alec: skipping meshSurf\n");
      //meshSurf(mesh);
    }
    meshBox(mesh,1);
    if ( !quiet || true)  meshInfo(mesh);

    /* read metric  */   
    if ( !loadSol(mesh,mesh->name,1) )
      bbfile(mesh);
    if ( !quiet && mesh->nbb )
      fprintf(stdout,"    Solutions  %8d\n",mesh->nbb);
  }
  while ( ++k < cv.nbm );
  cv.nbs = cv.nbm;

  ct = difftime(clock(),ct);
  fprintf(stdout,"  Input seconds:     %.2f\n",
          (double)ct/(double)CLOCKS_PER_SEC);

  return(cv.nbm);
}


int medit0_popen() {
  pMesh    mesh;
  char     data[128],*name;
  int      k,l,ret;
  clock_t  ct;

  /* default */
  /*fprintf(stdout," \n medit0() \n");*/
  fprintf(stdout," Loading data file(s)\n");
  ct = clock();


  /* enter number of mesh */
  if ( !cv.nbm ) { 
    fprintf(stdout,"  Number of mesh missing:. Please enter : "); 
    fflush(stdout); fflush(stdin);
    fgets(data,128,stdin);
    cv.nbm = atoi(data);
  }

  /* read mesh(es) */
  k = 0;
  do {
    // printf("mesh number %i\n",k+1);
    if ( !cv.mesh[k] ) {
#ifdef IGL
      cv.mesh[k] = static_cast<pMesh>(M_calloc(1,sizeof(Mesh),"medit0.mesh"));
#else
      cv.mesh[k] = M_calloc(1,sizeof(Mesh),"medit0.mesh");
#endif
      if ( !cv.mesh[k] )  return(0);
    }
    mesh = cv.mesh[k];
    mesh->typ = 0;

    //fgets(data,128,stdin);
    //name = data;
    //printf("data=%s\n",data);
    //name = "toto.dat";
    //strcpy(cv.mesh[k]->name,name);

    if(dpopenbin) 
      ret  = loadMesh_popen_bin(mesh);
    else
      ret  = loadMesh_popen(mesh);

    /* compute mesh box */  
    if ( (mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq) )  
      meshSurf(mesh);
    meshBox(mesh,1);
    if ( !quiet || true)  meshInfo(mesh);

    /*  /\* read metric *\/     // a changer lecture .sol et .bb */
    /*    if ( !loadSol_popen(mesh,mesh->name,1) ) */
    /*       bbfile_popen(mesh); */
    /*     if ( !quiet && mesh->nbb ) */
    /*       fprintf(stdout,"    Solutions  %8d\n",mesh->nbb); */
    if( dpopensol ){
      if(dpopenbin) 
	loadSol_popen_bin(mesh,mesh->name,1);
      else
	loadSol_popen(mesh,mesh->name,1);
    }
    if ( !quiet && mesh->nbb )
      fprintf(stdout,"    Solutions  %8d\n",mesh->nbb);
  }
  while ( ++k < cv.nbm );
  cv.nbs = cv.nbm;

  ct = difftime(clock(),ct);
  fprintf(stdout,"  Input seconds:     %.2f\n",
          (double)ct/(double)CLOCKS_PER_SEC);

  return(cv.nbm);
}

int medit1() {
  pScene   scene;
  pMesh    mesh;
  int      k;
  clock_t  ct;

  /* create grafix */
  fprintf(stdout,"\n medit1() \n");
  fprintf(stdout,"\n Building scene(s)\n");
  ct = clock();
  for (k=0; k<cv.nbs; k++) {
    if ( !cv.scene[k] ) {
      cv.scene[k] = (pScene)M_calloc(1,sizeof(Scene),"medit1.scene");
      if ( !cv.scene[k] )  return(0);
    }
    scene = cv.scene[k];
    if ( !cv.mesh[k] ) {
      cv.mesh[k] = (pMesh)M_calloc(1,sizeof(Mesh),"medit1.mesh");
      if ( !cv.mesh[k] )  return(0);
    }
    mesh  = cv.mesh[k];

    fprintf(stdout,"  Creating scene %d\n",k+1);
    parsop(scene,mesh);
    meshRef(scene,mesh);
    matSort(scene);
    
    if ( option == ISOSURF ) {
	  if ( !mesh->nbb ) return(0);
      setupPalette(scene,mesh);
	  tetraIsoPOVray(scene,mesh);
	}
    else if ( !createScene(scene,k) ) {
      fprintf(stderr,"  ## Unable to create scene\n");
      return(0);
    }
  }
  ct = difftime(clock(),ct);
  fprintf(stdout,"  Scene seconds:     %.2f\n",(double)ct/(double)CLOCKS_PER_SEC);

  return(1);
}


int main(int argc,char *argv[]) {
  int    type;
  char   pwd[1024];

#ifdef ppc
  if ( !getwd(pwd) )  exit(2);
#endif

  fprintf(stdout,"  -- Medit,  Release %s (%s)\n",ME_VER,ME_REL);
  fprintf(stdout,"     %s.\n",ME_CPY);
  fprintf(stdout,"     compiled: %s.\n\n",COMPIL);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* default values */
  option     = STANDARD;
  saveimg    = GL_FALSE;
  imgtype    = P6;
  animate    = FALSE;
  morphing   = FALSE;
  fullscreen = FALSE;
  animdep    = 0;
  animfin    = 0;
  ddebug     = FALSE;
  quiet      = TRUE;//FALSE;
  stereoMode = 0;
  cv.nbm = cv.nbs = 0;

  /* default value for popen */
  dpopen = FALSE;
  dpopenbin = FALSE;
  dpopensol = FALSE;

  /* init grafix */
  parsar(argc,argv);
  //printf("fin de parsar");

  if ( option == ISOSURF ) {
    fprintf(stdout,"ISOSURF");
    if ( !medit0() )  exit(1);
    if ( !medit1() )  exit(1);
    return(0);
  }

  glutInit(&argc,argv);
  //printf("fin de glutInit");
#ifdef ppc
  chdir(pwd);
#endif

  chrono(ON,&ctim[0]);
  if ( stereoMode == MONO )
    type = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH;
  else
    type = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO;
  glutInitDisplayMode(type);
#ifdef IGL
  glutInitDisplayString( "rgba depth double samples>=8 ");
#endif
  if ( infogl )  grInfo();

  /* call animate or normal mode */
  if( dpopen == FALSE ){
    
    switch (option) {
    case STANDARD:

    case SCHNAUZER:
      if ( !medit0() )  exit(1);
      if ( !medit1() )  exit(1);
      break;
    case SEQUENCE:
      if ( !animat() )  exit(1);
      break;
    case SEQUENCE+PARTICLE:
      
      if ( !animat() )  exit(1);
      break;
    case MORPHING:
      
      if ( !medit0() )  exit(1);
      if ( !modeMorphing() )  exit(1);
      morphing = GL_FALSE;
      break;
    default:
      fprintf(stderr,"  ## Unrecognized option %d\n",option);
      exit(1);
      break;
    }
  }
  else{
   switch (option) {
   case STANDARD:
     
   case SCHNAUZER:
     
     if ( !medit0_popen() )  exit(1);
     if ( !medit1() )  exit(1);
     break;
  
  case SEQUENCE:
    
    if ( !animat() )  exit(1);
    break;
  case SEQUENCE+PARTICLE:
    
    if ( !animat() )  exit(1);
    break;
  case MORPHING:
    
    if ( !medit0_popen() )  exit(1);
    if ( !modeMorphing() )  exit(1);
    morphing = GL_FALSE;
    break;
  default:
    fprintf(stderr,"  ## Unrecognized option %d\n",option);
    exit(1);
    break;
  }
  }

  /* main grafix loop */
  fprintf(stdout,"\n Rendering scene(s)\n");
  glGetBooleanv(GL_STEREO,&hasStereo);
  glutMainLoop();

  return(0);
}
