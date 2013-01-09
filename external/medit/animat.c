#include "medit.h"
#include "extern.h"
#include "sproto.h"

extern void resetLists(pScene ,pMesh );


static int getmesh(pMesh mesh,int range) {
  int     k,ret;
  char   *ptr,data[256];
  static char base[256];

  /* allocate mesh structure */
  if ( ddebug ) fprintf(stdout,"getmesh: create mesh %d\n",range);

  if ( mesh->tria )  M_free(mesh->tria);
  if ( mesh->quad )  M_free(mesh->quad);
  if ( mesh->edge )  M_free(mesh->edge);
  if ( mesh->tetra ) M_free(mesh->tetra);
  if ( mesh->hexa )  M_free(mesh->hexa);
  if ( mesh->adja )  M_free(mesh->adja);
  if ( mesh->voy  )  M_free(mesh->voy);
  if ( mesh->point ) M_free(mesh->point);
  if ( mesh->extra ) {
    if ( mesh->extra->iv )  M_free(mesh->extra->nv);
    if ( mesh->extra->it )  M_free(mesh->extra->nt);
    if ( mesh->extra->iq )  M_free(mesh->extra->nq);
    if ( mesh->extra->n )   M_free(mesh->extra->n);
    M_free(mesh->extra);
#ifdef IGL
    mesh->extra = NULL;
#else
    mesh->extra = (void*)0;
#endif
  }
  if ( mesh->sol && mesh->nbb ) {
    if ( (mesh->dim==2 && mesh->nfield==3) || (mesh->dim==3 && mesh->nfield==6) )
      for (k=1; k<=mesh->nbb; k++)
        free(mesh->sol[k].m);
    M_free(mesh->sol);
#ifdef IGL
    mesh->sol = NULL;
#else
    mesh->sol = (void*)0;
#endif
  }
#ifdef IGL
  mesh->point = NULL;
#else
  mesh->point = (void*)0;
#endif
#ifdef IGL
  mesh->tria  = NULL;
#else
  mesh->tria  = (void*)0;
#endif
#ifdef IGL
  mesh->quad  = NULL;
#else
  mesh->quad  = (void*)0;
#endif
#ifdef IGL
  mesh->edge  = NULL;
#else
  mesh->edge  = (void*)0;
#endif
#ifdef IGL
  mesh->tetra = NULL;
#else
  mesh->tetra = (void*)0;
#endif
#ifdef IGL
  mesh->hexa  = NULL;
#else
  mesh->hexa  = (void*)0;
#endif
#ifdef IGL
  mesh->adja  = NULL;
#else
  mesh->adja  = (void*)0;
#endif
#ifdef IGL
  mesh->voy   = NULL;
#else
  mesh->voy   = (void*)0;
#endif
  mesh->np = mesh->nt = mesh->nq = mesh->ne = 0;
  mesh->ntet = mesh->nhex = mesh->nbb = 0;

  if ( animdep == range ) {
    sprintf(data,".%d",range);
    ptr = (char *)strstr(mesh->name,data);
    if ( ptr )  *ptr = '\0';
    strcpy(base,mesh->name);
  }

  /* adjust file name */
  sprintf(mesh->name,"%s.%d",base,range);
  if ( ddebug )  fprintf(stdout,"load %s\n",mesh->name);

  /* read mesh */
  if ( quiet ) {  
    fprintf(stdout,"  Loading %s     \r",mesh->name);
    fflush(stdout);
  }
  mesh->typ = 0;
  ret = loadMesh(mesh);
  if ( ret < 1 ) {
    mesh->typ = 1;
    ret = inmsh2(mesh);
    if ( !ret ) {
      mesh->typ = 2;
      ret = loadGIS(mesh);
      if ( !ret )  return(0);
    }
  }
  if ( (mesh->ntet && !mesh->nt) || (mesh->nhex && !mesh->nq) )  
    meshSurf(mesh);
  meshBox(mesh,animdep==range);
  if ( !quiet )  meshInfo(mesh);

  /* read metric */
  if ( !loadSol(mesh,mesh->name,1) )
    bbfile(mesh);
  if ( !quiet && mesh->nbb )
    fprintf(stdout,"    Solutions  %8d\n",mesh->nbb);

  return(1);
}


int loadNextMesh(pMesh mesh,int k,int parse) {
  pScene    sc;
  int       is;

  if ( !getmesh(mesh,k) )  return(0);
  if ( ddebug ) printf("loadNextMesh: create %d window(s)\n",cv.nbs);

  /* compute mesh box */
  is = currentScene();
  if ( !cv.scene[is] ) {
    cv.scene[is] = (pScene)M_calloc(1,sizeof(Scene),"loadNextMesh.scene");
    if ( !cv.scene[is] )  return(0);
  }
  sc = cv.scene[is];
  if ( parse ) {
    parsop(sc,mesh);
  }
  setupPalette(sc,mesh);
  meshRef(sc,mesh);
  matSort(sc);

  return(1);
}


int playAnim(pScene sc,pMesh mesh,int deb,int fin) {
  int     k;
  char   *ptr,data[256],base[256];

  /* get basename */
  sprintf(data,".%d",deb);
  ptr = (char *)strstr(mesh->name,data);
  if ( ptr )  *ptr = '\0';
  strcpy(base,mesh->name);

  if ( saveimg ) {
    glDrawBuffer(GL_BACK_LEFT);
  }

  /* start animation */
  for (k=deb; k<=fin; k++) {
    strcpy(mesh->name,base);

    resetLists(sc,mesh);
    if ( !loadNextMesh(mesh,k,0) )  return(0);
    doLists(sc,mesh);
    sc->glist = geomList(sc,mesh);
    if ( sc->mode & S_MAP ) doMapLists(sc,mesh,1);
    if ( sc->isotyp )       doIsoLists(sc,mesh,1);

    if ( !saveimg ) {
      sprintf(data,"Medit - [%s] #%d",mesh->name,sc->idwin);
      glutSetWindowTitle(data);
    }

    glClearColor(sc->par.back[0],sc->par.back[1],sc->par.back[2],sc->par.back[3]);
    if ( sc->type & S_SCISSOR )
      scissorScene();
    else
      redrawScene();
  }
  fprintf(stdout,"\n  Done.\n");
  strcpy(mesh->name,base);

  if ( saveimg ) {
    glDrawBuffer(GL_FRONT | GL_BACK);
    if ( saveimg ) {
      sprintf(data,"Medit - [%s] #%d",mesh->name,sc->idwin);
      glutSetWindowTitle(data);
    }
    glutPostRedisplay();
    saveimg = 0;
  }

  return(1);
}


int animParticle(pScene sc,pMesh mesh) {
  int     cur;
  char   *ptr,data[256],base[256];

  /* get basename */
  strcpy(base,mesh->name);
  ptr = (char *)strstr(base,".");
  sscanf(ptr,".%d",&cur);
  cur++;
  if ( cur > animfin )  return(0);
  if ( ptr )  ptr[0] = '\0';
  sprintf(data,"%s.%.3d",base,cur);
  strcpy(mesh->name,base);

  resetLists(sc,mesh);
  if ( !loadNextMesh(mesh,cur,1) )  return(0);

  doLists(sc,mesh);
  sc->glist = geomList(sc,mesh);
  if ( sc->mode & S_MAP ) doMapLists(sc,mesh,1);
  if ( sc->isotyp )       doIsoLists(sc,mesh,2);

  if ( !advectParticle(sc,mesh) )  return(0);

  return(1);
}


int animat() {
  pMesh    mesh;
  pScene   sc;
  char     data[128],*name;

  /* default */
  if ( ddebug )  printf("animat: read file(s)\n");

  /* enter basename */
  if ( !cv.nbm ) {
    fprintf(stdout,"  File name(s) missing. Please enter : ");
    fflush(stdout); fflush(stdin);
    fgets(data,120,stdin);
    if ( !strlen(data) ) {
      fprintf(stdout,"  ## Error\n");
      return(0);
    }
    fprintf(stdout,"  Enter range [start,end] :");
    fflush(stdout); fflush(stdin);
    fscanf(stdin,"%d %d",&animdep,&animfin);

    /* parse file name(s) */
    name = strtok(data," \n");
    while( name ) {
      if ( !cv.mesh[cv.nbm] ) {
        cv.mesh[cv.nbm] = (pMesh)M_calloc(1,sizeof(Mesh),"animat.mesh");
        if ( !cv.mesh[cv.nbm] )  return(0);
      }
      /*(cv.mesh[cv.nbm])->name = calloc(strlen(name)+15,sizeof(char));*/
      strcpy(cv.mesh[cv.nbm]->name,name);
      break;
    }
    if ( !cv.nbm ) return(0);
  }
  cv.nbs = cv.nbm;

  /* read initial mesh */
  mesh = cv.mesh[0];
  if ( !loadNextMesh(mesh,animdep,1) )  return(0);

  /* create grafix */
  sc = cv.scene[0];
  if ( !createScene(sc,0) ) {
    fprintf(stderr,"  ## Unable to create scene\n");
    return(0);
  }

  quiet = 1;
  return(1);
}

