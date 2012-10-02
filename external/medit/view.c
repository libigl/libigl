#include "medit.h"
#include "sproto.h"
#include "extern.h"

static pTransform  cview   = 0;
static pCamera     ccam    = 0;
static pPersp      cpersp  = 0;
static int         cscene  = 0;
static ubyte       curview = 0;


void copyView(pTransform view,pCamera cam,pPersp persp) {
  cscene = currentScene();

  if ( !cview ) {
    cview = (pTransform)calloc(1,sizeof(struct transform));
    if ( !cview )  exit(2);
  }
  cview = (pTransform)memcpy(cview,view,sizeof(struct transform));
  if ( !cview )  exit(2);
  
  if ( !ccam ) {
    ccam = (pCamera)calloc(1,sizeof(struct camera));
    if ( !ccam )  exit(2);
  }
  ccam = (pCamera)memcpy(ccam,cam,sizeof(struct camera));
  if ( !ccam )  exit(2);
  
  if ( !cpersp ) {
    cpersp = (pPersp)calloc(1,sizeof(struct sperspective));
    if ( !cpersp )  exit(2);
  }
  cpersp = (pPersp)memcpy(cpersp,persp,sizeof(struct sperspective));
  if ( !cpersp )  exit(2);

  curview = 1;
}

int pasteView(pTransform view,pCamera cam,pPersp persp) {
  if ( !curview && !ccam && !persp )   return(0);
 
  view  = (pTransform)memcpy(view,cview,sizeof(struct transform));
  cam   = (pCamera)memcpy(cam,ccam,sizeof(struct camera));
  persp = memcpy(persp,cpersp,sizeof(struct sperspective));
  
  if ( !view || !cam || !persp )  exit(2);
  curview = 0;

  return(1);
}

/* link scene 1 (slave) to scene 2 (master) */
int linkView(pScene slave) {
  pScene master;
  int    idw = currentScene();

  /* default */
  if ( ddebug ) printf("link view\n");

  if ( !curview || idw == cscene )   return(0);

  /* link 2 scenes */
  master         = cv.scene[cscene];
  master->slave  = idw;
  master->master = -1;
  memcpy(slave->view,master->view,sizeof(struct transform));
  memcpy(slave->camera,master->camera,sizeof(struct camera));
  memcpy(slave->persp,master->persp,sizeof(struct sperspective));
  memcpy(slave->clip->eqn,master->clip->eqn,4*sizeof(float));
  memcpy(slave->clip->cliptr,master->clip->cliptr,sizeof(struct transform));
  slave->slave  = -1;
  slave->master = cscene;
  curview       = 0;

  return(1);
}

void unlinkView(pScene sc1) {
  pScene  sc2;

  if ( sc1->master == -1 )  return;

  /* copy master view */
  sc2 = cv.scene[sc1->master];
  sc1->view = (pTransform)createTransform();
  memcpy(sc1->view,sc2->view,sizeof(struct transform));

  /* remove link */
  sc1->master = -1;
  sc2->slave  = -1;
}
