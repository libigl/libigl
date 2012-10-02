#include "medit.h"
#include "extern.h"
#include "sproto.h"


void resetTransform(pTransform tr) {
  static float  itransf[16] = { 1., 0., 0., 0.,  0., 1., 0., 0.,
				0., 0., 1., 0.,  0., 0., 0., 1.};

  tr->pos[0] = tr->pos[1] = tr->pos[2] = 0.0f;
  tr->angle  = 0.0f;
  tr->panx   = tr->pany  = 0.0f;
  tr->opanx  = tr->opany = 0.0f;
  tr->mstate = 1;  
  tr->manim  = GL_FALSE;

  memcpy(tr->matrix,itransf,16*sizeof(float));
  memcpy(tr->rot,itransf,16*sizeof(float));
  memcpy(tr->tra,itransf,16*sizeof(float));
}

pTransform createTransform() {
  pTransform   tr;

  /* default */
  if ( ddebug) printf("create transformation\n");

  tr = (pTransform)M_calloc(1,sizeof(struct transform),"transform") ;
  assert(tr);

  /* set default values */
  resetTransform(tr);
  tr->mbutton = 0;

  return(tr);
}
