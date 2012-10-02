#include "medit.h"
#include "extern.h"
#include "sproto.h"


int zaldy2(pMesh mesh) {
  pSolution  ps;
  int        k,nbf;
  
  /* memory alloc. */
  mesh->sol = (pSolution)M_calloc(mesh->nbb+1,sizeof(struct solu),"zaldy2");
  assert(mesh->sol);

  if ( mesh->nfield == 1 )  return(1);
  if ( mesh->nfield == mesh->dim )
    nbf = mesh->dim;      /* vector field */
  else 
    nbf = mesh->dim*(mesh->dim+1)/2;  /* d*d matrix */        

  for (k=1; k<=mesh->nbb; k++) {
    ps = &mesh->sol[k];
    ps->m = (float*)malloc(nbf*sizeof(float));
  }

  return(1);
}
