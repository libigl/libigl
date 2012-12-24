#include "svd.h"

#include <cstdlib>
#include <Accelerate/Accelerate.h>
#include <cstdio>

bool igl::svd3x3(double * a, double * u, double * s, double * vt)
{
  /* Locals */
  int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;
  double wkopt;
  double* work;
  /* Local arrays */
  /* iwork dimension should be at least 8*min(m,n) */
  int iwork[8*3];
  //double s[3], u[3*3], vt[3*3];
  //double a[3*3] = {8,3,4,1,5,9,6,7,2};
  /* Query and allocate the optimal workspace */
  lwork = -1;
  dgesdd_( 
    "Singular vectors",
    &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Compute SVD */
  dgesdd_(
    "Singular vectors",
    &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );
  /* Check for convergence */
  if( info > 0 )
  {
    printf("The algorithm computing SVD failed to converge.\n" );
    return false;
  }
  /* Free workspace */
  free( (void*)work );
  return true;
}
