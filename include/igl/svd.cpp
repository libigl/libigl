// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//// This only works on MAC ...
//#ifdef __APPLE__
//#include "svd.h"
//#include <Accelerate/Accelerate.h>
//#include <cstdlib>
//#include <cstdio>
//
//bool igl::svd3x3(double * a, double * u, double * s, double * vt)
//{
//  /* Locals */
//  int m = 3, n = 3, lda = 3, ldu = 3, ldvt = 3, info, lwork;
//  double wkopt;
//  double* work;
//  /* Local arrays */
//  /* iwork dimension should be at least 8*min(m,n) */
//  int iwork[8*3];
//  //double s[3], u[3*3], vt[3*3];
//  //double a[3*3] = {8,3,4,1,5,9,6,7,2};
//  /* Query and allocate the optimal workspace */
//  lwork = -1;
//  dgesdd_( 
//    "Singular vectors",
//    &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, iwork, &info);
//  lwork = (int)wkopt;
//  work = (double*)malloc( lwork*sizeof(double) );
//  /* Compute SVD */
//  dgesdd_(
//    "Singular vectors",
//    &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info );
//  /* Check for convergence */
//  if( info > 0 )
//  {
//    printf("The algorithm computing SVD failed to converge.\n" );
//    return false;
//  }
//  /* Free workspace */
//  free( (void*)work );
//  return true;
//}
//#endif
