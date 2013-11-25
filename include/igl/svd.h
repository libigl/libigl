// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//THIS DEPENDS ON BLAS. WHY? WHY NOT EIGEN?

//#ifndef IGL_SVD_H
//#define IGL_SVD_H
//namespace igl
//{
//  // Compute 3x3 SVD using lapack's dgesdd_ function
//  //
//  // Input:
//  //   a  pointer to 3x3 matrix in COLUMN MAJOR order
//  // Outputs:
//  //   u  pointer to 3x3 matrix in COLUMN MAJOR order
//  //   s  pointer to 3 vector 
//  //   vt  pointer to 3x3 matrix in COLUMN MAJOR order, or think of this as v in
//  //     row-major order
//  // Returns true on success, false on failure
//  // 
//  // Known bugs: This only compiles on Mac and depends on Lapack rather than eigen
//  bool svd3x3(double * a, double * u, double * s, double * vt);
//};
//
//#ifdef IGL_HEADER_ONLY
//#include "svd.cpp"
//#endif
//
//#endif
