// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "procrustes.h"
#include "polar_svd.h"
#include "polar_dec.h"
 
template <typename DerivedV, typename Scalar, int DIM, int TType>
IGL_INLINE void igl::procrustes(
    const Eigen::PlainObjectBase<DerivedV>& X,
    const Eigen::PlainObjectBase<DerivedV>& Y,
    bool includeScaling,
    bool includeReflections,
    Eigen::Transform<Scalar,DIM,TType>& T)
{
   using namespace Eigen;
   assert (X.rows() == Y.rows() && "Same number of points");
   assert(X.cols() == Y.cols() && "Points have same dimensions");

   // Center data
   const VectorXd Xmean = X.colwise().mean();      
   const VectorXd Ymean = Y.colwise().mean();      
   MatrixXd XC = X.rowwise() - Xmean.transpose();
   MatrixXd YC = Y.rowwise() - Ymean.transpose();


   // Determine scale
   double scale = 1.;
   if (includeScaling)
   {
       double scaleX = XC.norm() / XC.rows();
       double scaleY = YC.norm() / YC.rows();
       scale = scaleY/scaleX;
       XC *= scale;
       assert (abs(XC.norm() / XC.rows() - scaleY) < 1e-8);
   }


   // Rotation 
   MatrixXd S = XC.transpose() * YC; 
   Matrix<Scalar,DIM,DIM> Rm,Tm;
   if (includeReflections)
     polar_dec(S,Rm,Tm);
   else
     polar_svd(S,Rm,Tm);

   // Translation
   Matrix<Scalar,DIM,1> t = Ymean - scale*Rm*Xmean;

   // Combine
   T = Translation<Scalar,DIM>(t) * Rm.transpose() * Scaling(scale);
}