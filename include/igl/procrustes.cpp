// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Stefan Brugger <stefanbrugger@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "procrustes.h"
 
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
   VectorXd Xmean = X.colwise().mean();      
   VectorXd Ymean = Y.colwise().mean();      
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
   MatrixXd M = XC.transpose() * YC; 
   JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

   MatrixXd sigma;
   sigma.setIdentity(DIM,DIM);
   if (!includeReflections && (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0)
       sigma(DIM-1,DIM-1) = -1.;

   Matrix<Scalar,DIM,DIM> R = svd.matrixU() * sigma * svd.matrixV().transpose();
   assert(abs(R.determinant() - 1) < 1e-10);


   // Translation
   Matrix<Scalar,DIM,1> t = Ymean - scale*R*Xmean;
     

   // Combine
   T = Translation<Scalar,DIM>(t) * R.transpose() * Scaling(scale);
}