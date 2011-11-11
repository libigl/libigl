//
//  moveFV.h
//  Preview3D
//
//  Created by Olga Diamanti on 11/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef Preview3D_moveFV_h
#define Preview3D_moveFV_h

namespace igl 
{
  // moveFV 
  // Move a scalar field defined on faces to vertices by averaging
  //
  // Input:
  // V,F: mesh
  // S: scalar field defined on faces, Fx1
  // 
  // Output:
  // SV: scalar field defined on vertices
  void moveFV(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &S, Eigen::MatrixXd &SV);
}

inline void igl::moveFV(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &S, Eigen::MatrixXd &SV)
{
  
  SV = Eigen::MatrixXd::Zero(V.rows(),S.cols());
  Eigen::VectorXd COUNT = Eigen::VectorXd::Zero(V.rows());
  for (int i = 0; i <F.rows(); ++i)
  {
    for (int j = 0; j<F.cols(); ++j)
    {
      SV.row(F(i,j)) += S.row(i);
      COUNT[F(i,j)] ++;
    }
  }
  for (int i = 0; i <V.rows(); ++i)
    SV.row(i) /= COUNT[i];
  
};

#endif
