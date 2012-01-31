//
//  moveFV.h
//  Preview3D
//
//  Created by Olga Diamanti on 11/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef IGL_MOVEFV_H
#define IGL_MOVEFV_H

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
  template <typename T>
  inline void moveFV(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
              const Eigen::MatrixXi &F,
              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &S,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &SV);
}

// Implementation

template <typename T>
inline void igl::moveFV(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V,
            const Eigen::MatrixXi &F,
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &S,
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &SV)
{
  
  SV = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(V.rows(),S.cols());
  Eigen::Matrix<T, Eigen::Dynamic, 1> COUNT = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(V.rows());
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
