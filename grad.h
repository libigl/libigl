//
//  grad.h
//  Preview3D
//
//  Created by Olga Diamanti on 11/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef Preview3D_grad_h
#define Preview3D_grad_h

#include <Eigen/Core>

namespace igl {
  // GRAD
  // G = grad(V,F,X)
  //
  // Compute the numerical gradient at every face of a triangle mesh.
  //
  // Inputs:
  //   V  #vertices by 3 list of mesh vertex positions
  //   F  #faces by 3 list of mesh face indices
  //   X  # vertices list of scalar function values
  // Outputs:
  //   G  #faces by 3 list of gradient values
  //
  
  // Gradient of a scalar function defined on piecewise linear elements (mesh)
  // is constant on each triangle i,j,k:
  // grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  // where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  // i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
  // 90 degrees
  //
  void grad(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &X, Eigen::MatrixXd &G );
}

inline void igl::grad(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &X, Eigen::MatrixXd &G)
{
  G = Eigen::MatrixXd::Zero(F.rows(),3);
  for (int i = 0; i<F.rows(); ++i)
  {
    // renaming indices of vertices of triangles for convenience
    int i1 = F(i,0);
    int i2 = F(i,1);
    int i3 = F(i,2);
    
    // #F x 3 matrices of triangle edge vectors, named after opposite vertices
    Eigen::RowVector3d v32 = V.row(i3) - V.row(i2);
    Eigen::RowVector3d v13 = V.row(i1) - V.row(i3);
    Eigen::RowVector3d v21 = V.row(i2) - V.row(i1);
    
    // area of parallelogram is twice area of triangle
    // area of parallelogram is || v1 x v2 || 
    Eigen::RowVector3d n  = v32.cross(v13); 
    
    // This does correct l2 norm of rows, so that it contains #F list of twice
    // triangle areas
    double dblA = std::sqrt(n.dot(n));
    
    // now normalize normals to get unit normals
    Eigen::RowVector3d u = n / dblA;
    
    // rotate each vector 90 degrees around normal
    double norm21 = std::sqrt(v21.dot(v21));
    double norm13 = std::sqrt(v13.dot(v13));
    Eigen::RowVector3d eperp21 = u.cross(v21);
    eperp21 = eperp21 / std::sqrt(eperp21.dot(eperp21));
    eperp21 *= norm21;
    Eigen::RowVector3d eperp13 = u.cross(v13);
    eperp13 = eperp13 / std::sqrt(eperp13.dot(eperp13));
    eperp13 *= norm13;
    
    G.row(i) = ((X[i2] -X[i1]) *eperp13 + (X[i3] - X[i1]) *eperp21) / dblA;
  };
}
  
  
#endif
