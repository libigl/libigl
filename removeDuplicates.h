//
//  removeDuplicates.h
//  Preview3D
//
//  Created by Olga Diamanti on 17/11/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef RemoveDuplicates_h
#define RemoveDuplicates_h

#include <Eigen/Core>
namespace igl 
{
  // [ NV, NF ] = removeDuplicates( V,F,epsilon )
  // Merge the duplicate vertices from V, fixing the topology accordingly
  //
  // Input:
  // V,F: mesh description
  // epsilon: minimal distance to consider two vertices identical
  //
  // Output:
  // NV, NF: new mesh without duplicate vertices
  
  void removeDuplicates(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I, const double epsilon);
  
}

// Implementation
inline void igl::removeDuplicates(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I, const double epsilon = 2.2204e-15)
{
  assert (V.cols() == 3);
  assert (F.cols() == 3);
  
  //// build collapse map
  int n = V.rows();
  
  I = Eigen::VectorXi (n);
  I[0] = 0;
  
  bool *VISITED = new bool[n];
  for (int i =0; i <n; ++i)
    VISITED[i] = false;
  
  NV = Eigen::MatrixXd(n,3);
  int count = 0;
  Eigen::VectorXd d(n);
  for (int i =0; i <n; ++i)
  {
    if(!VISITED[i])
    {
      NV.row(count) = V.row(i);
      I[i] = count;
      VISITED[i] = true;
      for (int j = i+1; j <n; ++j)
      {
        if((V.row(j) - V.row(i)).norm() < epsilon)
        {
          VISITED[j] = true;
          I[j] = count;
        }
      }
      count ++;
    }
  }
  
  NV.conservativeResize	(	count , Eigen::NoChange );

  count = 0;
  NF = Eigen::MatrixXi(F.rows(),3);
  for (int i =0; i <F.rows(); ++i)
  {
    int v0 = I[F(i,0)];
    int v1 = I[F(i,1)];
    int v2 = I[F(i,2)];
    if ( (v0 != v1) && (v2 != v1) && (v0 != v2) )
      NF.row(count++) << v0, v1, v2;
  }
  NF.conservativeResize	(	count , Eigen::NoChange );
  
  delete [] VISITED;
}






#endif
