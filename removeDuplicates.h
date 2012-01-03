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
  
  template <typename T>
  inline void removeDuplicates(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V, const Eigen::MatrixXi &F, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I, const double epsilon = 2.2204e-15);
  
}

// Implementation
template <typename T>
inline void igl::removeDuplicates(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V, const Eigen::MatrixXi &F, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &NV, Eigen::MatrixXi &NF, Eigen::VectorXi &I, const double epsilon)
{
  //// build collapse map
  int n = V.rows();
  
  I = Eigen::VectorXi (n);
  I[0] = 0;
  
  bool *VISITED = new bool[n];
  for (int i =0; i <n; ++i)
    VISITED[i] = false;
  
  NV = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(n,V.cols());
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
  std::vector<int> face;
  NF = Eigen::MatrixXi(F.rows(),F.cols());
  for (int i =0; i <F.rows(); ++i)
  {
    face.clear();
    for (int j = 0; j< F.cols(); ++j)
      if(std::find(face.begin(), face.end(), I[F(i,j)]) == face.end())
         face.push_back(I[F(i,j)]);
    if (face.size() == F.cols())
    {
      for (int j = 0; j< F.cols(); ++j)
        NF(count,j) = face[j];
      count ++;
    }
  }
  NF.conservativeResize	(	count , Eigen::NoChange );
  
  delete [] VISITED;
}






#endif
