#ifndef IGL_COTMATRIX_H
#define IGL_COTMATRIX_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

// History:
//  Used const references rather than copying the entire mesh 
//    Alec 9 October 2011
//  removed cotan (uniform weights) optional parameter it was building a buggy
//    half of the uniform laplacian, please see adjacency_matrix istead 
//    Alec 9 October 2011

namespace igl 
{
  // Constructs the cotangent stiffness matrix (discrete laplacian) for a given
  // mesh (V,F).
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by 3 list of mesh faces (must be triangles)
  // Outputs: 
  //   L  #V by #V cotangent matrix, each row i corresponding to V(i,:)
  //
  // See also: adjacency_matrix
  inline void cotmatrix(
    const Eigen::MatrixXd & V, 
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double>& L);
  // Helper function that computes the cotangent weights for each corner of a
  // given triangle
  // Inputs:
  //   v1  position of corner #1 of triangle
  //   v2  position of corner #2 of triangle
  //   v3  position of corner #3 of triangle
  // Outputs:
  //   cot1  cotangent of angle at corner #1
  //   cot2 cotangent of angle at corner #2
  //   cot3  cotangent of angle at corner #3
  //   
  inline void computeCotWeights(
    const Eigen::Vector3d& v1, 
    const Eigen::Vector3d& v2, 
    const Eigen::Vector3d& v3, 
    double& cot1, 
    double& cot2, 
    double& cot3);
}

// Implementation

// For error printing
#include <cstdio>

inline void igl::computeCotWeights(
  const Eigen::Vector3d& v1, 
  const Eigen::Vector3d& v2, 
  const Eigen::Vector3d& v3, 
  double& cot1, 
  double& cot2, 
  double& cot3)
{
  Eigen::Vector3d v12 = v2-v1;
  Eigen::Vector3d v13 = v3-v1;
  Eigen::Vector3d v23 = v3-v2;
  
  double halfArea = (v12.cross(v13)).norm();//squaredNorm();
  
  //improve numerical stability
  const double cotTolerance = 1e-10;
  if(halfArea < cotTolerance)
  {
    fprintf(stderr,"Cot weights are close to singular!\n");
    halfArea = cotTolerance;
  }
  
  cot1 = (v12.dot(v13)) / halfArea /2;
  cot2 =  (v23.dot(-v12)) / halfArea /2;
  cot3 = (-v23.dot(-v13)) / halfArea /2;
}

inline void igl::cotmatrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  Eigen::SparseMatrix<double>& L)
{
  // Assumes vertices are 3D or 2D
  assert((V.cols() == 3) || (V.cols() == 2));
  int dim = V.cols();
  // Assumes faces are triangles
  assert(F.cols() == 3);

  Eigen::DynamicSparseMatrix<double, Eigen::RowMajor> dyn_L (V.rows(), V.rows());
  // This is important! it could decrease the comptuation time by a factor of 2
  // Laplacian for a closed 2d manifold mesh will have on average 7 entries per
  // row
  dyn_L.reserve(7*V.rows());
  
  // Loop over triangles
  for (unsigned i = 0; i < F.rows(); i++)
  {
    // Corner indices of this triangle
    int vi1 = F(i,0);
    int vi2 = F(i,1);
    int vi3 = F(i,2);
    // Grab corner positions of this triangle
    Eigen::Vector3d v1(V(vi1,0), V(vi1,1), (dim==2?0:V(vi1,2)));
    Eigen::Vector3d v2(V(vi2,0), V(vi2,1), (dim==2?0:V(vi2,2)));
    Eigen::Vector3d v3(V(vi3,0), V(vi3,1), (dim==2?0:V(vi3,2)));
    // Compute cotangent of angles at each corner
    double cot1, cot2, cot3;
    computeCotWeights(v1, v2, v3, cot1, cot2, cot3);
    // Throw each corner's cotangent at opposite edge (in both directions)
    dyn_L.coeffRef(vi1, vi2) += cot3;
    dyn_L.coeffRef(vi2, vi1) += cot3;
    dyn_L.coeffRef(vi2, vi3) += cot1;
    dyn_L.coeffRef(vi3, vi2) += cot1;
    dyn_L.coeffRef(vi3, vi1) += cot2;
    dyn_L.coeffRef(vi1, vi3) += cot2;
  }

  for (int k=0; k < dyn_L.outerSize(); ++k)
  {
    double tmp = 0.0f;
    for(Eigen::DynamicSparseMatrix<double, Eigen::RowMajor>::InnerIterator it (dyn_L, k); it; ++it)
    {
      tmp += it.value();
    }
    dyn_L.coeffRef(k,k) = -tmp;
  }
  L = Eigen::SparseMatrix<double>(dyn_L);
}

#endif
