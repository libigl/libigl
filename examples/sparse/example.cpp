#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_IM_MAD_AS_HELL_AND_IM_NOT_GOING_TO_TAKE_IT_ANYMORE
#include <Eigen/Sparse>
#include <Eigen/SparseExtra>
using namespace Eigen;

#include <cstdio>
#include <iostream>
using namespace std;

#include "readOBJ.h"
#include "print_ijv.h"
#include "cotmatrix.h"
#include "get_seconds.h"
#include "EPS.h"
using namespace igl;

#include "SparseSolver.h"

// Build cotangent matrix by looping over triangles filling in eigen's dynamic
// sparse matrix then casting to sparse matrix
inline void cotmatrix_dyncast(
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
  
  // Loop over triangles
  for (unsigned i = 0; i < F.rows(); i++)
  {
    // Corner indices of this triangle
    int vi1 = F(i,0);
    int vi2 = F(i,1);
    int vi3 = F(i,2);
    // Grab corner positions of this triangle
    Eigen::Vector3d v1(V(vi1,0), V(vi1,1), V(vi1,2));
    Eigen::Vector3d v2(V(vi2,0), V(vi2,1), V(vi2,2));
    Eigen::Vector3d v3(V(vi3,0), V(vi3,1), V(vi3,2));
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
    //dyn_L.coeffRef(vi1, vi1) -= cot3;
    //dyn_L.coeffRef(vi2, vi2) -= cot3;
    //dyn_L.coeffRef(vi2, vi2) -= cot1;
    //dyn_L.coeffRef(vi3, vi3) -= cot1;
    //dyn_L.coeffRef(vi3, vi3) -= cot2;
    //dyn_L.coeffRef(vi1, vi1) -= cot2;
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

inline bool safe_AddToIJV(SparseSolver & S, int i, int j, double v)
{
  // Only add if not within double precision of zero
  if( v>DOUBLE_EPS|| v< -DOUBLE_EPS)
  {
    S.AddToIJV(i,j,v);
    return true;
  }
  //else
  return false;
}

inline void cotmatrix_sparsesolver(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  Eigen::SparseMatrix<double>& L)
{
  SparseSolver ssL(V.rows(),V.rows());

  // loop over triangles
  // Loop over triangles
  for (int i = 0; i < F.rows(); i++)
  {
    // Corner indices of this triangle
    int vi1 = F(i,0);
    int vi2 = F(i,1);
    int vi3 = F(i,2);
    // Grab corner positions of this triangle
    Eigen::Vector3d v1(V(vi1,0), V(vi1,1), V(vi1,2));
    Eigen::Vector3d v2(V(vi2,0), V(vi2,1), V(vi2,2));
    Eigen::Vector3d v3(V(vi3,0), V(vi3,1), V(vi3,2));
    // Compute cotangent of angles at each corner
    double cot1, cot2, cot3;
    computeCotWeights(v1, v2, v3, cot1, cot2, cot3);
    // Throw each corner's cotangent at opposite edge (in both directions)
    safe_AddToIJV(ssL,vi1,vi2,cot3);
    safe_AddToIJV(ssL,vi2,vi1,cot3);
    safe_AddToIJV(ssL,vi2,vi3,cot1);
    safe_AddToIJV(ssL,vi3,vi2,cot1);
    safe_AddToIJV(ssL,vi3,vi1,cot2);
    safe_AddToIJV(ssL,vi1,vi3,cot2);

    safe_AddToIJV(ssL,vi1,vi1,-cot2);
    safe_AddToIJV(ssL,vi1,vi1,-cot3);
    safe_AddToIJV(ssL,vi2,vi2,-cot3);
    safe_AddToIJV(ssL,vi2,vi2,-cot1);
    safe_AddToIJV(ssL,vi3,vi3,-cot1);
    safe_AddToIJV(ssL,vi3,vi3,-cot2);
  }
}

inline void cotmatrix_adjacencylist(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  Eigen::SparseMatrix<double>& L)
{
  // adjacency list
  std::vector<std::vector<int> > A;
  //adjacency_list(F,A);
}

int main(int argc, char * argv[])
{
  if(argc < 2)
  {
    printf("Usage:\n  ./example mesh.obj\n");
    return 1;
  }
  // Read in a triangle mesh
  MatrixXd V;
  MatrixXi F;
  readOBJ(argv[1],V,F);
  // Should be 3D triangle mesh
  assert(V.cols() == 3);
  assert(F.cols() == 3);
  // Print info about the mesh
  printf("#vertices: %d\n",(int)V.rows());
  printf("#faces: %d\n",(int)F.rows());
  printf("min face index: %d\n",(int)F.minCoeff());
  printf("max face index: %d\n",(int)F.maxCoeff());

  double before;
  int trials;
  int max_trials = 5;
  SparseMatrix<double> L;

  printf("cotmatrix_dyncast:\n  ");
  before = get_seconds();
  for(trials = 0;trials<max_trials;trials++)
  {
    cotmatrix_dyncast(V,F,L);
  }
  printf("%g\n",(get_seconds()-before)/(double)trials);

  printf("cotmatrix_sparsesolver:\n  ");
  before = get_seconds();
  for(trials = 0;trials<max_trials;trials++)
  {
    cotmatrix_sparsesolver(V,F,L);
  }
  printf("%g\n",(get_seconds()-before)/(double)trials);

  return 0;
}
