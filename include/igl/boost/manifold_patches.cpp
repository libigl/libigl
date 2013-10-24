#include "manifold_patches.h"
#include "components.h"
#include <igl/sort.h>
#include <igl/unique.h>
#include <vector>

template <typename DerivedF, typename DerivedC, typename AScalar>
void igl::manifold_patches(
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::SparseMatrix<AScalar> & A)
{
  using namespace Eigen;
  using namespace std;

  // simplex size
  assert(F.cols() == 3);

  // List of all "half"-edges: 3*#F by 2
  Matrix<typename DerivedF::Scalar, Dynamic, 2> allE,sortallE,uE;
  VectorXi IC,_;
  allE.block(0*F.rows(),0,F.rows(),0) = F.col(1);
  allE.block(0*F.rows(),0,F.rows(),1) = F.col(2);
  allE.block(1*F.rows(),0,F.rows(),0) = F.col(2);
  allE.block(1*F.rows(),0,F.rows(),1) = F.col(0);
  allE.block(2*F.rows(),0,F.rows(),0) = F.col(0);
  allE.block(2*F.rows(),0,F.rows(),1) = F.col(1);
  // Sort each row
  sort(allE,2,true,sortallE,_);
  //IC(i) tells us where to find sortallE(i,:) in uE: 
  // so that sortallE(i,:) = uE(IC(i),:)
  unique_rows(sortallE,uE,_,IC);
  // uE2F(e,f) = 1 means face f is adjacent to unique edge e
  vector<Triplet<AScalar> > uE2Fijv(IC.rows());
  for(int e = 0;e<IC.rows();e++)
  {
    uE2Fijv[e] = Triplet<AScalar>(IC(e),e%F.rows(),1);
  }
  SparseMatrix<AScalar> uE2F(uE.rows(),F.rows());
  uE2F.setFromTriplets(uE2Fijv.begin(),uE2Fijv.end());
  // kill non-manifold edges
  for(int j=0; j<uE2F.outerSize();j++)
  {
    // Iterate over inside
    for(typename SparseMatrix<AScalar>::InnerIterator it (uE2F,j); it; ++it)
    {
      if(it.value() > 2)
      {
        uE2F.coeffRef(it.row(),it.col()) = 0;
      }
    }
  }
  // Face-face Adjacency matrix
  SparseMatrix<AScalar> uE2FT;
  uE2FT = uE2F.transpose().eval();
  A = uE2FT*uE2F;
  // All ones
  for(int j=0; j<A.outerSize();j++)
  {
    // Iterate over inside
    for(typename SparseMatrix<AScalar>::InnerIterator it (A,j); it; ++it)
    {
      if(it.value() > 1)
      {
        A.coeffRef(it.row(),it.col()) = 1;
      }
    }
  }
  //% Connected components are patches
  //%C = components(A); % alternative to graphconncomp from matlab_bgl
  //[~,C] = graphconncomp(A);
  // graph connected components using boost
  components(A,C);

}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::manifold_patches<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::SparseMatrix<int, 0, int>&);
#endif
