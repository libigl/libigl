#ifndef IGL_CROUZEIX_RAVIART_COTMATRIX
#define IGL_CROUZEIX_RAVIART_COTMATRIX
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // CROUZEIX_RAVIART_COTMATRIX Compute the Crouzeix-Raviart cotangent
  // stiffness matrix.
  //
  // See for example "Discrete Quadratic Curvature Energies" [Wardetzky, Bergou,
  //
  // Harmon, Zorin, Grinspun 2007]
  // Inputs:
  //   V  #V by dim list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Outputs:
  //   L  #E by #E edge-based diagonal cotangent matrix
  //   E  #E by 2 list of edges
  //   EMAP  #F*3 list of indices mapping allE to E
  //
  template <typename DerivedV, typename DerivedF, typename LT, typename DerivedE, typename DerivedEMAP>
  void crouzeix_raviart_cotmatrix(
      const Eigen::MatrixBase<DerivedV> & V, 
      const Eigen::MatrixBase<DerivedF> & F, 
      Eigen::SparseMatrix<LT> & L,
      Eigen::PlainObjectBase<DerivedE> & E,
      Eigen::PlainObjectBase<DerivedEMAP> & EMAP);
  // wrapper if E and EMAP are already computed (better match!)
  template <typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedEMAP, typename LT>
  void crouzeix_raviart_cotmatrix(
      const Eigen::MatrixBase<DerivedV> & V, 
      const Eigen::MatrixBase<DerivedF> & F, 
      const Eigen::MatrixBase<DerivedE> & E,
      const Eigen::MatrixBase<DerivedEMAP> & EMAP,
      Eigen::SparseMatrix<LT> & L);
}
#ifndef IGL_STATIC_LIBRARY
#  include "crouzeix_raviart_cotmatrix.cpp"
#endif
#endif
