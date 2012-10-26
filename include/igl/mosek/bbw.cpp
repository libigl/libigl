#define VERBOSE
#include "bbw.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/speye.h>
#include <igl/slice_into.h>
#include <igl/normalize_row_sums.h>
#include <igl/verbose.h>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

#include <iostream>
#include <cstdio>


igl::BBWData::BBWData():
  partition_unity(false)
{}

void igl::BBWData::print()
{
  using namespace std;
  cout<<"partition_unity: "<<partition_unity<<endl;
  cout<<"W0=["<<endl<<W0<<endl<<"];"<<endl;
}


template <
  typename DerivedV, 
  typename DerivedEle, 
  typename Derivedb,
  typename Derivedbc, 
  typename DerivedW>
IGL_INLINE bool igl::bbw(
  const Eigen::MatrixBase<DerivedV> & V, 
  const Eigen::MatrixBase<DerivedEle> & Ele, 
  const Eigen::MatrixBase<Derivedb> & b, 
  const Eigen::MatrixBase<Derivedbc> & bc, 
  igl::BBWData & data,
  Eigen::MatrixBase<DerivedW> & W
  )
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;

  // number of domain vertices
  int n = V.rows();
  // number of handles
  int m = bc.cols();

  SparseMatrix<typename DerivedW::Scalar> L;
  cotmatrix(V,Ele,L);
  MassMatrixType mmtype = MASSMATRIX_VORONOI;
  if(Ele.cols() == 4)
  {
    mmtype = MASSMATRIX_BARYCENTRIC;
  }
  SparseMatrix<typename DerivedW::Scalar> M;
  SparseMatrix<typename DerivedW::Scalar> Mi;
  massmatrix(V,Ele,mmtype,M);

  invert_diag(M,Mi);

  // Biharmonic operator
  SparseMatrix<typename DerivedW::Scalar> Q = L.transpose() * Mi * L;

  W.derived().resize(n,m);
  if(data.partition_unity)
  {
    // Not yet implemented
    assert(false);
  }else
  {
    // No linear terms
    VectorXd c = VectorXd::Zero(n);
    // No linear constraints
    SparseMatrix<typename DerivedW::Scalar> A(0,n);
    VectorXd uc(0,1);
    VectorXd lc(0,1);
    // Upper and lower box constraints (Constant bounds)
    VectorXd ux = VectorXd::Ones(n);
    VectorXd lx = VectorXd::Zero(n);
    // Loop over handles
    for(int i = 0;i<m;i++)
    {
      verbose("\n^%s: Computing weight for handle %d out of %d.\n\n",
        __FUNCTION__,i+1,m);
      // impose boundary conditions
      VectorXd bci = bc.col(i);
      slice_into(bci,b,ux);
      slice_into(bci,b,lx);
      VectorXd Wi;
      bool r = igl::mosek_quadprog(Q,c,0,A,lc,uc,lx,ux,data.mosek_data,Wi);
      if(!r)
      {
        return false;
      }
      W.col(i) = Wi;
    }
    // Need to normalize
    igl::normalize_row_sums(W,W); 
  }

  return true;
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template bool igl::bbw<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::BBWData&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
