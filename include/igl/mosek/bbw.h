#ifndef BBW_H
#define BBW_H
#include "../igl_inline.h"

#include <Eigen/Dense>
#include "mosek_quadprog.h"

namespace igl
{
  // Container for BBW computation related data and flags
  class BBWData
  {
    public:
      // Enforce partition of unity during optimization (optimize all weight
      // simultaneously)
      bool partition_unity;
      // TODO: Is it safe if this is a reference?
      // Initial guess
      Eigen::MatrixXd W0;
      // TODO: Mosek options
      igl::MosekData mosek_data;
      // TODO: Active set options
    public:
      BBWData();
      // Print current state of object
      void print();
  };

  // Compute Bounded Biharmonic Weights on a given domain (V,Ele) with a given
  // set of boundary conditions
  //
  // Templates
  //   DerivedV  derived type of eigen matrix for V (e.g. MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. MatrixXi)
  //   Derivedb  derived type of eigen matrix for b (e.g. VectorXi)
  //   Derivedbc  derived type of eigen matrix for bc (e.g. MatrixXd)
  //   DerivedW  derived type of eigen matrix for W (e.g. MatrixXd)
  // Inputs:
  //   V  #V by dim vertex positions
  //   Ele  #Elements by simplex-size list of element indices
  //   b  #b boundary indices into V
  //   bc #b by #W list of boundary values
  //   data  object containing options, intial guess --> solution and results
  // Outputs:
  //   W  #V by #W list of weights
  // Returns true on success, false on failure
  template <
    typename DerivedV, 
    typename DerivedEle, 
    typename Derivedb,
    typename Derivedbc, 
    typename DerivedW>
  IGL_INLINE bool bbw(
    const Eigen::MatrixBase<DerivedV> & V, 
    const Eigen::MatrixBase<DerivedEle> & Ele, 
    const Eigen::MatrixBase<Derivedb> & b, 
    const Eigen::MatrixBase<Derivedbc> & bc, 
    BBWData & data,
    Eigen::MatrixBase<DerivedW> & W);
}
  
#ifdef IGL_HEADER_ONLY
#  include "bbw.cpp"
#endif

#endif
