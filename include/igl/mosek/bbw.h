#ifndef BBW_H
#define BBW_H
#include "../igl_inline.h"

#include <Eigen/Dense>
#include "mosek_quadprog.h"
#include <igl/active_set.h>

namespace igl
{
  enum QPSolver
  {
    QP_SOLVER_IGL_ACTIVE_SET = 0,
    QP_SOLVER_MOSEK = 1,
    NUM_QP_SOLVERS = 2
  };
  const char * const QPSolverNames[NUM_QP_SOLVERS] =
  {
    "QP_SOLVER_IGL_ACTIVE_SET",
    "QP_SOLVER_MOSEK"
  };
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
      igl::active_set_params active_set_params;
      // Which solver
      QPSolver qp_solver;
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
    const Eigen::PlainObjectBase<DerivedV> & V, 
    const Eigen::PlainObjectBase<DerivedEle> & Ele, 
    const Eigen::PlainObjectBase<Derivedb> & b, 
    const Eigen::PlainObjectBase<Derivedbc> & bc, 
    BBWData & data,
    Eigen::PlainObjectBase<DerivedW> & W);
}
  
#ifdef IGL_HEADER_ONLY
#  include "bbw.cpp"
#endif

#endif
