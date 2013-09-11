#ifndef ACTIVE_SET_H
#define ACTIVE_SET_H

#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
  enum SolverStatus
  {
    // Good
    SOLVER_STATUS_CONVERGED = 0,
    // OK
    SOLVER_STATUS_MAX_ITER = 1,
    // Bad
    SOLVER_STATUS_ERROR = 2,
    NUM_SOLVER_STATUSES = 3,
  };
  struct active_set_params;
  // ACTIVE_SET Minimize quadratic energy Z'*A*Z + Z'*B + C with constraints
  // that Z(known) = Y, optionally also subject to the constraints Aeq*Z = Beq,
  // and further optionally subject to the linear inequality constraints that
  // Aieq*Z <= Bieq and constant inequality constraints lx <= x <= ux
  //
  // Templates:
  // Inputs:
  //   A  n by n matrix of quadratic coefficients
  //   B  n by 1 column of linear coefficients
  //   known  list of indices to known rows in Z
  //   Y  list of fixed values corresponding to known rows in Z
  //   Aeq  meq by n list of linear equality constraint coefficients
  //   Beq  meq by 1 list of linear equality constraint constant values
  //   Aieq  mieq by n list of linear equality constraint coefficients
  //   Bieq  mieq by 1 list of linear equality constraint constant values
  //   lx  n by 1 list of lower bounds [] implies -Inf
  //   ux  n by 1 list of upper bounds [] implies Inf
  //   params  struct of additional parameters (see below)
  // Outputs:
  //   Z  n by 1 list of solution values
  // Returns true on success, false on error
  //
  // Benchmark: For a harmonic solve on a mesh with 325K facets, matlab 2.2
  // secs, igl/min_quad_with_fixed.h 7.1 secs
  //

  template <
    typename AT, 
    typename DerivedB,
    typename Derivedknown, 
    typename DerivedY,
    typename AeqT,
    typename DerivedBeq,
    typename AieqT,
    typename DerivedBieq,
    typename Derivedlx,
    typename Derivedux,
    typename DerivedZ
    >
  IGL_INLINE igl::SolverStatus active_set(
    const Eigen::SparseMatrix<AT>& A,
    const Eigen::PlainObjectBase<DerivedB> & B,
    const Eigen::PlainObjectBase<Derivedknown> & known,
    const Eigen::PlainObjectBase<DerivedY> & Y,
    const Eigen::SparseMatrix<AeqT>& Aeq,
    const Eigen::PlainObjectBase<DerivedBeq> & Beq,
    const Eigen::SparseMatrix<AieqT>& Aieq,
    const Eigen::PlainObjectBase<DerivedBieq> & Bieq,
    const Eigen::PlainObjectBase<Derivedlx> & lx,
    const Eigen::PlainObjectBase<Derivedux> & ux,
    const igl::active_set_params & params,
    Eigen::PlainObjectBase<DerivedZ> & Z
    );
};

#include "EPS.h"
struct igl::active_set_params
{
  // Input parameters for active_set:
  //   Auu_pd  whethter Auu is positive definite {false}
  //   max_iter  Maximum number of iterations ({0} = Infinity)
  //   inactive_threshold  Threshold on Lagrange multiplier values to determine
  //     whether to keep constraints active {EPS}
  bool Auu_pd;
  int max_iter;
  double inactive_threshold;
  active_set_params():
    Auu_pd(false),
    max_iter(-1),
    inactive_threshold(igl::DOUBLE_EPS)
    {};
};

#ifdef IGL_HEADER_ONLY
#  include "active_set.cpp"
#endif

#endif
