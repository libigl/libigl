#include "active_set.h"
#include "min_quad_with_fixed.h"
#include "slice.h"

#include <iostream>
#include <limits>
#include <algorithm>

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
IGL_INLINE igl::SolverStatus igl::active_set(
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
  )
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  // Linear inequality constraints are not supported yet
  assert(Aieq.size() == 0);
  assert(Bieq.size() == 0);
  SolverStatus ret = SOLVER_STATUS_ERROR;
  const int n = A.rows();
  assert(n == A.cols());
  // Discard const qualifiers
  //if(B.size() == 0)
  //{
  //  B = Eigen::PlainObjectBase<DerivedB>::Zero(n,1);
  //}
  assert(n == B.rows());
  assert(B.cols() == 1);
  assert(Y.cols() == 1);
  assert((Aeq.size() == 0 && Beq.size() == 0) || Aeq.cols() == n);
  assert((Aeq.size() == 0 && Beq.size() == 0) || Aeq.rows() == Beq.rows());
  assert((Aeq.size() == 0 && Beq.size() == 0) || Beq.cols() == 1);
  assert((Aieq.size() == 0 && Bieq.size() == 0) || Aieq.cols() == n);
  assert((Aieq.size() == 0 && Bieq.size() == 0) || Aieq.rows() == Bieq.rows());
  assert((Aieq.size() == 0 && Bieq.size() == 0) || Bieq.cols() == 1);
  // Discard const qualifiers
  //if(lx.size() == 0)
  //{
  //  lx = Eigen::PlainObjectBase<Derivedlx>::Constant(
  //    n,1,numeric_limits<typename Derivedlx::Scalar>::min());
  //}
  //if(ux.size() == 0)
  //{
  //  ux = Eigen::PlainObjectBase<Derivedux>::Constant(
  //    n,1,numeric_limits<typename Derivedux::Scalar>::max());
  //}
  assert(lx.rows() == n);
  assert(ux.rows() == n);
  assert(ux.cols() == 1);
  assert(lx.cols() == 1);
  assert((ux.array()-lx.array()).minCoeff() > 0);
  if(Z.size() != 0)
  {
    // Initial guess should have correct size
    assert(Z.rows() == n);
    assert(Z.cols() == 1);
  }

  // Initialize active sets
  Matrix<bool,Dynamic,1> as_lx(n,1);
  Matrix<bool,Dynamic,1> as_ux(n,1);
  Matrix<bool,Dynamic,1> as_ieq(Aieq.rows(),1);

  int iter = 0;
  while(true)
  {
    // FIND BREACHES OF CONSTRAINTS
    if(Z.size() > 0)
    {
      for(int z = 0;z < n;z++)
      {
        if(Z(z) < lx(z))
        {
          as_lx(z) = true;
        }
        if(Z(z) > ux(z))
        {
          as_ux(z) = true;
        }
      }
    }

    const int as_lx_count = count(as_lx.data(),as_lx.data()+n,true);
    const int as_ux_count = count(as_ux.data(),as_ux.data()+n,true);
    // PREPARE FIXED VALUES
    Eigen::PlainObjectBase<Derivedknown> known_i;
    known_i.resize((int)known.size() + as_lx_count + as_ux_count,1);
    Eigen::PlainObjectBase<DerivedY> Y_i;
    Y_i.resize((int)known.size() + as_lx_count + as_ux_count,1);
    {
      known_i.block(0,0,known.rows(),known.cols()) = known;
      Y_i.block(0,0,Y.rows(),Y.cols()) = Y;
      int k = known.size();
      // Then all lx
      for(int z = 0;z < n;z++)
      {
        if(as_lx(z))
        {
          known_i(k) = z;
          Y_i(k) = lx(z);
          k++;
        }
      }
      // Finally all ux
      for(int z = 0;z < n;z++)
      {
        if(as_ux(z))
        {
          known_i(k) = z;
          Y_i(k) = ux(z);
          k++;
        }
      }
      assert(k==Y_i.size());
      assert(k==known_i.size());
    }

    min_quad_with_fixed_data<AT> data;
    if(!min_quad_with_fixed_precompute(A,known_i,Aeq,params.Auu_pd,data))
    {
      cerr<<"Error: min_quad_with_fixed precomputation failed."<<endl;
      ret = SOLVER_STATUS_ERROR;
      break;
    }
    if(!min_quad_with_fixed_solve(data,B,Y_i,Beq,Z))
    {
      cerr<<"Error: min_quad_with_fixed solve failed."<<endl;
      ret = SOLVER_STATUS_ERROR;
      break;
    }

    // Compute Lagrange multiplier values for known_i
    // This needs to be adjusted slightly if A is not symmetric
    assert(data.Auu_sym);
    SparseMatrix<AT> Ak;
    // Slow
    slice(A,known_i,1,Ak);
    Eigen::PlainObjectBase<DerivedB> Bk;
    slice(B,known_i,Bk);
    MatrixXd Lambda_known_i = -(Ak*Z + 0.5*Bk);
    // reverse the lambda values for lx
    Lambda_known_i.block(known.size(),0,as_lx_count,1) = 
      (-1*Lambda_known_i.block(known.size(),0,as_lx_count,1)).eval();
    
    // Remove from active set
    for(int z = 0;z<as_lx_count;z++)
    {
      if(Lambda_known_i(known.size() + z) < params.inactive_threshold)
      {
        as_lx(known_i(z)) = false;
      }
    }
    for(int z = 0;z<as_ux_count;z++)
    {
      if(Lambda_known_i(known.size() + as_lx_count + z) < 
        params.inactive_threshold)
      {
        as_ux(known_i(z)) = false;
      }
    }

    iter++;
    if(params.max_iter>0 && iter>=params.max_iter)
    {
      ret = SOLVER_STATUS_MAX_ITER;
      break;
    }
  }

finish:
  return ret;
}


#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template igl::SolverStatus igl::active_set<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, igl::active_set_params const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
