#include "active_set.h"
#include "min_quad_with_fixed.h"
#include "slice.h"
#include "cat.h"
#include "matlab_format.h"

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
  assert(known.cols() == 1);
  // Number of knowns
  const int nk = known.size();

  // Initialize active sets
  typedef bool BOOL;
#define TRUE true
#define FALSE false
  Matrix<BOOL,Dynamic,1> as_lx = Matrix<BOOL,Dynamic,1>::Constant(n,1,FALSE);
  Matrix<BOOL,Dynamic,1> as_ux = Matrix<BOOL,Dynamic,1>::Constant(n,1,FALSE);
  Matrix<BOOL,Dynamic,1> as_ieq(Aieq.rows(),1);

  // Keep track of previous Z for comparison
  PlainObjectBase<DerivedZ> old_Z;
  old_Z = PlainObjectBase<DerivedZ>::Constant(
      n,1,numeric_limits<typename DerivedZ::Scalar>::max());

  int iter = 0;
  while(true)
  {
    // FIND BREACHES OF CONSTRAINTS
    int new_as_lx = 0;
    int new_as_ux = 0;
    int new_as_ieq = 0;
    if(Z.size() > 0)
    {
      for(int z = 0;z < n;z++)
      {
        if(Z(z) < lx(z))
        {
          new_as_lx += (as_lx(z)?0:1);
          //new_as_lx++;
          as_lx(z) = TRUE;
        }
        if(Z(z) > ux(z))
        {
          new_as_ux += (as_ux(z)?0:1);
          //new_as_ux++;
          as_ux(z) = TRUE;
        }
      }
      PlainObjectBase<DerivedZ> AieqZ;
      AieqZ = Aieq*Z;
      for(int a = 0;a<Aieq.rows();a++)
      {
        if(AieqZ(a) > Bieq(a))
        {
          new_as_ieq += (as_ieq(a)?0:1);
          as_ieq(a) = TRUE;
        }
      }
      //cout<<"new_as_lx: "<<new_as_lx<<endl;
      //cout<<"new_as_ux: "<<new_as_ux<<endl;
      const double diff = (Z-old_Z).squaredNorm();
      //cout<<"diff: "<<diff<<endl;
      if(diff < params.solution_diff_threshold)
      {
        ret = SOLVER_STATUS_CONVERGED;
        break;
      }
      old_Z = Z;
    }

    const int as_lx_count = count(as_lx.data(),as_lx.data()+n,TRUE);
    const int as_ux_count = count(as_ux.data(),as_ux.data()+n,TRUE);
    const int as_ieq_count = count(as_ieq.data(),as_ieq.data()+n,TRUE);

    // PREPARE FIXED VALUES
    PlainObjectBase<Derivedknown> known_i;
    known_i.resize(nk + as_lx_count + as_ux_count,1);
    PlainObjectBase<DerivedY> Y_i;
    Y_i.resize(nk + as_lx_count + as_ux_count,1);
    {
      known_i.block(0,0,known.rows(),known.cols()) = known;
      Y_i.block(0,0,Y.rows(),Y.cols()) = Y;
      int k = nk;
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
    //cout<<matlab_format((known_i.array()+1).eval(),"known_i")<<endl;
    // PREPARE EQUALITY CONSTRAINTS
    VectorXi as_ieq_list(as_ieq_count,1);
    // Gather active constraints and resp. rhss
    PlainObjectBase<DerivedBeq> Beq_i;
    Beq_i.resize(Beq.rows()+as_ieq_count,1);
    {
      int k =0;
      for(int a=0;a<as_ieq.size();a++)
      {
        if(a)
        {
          as_ieq_list(k)=a;
          Beq_i(Beq.rows()+k,1) = Bieq(k,1);
          k++;
        }
      }
      assert(k == as_ieq_count);
    }
    // extract active constraint rows
    SparseMatrix<AeqT> Aeq_i,Aieq_i;
    slice(Aieq,as_ieq_list,1,Aieq_i);
    // Append to equality constraints
    cat(1,Aeq,Aieq_i,Aeq_i);


    min_quad_with_fixed_data<AT> data;
    if(!min_quad_with_fixed_precompute(A,known_i,Aeq_i,params.Auu_pd,data))
    {
      cerr<<"Error: min_quad_with_fixed precomputation failed."<<endl;
      ret = SOLVER_STATUS_ERROR;
      break;
    }
    Eigen::PlainObjectBase<DerivedZ> sol;
    if(!min_quad_with_fixed_solve(data,B,Y_i,Beq_i,Z,sol))
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
    Lambda_known_i.block(nk,0,as_lx_count,1) = 
      (-1*Lambda_known_i.block(nk,0,as_lx_count,1)).eval();

    // Extract Lagrange multipliers for Aieq_i (always at back of sol)
    VectorXd Lambda_Aieq_i(Aieq_i.rows(),1);
    for(int l = 0;l<Aieq_i.rows();l++)
    {
      Lambda_Aieq_i(Aieq_i.rows()-1-l) = sol(sol.rows()-1-l);
    }
    
    // Remove from active set
    for(int l = 0;l<as_lx_count;l++)
    {
      if(Lambda_known_i(nk + l) < params.inactive_threshold)
      {
        as_lx(known_i(nk + l)) = FALSE;
      }
    }
    for(int u = 0;u<as_ux_count;u++)
    {
      if(Lambda_known_i(nk + as_lx_count + u) < 
        params.inactive_threshold)
      {
        as_ux(known_i(nk + as_lx_count + u)) = FALSE;
      }
    }
    for(int a = 0;a<as_ieq_count;a++)
    {
      if(Lambda_Aieq_i(a) < params.inactive_threshold)
      {
        as_ieq(as_ieq_list(a)) = FALSE;
      }
    }

    iter++;
    //cout<<iter<<endl;
    if(params.max_iter>0 && iter>=params.max_iter)
    {
      ret = SOLVER_STATUS_MAX_ITER;
      break;
    }

  }

  return ret;
}


#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template igl::SolverStatus igl::active_set<double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, double, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::SparseMatrix<double, 0, int> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, igl::active_set_params const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
