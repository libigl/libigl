#include "quadprog.h"
#include "min_quad_with_fixed.h"
//#include <igl/matlab_format.h>
//#include <igl/eigen_format.h>
#include <iostream>

template <typename Scalar, int n, int ni>
IGL_INLINE Eigen::Matrix<Scalar,n,1> igl::quadprog(
  const Eigen::Matrix<Scalar,n,n> & H,
  const Eigen::Matrix<Scalar,n,1> & f,
  const Eigen::Matrix<Scalar,ni,n> & Ai,
  const Eigen::Matrix<Scalar,ni,1> & lbi,
  const Eigen::Matrix<Scalar,ni,1> & ubi,
  const Eigen::Matrix<Scalar,n,1> & lb,
  const Eigen::Matrix<Scalar,n,1> & ub)
{
  const auto dyn_n = n == Eigen::Dynamic ? H.rows() : n;
  const auto dyn_ni = ni == Eigen::Dynamic ? Ai.rows() : ni;
  // min_x ½ xᵀ H x + xᵀ f   
  // subject to lbi ≤ Ai x ≤ ubi, lb≤x≤ub
  //
  // min_x,xi ½ xᵀ H x + xᵀ f
  // subject to Ai x - xi = 0, lbi≤xi≤ubi, lb≤x≤ub 
  //
  // min_z ½ zᵀ [H 0;0 0] z + zᵀ [f;0]
  // subject to [Ai -I] z = 0, [lb;lbi]≤z≤[ub;ubi]
  const int nn = n == Eigen::Dynamic ? Eigen::Dynamic : n+ni;
  const auto dyn_nn = nn == Eigen::Dynamic ? dyn_n+dyn_ni : nn;
  const auto make_HH = [&]()
  {
    Eigen::Matrix<Scalar,nn,nn> HH =
      Eigen::Matrix<Scalar,nn,nn>::Zero(dyn_nn,dyn_nn);
    HH.topLeftCorner(dyn_n,dyn_n) = H;
    return HH;
  };
  const Eigen::Matrix<Scalar,nn,nn> HH = make_HH();
  const auto make_ff = [&]()
  {
    Eigen::Matrix<Scalar,nn,1> ff = Eigen::Matrix<Scalar,nn,1>::Zero(dyn_nn,1);
    ff.head(dyn_n) = f;
    return ff;
  };
  const Eigen::Matrix<Scalar,nn,1> ff = make_ff();
  const auto make_lblb = [&]()
  {
    Eigen::Matrix<Scalar,nn,1> lblb =
      Eigen::Matrix<Scalar,nn,1>::Zero(dyn_nn,1);
    lblb.head(dyn_n) = lb;
    lblb.tail(dyn_ni) = lbi;
    return lblb;
  };
  const Eigen::Matrix<Scalar,nn,1> lblb = make_lblb();
  const auto make_ubub = [&]()
  {
    Eigen::Matrix<Scalar,nn,1> ubub =
      Eigen::Matrix<Scalar,nn,1>::Zero(dyn_nn,1);
    ubub.head(dyn_n) = ub;
    ubub.tail(dyn_ni) = ubi;
    return ubub;
  };
  const Eigen::Matrix<Scalar,nn,1> ubub = make_ubub();
  const auto make_AA = [&]()
  {
    Eigen::Matrix<Scalar,ni,nn> AA =
      Eigen::Matrix<Scalar,ni,nn>::Zero(dyn_ni,dyn_nn);
    AA.leftCols(dyn_n) = Ai;
    AA.rightCols(dyn_ni) =
      -Eigen::Matrix<Scalar,ni,ni>::Identity(dyn_ni,dyn_ni);
    return AA;
  };
  const Eigen::Matrix<Scalar,ni,nn> AA = make_AA();
  const Eigen::Matrix<Scalar,ni,1> bb = 
   Eigen::Matrix<Scalar,ni,1>::Zero(dyn_ni,1);
  // This ends make a system that looks like [H 0 Aᵀ;0 0 -I;A I 0] which is rank
  // deficient.
  const Eigen::Matrix<Scalar,nn,1> xx = 
    igl::quadprog(HH,ff,AA,bb,lblb,ubub);
  return xx.head(dyn_n);
}


template <typename Scalar, int n, int m>
IGL_INLINE Eigen::Matrix<Scalar,n,1> igl::quadprog(
  const Eigen::Matrix<Scalar,n,n> & H,
  const Eigen::Matrix<Scalar,n,1> & f,
  const Eigen::Matrix<Scalar,m,n> & A,
  const Eigen::Matrix<Scalar,m,1> & b,
  const Eigen::Matrix<Scalar,n,1> & lb,
  const Eigen::Matrix<Scalar,n,1> & ub)
{
  //std::cout<<igl::matlab_format(H,"H")<<std::endl;
  //std::cout<<igl::matlab_format(f,"f")<<std::endl;
  //std::cout<<igl::matlab_format(A,"A")<<std::endl;
  //std::cout<<igl::matlab_format(b,"b")<<std::endl;
  //std::cout<<igl::matlab_format(lb,"lb")<<std::endl;
  //std::cout<<igl::matlab_format(ub,"ub")<<std::endl;
  typedef Eigen::Matrix<Scalar,n,1> VectorSn;
  typedef Eigen::Array<bool,n,1>    Arraybn;
  assert( (lb.array() < ub.array() ).all() );
  const int dyn_n = n == Eigen::Dynamic ? H.rows() : n;
  VectorSn x(dyn_n);
  VectorSn bc = VectorSn::Constant(dyn_n,1,-1e26);
  Arraybn k = Arraybn::Constant(dyn_n,1,false);
  Eigen::Index iter;
  // n³ is probably way too conservative. 
  for(iter = 0;iter<dyn_n*dyn_n*dyn_n;iter++)
  {
    // For dual contouring 99% of the time the active set is empty.
    // Optimize for this common case.
    x = min_quad_with_fixed(H,f,k,bc,A,b);
    //std::cout<<igl::matlab_format(x,"x")<<std::endl;
    // constraint violations 
    VectorSn vl = lb-x;
    VectorSn vu = x-ub;

    // try to add/remove constraints
    Eigen::Index best_add = -1; Scalar worst_offense = 0;
    bool add_lower;
    Eigen::Index best_remove = -1; Scalar worst_lambda = 0;
    for(Eigen::Index i = 0;i<dyn_n;i++)
    {
      if(vl(i)>worst_offense)
      {
        best_add = i;
        add_lower = true;
        worst_offense = vl(i);
      }
      if(vu(i)>worst_offense)
      {
        best_add = i;
        add_lower = false;
        worst_offense = vu(i);
      }
      // bias toward adding constraints
      if(best_add<0 && k(i))
      {
        const Scalar sign = bc(i)==ub(i)?1:-1;
        const Scalar lambda_i = sign * (H.row(i)*x+f(i));
        //printf("  considering k(%d) (λ = %g)\n",i,lambda_i);
        if(lambda_i > worst_lambda)
        {
          //printf("  removing k(%d) (λ = %g)\n",i,lambda_i);
          best_remove = i;
          worst_lambda = lambda_i;
        }
      }
    }
    // bias toward adding constraints
    if(best_add >= 0)
    {
      const auto i = best_add;
      assert(!k(i));
      //add_lower ? printf("  adding lb(%d)\n",i) : printf("  adding lb(%d)\n",i);
      bc(i) = add_lower ? lb(i) : ub(i);
      k(i) = true;
    }else if(best_remove >= 0)
    {
      const auto i = best_remove;
      assert(k(i));
      //printf("  removing k(%d) (λ = %g)\n",i,worst_lambda);
      k(i) = false;
    }else /*if(best_add < 0 && best_remove < 0)*/
    {
      return x;
    }
  }
  // Should never happen.
  assert(false && "quadprog failed after too many iterations");
  //std::cout<<igl::eigen_format(H,"H")<<std::endl;
  //std::cout<<igl::eigen_format(f,"f")<<std::endl;
  //std::cout<<igl::eigen_format(lb,"lb")<<std::endl;
  //std::cout<<igl::eigen_format(ub,"ub")<<std::endl;
  return VectorSn::Zero(dyn_n);
}

template <typename Scalar, int n>
IGL_INLINE Eigen::Matrix<Scalar,n,1> igl::quadprog(
  const Eigen::Matrix<Scalar,n,n> & H,
  const Eigen::Matrix<Scalar,n,1> & f,
  const Eigen::Matrix<Scalar,n,1> & lb,
  const Eigen::Matrix<Scalar,n,1> & ub)
{
  const int m = n == Eigen::Dynamic ? Eigen::Dynamic : 0;
  return quadprog(
    H,f,
    Eigen::Matrix<Scalar,m,n>(0,H.cols()),
    Eigen::Matrix<Scalar,m,1>(0,1),
    lb,ub);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> igl::quadprog<double, -1, -1>(Eigen::Matrix<double, -1, -1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)1) : ((((-1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, -1> const&, Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> const&, Eigen::Matrix<double, -1, -1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)1) : ((((-1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, -1> const&, Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> const&, Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> const&, Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> const&, Eigen::Matrix<double, -1, 1, ((Eigen::StorageOptions)0) | ((((-1) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((-1) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), -1, 1> const&);
// generated by autoexplicit.sh
template Eigen::Matrix<double, 6, 1, ((Eigen::StorageOptions)0) | ((((6) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 6, 1> igl::quadprog<double, 6, 2>(Eigen::Matrix<double, 6, 6, ((Eigen::StorageOptions)0) | ((((6) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)1) : ((((6) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 6, 6> const&, Eigen::Matrix<double, 6, 1, ((Eigen::StorageOptions)0) | ((((6) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 6, 1> const&, Eigen::Matrix<double, 2, 6, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)1) : ((((6) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 6> const&, Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> const&, Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> const&, Eigen::Matrix<double, 6, 1, ((Eigen::StorageOptions)0) | ((((6) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 6, 1> const&, Eigen::Matrix<double, 6, 1, ((Eigen::StorageOptions)0) | ((((6) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((6) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 6, 1> const&);
// generated by autoexplicit.sh
template Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> igl::quadprog<double, 2>(Eigen::Matrix<double, 2, 2, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)1) : ((((2) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 2> const&, Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> const&, Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> const&, Eigen::Matrix<double, 2, 1, ((Eigen::StorageOptions)0) | ((((2) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((2) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 2, 1> const&);
// generated by autoexplicit.sh
template Eigen::Matrix<double, 3, 1, ((Eigen::StorageOptions)0) | ((((3) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 3, 1> igl::quadprog<double, 3>(Eigen::Matrix<double, 3, 3, ((Eigen::StorageOptions)0) | ((((3) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)1) : ((((3) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 3, 3> const&, Eigen::Matrix<double, 3, 1, ((Eigen::StorageOptions)0) | ((((3) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 3, 1> const&, Eigen::Matrix<double, 3, 1, ((Eigen::StorageOptions)0) | ((((3) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 3, 1> const&, Eigen::Matrix<double, 3, 1, ((Eigen::StorageOptions)0) | ((((3) == (1)) && ((1) != (1))) ? ((Eigen::StorageOptions)1) : ((((1) == (1)) && ((3) != (1))) ? ((Eigen::StorageOptions)0) : ((Eigen::StorageOptions)0))), 3, 1> const&);
#endif
