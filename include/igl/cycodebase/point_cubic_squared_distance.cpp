#include "point_cubic_squared_distance.h"
#include "../parallel_for.h"
#include "../cubic.h"
#include <cyPolynomial.h>
#include <limits>

template <
  typename DerivedQ,
  typename DerivedC,
  typename DerivedsqrD,
  typename DerivedS,
  typename DerivedK>
void igl::cycodebase::point_cubic_squared_distance(
  const Eigen::MatrixBase<DerivedQ>& Q,
  const Eigen::MatrixBase<DerivedC>& C,
  Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
  Eigen::PlainObjectBase<DerivedS>& S,
  Eigen::PlainObjectBase<DerivedK>& K)
{
  using Scalar = typename DerivedQ::Scalar;
  constexpr int ColsAtCompileTime = DerivedQ::ColsAtCompileTime;
  using cy::PolynomialRoots;
  const int dim = C.cols();
  assert(C.rows() == 4 && "C should be 4 x dim.");
  assert(Q.cols() == dim && "Q and C should have the same number of columns.");
  // min_t ‖ Q - C(t) ‖^2
  //
  // Necessary condition for a minimum:
  //
  // d/dt ‖ Q - C(t) ‖^2 = 0
  // => d/dt (Q - C(t)) · (Q - C(t)) = 0
  // => -2 (dC/dt) · (Q - C(t)) = 0
  // => (dC/dt) · (C(t) - Q) = f(t) = 0
  //
  const auto C0 = C.row(0);
  const auto C1 = C.row(1);
  const auto C2 = C.row(2);
  const auto C3 = C.row(3);
  using RowVectorSD = Eigen::RowVector<Scalar,ColsAtCompileTime>;
  using MatrixS4D = Eigen::Matrix<Scalar,4,ColsAtCompileTime>;
  using MatrixS3D = Eigen::Matrix<Scalar,3,ColsAtCompileTime>;
  // monomial coefficients for C(t)
  MatrixS4D M(4,C.cols());
  M <<
    C.row(0).eval(),
    3 * (C1 - C0),
    3 * (C0 - 2 * C1 + C2),
    C3 - C0 + 3 * (C1 - C2);
  // Monomial coefficients for dC/dt
  MatrixS3D D(3,C.cols());
  D << M.row(1),
     2 * M.row(2),
     3 * M.row(3);
  Eigen::RowVector<Scalar,6> B(0,0,0,0,0,0);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      B(i + j) += D.row(i).dot(M.row(j));
    }
  }
  return point_cubic_squared_distance(
    Q, C, D, B, sqrD, S, K);
}

template <
  typename DerivedQ,
  typename DerivedC,
  typename DerivedD,
  typename DerivedB,
  typename DerivedsqrD,
  typename DerivedS,
  typename DerivedK>
void igl::cycodebase::point_cubic_squared_distance(
  const Eigen::MatrixBase<DerivedQ>& Q,
  const Eigen::MatrixBase<DerivedC>& C,
  const Eigen::MatrixBase<DerivedD>& D,
  const Eigen::MatrixBase<DerivedB>& B,
  Eigen::PlainObjectBase<DerivedsqrD>& sqrD,
  Eigen::PlainObjectBase<DerivedS>& S,
  Eigen::PlainObjectBase<DerivedK>& K)
{
  using Scalar = typename DerivedQ::Scalar;
  constexpr int ColsAtCompileTime = DerivedQ::ColsAtCompileTime;
  using cy::PolynomialRoots;
  const int dim = C.cols();
  using RowVectorSD = Eigen::RowVector<Scalar,ColsAtCompileTime>;
  sqrD.setConstant(Q.rows(),1,std::numeric_limits<Scalar>::infinity());
  S.resize(Q.rows());
  K.resize(Q.rows(),dim);
  const int n = Q.rows();
  for(int i = 0;i<n;i++)
  igl::parallel_for(n,[&](const int i)
  {
    // Fill in coefficients using precomputed data and Q
    // => (dC/dt) · (C(t) - Q) = f(t) = 0
    //
    // C(t) = (1-t)^3 C0 + 3(1-t)^2 t C1 + 3(1-t) t^2 C2 + t^3 C3
    // C(t) = C0 + t⋅3(C1 - C0) + t^2⋅3(C0 - 2C1 + C2) + t^3⋅(C3 - C0 + 3(C1 - C2))
    // dC/dt =       3(C1 - C0) +   t⋅6(C0 - 2C1 + C2) + t^2⋅3(C3 - C0 + 3(C1 - C2))
    Eigen::RowVector<Scalar,6> coef = B;
    for (int j = 0; j < 3; ++j)
    {
      coef(j) -= D.row(j).dot(Q.row(i));
    }

    // f is a quintic polynomial:
    constexpr int N = 5;
    Scalar r[N];
    int nr = PolynomialRoots<N>(r, coef.data(),Scalar(0),Scalar(1));
    for(int j = 0;j<nr+2;j++)
    {
      Scalar t;
      if(j==nr)
      {
        t = Scalar(0);
      }
      else if(j==nr+1)
      {
        t = Scalar(1);
      }else
      {
        t = r[j];
      }
      RowVectorSD Ct;
      igl::cubic(C,t,Ct);
      const Scalar sqrD_j = (Ct - Q.row(i)).squaredNorm();
      if(sqrD_j < sqrD(i))
      {
        sqrD(i) = sqrD_j;
        S(i) = t;
        K.row(i) = Ct;
      }
    }
  }
  ,1000);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::cycodebase::point_cubic_squared_distance<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>>(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&);
#endif
