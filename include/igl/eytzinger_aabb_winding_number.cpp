#include "eytzinger_aabb_winding_number.h"
template <
  typename Derivedp,
  typename DerivedV,
  typename DerivedE,
  typename DerivedB,
  typename Derivedleaf,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE void igl::eytzinger_aabb_winding_number(
  const Eigen::MatrixBase<Derivedp> & p,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedB> & B1,
  const Eigen::MatrixBase<DerivedB> & B2,
  const Eigen::MatrixBase<Derivedleaf> & leaf,
  const Eigen::MatrixBase<DerivedI> & I,
  const Eigen::MatrixBase<DerivedC> & C,
  typename DerivedV::Scalar & wn)
{
  const auto signed_angle = [&V,&p]( const int e0, const int e1)->typename DerivedV::Scalar
  {
    const Eigen::Matrix<typename DerivedV::Scalar,1,2> v0 = V.row(e0) - p;
    const Eigen::Matrix<typename DerivedV::Scalar,1,2> v1 = V.row(e1) - p;
    const typename DerivedV::Scalar angle = std::atan2(
      v0(0)*v1(1) - v0(1)*v1(0),
      v0(0)*v1(0) + v0(1)*v1(1));
    return angle;
  };

  wn = 0;
  if(leaf.size() == 0) { wn = 0; return; }
  // I don't think stack or queue matters
  std::vector<int> stack;
  stack.push_back(0);
  while(!stack.empty())
  {
    int r = stack.back();
    stack.pop_back();
    if(leaf(r) >= 0)
    {
      wn += signed_angle(E(leaf(r),0),E(leaf(r),1));
      continue;
    }
    // otherwise is p outside the box B1.row(r), B2.row(r)?
    if( (p(0) < B1(r,0)) || (p(0) > B2(r,0)) ||
        (p(1) < B1(r,1)) || (p(1) > B2(r,1)) )
    {
      // consider the pairs in the sequence at I(C(r)) to I(C(r)+1)
      assert(C(r+1) - C(r) % 2 == 0);
      for(int i = C(r); i < C(r+1); i+=2)
      {
        wn += signed_angle(I(i),I(i+1));
      }
      continue;
    }
    // else traverse children
    const int left_r = 2*r + 1;
    const int right_r = 2*r + 2;
    stack.push_back(left_r);
    stack.push_back(right_r);
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::eytzinger_aabb_winding_number<Eigen::Matrix<double, 1, 2, 1, 1, 2>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 2, 1, -1, 2>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>>(Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2>> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 2, 1, -1, 2>> const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 2, 1, -1, 2>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1>> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar&);
#endif

