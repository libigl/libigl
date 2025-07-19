#ifndef IGL_EYTZINGER_AABB_SDF_H
#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// 
  /// @param[in] B1  #B by dim list of minimum corners of the Eytzinger AABBs
  /// @param[in] B2  #B by dim list of maximum corners of the Eytzinger AABBs
  /// @param[in] leaf #B list of leaf indices, -1 indicates internal node, -2
  /// indicates empty node
  /// @param[in] primitive  function handle that takes as input a primitive id
  /// and returns the SDF to that primitive `primitive(i)`
  /// @param[out] f  SDF value at query point x
  ///
  template <
    typename Derivedp,
    typename Func,
    typename DerivedB,
    typename Derivedleaf
  >
  IGL_INLINE void eytzinger_aabb_sdf(
    const Eigen::MatrixBase<Derivedp> & p,
    const Func & primitive,
    const Eigen::MatrixBase<DerivedB> & B1,
    const Eigen::MatrixBase<DerivedB> & B2,
    const Eigen::MatrixBase<Derivedleaf> & leaf,
    typename Derivedp::Scalar & f);
  template <
    typename DerivedP,
    typename Func,
    typename DerivedB,
    typename Derivedleaf,
    typename DerivedS
  >
  IGL_INLINE void eytzinger_aabb_sdf(
    const Eigen::MatrixBase<DerivedP> & P,
    const Func & primitive,
    const Eigen::MatrixBase<DerivedB> & B1,
    const Eigen::MatrixBase<DerivedB> & B2,
    const Eigen::MatrixBase<Derivedleaf> & leaf,
    Eigen::PlainObjectBase<DerivedS> & S);
}

#ifndef IGL_STATIC_LIBRARY
#  include "eytzinger_aabb_sdf.cpp"
#endif
#endif 

