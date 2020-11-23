#ifndef IGL_SCREEN_SPACE_SELECTION_H
#define IGL_SCREEN_SPACE_SELECTION_H

#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <vector>
// Forward declaration
namespace igl { template <typename DerivedV, int DIM> class AABB; }

namespace igl
{
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedM,
    typename DerivedN,
    typename DerivedO,
    typename Ltype,
    typename DerivedW,
    typename Deriveda>
  IGL_INLINE void screen_space_selection(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const igl::AABB<DerivedV, 3> & tree,
    const Eigen::MatrixBase<DerivedM>& model,
    const Eigen::MatrixBase<DerivedN>& proj,
    const Eigen::MatrixBase<DerivedO>& viewport,
    const std::vector<Eigen::Matrix<Ltype,1,2> > & L,
    Eigen::PlainObjectBase<DerivedW> & W,
    Eigen::PlainObjectBase<Deriveda> & and_visible);
  template <
    typename DerivedV,
    typename DerivedM,
    typename DerivedN,
    typename DerivedO,
    typename Ltype,
    typename DerivedW>
  IGL_INLINE void screen_space_selection(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedM>& model,
    const Eigen::MatrixBase<DerivedN>& proj,
    const Eigen::MatrixBase<DerivedO>& viewport,
    const std::vector<Eigen::Matrix<Ltype,1,2> > & L,
    Eigen::PlainObjectBase<DerivedW> & W);
  template <
    typename DerivedV,
    typename DerivedM,
    typename DerivedN,
    typename DerivedO,
    typename DerivedP,
    typename DerivedE,
    typename DerivedW>
  IGL_INLINE void screen_space_selection(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedM>& model,
    const Eigen::MatrixBase<DerivedN>& proj,
    const Eigen::MatrixBase<DerivedO>& viewport,
    const Eigen::MatrixBase<DerivedP> & P,
    const Eigen::MatrixBase<DerivedE> & E,
    Eigen::PlainObjectBase<DerivedW> & W);
}

#ifndef IGL_STATIC_LIBRARY
#include "screen_space_selection.cpp"
#endif
  
#endif
