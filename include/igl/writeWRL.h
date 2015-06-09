#ifndef IGL_WRITE_WRL_H
#define IGL_WRITE_WRL_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Write mesh to a .wrl file
  //
  // Inputs:
  //   str  path to .wrl file
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  // Returns true iff succes
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeWRL(
    const std::string & str,
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F);
}
#ifndef IGL_STATIC_LIBRARY
#include "writeWRL.cpp"
#endif
#endif
