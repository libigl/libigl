#ifndef IGL_WRITEOBJ_H
#define IGL_WRITEOBJ_H
#include "igl_inline.h"
// History:
//  return type changed from void to bool  Alec 20 Sept 2011

#include <Eigen/Core>
#include <string>

namespace igl 
{
  // Write a mesh in an ascii obj file
  // Inputs:
  //   str  path to outputfile
  //   V  eigen double matrix #V by 3 (mesh vertices)
  //   F  eigen int matrix #F by 3 (mesh indices)
  // Returns true on success, false on error
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool writeOBJ(
    const std::string str,
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F);
  
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool writeOBJ(
    const std::string str,
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const Eigen::PlainObjectBase<DerivedV>& CN,
    const Eigen::PlainObjectBase<DerivedF>& FN,
    const Eigen::PlainObjectBase<DerivedT>& TC,
    const Eigen::PlainObjectBase<DerivedF>& FTC);

}

#ifdef IGL_HEADER_ONLY
#  include "writeOBJ.cpp"
#endif

#endif
