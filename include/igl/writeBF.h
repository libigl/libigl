#ifndef IGL_WRITEBF_H
#define IGL_WRITEBF_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Write a bones forest to a file
  //
  // Input:
  //   file_name  path to .bf bones tree file
  //   WI  #B list of unique weight indices
  //   P  #B list of parent indices into B, -1 for roots
  //   O  #B list of tip offsets
  // Returns true on success, false on errors
  template < 
    typename DerivedWI,
    typename DerivedP,
    typename DerivedO>
  IGL_INLINE bool writeBF(
    const std::string & filename,
    const Eigen::PlainObjectBase<DerivedWI> & WI,
    const Eigen::PlainObjectBase<DerivedP> & P,
    const Eigen::PlainObjectBase<DerivedO> & O);
}

#ifndef IGL_STATIC_LIBRARY
#  include "writeBF.cpp"
#endif
#endif
