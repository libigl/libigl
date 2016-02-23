#ifndef IGL_READBF_H
#define IGL_READBF_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <string>
namespace igl
{
  // Read a bones forest from a file, returns a list of bone roots
  // Input:
  //   file_name  path to .bf bones tree file
  // Output:
  //   WI  #B list of unique weight indices
  //   P  #B list of parent indices into B, -1 for roots
  //   C  #B list of tip positions
  // Returns true on success, false on errors
  template < 
    typename DerivedWI,
    typename DerivedP,
    typename DerivedC>
  IGL_INLINE bool readBF(
    const std::string & filename,
    Eigen::PlainObjectBase<DerivedWI> & WI,
    Eigen::PlainObjectBase<DerivedP> & P,
    Eigen::PlainObjectBase<DerivedC> & C);
}

#ifndef IGL_STATIC_LIBRARY
#  include "readBF.cpp"
#endif
#endif
