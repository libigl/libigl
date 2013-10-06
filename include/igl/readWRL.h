#ifndef IGL_READWRL_H
#define IGL_READWRL_H
#include "igl_inline.h"

#include <string>
#include <vector>

namespace igl 
{
  // Read a mesh from an ascii wrl file, filling in vertex positions and face
  // indices of the first model. Mesh may have faces of any number of degree
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .wrl file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   F  #F list of face indices into vertex positions
  // Returns true on success, false on errors
  template <typename Scalar, typename Index>
  IGL_INLINE bool readWRL(
    const std::string wrl_file_name, 
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Index > > & F);

}

#ifdef IGL_HEADER_ONLY
#  include "readWRL.cpp"
#endif

#endif

