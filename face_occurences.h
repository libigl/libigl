#ifndef IGL_FACE_OCCURENCES
#define IGL_FACE_OCCURENCES
#include "igl_inline.h"

#include <vector>
namespace igl
{
  // Count the occruances of each face (row) in a list of face indices
  // (irrespecitive of order)
  // Inputs:
  //   F  #F by simplex-size
  // Outputs
  //   C  #F list of counts
  // Known bug: triangles/tets only (where ignoring order still gives simplex)
  template <typename IntegerF, typename IntegerC>
  IGL_INLINE void face_occurences(
    const std::vector<std::vector<IntegerF> > & F,
    std::vector<IntegerC> & C);
}

#ifdef IGL_HEADER_ONLY
#  include "face_occurences.cpp"
#endif

#endif


