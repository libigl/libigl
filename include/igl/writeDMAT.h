#ifndef IGL_WRITEDMAT_H
#define IGL_WRITEDMAT_H
#include "igl_inline.h"
// See writeDMAT.h for a description of the .dmat file type
#include <string>
#include <vector>
namespace igl
{
  // Write a matrix using ascii dmat file type
  //
  // Template:
  //   Mat  matrix type that supports .rows(), .cols(), operator(i,j)
  // Inputs:
  //   file_name  path to .dmat file
  //   W  eigen matrix containing to-be-written coefficients
  // Returns true on success, false on error
  //
  template <class Mat>
  IGL_INLINE bool writeDMAT(const std::string file_name, const Mat & W);
  template <typename Scalar>
  IGL_INLINE bool writeDMAT(
    const std::string file_name, 
    const std::vector<std::vector<Scalar> > W);
}

#ifdef IGL_HEADER_ONLY
#  include "writeDMAT.cpp"
#endif

#endif
