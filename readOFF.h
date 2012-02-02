//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#ifndef IGL_READOFF_H
#define IGL_READOFF_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>

namespace igl 
{
  
  // Read a mesh from an ascii obj file, filling in vertex positions, normals
  // and texture coordinates. Mesh may have faces of any number of degree
  //
  // Templates:
  //   Scalar  type for positions and vectors (will be read as double and cast
  //     to Scalar)
  //   Index  type for indices (will be read as int and cast to Index)
  // Inputs:
  //  str  path to .obj file
  // Outputs:
  //   V  double matrix of vertex positions  #V by 3
  //   F  #F list of face indices into vertex positions
  //   TC  double matrix of texture coordinats #TC by 2
  //   FTC  #F list of face indices into vertex texture coordinates
  //   N  double matrix of corner normals #N by 3
  //   FN  #F list of face indices into vertex normals
  // Returns true on success, false on errors
  template <typename Scalar, typename Index>
  IGL_INLINE bool readOFF(
                          const std::string off_file_name, 
                          std::vector<std::vector<Scalar > > & V,
                          std::vector<std::vector<Index > > & F);
  
  
  // read mesh from a ascii off file
  // Inputs:
  //   str  path to .off file
  // Outputs:
  //   V  eigen double matrix #V by 3
  //   F  eigen int matrix #F by 3
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool readOFF(
                          const std::string str,
                          Eigen::PlainObjectBase<DerivedV>& V,
                          Eigen::PlainObjectBase<DerivedF>& F);
}

#ifdef IGL_HEADER_ONLY
#  include "readOFF.cpp"
#endif

#endif
