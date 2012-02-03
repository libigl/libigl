//
//  IGL Lib - Simple C++ mesh library 
//
//  Copyright 2011, Daniele Panozzo. All rights reserved.

// History:
//  return type changed from void to bool  Alec 18 Sept 2011
//  added pure vector of vectors version that has much more support Alec 31 Oct
//    2011

#ifndef IGL_READOBJ_H
#define IGL_READOBJ_H
#include "igl_inline.h"

#include <Eigen/Core>
#include <string>
#include <vector>

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
  IGL_INLINE bool readOBJ(
    const std::string obj_file_name, 
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Scalar > > & TC,
    std::vector<std::vector<Scalar > > & N,
    std::vector<std::vector<Index > > & F,
    std::vector<std::vector<Index > > & FTC,
    std::vector<std::vector<Index > > & FN);

  //! Read a mesh from an ascii obj file
  // Inputs:
  //   str  path to .obj file
  // Outputs:
  //   V  eigen matrix #V by 3
  //   F  eigen matrix #F by 3
  //
  // KNOWN BUG: This only knows how to read *triangle* meshes. It will probably
  // crash or give garbage on anything else.
  //
  // KNOWN BUG: This only knows how to face lines without normal or texture
  // indices. It will probably crash or give garbage on anything else.
  template <typename DerivedV, typename DerivedF, typename DerivedT>
  IGL_INLINE bool readOBJ(
                          const std::string str,
                          Eigen::PlainObjectBase<DerivedV>& V,
                          Eigen::PlainObjectBase<DerivedF>& F,
                          Eigen::PlainObjectBase<DerivedV>& CN,
                          Eigen::PlainObjectBase<DerivedF>& FN,
                          Eigen::PlainObjectBase<DerivedT>& TC,
                          Eigen::PlainObjectBase<DerivedF>& FTC);

  //! Read a mesh from an ascii obj file
  // Inputs:
  //   str  path to .obj file
  // Outputs:
  //   V  eigen matrix #V by 3
  //   F  eigen matrix #F by 3
  //
  // KNOWN BUG: This only knows how to read *triangle* meshes. It will probably
  // crash or give garbage on anything else.
  //
  // KNOWN BUG: This only knows how to face lines without normal or texture
  // indices. It will probably crash or give garbage on anything else.
  template <typename DerivedV, typename DerivedF>
  IGL_INLINE bool readOBJ(
                          const std::string str,
                          Eigen::PlainObjectBase<DerivedV>& V,
                          Eigen::PlainObjectBase<DerivedF>& F);

}

#ifdef IGL_HEADER_ONLY
#  include "readOBJ.cpp"
#endif

#endif
