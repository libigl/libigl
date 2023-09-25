// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_TETGEN_TETGENIO_TO_TETMESH_H
#define IGL_COPYLEFT_TETGEN_TETGENIO_TO_TETMESH_H
#include "../../igl_inline.h"

#ifndef TETLIBRARY
#define TETLIBRARY 
#endif
#include "tetgen.h" // Defined tetgenio, REAL
#include <vector>
#include <unordered_map>
#include <Eigen/Core>
namespace igl
{
  namespace copyleft
  {
    namespace tetgen
    {
      template <
        typename DerivedV, 
        typename DerivedT,
        typename DerivedF,
        typename DerivedTM,
        typename DerivedR,
        typename DerivedN,
        typename DerivedPT,
        typename DerivedFT>
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedT>& T,
        Eigen::PlainObjectBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedTM>& TM,
        Eigen::PlainObjectBase<DerivedR>& R,
        Eigen::PlainObjectBase<DerivedN>& N,
        Eigen::PlainObjectBase<DerivedPT>& PT,
        Eigen::PlainObjectBase<DerivedFT>& FT,
        int & num_regions);
    }
  }
}


#ifndef IGL_STATIC_LIBRARY
#  include "tetgenio_to_tetmesh.cpp"
#endif

#endif
