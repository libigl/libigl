// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef IGL_RELABEL_SMALL_IMMERSED_CELLS
#define IGL_RELABEL_SMALL_IMMERSED_CELLS

#include "../../igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      template<
        typename DerivedV,
        typename DerivedF,
        typename DerivedP,
        typename DerivedC,
        typename FT,
        typename DerivedW>
      IGL_INLINE void relabel_small_immersed_cells(
          const Eigen::PlainObjectBase<DerivedV>& V,
          const Eigen::PlainObjectBase<DerivedF>& F,
          const size_t num_patches,
          const Eigen::PlainObjectBase<DerivedP>& P,
          const size_t num_cells,
          const Eigen::PlainObjectBase<DerivedC>& C,
          const FT vol_threashold,
          Eigen::PlainObjectBase<DerivedW>& W);
    }
  }
}

#ifndef IGL_STATIC_LIBRARY
#  include "relabel_small_immersed_cells.cpp"
#endif
#endif
