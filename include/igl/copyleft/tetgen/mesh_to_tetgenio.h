// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_TETGEN_MESH_TO_TETGENIO_H
#define IGL_COPYLEFT_TETGEN_MESH_TO_TETGENIO_H
#include "../../igl_inline.h"

#ifndef TETLIBRARY
#  define TETLIBRARY 
#endif
#include "tetgen.h" // Defined tetgenio, REAL
#include <vector>
#include <Eigen/Core>

namespace igl
{
  namespace copyleft
  {
    namespace tetgen
    {
      /// Load a vertex list and face list into a tetgenio object
      ///
      /// @param[in] V  #V by 3 vertex position list
      /// @param[in] F  #F list of polygon face indices into V (0-indexed)
      /// @param[out] in  tetgenio input object
      ///  @param[out] H  #H list of seed point inside each hole
      ///  @param[out] R  #R list of seed point inside each region	
      /// @return true on success, false on error
      IGL_INLINE bool mesh_to_tetgenio(
        const std::vector<std::vector<REAL> > & V,
	const std::vector<std::vector<int> > & F,
	const std::vector<std::vector<REAL > > & H, 
	const std::vector<std::vector<REAL > > & R, 
	tetgenio & in);	
      IGL_INLINE bool mesh_to_tetgenio(
        const std::vector<std::vector<REAL > > & V, 
        const std::vector<std::vector<int> > & F, 
        tetgenio & in);
      /// \overload
      template <typename DerivedV, typename DerivedF>
      IGL_INLINE bool mesh_to_tetgenio(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        tetgenio & in);
      /// \overload
      template <typename DerivedV, typename DerivedF, typename DerivedH, typename DerivedR>
      IGL_INLINE bool mesh_to_tetgenio(
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	const Eigen::PlainObjectBase<DerivedH>& H, 
	const Eigen::PlainObjectBase<DerivedR>& R, 
	tetgenio& in);	
    }
  }
}


#ifndef IGL_STATIC_LIBRARY
#  include "mesh_to_tetgenio.cpp"
#endif

#endif
