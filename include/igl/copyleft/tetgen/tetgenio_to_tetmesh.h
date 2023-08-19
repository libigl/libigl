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
      /// Extract a tetrahedral mesh from a tetgenio object
      ///
      /// @param[in] out tetgenio output object
      /// @param[out] V  #V by 3 vertex position list
      /// @param[out] T  #T by 4 list of tetrahedra indices into V
      /// @param[out] F  #F by 3 list of marked facets
      /// @param[out] R  #T list of region IDs for tetrahedra
      /// @param[out] N  #T by 2 list of neighbors for each tetrahedron
      /// @param[out] PT #V list of incident tetrahedron for each vertex
      /// @param[out] FT #F by 2 list of tetrahedra sharing each face 
      /// @param[out] nR number of regions in output mesh
      /// @return true on success, false on error
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
	std::vector<std::vector<REAL > > & V,
	std::vector<std::vector<int> > & T,
        std::vector<std::vector<int> > & F, 
	std::vector<std::vector<REAL> > & R,// region marks for tetrahedrons
	std::vector<std::vector<int > > &N, // neighborlist per tet
	std::vector<std::vector<int > >	&PT, // Point to tet list per point
	std::vector<std::vector<int > > &FT, // face to tet list
	size_t & nR); // number of regions    
      /// \overload
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
        std::vector<std::vector<REAL > > & V, 
        std::vector<std::vector<int> > & T,
        std::vector<std::vector<int> > & F);
      /// \overload
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
        std::vector<std::vector<REAL > > & V, 
        std::vector<std::vector<int> > & T);
      /// \overload
      template <typename DerivedV, typename DerivedT, typename DerivedF>
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedT>& T,
        Eigen::PlainObjectBase<DerivedF>& F);
      /// \overload
      template <typename DerivedV, typename DerivedT>
      IGL_INLINE bool tetgenio_to_tetmesh(
        const tetgenio & out,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedT>& T);

  }
  }
}


#ifndef IGL_STATIC_LIBRARY
#  include "tetgenio_to_tetmesh.cpp"
#endif

#endif
