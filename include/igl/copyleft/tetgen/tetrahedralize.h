// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_TETGEN_TETRAHEDRALIZE_H
#define IGL_COPYLEFT_TETGEN_TETRAHEDRALIZE_H
#include "../../igl_inline.h"

#include <vector>
#include <string>
#include <Eigen/Core>
#ifndef TETLIBRARY
#define TETLIBRARY 
#endif
#include <tetgen.h> // Defined REAL

namespace igl
{
  namespace copyleft
  {
    namespace tetgen
    {
      /// Mesh the interior of a surface mesh (V,F) using tetgen
      ///
      /// @param[in] V  #V by 3 vertex position list
      /// @param[in] F  #F list of polygon face indices into V (0-indexed)
      /// @param[in] H  #H by 3 list of seed points inside holes
      /// @param[in] R  #R by 5 list of region attributes            
      /// @param[in] switches  string of tetgen options (See tetgen documentation) e.g.
      ///     "pq1.414a0.01" tries to mesh the interior of a given surface with
      ///       quality and area constraints
      ///     "" will mesh the convex hull constrained to pass through V (ignores F)
      /// @param[out] TV  #TV by 3 vertex position list
      /// @param[out] TT  #TT by 4 list of tet face indices
      /// @param[out] TF  #TF by 3 list of triangle face indices
      /// @param[out] TR  #TT list of region ID for each tetrahedron      
      /// @param[out] TN  #TT by 4 list of indices neighbors for each tetrahedron
      /// @param[out] PT  #TV list of incident tetrahedron for a vertex
      /// @param[out] FT  #TF by 2 list of tetrahedrons sharing a triface      
      /// @param[out] numRegions Number of regions in output mesh
      /// @return status:
      ///   0 success
      ///   1 tetgen threw exception
      ///   2 tetgen did not crash but could not create any tets (probably there are
      ///     holes, duplicate faces etc.)
      ///   -1 other error
      IGL_INLINE int tetrahedralize(
        const std::vector<std::vector<REAL> > &V, 
        const std::vector<std::vector<int> >  &F, 
        const std::vector<std::vector<REAL> > &H, 
        const std::vector<std::vector<REAL> > &R, 
        const std::string switches, 
        std::vector<std::vector<REAL > > & TV,
        std::vector<std::vector<int > >  & TT,
        std::vector<std::vector<int > >  & TF,
        std::vector<std::vector<REAL > > &TR,  
        std::vector<std::vector<int > > &TN, 
        std::vector<std::vector<int > > &PT, 
        std::vector<std::vector<int > > &FT, 
        size_t & numRegions);           
      /// \overload
      IGL_INLINE int tetrahedralize(
        const std::vector<std::vector<REAL > > & V, 
        const std::vector<std::vector<int> > & F, 
        const std::string switches,
        std::vector<std::vector<REAL > > & TV, 
        std::vector<std::vector<int > > & TT, 
        std::vector<std::vector<int> > & TF);
      /// \overload
      template <
        typename DerivedV, 
        typename DerivedF, 
        typename DerivedTV, 
        typename DerivedTT, 
        typename DerivedTF>
      IGL_INLINE int tetrahedralize(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const std::string switches,
        Eigen::PlainObjectBase<DerivedTV>& TV,
        Eigen::PlainObjectBase<DerivedTT>& TT,
        Eigen::PlainObjectBase<DerivedTF>& TF);
      /// \overload
      IGL_INLINE int tetrahedralize(
        const std::vector<std::vector<REAL > > & V, 
        const std::vector<std::vector<int> > & F, 
        const std::vector<int> & VM,
        const std::vector<int> & FM,
        const std::string switches,
        std::vector<std::vector<REAL > > & TV, 
        std::vector<std::vector<int > > & TT, 
        std::vector<std::vector<int> > & TF,
        std::vector<int> & TM);
      /// \overload
      template <
        typename DerivedV, 
        typename DerivedF, 
        typename DerivedVM,
        typename DerivedFM,
        typename DerivedTV, 
        typename DerivedTT, 
        typename DerivedTF, 
        typename DerivedTM>
      IGL_INLINE int tetrahedralize(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const Eigen::MatrixBase<DerivedVM>& VM,
        const Eigen::MatrixBase<DerivedFM>& FM,
        const std::string switches,
        Eigen::PlainObjectBase<DerivedTV>& TV,
        Eigen::PlainObjectBase<DerivedTT>& TT,
        Eigen::PlainObjectBase<DerivedTF>& TF,
        Eigen::PlainObjectBase<DerivedTM>& TM);
      /// \overload
      template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedH,
        typename DerivedR,
        typename DerivedTV,
        typename DerivedTT,
        typename DerivedTF,
        typename DerivedTR>      
      IGL_INLINE int tetrahedralize(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const Eigen::MatrixBase<DerivedH>& H,
        const Eigen::MatrixBase<DerivedR>& R,
        const std::string switches,
        Eigen::PlainObjectBase<DerivedTV>& TV,
        Eigen::PlainObjectBase<DerivedTT>& TT,
        Eigen::PlainObjectBase<DerivedTF>& TF,
        Eigen::PlainObjectBase<DerivedTR>& TR, 
        Eigen::PlainObjectBase<DerivedTT>& TN, 
        Eigen::PlainObjectBase<DerivedTT>& PT, 
        Eigen::PlainObjectBase<DerivedTT>& FT, 
        size_t & numRegions);            
   }
  }
}


#ifndef IGL_STATIC_LIBRARY
#  include "tetrahedralize.cpp"
#endif

#endif

