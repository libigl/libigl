// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#include "fast_find_intersections.h"
#include "AABB.h"
#include "tri_tri_intersect.h"



template <
  typename DerivedV1,
  typename DerivedF1,
  typename DerivedV2,
  typename DerivedF2,
  typename DerivedI,
  typename DerivedE>
IGL_INLINE void igl::fast_find_intersections(
  const Eigen::MatrixBase<DerivedV1>& V1,
  const Eigen::MatrixBase<DerivedF1>& F1,
  const Eigen::MatrixBase<DerivedV2>& V2,
  const Eigen::MatrixBase<DerivedF2>& F2,
        Eigen::PlainObjectBase<DerivedI>& intersect_pairs,
        Eigen::PlainObjectBase<DerivedE>& edges )
{
  using AABBTree=igl::AABB<DerivedV1,3>;
  AABBTree tree;
  tree.init(V1,F1);

  fast_find_intersections(tree, V1, F1, V2, F2, intersect_pairs, edges);
}


template <
  typename DerivedV1,
  typename DerivedF1,
  typename DerivedV2,
  typename DerivedF2,
  typename DerivedI,
  typename DerivedE>
IGL_INLINE void igl::fast_find_intersections(
  const AABB<DerivedV1,3>           & tree,
  const Eigen::MatrixBase<DerivedV1>& V1,
  const Eigen::MatrixBase<DerivedF1>& F1,
  const Eigen::MatrixBase<DerivedV2>& V2,
  const Eigen::MatrixBase<DerivedF2>& F2,
        Eigen::PlainObjectBase<DerivedI>& intersect_pairs,
        Eigen::PlainObjectBase<DerivedE>& edges )
{
    using AABBTree=igl::AABB<DerivedV1,3>;
    using Scalar=typename DerivedV1::Scalar;
    using BBOX=Eigen::AlignedBox<Scalar,3>;
    using VERTEX=Eigen::Matrix<typename DerivedE::Scalar,1,3,Eigen::RowMajor>;

    std::vector<VERTEX> _edges;
    std::vector<int>    _intersect_pairs;

    for(int i=0; i<F2.rows(); ++i)
    {
      BBOX tri_box;

      for(int j=0;j<3;++j)
        tri_box.extend( V2.row( F2(i,j) ).transpose() );
      
      // find leaf nodes containing intersecting tri_box
      // need to specify exact type instead of auto to allow recursion
      std::function<void(const AABBTree &,int)> check_intersect = 
        [&](const AABBTree &t,int d) -> void
      {
        if(t.m_primitive != -1) //check for the actual intersection //t.is_leaf()
        {
          bool coplanar=false;
          VERTEX edge1,edge2;

          if(igl::tri_tri_intersection_test_3d(
            V2.row(F2(i,0)),            V2.row(F2(i,1)),            V2.row(F2(i,2)),
            V1.row(F1(t.m_primitive,0)),V1.row(F1(t.m_primitive,1)),V1.row(F1(t.m_primitive,2)),
            coplanar,
            edge1,edge2))
          {
            if(!coplanar)
            {
              _intersect_pairs.push_back(t.m_primitive);
              _intersect_pairs.push_back(i);

              _edges.push_back(edge1);
              _edges.push_back(edge2);
            }
          }
        } else {
          if(t.m_box.intersects( tri_box )) 
          { // need to check all branches 
            check_intersect(*t.m_left, d+1);
            check_intersect(*t.m_right,d+1);
          }
        }
      };

      // run actual search
      check_intersect(tree, 0);
    }
    edges.resize(_edges.size(), 3);

    for(int i=0; i!=_edges.size(); ++i)
    {
      edges.row(i) = _edges[i];
    }

    intersect_pairs.resize(_intersect_pairs.size()/2,2);
    for(int i=0; i!=_intersect_pairs.size(); ++i)
    {
      intersect_pairs(i/2, i%2) = _intersect_pairs[i];
    }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::fast_find_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template void igl::fast_find_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);

#endif
