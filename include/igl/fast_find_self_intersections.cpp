// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#include "fast_find_self_intersections.h"
#include "AABB.h"
#include "tri_tri_intersect.h"

namespace igl{
namespace internal {

    // helper function, to check if two faces have shared vertices
    template <typename Derived>
    bool adjacent_faces(const Eigen::MatrixBase<Derived>& A,
                        const Eigen::MatrixBase<Derived>& B)
    {
        for(int i=0;i<3;++i)
          for(int j=0;j<3;j++)
          {
            if(A(i)==B(j))
              return true;
          }
        return false;
    }

}
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedI>
IGL_INLINE bool igl::fast_find_self_intersections(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedI>& intersect)
{
    using Scalar=typename DerivedV::Scalar;
    using BBOX=Eigen::AlignedBox<Scalar,3>;
    using AABBTree=igl::AABB<DerivedV,3>;
    AABBTree tree;

    tree.init(V,F);
    bool _intersects=false;

    intersect.resize(F.rows(),1);
    intersect.setConstant(0);

    for(int i=0; i<F.rows(); ++i)
    {
      if( intersect(i) ) continue;

      BBOX tri_box;

      for(int j=0;j<3;++j)
        tri_box.extend( V.row( F(i,j) ).transpose() );
      
      // find leaf nodes containing intersecting tri_box
      // need to declare full type to enable recursion
      std::function<bool(const AABBTree &,int)> check_intersect = 
        [&](const AABBTree &t,int d) -> bool
      {
        if(t.m_primitive != -1) //check for the actual intersection (is_leaf)
        {
          if(t.m_primitive==i) //itself
              return false;
          if(igl::internal::adjacent_faces(F.row(i), F.row(t.m_primitive)) )
              return false;

          bool coplanar=false;
          Eigen::Matrix<Scalar,1,3,Eigen::RowMajor> edge1,edge2;

          if(igl::tri_tri_intersection_test_3d(
            V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),
            V.row(F(t.m_primitive,0)),V.row(F(t.m_primitive,1)),V.row(F(t.m_primitive,2)),
            coplanar,
            edge1,edge2))
          {
            if(!coplanar)
            {
              intersect(i)=1;
              intersect(t.m_primitive)=1;
              return true;
            }
          }
        } else {
          if(t.m_box.intersects(tri_box)) {
            // need to check both subtrees
            bool r1=check_intersect(*t.m_left ,d+1);
            bool r2=check_intersect(*t.m_right,d+1);
            return r1||r2;
          }
        }
        return false;
      };

      bool r=check_intersect(tree,0);
      _intersects = _intersects || r;
    }
    return _intersects;
}


template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedI,
  typename DerivedE>
IGL_INLINE bool igl::fast_find_self_intersections(
  const Eigen::MatrixBase<DerivedV>& V,
  const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedI>& intersect,
        Eigen::PlainObjectBase<DerivedE>& edges )
{
    using Scalar=typename DerivedV::Scalar;
    using BBOX=Eigen::AlignedBox<Scalar,3>;
    using AABBTree=igl::AABB<DerivedV,3>;
    using EDGE=Eigen::Matrix<typename DerivedE::Scalar,1,3,Eigen::RowMajor>;

    std::vector<EDGE> _edges;
    AABBTree tree;

    tree.init(V,F);
    bool _intersects=false;

    intersect.resize(F.rows(),1);
    intersect.setConstant(0);

    for(int i=0; i<F.rows(); ++i)
    {
      if( intersect(i) ) continue;

      BBOX tri_box;

      for(int j=0;j<3;++j)
        tri_box.extend( V.row( F(i,j) ).transpose() );
      
      // find leaf nodes containing intersecting tri_box
      std::function<bool(const AABBTree &,int)> check_intersect = 
        [&](const AABBTree &t,int d) -> bool
      {
        if(t.m_primitive != -1) //check for the actual intersection //t.is_leaf()
        {
          if(t.m_primitive==i) //itself
              return false;
          if(igl::internal::adjacent_faces(F.row(i), F.row(t.m_primitive)) )
              return false;

          bool coplanar=false;
          EDGE edge1,edge2;

          if(igl::tri_tri_intersection_test_3d(
            V.row(F(i,0)),            V.row(F(i,1)),            V.row(F(i,2)),
            V.row(F(t.m_primitive,0)),V.row(F(t.m_primitive,1)),V.row(F(t.m_primitive,2)),
            coplanar,
            edge1,edge2))
          {
            if(!coplanar)
            {
              intersect(i)=1;
              intersect(t.m_primitive)=1;
              _edges.push_back(edge1);
              _edges.push_back(edge2);
              return true;
            }
          }
        } else {
          if(t.m_box.intersects( tri_box ))
          {
            bool r1=check_intersect(*t.m_left, d+1);
            bool r2=check_intersect(*t.m_right,d+1);
            return r1||r2;
          }
        }
        return false;
      };

      // run actual search
      bool r=check_intersect(tree, 0);
      _intersects = _intersects || r;
    }

    edges.resize(_edges.size(),3);
    int i=0;

    for(auto e=std::begin(_edges);e!=std::end(_edges);++e,++i)
    {
      edges.row(i) = *e;
    }
    return _intersects;
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::fast_find_self_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::fast_find_self_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
