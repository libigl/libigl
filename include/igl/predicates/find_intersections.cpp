// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Vladimir S. FONOV <vladimir.fonov@gmail.com>
// Copyright (C) 2024 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/
#include "find_intersections.h"

#include "../AABB.h"
#include "../triangle_triangle_intersect_shared_edge.h"
#include "../triangle_triangle_intersect_shared_vertex.h"
#include "triangle_triangle_intersect.h"
#include <stdio.h>

template <
  typename DerivedV1,
  typename DerivedF1,
  typename DerivedV2,
  typename DerivedF2,
  typename DerivedIF>
IGL_INLINE bool igl::predicates::find_intersections(
  const igl::AABB<DerivedV1,3> & tree,
  const Eigen::MatrixBase<DerivedV1> & V1,
  const Eigen::MatrixBase<DerivedF1> & F1,
  const Eigen::MatrixBase<DerivedV2> & V2,
  const Eigen::MatrixBase<DerivedF2> & F2,
  const bool first_only,
  Eigen::PlainObjectBase<DerivedIF> & IF)
{
  const bool detect_only = true;
  constexpr bool stinker = false;
  using AABBTree=igl::AABB<DerivedV1,3>;
  using Scalar=typename DerivedV1::Scalar;
  static_assert(
    std::is_same<Scalar,typename DerivedV2::Scalar>::value,
    "V1 and V2 must have the same scalar type");
  using RowVector3S=Eigen::Matrix<Scalar,1,3>;

  // Determine if V1,F1 and V2,F2 point to the same data
  const bool self_test = (&V1 == &V2) && (&F1 == &F2);
  if(stinker){ printf("%s\n",self_test?"🍎&(V1,F1) == 🍎&(V2,F2)":"🍎≠🍊"); }

  int num_if = 0;
  const auto append_intersection = [&IF,&num_if]( const int f1, const int f2)
  {
    if(num_if >= IF.rows())
    {
      IF.conservativeResize(2*IF.rows()+1,2);
    }
    IF.row(num_if++) << f1,f2;
  };

  // Returns corner in ith face opposite of shared edge; -1 otherwise
  const auto shared_edge = [&F1](const int f, const int g)->int
  {
    for(int c = 0;c<3;c++)
    {
      int s = F1(f,(c+1)%3);
      int d = F1(f,(c+2)%3);
      for(int e = 0;e<3;e++)
      {
        // Find in opposite direction on jth face
        if(F1(g,e) == d && F1(g,(e+1)%3) == s)
        {
          return c;
        }
      }
    }
    return -1;
  };
  // Returns corner of shared vertex in ith face; -1 otherwise
  const auto shared_vertex = [&F1](const int f, const int g, int & sf, int & sg)->bool
  {
    for(sf = 0;sf<3;sf++)
    {
      for(sg = 0;sg<3;sg++)
      {
        if(F1(g,sg) == F1(f,sf))
        {
          return true;
        }
      }
    }
    return false;
  };

  RowVector3S dummy;
  for(int f2 = 0; f2<F2.rows(); ++f2)
  {
    if(stinker){ printf("f2: %d\n",f2); }
    Eigen::AlignedBox<Scalar,3> box;
    box.extend( V2.row( F2(f2,0) ).transpose() );
    box.extend( V2.row( F2(f2,1) ).transpose() );
    box.extend( V2.row( F2(f2,2) ).transpose() );
    std::vector<const AABBTree*> candidates;
    tree.append_intersecting_leaves(box, candidates);
    for(const auto * candidate : candidates)
    {
      const int f1 = candidate->m_primitive;
      if(stinker){ printf("  f1: %d\n",f1); }
      bool found_intersection = false;
      bool yes_shared_verted = false;
      bool yes_shared_edge = false;
      if(self_test)
      {
        // Skip self-test and direction f2>=f1 (assumes by symmetry we'll find
        // the other direction since we're testing all pairs)
        if(f1 >= f2)
        {
          if(stinker){ printf("    ⏭\n"); }
          continue;
        }
        const int c = shared_edge(f1,f2);
        yes_shared_edge = c != -1;
        if(yes_shared_edge)
        {
          if(stinker){ printf("    ⚠️  shared edge\n"); }
          if(stinker)
          {
            printf("    %d: %d %d %d\n",f1,F1(f1,0),F1(f1,1),F1(f1,2));
            printf("    %d: %d %d %d\n",f2,F1(f2,0),F1(f2,1),F1(f2,2));
            printf("   edge: %d %d\n",F1(f1,(c+1)%3),F1(f1,(c+2)%3));
          }
          found_intersection = igl::triangle_triangle_intersect_shared_edge(
            V1,F1,f1,c,V1.row(F1(f1,c)),f2,1e-8);
          if(found_intersection)
          {
            append_intersection(f1,f2);
          }
        }else
        {
          int sf,sg;
          yes_shared_verted = shared_vertex(f1,f2,sf,sg);
          if(yes_shared_verted)
          {
            if(stinker){ printf("    ⚠️  shared vertex\n"); }
            // Just to be sure. c≠sf
            const int c = (sf+1)%3;
            assert(F1(f1,sf) == F1(f2,sg));
            found_intersection = igl::triangle_triangle_intersect_shared_vertex(
              V1,F1,f1,sf,c,V1.row(F1(f1,c)),f2,sg,1e-14);
            if(found_intersection && detect_only)
            {
              append_intersection(f1,f2);
            }
          }
          
        }
      }
      if(
        !self_test || 
        (!yes_shared_verted && !yes_shared_edge) || 
        (yes_shared_verted && found_intersection && !detect_only))
      {
        bool coplanar = false;
        RowVector3S v1,v2;
        const bool tt_found_intersection = 
          triangle_triangle_intersect(
            V2.row(F2(f2,0)).template head<3>().eval(),
            V2.row(F2(f2,1)).template head<3>().eval(),
            V2.row(F2(f2,2)).template head<3>().eval(),
            V1.row(F1(f1,0)).template head<3>().eval(),
            V1.row(F1(f1,1)).template head<3>().eval(),
            V1.row(F1(f1,2)).template head<3>().eval());
        if(found_intersection && !tt_found_intersection)
        {
          // We failed to find the edge. Mark it as an intersection but don't
          // include edge.
          append_intersection(f1,f2);
        }else if(tt_found_intersection)
        {
          found_intersection = true;
          append_intersection(f1,f2);
        }
      }
      if(stinker) { printf("    %s\n",found_intersection? "☠️":"❌"); }
      if(num_if && first_only) { break; }
    }
    if(num_if && first_only) { break; }
  }
  IF.conservativeResize(num_if,2);
  return IF.rows();
}

template <
  typename DerivedV1,
  typename DerivedF1,
  typename DerivedV2,
  typename DerivedF2,
  typename DerivedIF>
IGL_INLINE bool igl::predicates::find_intersections(
  const Eigen::MatrixBase<DerivedV1> & V1,
  const Eigen::MatrixBase<DerivedF1> & F1,
  const Eigen::MatrixBase<DerivedV2> & V2,
  const Eigen::MatrixBase<DerivedF2> & F2,
  const bool first_only,
  Eigen::PlainObjectBase<DerivedIF> & IF)
{
  igl::AABB<DerivedV1,3> tree1;
  tree1.init(V1,F1);
  return find_intersections(tree1,V1,F1,V2,F2,first_only,IF);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::predicates::find_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
