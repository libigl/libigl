// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "intrinsic_delaunay_triangulation.h"
#include "is_intrinsic_delaunay.h"
#include "tan_half_angle.h"
#include "unique_edge_map.h"
#include "flip_edge.h"
#include "EPS.h"
#include "matlab_format.h"
#include <iostream>
#include <queue>
#include <map>

template <
  typename Derivedl_in,
  typename DerivedF_in,
  typename Derivedl,
  typename DerivedF>
IGL_INLINE void igl::intrinsic_delaunay_triangulation(
  const Eigen::MatrixBase<Derivedl_in> & l_in,
  const Eigen::MatrixBase<DerivedF_in> & F_in,
  Eigen::PlainObjectBase<Derivedl> & l,
  Eigen::PlainObjectBase<DerivedF> & F)
{
  // We're going to work in place
  l = l_in;
  F = F_in;

  typedef Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,2> MatrixX2I;
  typedef Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> VectorXI;
  MatrixX2I E,uE;
  VectorXI EMAP;
  std::vector<std::vector<typename DerivedF::Scalar> > uE2E;
  igl::unique_edge_map(F, E, uE, EMAP, uE2E);
  typedef typename DerivedF::Scalar Index;
  typedef typename Derivedl::Scalar Scalar;
  const Index num_faces = F.rows();

  // Does edge (a,b) exist in the edges of all faces incident on
  // existing unique edge uei.
  //
  // Inputs:
  //   a  1st end-point of query edge
  //   b  2nd end-point of query edge
  //   uei  index into uE/uE2E of unique edge
  //   uE2E  map from unique edges to half-edges (see unique_edge_map)
  //   E  #F*3 by 2 list of half-edges
  //
  std::vector<Index> face_queue;
  face_queue.reserve(32);
  std::vector<Index> pushed;
  // 32 is faster than 8
  pushed.reserve(32);
  const auto edge_exists_near = 
    [&](const Index & a,const Index & b,const Index & uei)->bool
    {
      face_queue.clear();
      pushed.clear();
      assert(a!=b);
      // Not handling case where (a,b) is edge of face incident on uei
      // since this can't happen for edge-flipping.
      assert(a!=uE(uei,0));
      assert(a!=uE(uei,1));
      assert(b!=uE(uei,0));
      assert(b!=uE(uei,1));
      // starting with the (2) faces incident on e, consider all faces
      // incident on edges containing either a or b.
      //
      // face_queue  Queue containing faces incident on exactly one of a/b
      // Using a vector seems mildly faster
      const Index f1 = uE2E[uei][0]%num_faces;
      const Index f2 = uE2E[uei][1]%num_faces;
      // map is faster than unordered_map here, and vector + brute force
      // is_member check is even faster
      face_queue.push_back(f1);
      pushed.push_back(f1);
      face_queue.push_back(f2);
      pushed.push_back(f2);
      while(!face_queue.empty())
      {
        const Index f = face_queue.back();
        face_queue.pop_back();
        // consider each edge of this face
        for(int c = 0;c<3;c++)
        {
          // Unique edge id
          const Index uec = EMAP(c*num_faces+f);
          const Index s = uE(uec,0);
          const Index d = uE(uec,1);
          const bool ona = s == a || d == a;
          const bool onb = s == b || d == b;
          // Is this the edge we're looking for?
          if(ona && onb)
          {
            return true;
          }
          // not incident on either?
          if(!ona && !onb)
          {
            continue;
          }
          // loop over all incident half-edges
          for(const auto & he : uE2E[uec])
          {
            // face of this he
            const Index fhe = he%num_faces;
            bool already_pushed = false;
            for(const auto & fp : pushed)
            {
              if(fp == fhe)
              {
                already_pushed = true;
                break;
              }
            }
            if(!already_pushed)
            {
              pushed.push_back(fhe);
              face_queue.push_back(fhe);
            }
          }
        }
      }
      return false;
    };

  // Vector is faster than queue...
  std::vector<Index> Q;
  Q.reserve(uE2E.size());
  for (size_t uei=0; uei<uE2E.size(); uei++) 
  {
    Q.push_back(uei);
  }

  while(!Q.empty())
  {
    const Index uei = Q.back();
    Q.pop_back();
    if (uE2E[uei].size() == 2) 
    {
      if(!is_intrinsic_delaunay(l,F,uE2E,uei)) 
      {
        // update l just before flipping edge
        //      .        //
        //     /|\       //
        //   a/ | \d     //
        //   /  e  \     //
        //  /   |   \    //
        // .----|-f--.   //
        //  \   |   /    //
        //   \  |  /     //
        //   b\α|δ/c     //
        //     \|/       //
        //      .        //
        // Annotated from flip_edge:
        // Edge to flip [v1,v2] --> [v3,v4]
        // Before:
        // F(f1,:) = [v1,v2,v4] // in some cyclic order
        // F(f2,:) = [v1,v3,v2] // in some cyclic order
        // After: 
        // F(f1,:) = [v1,v3,v4] // in *this* order 
        // F(f2,:) = [v2,v4,v3] // in *this* order
        //
        //          v1                 v1
        //          /|\                / \
        //        c/ | \b            c/f1 \b
        //     v3 /f2|f1\ v4  =>  v3 /__f__\ v4
        //        \  e  /            \ f2  /
        //        d\ | /a            d\   /a
        //          \|/                \ /
        //          v2                 v2
        //
        // Compute intrinsic length of oppposite edge
        assert(uE2E[uei].size() == 2 && "edge should have 2 incident faces");
        const Index f1 = uE2E[uei][0]%num_faces;
        const Index f2 = uE2E[uei][1]%num_faces;
        const Index c1 = uE2E[uei][0]/num_faces;
        const Index c2 = uE2E[uei][1]/num_faces;
        assert(c1 < 3);
        assert(c2 < 3);
        assert(f1 != f2);
        const Index v1 = F(f1, (c1+1)%3);
        const Index v2 = F(f1, (c1+2)%3);
        const Index v4 = F(f1, c1);
        const Index v3 = F(f2, c2);
        assert(F(f2, (c2+2)%3) == v1);
        assert(F(f2, (c2+1)%3) == v2);
        // From gptoolbox/mesh/flip_edge.m
        // "If edge-after-flip already  exists then this will create a non-manifold
        // edge"
        // Yes, this can happen: e.g., an edge of a tetrahedron."
        // "If two edges will be the same edge after flip then this will create a
        // non-manifold edge."
        // I dont' think this can happen if we flip one at a time. gptoolbox
        // flips in parallel.

        // Over 50% of the time is spent doing this check...
        bool flippable = !edge_exists_near(v3,v4,uei);
        if(flippable)
        {
          assert( std::abs(l(f1,c1)-l(f2,c2)) < igl::EPS<Scalar>() );
          const Scalar e = l(f1,c1);
          const Scalar a = l(f1,(c1+1)%3);
          const Scalar b = l(f1,(c1+2)%3);
          const Scalar c = l(f2,(c2+1)%3);
          const Scalar d = l(f2,(c2+2)%3);
          // tan(α/2)
          const Scalar tan_a_2= tan_half_angle(a,b,e);
          // tan(δ/2)
          const Scalar tan_d_2 = tan_half_angle(d,e,c);
          // tan((α+δ)/2)
          const Scalar tan_a_d_2 = (tan_a_2 + tan_d_2)/(1.0-tan_a_2*tan_d_2);
          // cos(α+δ)
          const Scalar cos_a_d = 
            (1.0 - tan_a_d_2*tan_a_d_2)/(1.0+tan_a_d_2*tan_a_d_2);
          const Scalar f = sqrt(b*b + c*c - 2.0*b*c*cos_a_d);
          l(f1,0) = f;
          l(f1,1) = b;
          l(f1,2) = c;
          l(f2,0) = f;
          l(f2,1) = d;
          l(f2,2) = a;
          flip_edge(F, E, uE, EMAP, uE2E, uei);
          // append neighbors to back
          const size_t e_24 = f1 + ((c1 + 1) % 3) * num_faces;
          const size_t e_41 = f1 + ((c1 + 2) % 3) * num_faces;
          const size_t e_13 = f2 + ((c2 + 1) % 3) * num_faces;
          const size_t e_32 = f2 + ((c2 + 2) % 3) * num_faces;
          const size_t ue_24 = EMAP(e_24);
          const size_t ue_41 = EMAP(e_41);
          const size_t ue_13 = EMAP(e_13);
          const size_t ue_32 = EMAP(e_32);
          Q.push_back(ue_24);
          Q.push_back(ue_41);
          Q.push_back(ue_13);
          Q.push_back(ue_32);
        }
      }
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::intrinsic_delaunay_triangulation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
