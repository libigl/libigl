// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Qingnan Zhou <qnzhou@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "delaunay_triangulation.h"
#include "../../unique_edge_map.h"
#include "../../flip_edge.h"
#include "ExactPredicate.h"
#include "lexicographic_triangulation.h"

#include <vector>
#include <sstream>

template<
  typename DerivedV,
  typename DerivedF>
IGL_INLINE void igl::copyleft::cgal::delaunay_triangulation(
    const Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedF>& F)
{
  typedef typename DerivedF::Scalar Index;
  typedef typename DerivedV::Scalar Scalar;
  typedef typename igl::copyleft::cgal::ExactPredicate<Scalar> Predicate;
  igl::copyleft::cgal::lexicographic_triangulation(V, F);
  const size_t num_faces = F.rows();
  assert(F.cols() == 3);

  Eigen::MatrixXi E;
  Eigen::MatrixXi uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<Index> > uE2E;
  igl::unique_edge_map(F, E, uE, EMAP, uE2E);

  auto is_delaunay = [&V,&F,&uE2E,num_faces](size_t uei) {
    auto& half_edges = uE2E[uei];
    if (half_edges.size() != 2) {
      throw "Cannot flip non-manifold or boundary edge";
    }

    const size_t f1 = half_edges[0] % num_faces;
    const size_t f2 = half_edges[1] % num_faces;
    const size_t c1 = half_edges[0] / num_faces;
    const size_t c2 = half_edges[1] / num_faces;
    assert(c1 < 3);
    assert(c2 < 3);
    assert(f1 != f2);
    const size_t v1 = F(f1, (c1+1)%3);
    const size_t v2 = F(f1, (c1+2)%3);
    const size_t v4 = F(f1, c1);
    const size_t v3 = F(f2, c2);
    const Scalar p1[] = {V(v1, 0), V(v1, 1)};
    const Scalar p2[] = {V(v2, 0), V(v2, 1)};
    const Scalar p3[] = {V(v3, 0), V(v3, 1)};
    const Scalar p4[] = {V(v4, 0), V(v4, 1)};
    auto orientation = Predicate::incircle(p1, p2, p4, p3);
    return orientation <= 0;
  };

  bool all_delaunay = false;
  while(!all_delaunay) {
    all_delaunay = true;
    for (size_t i=0; i<uE2E.size(); i++) {
      if (uE2E[i].size() == 2) {
        if (!is_delaunay(i)) {
          all_delaunay = false;
          flip_edge(F, E, uE, EMAP, uE2E, i);
        }
      }
    }
  }
}

