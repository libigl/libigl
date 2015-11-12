// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "outer_facet.h"
#include "outer_element.h"
#include "order_facets_around_edge.h"
#include <algorithm>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedI,
    typename IndexType
    >
IGL_INLINE void igl::copyleft::cgal::outer_facet(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::PlainObjectBase<DerivedF> & F,
        const Eigen::PlainObjectBase<DerivedI> & I,
        IndexType & f,
        bool & flipped) {

    // Algorithm:
    //
    //    1. Find an outer edge (s, d).
    //
    //    2. Order adjacent facets around this edge. Because the edge is an
    //    outer edge, there exists a plane passing through it such that all its
    //    adjacent facets lie on the same side. The implementation of
    //    order_facets_around_edge() will find a natural start facet such that
    //    The first and last facets according to this order are on the outside.
    //
    //    3. Because the vertex s is an outer vertex by construction (see
    //    implemnetation of outer_edge()). The first adjacent facet is facing
    //    outside (i.e. flipped=false) if it has positive X normal component.
    //    If it has zero normal component, it is facing outside if it contains
    //    directed edge (s, d).  

    //typedef typename DerivedV::Scalar Scalar;
    typedef typename DerivedV::Index Index;

    Index s,d;
    Eigen::Matrix<Index,Eigen::Dynamic,1> incident_faces;
    outer_edge(V, F, I, s, d, incident_faces);
    assert(incident_faces.size() > 0);

    auto convert_to_signed_index = [&](size_t fid) -> int{
        if ((F(fid, 0) == s && F(fid, 1) == d) ||
            (F(fid, 1) == s && F(fid, 2) == d) ||
            (F(fid, 2) == s && F(fid, 0) == d) ) {
            return int(fid+1) * -1;
        } else {
            return int(fid+1);
        }
    };

    auto signed_index_to_index = [&](int signed_id) -> size_t {
        return size_t(abs(signed_id) - 1);
    };

    std::vector<int> adj_faces(incident_faces.size());
    std::transform(incident_faces.data(),
            incident_faces.data() + incident_faces.size(),
            adj_faces.begin(),
            convert_to_signed_index);

    DerivedV pivot_point = V.row(s);
    pivot_point(0, 0) += 1.0;

    Eigen::VectorXi order;
    order_facets_around_edge(V, F, s, d, adj_faces, pivot_point, order);

    f = signed_index_to_index(adj_faces[order[0]]);
    flipped = adj_faces[order[0]] > 0;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
template void igl::copyleft::cgal::outer_facet<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::copyleft::cgal::outer_facet<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::copyleft::cgal::outer_facet<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::copyleft::cgal::outer_facet<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, unsigned long>(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, unsigned long&, bool&);
template void igl::copyleft::cgal::outer_facet<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
#endif
