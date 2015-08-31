// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exterior_element.h"
#include "vertex_triangle_adjacency.h"
#include <iostream>

template <
     typename DerivedV,
     typename DerivedF,
     typename DerivedI,
     typename IndexType,
     typename DerivedA
     >
IGL_INLINE void igl::exterior_vertex(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::PlainObjectBase<DerivedF> & F,
        const Eigen::PlainObjectBase<DerivedI> & I,
        IndexType & v_index,
        Eigen::PlainObjectBase<DerivedA> & A) {
    // Algorithm: 
    //    Find an exterior vertex (i.e. vertex reachable from infinity)
    //    Return the vertex with the largest X value.
    //    If there is a tie, pick the one with largest Y value.
    //    If there is still a tie, pick the one with the largest Z value.
    //    If there is still a tie, then there are duplicated vertices within the
    //    mesh, which violates the precondition.
    const size_t INVALID = std::numeric_limits<size_t>::max();
    const size_t num_selected_faces = I.rows();
    std::vector<size_t> candidate_faces;
    size_t exterior_vid = INVALID;
    typename DerivedV::Scalar exterior_val = 0;
    for (size_t i=0; i<num_selected_faces; i++) {
        size_t f = I[i];
        for (size_t j=0; j<3; j++) {
            auto v = F(f, j);
            auto vx = V(v, 0);
            if (exterior_vid == INVALID || vx > exterior_val) {
                exterior_val = vx;
                exterior_vid = v;
                candidate_faces = {f};
            } else if (v == exterior_vid) {
                candidate_faces.push_back(f);
            } else if (vx == exterior_val) {
                // Break tie.
                auto vy = V(v,1);
                auto vz = V(v, 2);
                auto exterior_y = V(exterior_vid, 1);
                auto exterior_z = V(exterior_vid, 2);
                assert(!(vy == exterior_y && vz == exterior_z));
                bool replace = (vy > exterior_y) ||
                    ((vy == exterior_y) && (vz > exterior_z));
                if (replace) {
                    exterior_val = vx;
                    exterior_vid = v;
                    candidate_faces = {f};
                }
            }
        }
    }

    assert(exterior_vid != INVALID);
    assert(candidate_faces.size() > 0);
    v_index = exterior_vid;
    A.resize(candidate_faces.size());
    std::copy(candidate_faces.begin(), candidate_faces.end(), A.data());
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedI,
    typename IndexType,
    typename DerivedA
    >
IGL_INLINE void igl::exterior_edge(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::PlainObjectBase<DerivedF> & F,
        const Eigen::PlainObjectBase<DerivedI> & I,
        IndexType & v1,
        IndexType & v2,
        Eigen::PlainObjectBase<DerivedA> & A) {
    // Algorithm:
    //    Find an exterior vertex first.
    //    Find the incident edge with largest slope when projected onto XY plane.
    //    If there is still a tie, break it using the projected splot onto ZX plane.
    //    If there is still a tie, then there are multiple overlapping edges,
    //    which violates the precondition.
    typedef typename DerivedV::Scalar Scalar;
    typedef typename DerivedV::Index Index;
    typedef typename Eigen::Matrix<Scalar, 3, 1> ScalarArray3;
    typedef typename Eigen::Matrix<typename DerivedF::Scalar, 3, 1> IndexArray3;
    const size_t INVALID = std::numeric_limits<size_t>::max();

    Index exterior_vid;
    DerivedI candidate_faces;
    exterior_vertex(V, F, I, exterior_vid, candidate_faces);
    const ScalarArray3& exterior_v = V.row(exterior_vid);
    assert(candidate_faces.size() > 0);

    auto get_vertex_index = [&](const IndexArray3& f,
            Index vid) {
        if (f[0] == vid) return 0;
        if (f[1] == vid) return 1;
        if (f[2] == vid) return 2;
        assert(false);
    };

    Scalar exterior_slope_YX = 0;
    Scalar exterior_slope_ZX = 0;
    size_t exterior_opp_vid = INVALID;
    bool infinite_slope_detected = false;
    std::vector<Index> incident_faces;
    auto check_and_update_exterior_edge = [&](Index opp_vid, Index fid) {
        if (opp_vid == exterior_opp_vid) {
            incident_faces.push_back(fid);
            return;
        }

        const ScalarArray3 opp_v = V.row(opp_vid);
        if (!infinite_slope_detected && exterior_v[0] != opp_v[0]) {
            // Finite slope
            const ScalarArray3 diff = opp_v - exterior_v;
            const Scalar slope_YX = diff[1] / diff[0];
            const Scalar slope_ZX = diff[2] / diff[0];
            if (exterior_opp_vid == INVALID ||
                    slope_YX > exterior_slope_YX ||
                    (slope_YX == exterior_slope_YX &&
                     slope_ZX > exterior_slope_ZX)) {
                exterior_opp_vid = opp_vid;
                exterior_slope_YX = slope_YX;
                exterior_slope_ZX = slope_ZX;
                incident_faces = {fid};
            }
        } else if (!infinite_slope_detected) {
            // Infinite slope
            exterior_opp_vid = opp_vid;
            infinite_slope_detected = true;
            incident_faces = {fid};
        }
    };

    const auto num_candidate_faces = candidate_faces.size();
    for (size_t i=0; i<num_candidate_faces; i++) {
        const Index fid = candidate_faces[i];
        const IndexArray3& f = F.row(fid);
        size_t id = get_vertex_index(f, exterior_vid);
        Index next_vid = f[(id+1)%3];
        Index prev_vid = f[(id+2)%3];
        check_and_update_exterior_edge(next_vid, fid);
        check_and_update_exterior_edge(prev_vid, fid);
    }

    v1 = exterior_vid;
    v2 = exterior_opp_vid;
    A.resize(incident_faces.size());
    std::copy(incident_faces.begin(), incident_faces.end(), A.data());
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename DerivedI,
    typename IndexType
    >
IGL_INLINE void igl::exterior_facet(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::PlainObjectBase<DerivedF> & F,
        const Eigen::PlainObjectBase<DerivedN> & N,
        const Eigen::PlainObjectBase<DerivedI> & I,
        IndexType & f,
        bool & flipped) {
    // Algorithm:
    //    Compute the exterior edge.
    //    Find the facet with the largest absolute X normal component.
    //    If there is a tie, keep the one with positive X component.
    //    If there is still a tie, pick the face with the larger index.
    typedef typename DerivedV::Scalar Scalar;
    typedef typename DerivedV::Index Index;
    const size_t INVALID = std::numeric_limits<size_t>::max();

    Index v1,v2;
    DerivedI incident_faces;
    exterior_edge(V, F, I, v1, v2, incident_faces);
    assert(incident_faces.size() > 0);

    auto generic_fabs = [&](const Scalar& val) -> const Scalar {
        if (val >= 0) return val;
        else return -val;
    };

    Scalar max_nx = 0;
    size_t exterior_fid = INVALID;
    const size_t num_incident_faces = incident_faces.size();
    for (size_t i=0; i<num_incident_faces; i++) {
        const auto& fid = incident_faces[i];
        const Scalar nx = N(fid, 0);
        if (exterior_fid == INVALID) {
            max_nx = nx;
            exterior_fid = fid;
        } else {
            if (generic_fabs(nx) > generic_fabs(max_nx)) {
                max_nx = nx;
                exterior_fid = fid;
            } else if (nx == -max_nx && nx > 0) {
                max_nx = nx;
                exterior_fid = fid;
            } else if (nx == max_nx) {
                if ((max_nx >= 0 && exterior_fid < fid) ||
                    (max_nx <  0 && exterior_fid > fid)) {
                    max_nx = nx;
                    exterior_fid = fid;
                }
            }
        }
    }

    assert(exterior_fid != INVALID);
    f = exterior_fid;
    flipped = max_nx < 0;
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::exterior_facet<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::exterior_facet<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 1, -1, -1>, Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 1, -1, -1>, unsigned long>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 1, -1, -1> > const&, unsigned long&, bool&);
template void igl::exterior_facet<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::exterior_facet<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
#endif
