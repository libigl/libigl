// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#include "remesh_intersections.h"
#include "assign_scalar.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>

template <
typename DerivedV,
typename DerivedF,
typename Kernel,
typename DerivedVV,
typename DerivedFF,
typename DerivedJ,
typename DerivedIM>
IGL_INLINE void igl::copyleft::cgal::remesh_intersections(
        const Eigen::PlainObjectBase<DerivedV> & V,
        const Eigen::PlainObjectBase<DerivedF> & F,
        const std::vector<CGAL::Triangle_3<Kernel> > & T,
        const std::map<
        typename DerivedF::Index,
        std::pair<typename DerivedF::Index,
        std::vector<CGAL::Object> > > & offending,
        const std::map<
        std::pair<typename DerivedF::Index,typename DerivedF::Index>,
        std::vector<typename DerivedF::Index> > & /*edge2faces*/,
        Eigen::PlainObjectBase<DerivedVV> & VV,
        Eigen::PlainObjectBase<DerivedFF> & FF,
        Eigen::PlainObjectBase<DerivedJ> & J,
        Eigen::PlainObjectBase<DerivedIM> & IM) {

    typedef CGAL::Point_3<Kernel>    Point_3;
    typedef CGAL::Segment_3<Kernel>  Segment_3; 
    typedef CGAL::Triangle_3<Kernel> Triangle_3; 
    typedef CGAL::Plane_3<Kernel>    Plane_3;
    //typedef CGAL::Point_2<Kernel>    Point_2;
    //typedef CGAL::Segment_2<Kernel>  Segment_2; 
    //typedef CGAL::Triangle_2<Kernel> Triangle_2; 
    typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
    typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
    typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
    typedef CGAL::Exact_intersections_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> 
        CDT_2;
    typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;

    typedef typename DerivedF::Index Index;
    typedef std::pair<Index, Index> Edge;
    struct EdgeHash {
        size_t operator()(const Edge& e) const {
            return (e.first * 805306457) ^ (e.second * 201326611);
        }
    };
    typedef std::unordered_map<Edge, std::vector<Index>, EdgeHash > EdgeMap;

    auto normalize_plane_coeff = [](const Plane_3& P) {
        std::vector<typename Kernel::FT> coeffs = {
            P.a(), P.b(), P.c(), P.d()
        };
        const auto max_itr = std::max_element(coeffs.begin(), coeffs.end());
        const auto min_itr = std::min_element(coeffs.begin(), coeffs.end());
        typename Kernel::FT max_coeff;
        if (*max_itr < -1 * *min_itr) {
            max_coeff = *min_itr;
        } else {
            max_coeff = *max_itr;
        }
        std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                [&](const typename Kernel::FT& val)
                {return val / max_coeff; } );
        return coeffs;
    };

    auto plane_comp = [&](const Plane_3& p1, const Plane_3& p2) {
        const auto p1_coeffs = normalize_plane_coeff(p1);
        const auto p2_coeffs = normalize_plane_coeff(p2);
        if (p1_coeffs[0] != p2_coeffs[0])
            return p1_coeffs[0] < p2_coeffs[0];
        if (p1_coeffs[1] != p2_coeffs[1])
            return p1_coeffs[1] < p2_coeffs[1];
        if (p1_coeffs[2] != p2_coeffs[2])
            return p1_coeffs[2] < p2_coeffs[2];
        if (p1_coeffs[3] != p2_coeffs[3])
            return p1_coeffs[3] < p2_coeffs[3];
        return false;
    };
    std::map<Plane_3, std::vector<Index>, decltype(plane_comp)>
        unique_planes(plane_comp);

    const size_t num_faces = F.rows();
    const size_t num_base_vertices = V.rows();
    assert(num_faces == T.size());
    std::vector<bool> is_offending(num_faces, false);
    for (const auto itr : offending) {
        const auto& fid = itr.first;
        is_offending[fid] = true;

        Plane_3 key = T[fid].supporting_plane();
        assert(!key.is_degenerate());
        const auto jtr = unique_planes.find(key);
        if (jtr == unique_planes.end()) {
            unique_planes.insert({key, {fid}});
        } else {
            jtr->second.push_back(fid);
        }
    }

    std::vector<std::vector<Index> > resolved_faces;
    std::vector<Index> source_faces;
    std::vector<Point_3> new_vertices;
    EdgeMap edge_vertices;

    /**
     * Run constraint Delaunay triangulation on the plane.
     */
    auto run_delaunay_triangulation = [&](const Plane_3& P,
            const std::vector<Index>& involved_faces,
            std::vector<Point_3>& vertices,
            std::vector<std::vector<Index> >& faces) {
        CDT_plus_2 cdt;
        for (const auto& fid : involved_faces) {
            const auto itr = offending.find(fid);
            const auto& triangle = T[fid];
            cdt.insert_constraint(P.to_2d(triangle[0]), P.to_2d(triangle[1]));
            cdt.insert_constraint(P.to_2d(triangle[1]), P.to_2d(triangle[2]));
            cdt.insert_constraint(P.to_2d(triangle[2]), P.to_2d(triangle[0]));

            if (itr == offending.end()) continue;
            for (const auto& obj : itr->second.second) {
                if(const Segment_3 *iseg = CGAL::object_cast<Segment_3 >(&obj)) {
                    // Add segment constraint
                    cdt.insert_constraint(
                            P.to_2d(iseg->vertex(0)),P.to_2d(iseg->vertex(1)));
                }else if(const Point_3 *ipoint = CGAL::object_cast<Point_3 >(&obj)) {
                    // Add point
                    cdt.insert(P.to_2d(*ipoint));
                } else if(const Triangle_3 *itri = CGAL::object_cast<Triangle_3 >(&obj)) {
                    // Add 3 segment constraints
                    cdt.insert_constraint(
                            P.to_2d(itri->vertex(0)),P.to_2d(itri->vertex(1)));
                    cdt.insert_constraint(
                            P.to_2d(itri->vertex(1)),P.to_2d(itri->vertex(2)));
                    cdt.insert_constraint(
                            P.to_2d(itri->vertex(2)),P.to_2d(itri->vertex(0)));
                } else if(const std::vector<Point_3 > *polyp = 
                        CGAL::object_cast< std::vector<Point_3 > >(&obj)) {
                    //cerr<<REDRUM("Poly...")<<endl;
                    const std::vector<Point_3 > & poly = *polyp;
                    const Index m = poly.size();
                    assert(m>=2);
                    for(Index p = 0;p<m;p++)
                    {
                        const Index np = (p+1)%m;
                        cdt.insert_constraint(P.to_2d(poly[p]),P.to_2d(poly[np]));
                    }
                }else {
                    throw std::runtime_error("Unknown intersection object!");
                }
            }
        }
        std::map<typename CDT_plus_2::Vertex_handle,Index> v2i;
        size_t count=0;
        for (auto itr = cdt.finite_vertices_begin();
                itr != cdt.finite_vertices_end(); itr++) {
            vertices.push_back(P.to_3d(itr->point()));
            v2i[itr] = count;
            count++;
        }
        for (auto itr = cdt.finite_faces_begin();
                itr != cdt.finite_faces_end(); itr++) {
            faces.push_back( {
                    v2i[itr->vertex(0)],
                    v2i[itr->vertex(1)],
                    v2i[itr->vertex(2)] });
        }
    };

    /**
     * Given p on triangle indexed by ori_f, determine the index of p.
     */
    auto determine_point_index = [&](
            const Point_3& p, const size_t ori_f) -> Index {
        const auto& triangle = T[ori_f];
        const auto& f = F.row(ori_f).eval();

        // Check if p is one of the triangle corners.
        for (size_t i=0; i<3; i++) {
            if (p == triangle[i]) return f[i];
        }

        // Check if p is on one of the edges.
        for (size_t i=0; i<3; i++) {
            const Point_3 curr_corner = triangle[i];
            const Point_3 next_corner = triangle[(i+1)%3];
            const Segment_3 edge(curr_corner, next_corner);
            if (edge.has_on(p)) {
                const Index curr = f[i];
                const Index next = f[(i+1)%3];
                Edge key;
                key.first = curr<next?curr:next;
                key.second = curr<next?next:curr;
                auto itr = edge_vertices.find(key);
                if (itr == edge_vertices.end()) {
                    const Index index =
                        num_base_vertices + new_vertices.size();
                    edge_vertices.insert({key, {index}});
                    new_vertices.push_back(p);
                    return index;
                } else {
                    for (const auto vid : itr->second) {
                        if (p == new_vertices[vid - num_base_vertices]) {
                            return vid;
                        }
                    }
                    const size_t index = num_base_vertices + new_vertices.size();
                    new_vertices.push_back(p);
                    itr->second.push_back(index);
                    return index;
                }
            }
        }

        // p must be in the middle of the triangle.
        const size_t index = num_base_vertices + new_vertices.size();
        new_vertices.push_back(p);
        return index;
    };

    /**
     * Determine the vertex indices for each corner of each output triangle.
     */
    auto post_triangulation_process = [&](
            const std::vector<Point_3>& vertices,
            const std::vector<std::vector<Index> >& faces,
            const std::vector<Index>& involved_faces) {
        for (const auto& f : faces) {
            const Point_3& v0 = vertices[f[0]];
            const Point_3& v1 = vertices[f[1]];
            const Point_3& v2 = vertices[f[2]];
            Point_3 center(
                    (v0[0] + v1[0] + v2[0]) / 3.0,
                    (v0[1] + v1[1] + v2[1]) / 3.0,
                    (v0[2] + v1[2] + v2[2]) / 3.0);
            for (const auto& ori_f : involved_faces) {
                const auto& triangle = T[ori_f];
                const Plane_3 P = triangle.supporting_plane();
                if (triangle.has_on(center)) {
                    std::vector<Index> corners(3);
                    corners[0] = determine_point_index(v0, ori_f);
                    corners[1] = determine_point_index(v1, ori_f);
                    corners[2] = determine_point_index(v2, ori_f);
                    if (CGAL::orientation(
                                P.to_2d(v0), P.to_2d(v1), P.to_2d(v2))
                            == CGAL::RIGHT_TURN) {
                        std::swap(corners[0], corners[1]);
                    }
                    resolved_faces.emplace_back(corners);
                    source_faces.push_back(ori_f);
                }
            }
        }
    };

    // Process un-touched faces.
    for (size_t i=0; i<num_faces; i++) {
        if (!is_offending[i] && !T[i].is_degenerate()) {
            resolved_faces.push_back(
                    { F(i,0), F(i,1), F(i,2) } );
            source_faces.push_back(i);
        }
    }

    // Process self-intersecting faces.
    for (const auto itr : unique_planes) {
        Plane_3 P = itr.first;
        const auto& involved_faces = itr.second;

        std::vector<Point_3> vertices;
        std::vector<std::vector<Index> > faces;
        run_delaunay_triangulation(P, involved_faces, vertices, faces);
        post_triangulation_process(vertices, faces, involved_faces);
    }

    // Output resolved mesh.
    const size_t num_out_vertices = new_vertices.size() + num_base_vertices;
    VV.resize(num_out_vertices, 3);
    for (size_t i=0; i<num_base_vertices; i++) {
        assign_scalar(V(i,0), VV(i,0));
        assign_scalar(V(i,1), VV(i,1));
        assign_scalar(V(i,2), VV(i,2));
    }
    for (size_t i=num_base_vertices; i<num_out_vertices; i++) {
        assign_scalar(new_vertices[i-num_base_vertices][0], VV(i,0));
        assign_scalar(new_vertices[i-num_base_vertices][1], VV(i,1));
        assign_scalar(new_vertices[i-num_base_vertices][2], VV(i,2));
    }

    const size_t num_out_faces = resolved_faces.size();
    FF.resize(num_out_faces, 3);
    for (size_t i=0; i<num_out_faces; i++) {
        FF(i,0) = resolved_faces[i][0];
        FF(i,1) = resolved_faces[i][1];
        FF(i,2) = resolved_faces[i][2];
    }

    J.resize(num_out_faces);
    std::copy(source_faces.begin(), source_faces.end(), J.data());

    // Extract unique vertex indices.
    IM.resize(VV.rows(),1);
    std::map<Point_3,Index> vv2i;
    // Safe to check for duplicates using double for original vertices: if
    // incoming reps are different then the points are unique.
    for(Index v = 0;v<VV.rows();v++) {
        typename Kernel::FT p0,p1,p2;
        assign_scalar(VV(v,0),p0);
        assign_scalar(VV(v,1),p1);
        assign_scalar(VV(v,2),p2);
        const Point_3 p(p0,p1,p2);
        if(vv2i.count(p)==0) {
            vv2i[p] = v;
        }
        assert(vv2i.count(p) == 1);
        IM(v) = vv2i[p];
    }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif

