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
#include "../../get_seconds.h"
#include "../../unique.h"

#include <vector>
#include <map>
#include <queue>
#include <unordered_map>
#include <iostream>

//#define REMESH_INTERSECTIONS_TIMING

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
        std::vector<
        std::pair<typename DerivedF::Index, CGAL::Object> > > & offending,
        const std::map<
        std::pair<typename DerivedF::Index,typename DerivedF::Index>,
        std::vector<typename DerivedF::Index> > & edge2faces,
        Eigen::PlainObjectBase<DerivedVV> & VV,
        Eigen::PlainObjectBase<DerivedFF> & FF,
        Eigen::PlainObjectBase<DerivedJ> & J,
        Eigen::PlainObjectBase<DerivedIM> & IM) {
    igl::copyleft::cgal::remesh_intersections(
            V, F, T, offending, edge2faces, false, VV, FF, J, IM);
}

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
        std::vector<
        std::pair<typename DerivedF::Index, CGAL::Object> > > & offending,
        const std::map<
        std::pair<typename DerivedF::Index,typename DerivedF::Index>,
        std::vector<typename DerivedF::Index> > & /*edge2faces*/,
        bool stitch_all,
        Eigen::PlainObjectBase<DerivedVV> & VV,
        Eigen::PlainObjectBase<DerivedFF> & FF,
        Eigen::PlainObjectBase<DerivedJ> & J,
        Eigen::PlainObjectBase<DerivedIM> & IM) {

#ifdef REMESH_INTERSECTIONS_TIMING
    const auto & tictoc = []() -> double
    {
      static double t_start = igl::get_seconds();
      double diff = igl::get_seconds()-t_start;
      t_start += diff;
      return diff;
    };
    const auto log_time = [&](const std::string& label) -> void {
      std::cout << "remesh_intersections." << label << ": "
          << tictoc() << std::endl;
    };
    tictoc();
#endif

    typedef CGAL::Point_3<Kernel>    Point_3;
    typedef CGAL::Segment_3<Kernel>  Segment_3; 
    typedef CGAL::Triangle_3<Kernel> Triangle_3; 
    typedef CGAL::Plane_3<Kernel>    Plane_3;
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

    const size_t num_faces = F.rows();
    const size_t num_base_vertices = V.rows();
    assert(num_faces == T.size());

    std::vector<bool> is_offending(num_faces, false);
    for (const auto itr : offending) {
        const auto& fid = itr.first;
        is_offending[fid] = true;
    }

    std::unordered_map<Index, std::vector<Index> > intersecting_and_coplanar;
    for (const auto itr : offending) {
        const auto& fi = itr.first;
        const auto P = T[fi].supporting_plane();
        assert(!P.is_degenerate());
        for (const auto jtr : itr.second) {
            const auto& fj = jtr.first;
            const auto& tj = T[fj];
            if (P.has_on(tj[0]) && P.has_on(tj[1]) && P.has_on(tj[2])) {
                auto loc = intersecting_and_coplanar.find(fi);
                if (loc == intersecting_and_coplanar.end()) {
                    intersecting_and_coplanar[fi] = {fj};
                } else {
                    loc->second.push_back(fj);
                }
            }
        }
    }
#ifdef REMESH_INTERSECTIONS_TIMING
    log_time("overlap_analysis");
#endif

    std::vector<std::vector<Index> > resolved_faces;
    std::vector<Index> source_faces;
    std::vector<Point_3> new_vertices;
    EdgeMap edge_vertices;

    /**
     * Run constraint Delaunay triangulation on the plane.
     */
    auto run_delaunay_triangulation = [&offending, &T](const Plane_3& P,
            const std::vector<Index>& involved_faces,
            std::vector<Point_3>& vertices,
            std::vector<std::vector<Index> >& faces) -> void {
        CDT_plus_2 cdt;
        for (const auto& fid : involved_faces) {
            const auto itr = offending.find(fid);
            const auto& triangle = T[fid];
            cdt.insert_constraint(P.to_2d(triangle[0]), P.to_2d(triangle[1]));
            cdt.insert_constraint(P.to_2d(triangle[1]), P.to_2d(triangle[2]));
            cdt.insert_constraint(P.to_2d(triangle[2]), P.to_2d(triangle[0]));

            if (itr == offending.end()) continue;
            for (const auto& index_obj : itr->second) {
                const auto& ofid = index_obj.first;
                const auto& obj = index_obj.second;
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
      if (stitch_all) {
        // No need to check if p shared by multiple triangles because all shared
        // vertices would be merged later on.
        const size_t index = num_base_vertices + new_vertices.size();
        new_vertices.push_back(p);
        return index;
      } else {
        // Stitching triangles according to input connectivity.
        // This step is potentially costly.
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
      }
    };

    /**
     * Determine the vertex indices for each corner of each output triangle.
     */
    auto post_triangulation_process = [&](
            const std::vector<Point_3>& vertices,
            const std::vector<std::vector<Index> >& faces,
            const std::vector<Index>& involved_faces) -> void {
      assert(involved_faces.size() > 0);
      for (const auto& f : faces) {
        const Point_3& v0 = vertices[f[0]];
        const Point_3& v1 = vertices[f[1]];
        const Point_3& v2 = vertices[f[2]];
        Point_3 center(
            (v0[0] + v1[0] + v2[0]) / 3.0,
            (v0[1] + v1[1] + v2[1]) / 3.0,
            (v0[2] + v1[2] + v2[2]) / 3.0);
        if (involved_faces.size() == 1) {
          // If only there is only one involved face, all sub-triangles must
          // belong to it and have the correct orientation.
          const auto& ori_f = involved_faces[0];
          std::vector<Index> corners(3);
          corners[0] = determine_point_index(v0, ori_f);
          corners[1] = determine_point_index(v1, ori_f);
          corners[2] = determine_point_index(v2, ori_f);
          resolved_faces.emplace_back(corners);
          source_faces.push_back(ori_f);
        } else {
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
    std::vector<bool> processed(num_faces, false);
    std::vector<std::pair<Plane_3, std::vector<Index> > > cdt_inputs;
    for (const auto itr : offending) {
        const auto fid = itr.first;
        if (processed[fid]) continue;
        processed[fid] = true;

        const auto loc = intersecting_and_coplanar.find(fid);
        std::vector<Index> involved_faces;
        if (loc == intersecting_and_coplanar.end()) {
            involved_faces.push_back(fid);
        } else {
            std::queue<Index> Q;
            Q.push(fid);
            while (!Q.empty()) {
                const auto index = Q.front();
                involved_faces.push_back(index);
                Q.pop();

                const auto overlapping_faces =
                    intersecting_and_coplanar.find(index);
                assert(overlapping_faces != intersecting_and_coplanar.end());

                for (const auto other_index : overlapping_faces->second) {
                    if (processed[other_index]) continue;
                    processed[other_index] = true;
                    Q.push(other_index);
                }
            }
        }

        Plane_3 P = T[fid].supporting_plane();
        cdt_inputs.emplace_back(P, involved_faces);
    }
#ifdef REMESH_INTERSECTIONS_TIMING
    log_time("preprocess");
#endif

    const size_t num_cdts = cdt_inputs.size();
    std::vector<std::vector<Point_3> > cdt_vertices(num_cdts);
    std::vector<std::vector<std::vector<Index> > > cdt_faces(num_cdts);

    const auto run_cdt = [&](const size_t first, const size_t last) {
      for (size_t i=first; i<last; i++) {
        auto& vertices = cdt_vertices[i];
        auto& faces = cdt_faces[i];
        const auto& P = cdt_inputs[i].first;
        const auto& involved_faces = cdt_inputs[i].second;
        run_delaunay_triangulation(P, involved_faces, vertices, faces);
      }
    };
    run_cdt(0, num_cdts);
#ifdef REMESH_INTERSECTIONS_TIMING
    log_time("cdt");
#endif

    for (size_t i=0; i<num_cdts; i++) {
      const auto& vertices = cdt_vertices[i];
      const auto& faces = cdt_faces[i];
      const auto& involved_faces = cdt_inputs[i].second;
      post_triangulation_process(vertices, faces, involved_faces);
    }
#ifdef REMESH_INTERSECTIONS_TIMING
    log_time("stitching");
#endif

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
    const size_t VV_size = VV.rows();
    IM.resize(VV_size,1);

    DerivedVV unique_vv;
    Eigen::VectorXi unique_to_vv, vv_to_unique;
    igl::unique_rows(VV, unique_vv, unique_to_vv, vv_to_unique);
    if (!stitch_all) {
      // Vertices with the same coordinates would be represented by one vertex.
      // The IM value of an vertex is the index of the representative vertex.
      for (Index v=0; v<VV_size; v++) {
        IM(v) = unique_to_vv[vv_to_unique[v]];
      }
    } else {
      // Screw IM and representative vertices.  Merge all vertices having the
      // same coordinates into a single vertex and set IM to identity map.
      VV = unique_vv;
      std::transform(FF.data(), FF.data() + FF.rows()*FF.cols(),
          FF.data(), [&vv_to_unique](const typename DerivedFF::Scalar& a)
          { return vv_to_unique[a]; });
      IM.setLinSpaced(unique_vv.rows(), 0, unique_vv.rows()-1);
    }
#ifdef REMESH_INTERSECTIONS_TIMING
    log_time("store_results");
#endif
}

#ifdef IGL_STATIC_LIBRARY
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index,CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int , -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1 > >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0,-1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1,3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1,3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1,3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::vector<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#undef IGL_STATIC_LIBRARY
#include "../../unique.cpp"
template void igl::unique_rows<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
