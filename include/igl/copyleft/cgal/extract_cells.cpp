// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#include "extract_cells.h"
#include "../../extract_manifold_patches.h"
#include "../../facet_components.h"
#include "../../triangle_triangle_adjacency.h"
#include "../../unique_edge_map.h"
#include "../../get_seconds.h"
#include "closest_facet.h"
#include "order_facets_around_edge.h"
#include "outer_facet.h"

#include <vector>
#include <queue>

//#define EXTRACT_CELLS_DEBUG

template<
typename DerivedV,
typename DerivedF,
typename DerivedP,
typename DeriveduE,
typename uE2EType,
typename DerivedEMAP,
typename DerivedC >
IGL_INLINE size_t igl::copyleft::cgal::extract_cells_single_component(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const std::vector<std::vector<uE2EType> >& uE2E,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        Eigen::PlainObjectBase<DerivedC>& cells) {
    //typedef typename DerivedF::Scalar Index;
    const size_t num_faces = F.rows();
    auto edge_index_to_face_index = [&](size_t index) {
        return index % num_faces;
    };
    auto is_consistent = [&](const size_t fid, const size_t s, const size_t d) {
        if ((size_t)F(fid, 0) == s && (size_t)F(fid, 1) == d) return false;
        if ((size_t)F(fid, 1) == s && (size_t)F(fid, 2) == d) return false;
        if ((size_t)F(fid, 2) == s && (size_t)F(fid, 0) == d) return false;

        if ((size_t)F(fid, 0) == d && (size_t)F(fid, 1) == s) return true;
        if ((size_t)F(fid, 1) == d && (size_t)F(fid, 2) == s) return true;
        if ((size_t)F(fid, 2) == d && (size_t)F(fid, 0) == s) return true;
        throw "Invalid face!";
        return false;
    };

    const size_t num_unique_edges = uE.rows();
    const size_t num_patches = P.maxCoeff() + 1;

    std::vector<std::vector<size_t> > patch_edge_adj(num_patches);
    std::vector<Eigen::VectorXi> orders(num_unique_edges);
    std::vector<std::vector<bool> > orientations(num_unique_edges);
    for (size_t i=0; i<num_unique_edges; i++) {
        const size_t s = uE(i,0);
        const size_t d = uE(i,1);
        const auto adj_faces = uE2E[i];
        if (adj_faces.size() > 2) {
            std::vector<int> signed_adj_faces;
            for (auto ei : adj_faces) {
                const size_t fid = edge_index_to_face_index(ei);
                bool cons = is_consistent(fid, s, d);
                signed_adj_faces.push_back((fid+1)*(cons ? 1:-1));
            }
            auto& order = orders[i];
            igl::copyleft::cgal::order_facets_around_edge(
                    V, F, s, d, signed_adj_faces, order);
            auto& orientation = orientations[i];
            orientation.resize(order.size());
            std::transform(order.data(), order.data() + order.size(),
                    orientation.begin(), [&](int index) { return
                    signed_adj_faces[index] > 0; });
            std::transform(order.data(), order.data() + order.size(),
                    order.data(), [&](int index) { return adj_faces[index]; });

            for (auto ei : adj_faces) {
                const size_t fid = edge_index_to_face_index(ei);
                patch_edge_adj[P[fid]].push_back(ei);
            }
        }
    }

    const int INVALID = std::numeric_limits<int>::max();
    cells.resize(num_patches, 2);
    cells.setConstant(INVALID);

    auto peel_cell_bd = [&](
            size_t seed_patch_id,
            short seed_patch_side,
            size_t cell_idx) {
        typedef std::pair<size_t, short> PatchSide;
        std::queue<PatchSide> Q;
        Q.emplace(seed_patch_id, seed_patch_side);
        cells(seed_patch_id, seed_patch_side) = cell_idx;
        while (!Q.empty()) {
            auto entry = Q.front();
            Q.pop();
            const size_t patch_id = entry.first;
            const short side = entry.second;
            const auto& adj_edges = patch_edge_adj[patch_id];
            for (const auto& ei : adj_edges) {
                const size_t uei = EMAP[ei];
                const auto& order = orders[uei];
                const auto& orientation = orientations[uei];
                const size_t edge_valance = order.size();
                size_t curr_i = 0;
                for (curr_i=0; curr_i < edge_valance; curr_i++) {
                    if ((size_t)order[curr_i] == ei) break;
                }
                assert(curr_i < edge_valance);

                const bool cons = orientation[curr_i];
                size_t next;
                if (side == 0) {
                    next = (cons ? (curr_i + 1) :
                            (curr_i + edge_valance - 1)) % edge_valance;
                } else {
                    next = (cons ? (curr_i + edge_valance - 1) :
                            (curr_i + 1)) % edge_valance;
                }
                const size_t next_ei = order[next];
                const bool next_cons = orientation[next];
                const size_t next_patch_id = P[next_ei % num_faces];
                const short next_patch_side = (next_cons != cons) ?
                    side:abs(side-1);
                if (cells(next_patch_id, next_patch_side) == INVALID) {
                    Q.emplace(next_patch_id, next_patch_side);
                    cells(next_patch_id, next_patch_side) = cell_idx;
                } else {
                    assert(
                      (size_t)cells(next_patch_id, next_patch_side) == 
                      cell_idx);
                }
            }
        }
    };

    size_t count=0;
    for (size_t i=0; i<num_patches; i++) {
        if (cells(i, 0) == INVALID) {
            peel_cell_bd(i, 0, count);
            count++;
        }
        if (cells(i, 1) == INVALID) {
            peel_cell_bd(i, 1, count);
            count++;
        }
    }
    return count;
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename DerivedE,
    typename DeriveduE,
    typename uE2EType,
    typename DerivedEMAP,
    typename DerivedC >
IGL_INLINE size_t igl::copyleft::cgal::extract_cells(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const Eigen::PlainObjectBase<DerivedE>& E,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const std::vector<std::vector<uE2EType> >& uE2E,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        Eigen::PlainObjectBase<DerivedC>& cells) {
#ifdef EXTRACT_CELLS_DEBUG
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  tictoc();
#endif
    const size_t num_faces = F.rows();
    typedef typename DerivedF::Scalar Index;
    const size_t num_patches = P.maxCoeff()+1;

    DerivedC raw_cells;
    const size_t num_raw_cells =
        igl::copyleft::cgal::extract_cells_single_component(
                V, F, P, uE, uE2E, EMAP, raw_cells);
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "Extract single component cells: " << tictoc() << std::endl;
#endif

    std::vector<std::vector<std::vector<Index > > > TT,_1;
    igl::triangle_triangle_adjacency(E, EMAP, uE2E, false, TT, _1);
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "face adj: " << tictoc() << std::endl;
#endif

    Eigen::VectorXi C, counts;
    igl::facet_components(TT, C, counts);
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "face comp: " << tictoc() << std::endl;
#endif

    const size_t num_components = counts.size();
    std::vector<std::vector<size_t> > components(num_components);
    for (size_t i=0; i<num_faces; i++) {
        components[C[i]].push_back(i);
    }
    std::vector<Eigen::VectorXi> Is(num_components);
    for (size_t i=0; i<num_components; i++) {
        Is[i].resize(components[i].size());
        std::copy(components[i].begin(), components[i].end(),
                Is[i].data());
    }

    Eigen::VectorXi outer_facets(num_components);
    Eigen::VectorXi outer_facet_orientation(num_components);
    Eigen::VectorXi outer_cells(num_components);
    for (size_t i=0; i<num_components; i++) {
        bool flipped;
        igl::copyleft::cgal::outer_facet(V, F, Is[i], outer_facets[i], flipped);
        outer_facet_orientation[i] = flipped?1:0;
        outer_cells[i] = raw_cells(P[outer_facets[i]], outer_facet_orientation[i]);
    }
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "Per comp outer facet: " << tictoc() << std::endl;
#endif

    auto get_triangle_center = [&](const size_t fid) {
        return ((V.row(F(fid, 0)) + V.row(F(fid, 1)) + V.row(F(fid, 2)))
                /3.0).eval();
    };
    std::vector<std::vector<size_t> > nested_cells(num_raw_cells);
    std::vector<std::vector<size_t> > ambient_cells(num_raw_cells);
    std::vector<std::vector<size_t> > ambient_comps(num_components);
    if (num_components > 1) {
        DerivedV bbox_min(num_components, 3);
        DerivedV bbox_max(num_components, 3);
        bbox_min.rowwise() = V.colwise().maxCoeff().eval();
        bbox_max.rowwise() = V.colwise().minCoeff().eval();
        for (size_t i=0; i<num_faces; i++) {
          const auto comp_id = C[i];
          const auto& f = F.row(i);
          for (size_t j=0; j<3; j++) {
            bbox_min(comp_id, 0) = std::min(bbox_min(comp_id, 0), V(f[j], 0));
            bbox_min(comp_id, 1) = std::min(bbox_min(comp_id, 1), V(f[j], 1));
            bbox_min(comp_id, 2) = std::min(bbox_min(comp_id, 2), V(f[j], 2));
            bbox_max(comp_id, 0) = std::max(bbox_max(comp_id, 0), V(f[j], 0));
            bbox_max(comp_id, 1) = std::max(bbox_max(comp_id, 1), V(f[j], 1));
            bbox_max(comp_id, 2) = std::max(bbox_max(comp_id, 2), V(f[j], 2));
          }
        }
        auto bbox_intersects = [&](size_t comp_i, size_t comp_j) {
          return !(
              bbox_max(comp_i,0) < bbox_min(comp_j,0) ||
              bbox_max(comp_i,1) < bbox_min(comp_j,1) ||
              bbox_max(comp_i,2) < bbox_min(comp_j,2) ||
              bbox_max(comp_j,0) < bbox_min(comp_i,0) ||
              bbox_max(comp_j,1) < bbox_min(comp_i,1) ||
              bbox_max(comp_j,2) < bbox_min(comp_i,2));
        };

        for (size_t i=0; i<num_components; i++) {
            std::vector<size_t> candidate_comps;
            candidate_comps.reserve(num_components);
            for (size_t j=0; j<num_components; j++) {
                if (i == j) continue;
                if (bbox_intersects(i,j)) candidate_comps.push_back(j);
            }

            const size_t num_candidate_comps = candidate_comps.size();
            if (num_candidate_comps == 0) continue;

            DerivedV queries(num_candidate_comps, 3);
            for (size_t j=0; j<num_candidate_comps; j++) {
                const size_t index = candidate_comps[j];
                queries.row(j) = get_triangle_center(outer_facets[index]);
            }

            const auto& I = Is[i];
            Eigen::VectorXi closest_facets, closest_facet_orientations;
            igl::copyleft::cgal::closest_facet(V, F, I, queries,
                    uE2E, EMAP, closest_facets, closest_facet_orientations);

            for (size_t j=0; j<num_candidate_comps; j++) {
                const size_t index = candidate_comps[j];
                const size_t closest_patch = P[closest_facets[j]];
                const size_t closest_patch_side = closest_facet_orientations[j]
                    ? 0:1;
                const size_t ambient_cell = raw_cells(closest_patch,
                        closest_patch_side);
                if (ambient_cell != (size_t)outer_cells[i]) {
                    nested_cells[ambient_cell].push_back(outer_cells[index]);
                    ambient_cells[outer_cells[index]].push_back(ambient_cell);
                    ambient_comps[index].push_back(i);
                }
            }
        }
    }
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "Determine nested relaitonship: " << tictoc() << std::endl;
#endif

    const size_t INVALID = std::numeric_limits<size_t>::max();
    const size_t INFINITE_CELL = num_raw_cells;
    std::vector<size_t> embedded_cells(num_raw_cells, INVALID);
    for (size_t i=0; i<num_components; i++) {
        const size_t outer_cell = outer_cells[i];
        const auto& ambient_comps_i = ambient_comps[i];
        const auto& ambient_cells_i = ambient_cells[outer_cell];
        const size_t num_ambient_comps = ambient_comps_i.size();
        assert(num_ambient_comps == ambient_cells_i.size());
        if (num_ambient_comps > 0) {
            size_t embedded_comp = INVALID;
            size_t embedded_cell = INVALID;
            for (size_t j=0; j<num_ambient_comps; j++) {
                if (ambient_comps[ambient_comps_i[j]].size() ==
                        num_ambient_comps-1) {
                    embedded_comp = ambient_comps_i[j];
                    embedded_cell = ambient_cells_i[j];
                    break;
                }
            }
            assert(embedded_comp != INVALID);
            assert(embedded_cell != INVALID);
            embedded_cells[outer_cell] = embedded_cell;
        } else {
            embedded_cells[outer_cell] = INFINITE_CELL;
        }
    }
    for (size_t i=0; i<num_patches; i++) {
        if (embedded_cells[raw_cells(i,0)] != INVALID) {
            raw_cells(i,0) = embedded_cells[raw_cells(i, 0)];
        }
        if (embedded_cells[raw_cells(i,1)] != INVALID) {
            raw_cells(i,1) = embedded_cells[raw_cells(i, 1)];
        }
    }

    size_t count = 0;
    std::vector<size_t> mapped_indices(num_raw_cells+1, INVALID);
    for (size_t i=0; i<num_patches; i++) {
        const size_t old_positive_cell_id = raw_cells(i, 0);
        const size_t old_negative_cell_id = raw_cells(i, 1);
        size_t positive_cell_id, negative_cell_id;
        if (mapped_indices[old_positive_cell_id] == INVALID) {
            mapped_indices[old_positive_cell_id] = count;
            positive_cell_id = count;
            count++;
        } else {
            positive_cell_id = mapped_indices[old_positive_cell_id];
        }
        if (mapped_indices[old_negative_cell_id] == INVALID) {
            mapped_indices[old_negative_cell_id] = count;
            negative_cell_id = count;
            count++;
        } else {
            negative_cell_id = mapped_indices[old_negative_cell_id];
        }
        raw_cells(i, 0) = positive_cell_id;
        raw_cells(i, 1) = negative_cell_id;
    }
    cells = raw_cells;
#ifdef EXTRACT_CELLS_DEBUG
    std::cout << "Finalize and output: " << tictoc() << std::endl;
#endif
    return count;
}

template<
typename DerivedV,
typename DerivedF,
typename DerivedC >
IGL_INLINE size_t igl::copyleft::cgal::extract_cells(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedC>& cells) {
    const size_t num_faces = F.rows();
    //typedef typename DerivedF::Scalar Index;

    Eigen::MatrixXi E, uE;
    Eigen::VectorXi EMAP;
    std::vector<std::vector<size_t> > uE2E;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);

    Eigen::VectorXi P;
    //const size_t num_patches = 
      igl::extract_manifold_patches(F, EMAP, uE2E, P);

    DerivedC per_patch_cells;
    const size_t num_cells =
        igl::copyleft::cgal::extract_cells(V, F, P, E, uE, uE2E, EMAP, per_patch_cells);

    cells.resize(num_faces, 2);
    for (size_t i=0; i<num_faces; i++) {
        cells.row(i) = per_patch_cells.row(P[i]);
    }
    return num_cells;
}

#ifdef IGL_STATIC_LIBRARY
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
template unsigned long igl::copyleft::cgal::extract_cells<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, unsigned long, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
