#include "propagate_winding_numbers.h"
#include "../extract_manifold_patches.h"
#include "../extract_non_manifold_edge_curves.h"
#include "../facet_components.h"
#include "../unique_edge_map.h"
#include "order_facets_around_edge.h"
#include "outer_facet.h"
#include "closest_facet.h"

#include <stdexcept>
#include <limits>
#include <vector>

namespace propagate_winding_numbers_helper {
    template<typename DerivedW >
    bool winding_number_assignment_is_consistent(
            const std::vector<Eigen::VectorXi>& orders,
            const std::vector<std::vector<bool> >& orientations,
            const Eigen::PlainObjectBase<DerivedW>& per_patch_winding_number) {
        const size_t num_edge_curves = orders.size();
        const size_t num_labels = per_patch_winding_number.cols() / 2;
        for (size_t i=0; i<num_edge_curves; i++) {
            const auto& order = orders[i];
            const auto& orientation = orientations[i];
            assert(order.size() == orientation.size());
            const size_t order_size = order.size();
            for (size_t j=0; j<order_size; j++) {
                const size_t curr = j;
                const size_t next = (j+1) % order_size;

                for (size_t k=0; k<num_labels; k++) {
                    // Retrieve the winding numbers of the current partition from
                    // the current patch and the next patch.  If the patches forms
                    // the boundary of a 3D volume, the winding number assignments
                    // should be consistent.
                    int curr_winding_number = per_patch_winding_number(
                            order[curr], k*2+(orientation[curr]? 0:1));
                    int next_winding_number = per_patch_winding_number(
                            order[next], k*2+(orientation[next]? 1:0));
                    if (curr_winding_number != next_winding_number) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DeriveduE,
    typename uE2EType,
    typename DerivedC,
    typename DerivedP,
    typename DerivedW >
IGL_INLINE bool igl::cgal::propagate_winding_numbers_single_component_patch_wise(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const std::vector<std::vector<uE2EType> >& uE2E,
        const Eigen::PlainObjectBase<DerivedC>& labels,
        const Eigen::PlainObjectBase<DerivedP>& P,
        const std::vector<std::vector<size_t> >& intersection_curves,
        Eigen::PlainObjectBase<DerivedW>& patch_W) {
    const size_t num_faces = F.rows();
    const size_t num_patches = P.maxCoeff() + 1;
    assert(labels.size() == num_patches);
    // Utility functions.
    auto edge_index_to_face_index = [&](size_t ei) { return ei % num_faces; };
    auto edge_index_to_corner_index = [&](size_t ei) { return ei / num_faces; };
    auto face_and_corner_index_to_edge_index = [&](size_t fi, size_t ci) {
        return ci*num_faces + fi;
    };
    auto get_edge_end_points = [&](size_t ei, size_t& s, size_t& d) {
        const size_t fi = edge_index_to_face_index(ei);
        const size_t ci = edge_index_to_corner_index(ei);
        s = F(fi, (ci+1)%3);
        d = F(fi, (ci+2)%3);
    };
    auto is_positively_orientated =
        [&](size_t fi, size_t s, size_t d) -> bool{
            const Eigen::Vector3i f = F.row(fi);
            if (f[0] == d && f[1] == s) return true;
            if (f[1] == d && f[2] == s) return true;
            if (f[2] == d && f[0] == s) return true;
            if (f[0] == s && f[1] == d) return false;
            if (f[1] == s && f[2] == d) return false;
            if (f[2] == s && f[0] == d) return false;
            throw std::runtime_error("Edge does not belong to face.");
            return false;
        };

    auto compute_signed_index = [&](size_t fi, size_t s, size_t d) -> int{
        return int(fi+1) * (is_positively_orientated(fi, s, d) ? 1:-1);
    };
    auto compute_unsigned_index = [&](int signed_idx) -> size_t {
        return abs(signed_idx) - 1;
    };

    // Order patches around each intersection curve.
    const size_t num_edge_curves = intersection_curves.size();
    std::vector<Eigen::VectorXi> orders(num_edge_curves);
    std::vector<std::vector<bool> > orientations(num_edge_curves);
    std::vector<std::vector<size_t> > patch_curve_adjacency(num_patches);
    for (size_t i=0; i<num_edge_curves; i++) {
        const auto& curve = intersection_curves[i];
        const size_t uei = curve[0];
        size_t s = uE(uei, 0);
        size_t d = uE(uei, 1);
        std::vector<int> adj_faces;
        for (size_t ei : uE2E[uei]) {
            const size_t fi = edge_index_to_face_index(ei);
            const size_t signed_fi =
                compute_signed_index(fi, s, d);
            adj_faces.push_back(signed_fi);
            patch_curve_adjacency[P[fi]].push_back(i);
        }

        auto& order = orders[i];
        igl::cgal::order_facets_around_edge(
                V, F, s, d, adj_faces, order);
        assert(order.minCoeff() == 0);
        assert(order.maxCoeff() == adj_faces.size() - 1);

        auto& orientation = orientations[i];
        orientation.resize(order.size());
        std::transform(order.data(), order.data()+order.size(),
                orientation.begin(), 
                [&](size_t index) { return adj_faces[index] > 0; });
        std::transform(order.data(), order.data()+order.size(),
                order.data(), [&](size_t index) {
                return P[compute_unsigned_index(adj_faces[index])];
                } );
    }

    // Propagate winding number from infinity.
    // Assuming infinity has winding number 0.
    const size_t num_labels = labels.maxCoeff() + 1;
    const int INVALID = std::numeric_limits<int>::max();
    patch_W.resize(num_patches, 2*num_labels);
    patch_W.setConstant(INVALID);

    size_t outer_facet_idx;
    bool outer_facet_is_flipped;
    Eigen::VectorXi face_indices(num_faces);
    face_indices.setLinSpaced(num_faces, 0, num_faces-1);
    igl::cgal::outer_facet(V, F, face_indices,
            outer_facet_idx, outer_facet_is_flipped);
    size_t outer_patch_idx = P[outer_facet_idx];
    size_t outer_patch_label = labels[outer_patch_idx];
    patch_W.row(outer_patch_idx).setZero();
    if (outer_facet_is_flipped) {
        patch_W(outer_patch_idx, outer_patch_label*2) = -1;
    } else {
        patch_W(outer_patch_idx, outer_patch_label*2+1) = 1;
    }

    auto winding_num_assigned = [&](size_t patch_idx) -> bool{
        return (patch_W.row(patch_idx).array() != INVALID).all();
    };

    std::queue<size_t> Q;
    Q.push(outer_patch_idx);
    while (!Q.empty()) {
        const size_t curr_patch_idx = Q.front();
        const size_t curr_patch_label = labels[curr_patch_idx];
        Q.pop();

        const auto& adj_curves = patch_curve_adjacency[curr_patch_idx];
        for (size_t curve_idx : adj_curves) {
            const auto& order = orders[curve_idx];
            const auto& orientation = orientations[curve_idx];
            const size_t num_adj_patches = order.size();
            assert(num_adj_patches == orientation.size());

            size_t curr_i = std::numeric_limits<size_t>::max();
            for (size_t i=0; i<num_adj_patches; i++) {
                if (order[i] == curr_patch_idx) {
                    curr_i = i;
                    break;
                }
            }
            assert(curr_i < num_adj_patches);
            const bool curr_ori = orientation[curr_i];

            const size_t next_i = curr_ori ? (curr_i+1) % num_adj_patches
                : (curr_i+num_adj_patches-1)%num_adj_patches;
            const size_t prev_i = curr_ori ?
                (curr_i+num_adj_patches-1)%num_adj_patches
                : (curr_i+1) % num_adj_patches;
            const size_t next_patch_idx = order[next_i];
            const size_t prev_patch_idx = order[prev_i];

            if (!winding_num_assigned(next_patch_idx)) {
                const bool next_ori = orientation[next_i];
                const bool next_cons = next_ori != curr_ori;
                const size_t next_patch_label = labels[next_patch_idx];
                for (size_t i=0; i<num_labels; i++) {
                    const int shared_winding_number =
                        patch_W(curr_patch_idx, i*2);

                    if (i == next_patch_label) {
                        // Truth table:
                        // curr_ori  next_ori  wind_# inc
                        // True      True      -1
                        // True      False      1
                        // False     True       1
                        // False     False     -1

                        patch_W(next_patch_idx, i*2+(next_cons ?0:1)) =
                            shared_winding_number;
                        patch_W(next_patch_idx, i*2+(next_cons ?1:0)) = 
                            shared_winding_number + (next_cons ? 1:-1);
                    } else {
                        patch_W(next_patch_idx, i*2  ) = shared_winding_number;
                        patch_W(next_patch_idx, i*2+1) = shared_winding_number;
                    }
                }
                Q.push(next_patch_idx);
            }
            if (!winding_num_assigned(prev_patch_idx)) {
                const bool prev_ori = orientation[prev_i];
                const bool prev_cons = prev_ori != curr_ori;
                const size_t prev_patch_label = labels[prev_patch_idx];

                for (size_t i=0; i<num_labels; i++) {
                    const int shared_winding_number =
                        patch_W(curr_patch_idx, i*2+1);

                    if (i == prev_patch_label) {
                        // Truth table:
                        // curr_ori  next_ori  wind_# inc
                        // True      True       1
                        // True      False     -1
                        // False     True      -1
                        // False     False      1

                        patch_W(prev_patch_idx, i*2+(prev_cons ?1:0)) =
                            shared_winding_number;
                        patch_W(prev_patch_idx, i*2+(prev_cons ?0:1)) =
                            shared_winding_number + (prev_cons ? -1:1);
                    } else {
                        patch_W(prev_patch_idx, i*2  ) = shared_winding_number; 
                        patch_W(prev_patch_idx, i*2+1) = shared_winding_number; 
                    }
                }
                Q.push(prev_patch_idx);
            }
        }
    }

    using namespace propagate_winding_numbers_helper;
    return winding_number_assignment_is_consistent(orders, orientations, patch_W);
}

template<
typename DerivedV,
typename DerivedF,
typename DerivedC,
typename DerivedW>
IGL_INLINE bool igl::cgal::propagate_winding_numbers_single_component(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedC>& labels,
        Eigen::PlainObjectBase<DerivedW>& W) {
    typedef typename DerivedF::Scalar Index;
    const size_t num_faces = F.rows();

    // Extract unique edges.
    std::vector<std::vector<size_t> > uE2E;
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> E, uE;
    Eigen::Matrix<Index, Eigen::Dynamic, 1> EMAP;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);

    // Extract manifold patches and intersection curves.
    Eigen::VectorXi P;
    std::vector<std::vector<size_t> > intersection_curves;
    size_t num_patches =
        igl::extract_manifold_patches(F, EMAP, uE2E, P);
    igl::extract_non_manifold_edge_curves(F, EMAP, uE2E,
            intersection_curves);
    assert(P.size() == num_faces);
    assert(P.maxCoeff() + 1 == num_patches);

    Eigen::VectorXi patch_labels(num_patches);
    const int INVALID = std::numeric_limits<int>::max();
    patch_labels.setConstant(INVALID);
    for (size_t i=0; i<num_faces; i++) {
        if (patch_labels[P[i]] == INVALID) {
            patch_labels[P[i]] = labels[i];
        } else {
            assert(patch_labels[P[i]] == labels[i]);
        }
    }
    assert((patch_labels.array() != INVALID).all());
    const size_t num_labels = patch_labels.maxCoeff()+1;

    Eigen::MatrixXi winding_numbers;
    bool is_consistent =
        igl::cgal::propagate_winding_numbers_single_component_patch_wise(
            V, F, uE, uE2E, patch_labels, P, intersection_curves, winding_numbers);
    assert(winding_numbers.rows() == num_patches);
    assert(winding_numbers.cols() == 2 * num_labels);

    W.resize(num_faces, 2*num_labels);
    for (size_t i=0; i<num_faces; i++) {
        W.row(i) = winding_numbers.row(P[i]);
    }

    return is_consistent;
}

template<
typename DerivedV,
typename DerivedF,
typename DerivedW>
IGL_INLINE bool igl::cgal::propagate_winding_numbers_single_component(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedW>& W) {
    const size_t num_faces = F.rows();
    Eigen::VectorXi labels(num_faces);
    labels.setZero();
    return igl::cgal::propagate_winding_numbers_single_component(V, F, labels, W);
}

template<
typename DerivedV,
typename DerivedF,
typename DerivedL,
typename DerivedW>
IGL_INLINE void igl::cgal::propagate_winding_numbers(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedL>& labels,
        Eigen::PlainObjectBase<DerivedW>& W) {
    typedef typename DerivedF::Scalar Index;
    const size_t num_faces = F.rows();

    // Extract unique edges.
    std::vector<std::vector<size_t> > uE2E;
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> E, uE;
    Eigen::Matrix<Index, Eigen::Dynamic, 1> EMAP;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);

    // Check to make sure there is no boundaries and no non-manifold edges with
    // odd number of adjacent faces.
    for (const auto& adj_faces : uE2E) {
        if (adj_faces.size() % 2 == 1) {
            std::stringstream err_msg;
            err_msg << "Input mesh contains odd number of faces "
                << "sharing a single edge" << std::endl;
            err_msg << "Indicating the input mesh does not represent a valid volume"
                << ", and winding number cannot be propagated." << std::endl;
            throw std::runtime_error(err_msg.str());
        }
    }

    // Gather connected components.
    std::vector<std::vector<std::vector<Index > > > TT,_1;
    triangle_triangle_adjacency(E,EMAP,uE2E,false,TT,_1);
    Eigen::VectorXi counts;
    Eigen::VectorXi C;
    igl::facet_components(TT,C,counts);

    const size_t num_components = counts.size();
    std::vector<std::vector<size_t> > components(num_components);
    for (size_t i=0; i<num_faces; i++) {
        components[C[i]].push_back(i);
    }
    std::vector<Eigen::MatrixXi> comp_faces(num_components);
    std::vector<Eigen::VectorXi> comp_labels(num_components);
    for (size_t i=0; i<num_components; i++) {
        const auto& comp = components[i];
        auto& faces = comp_faces[i];
        auto& c_labels = comp_labels[i];
        const size_t comp_size = comp.size();
        faces.resize(comp_size, 3);
        c_labels.resize(comp_size);
        for (size_t j=0; j<comp_size; j++) {
            faces.row(j) = F.row(comp[j]);
            c_labels[j] = labels[comp[j]];
        }
    }

    // Compute winding number for each component.
    const size_t num_labels = labels.maxCoeff()+1;
    W.resize(num_faces, 2*num_labels);
    W.setZero();
    std::vector<Eigen::MatrixXi> Ws(num_components);
    for (size_t i=0; i<num_components; i++) {
        bool is_consistent =
            propagate_winding_numbers_single_component(V, comp_faces[i], comp_labels[i], Ws[i]);
        const size_t num_labels_in_i = comp_labels[i].maxCoeff()+1;
        const size_t num_faces_in_comp = comp_faces[i].rows();
        assert(Ws[i].cols() == num_labels_in_i*2);
        assert(Ws[i].rows() == num_faces_in_comp);
        const auto& comp = components[i];
        for (size_t j=0; j<num_faces_in_comp; j++) {
            const size_t fid = comp[j];
            W.block(fid, 0, 1, num_labels_in_i*2) = Ws[i].row(j);
        }

        if (!is_consistent) {
            std::stringstream err_msg;
            err_msg << "Component " << i
                << " has inconsistant winding number assignment." << std::endl;
            throw std::runtime_error(err_msg.str());
        }
    }

    auto sample_component = [&](size_t cid) {
        const auto& f = comp_faces[cid].row(0).eval();
        return ((V.row(f[0]) + V.row(f[1]) + V.row(f[2])) / 3.0).eval();
    };

    std::vector<Eigen::MatrixXi> nested_Ws = Ws;
    Eigen::MatrixXi ambient_correction(num_components, 2*num_labels);
    ambient_correction.setZero();
    for (size_t i=0; i<num_components; i++) {
        DerivedV samples(num_components-1, 3);
        Eigen::VectorXi is_inside;
        auto index_without_i = [&](size_t index) {
            return index < i ? index:index-1;
        };
        for (size_t j=0; j<num_components; j++) {
            if (i == j) continue;
            samples.row(index_without_i(j)) = sample_component(j);
        }
        Eigen::VectorXi fids;
        Eigen::Matrix<bool, Eigen::Dynamic, 1> orientation;
        igl::cgal::closest_facet(V, comp_faces[i], samples,
                fids, orientation);

        const auto& comp = components[i];
        for (size_t j=0; j<num_components; j++) {
            if (i == j) continue;
            const size_t index = index_without_i(j);
            const size_t fid = fids(index, 0);
            const bool ori = orientation(index, 0);
            for (size_t k=0; k<num_labels; k++) {
                const int correction = W(comp[fid], k*2+(ori?0:1));
                ambient_correction(j, k*2  ) += correction;
                ambient_correction(j, k*2+1) += correction;
            }
        }
    }

    for (size_t i=0; i<num_components; i++) {
        const auto& comp = components[i];
        const auto& correction = ambient_correction.row(i).eval();
        for (const auto& fid : comp) {
            W.row(fid) += correction;
        }
    }
}

