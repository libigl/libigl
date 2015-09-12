// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "order_facets_around_edges.h"
#include "../sort_angles.h"
#include <type_traits>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2EType,
    typename uE2oEType,
    typename uE2CType >
IGL_INLINE
typename std::enable_if<!std::is_same<typename DerivedV::Scalar,
typename CGAL::Exact_predicates_exact_constructions_kernel::FT>::value, void>::type
igl::cgal::order_facets_around_edges(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedN>& N,
        const Eigen::PlainObjectBase<DerivedE>& E,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const std::vector<std::vector<uE2EType> >& uE2E,
        std::vector<std::vector<uE2oEType> >& uE2oE,
        std::vector<std::vector<uE2CType > >& uE2C ) {

    typedef Eigen::Matrix<typename DerivedN::Scalar, 3, 1> Vector3F;
    const typename DerivedV::Scalar EPS = 1e-12;
    const size_t num_faces = F.rows();
    const size_t num_undirected_edges = uE.rows();

    auto edge_index_to_face_index = [&](size_t ei) { return ei % num_faces; };
    auto edge_index_to_corner_index = [&](size_t ei) { return ei / num_faces; };

    uE2oE.resize(num_undirected_edges);
    uE2C.resize(num_undirected_edges);

    for(size_t ui = 0;ui<num_undirected_edges;ui++)
    {
        const auto& adj_edges = uE2E[ui];
        const size_t edge_valance = adj_edges.size();
        assert(edge_valance > 0);

        const auto ref_edge = adj_edges[0];
        const auto ref_face = edge_index_to_face_index(ref_edge);
        Vector3F ref_normal = N.row(ref_face);

        const auto ref_corner_o = edge_index_to_corner_index(ref_edge);
        const auto ref_corner_s = (ref_corner_o+1)%3;
        const auto ref_corner_d = (ref_corner_o+2)%3;

        const typename DerivedF::Scalar o = F(ref_face, ref_corner_o);
        const typename DerivedF::Scalar s = F(ref_face, ref_corner_s);
        const typename DerivedF::Scalar d = F(ref_face, ref_corner_d);

        Vector3F edge = V.row(d) - V.row(s);
        auto edge_len = edge.norm();
        bool degenerated = edge_len < EPS;
        if (degenerated) {
            if (edge_valance <= 2) {
                // There is only one way to order 2 or less faces.
                edge.setZero();
            } else {
                edge.setZero();
                Eigen::Matrix<typename DerivedN::Scalar, Eigen::Dynamic, 3>
                    normals(edge_valance, 3);
                for (size_t fei=0; fei<edge_valance; fei++) {
                    const auto fe = adj_edges[fei];
                    const auto f = edge_index_to_face_index(fe);
                    normals.row(fei) = N.row(f);
                }
                for (size_t i=0; i<edge_valance; i++) {
                    size_t j = (i+1) % edge_valance;
                    Vector3F ni = normals.row(i);
                    Vector3F nj = normals.row(j);
                    edge = ni.cross(nj);
                    edge_len = edge.norm();
                    if (edge_len >= EPS) {
                        edge.normalize();
                        break;
                    }
                }

                // Ensure edge direction are consistent with reference face.
                Vector3F in_face_vec = V.row(o) - V.row(s);
                if (edge.cross(in_face_vec).dot(ref_normal) < 0) {
                    edge *= -1;
                }

                if (edge.norm() < EPS) {
                    std::cerr << "=====================================" << std::endl;
                    std::cerr << "  ui: " << ui << std::endl;
                    std::cerr << "edge: " << ref_edge << std::endl;
                    std::cerr << "face: " << ref_face << std::endl;
                    std::cerr << "  vs: " << V.row(s) << std::endl;
                    std::cerr << "  vd: " << V.row(d) << std::endl;
                    std::cerr << "adj face normals: " << std::endl;
                    std::cerr << normals << std::endl;
                    std::cerr << "Very degenerated case detected:" << std::endl;
                    std::cerr << "Near zero edge surrounded by "
                        << edge_valance << " neearly colinear faces" <<
                        std::endl;
                    std::cerr << "=====================================" << std::endl;
                }
            }
        } else {
            edge.normalize();
        }

        Eigen::MatrixXd angle_data(edge_valance, 3);
        std::vector<bool> cons(edge_valance);

        for (size_t fei=0; fei<edge_valance; fei++) {
            const auto fe = adj_edges[fei];
            const auto f = edge_index_to_face_index(fe);
            const auto c = edge_index_to_corner_index(fe);
            cons[fei] = (d == F(f, (c+1)%3));
            assert( cons[fei] ||  (d == F(f,(c+2)%3)));
            assert(!cons[fei] || (s == F(f,(c+2)%3)));
            assert(!cons[fei] || (d == F(f,(c+1)%3)));
            Vector3F n = N.row(f);
            angle_data(fei, 0) = ref_normal.cross(n).dot(edge);
            angle_data(fei, 1) = ref_normal.dot(n);
            if (cons[fei]) {
                angle_data(fei, 0) *= -1;
                angle_data(fei, 1) *= -1;
            }
            angle_data(fei, 0) *= -1; // Sort clockwise.
            angle_data(fei, 2) = (cons[fei]?1.:-1.)*(f+1);
        }

        std::cout << "angle_data:" << std::endl;
        std::cout << angle_data << std::endl;
        Eigen::VectorXi order;
        igl::sort_angles(angle_data, order);
        std::cout << "order: " << order.transpose() << std::endl;

        auto& ordered_edges = uE2oE[ui];
        auto& consistency = uE2C[ui];

        ordered_edges.resize(edge_valance);
        consistency.resize(edge_valance);
        for (size_t fei=0; fei<edge_valance; fei++) {
            ordered_edges[fei] = adj_edges[order[fei]];
            consistency[fei] = cons[order[fei]];
        }
    }
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2EType,
    typename uE2oEType,
    typename uE2CType >
IGL_INLINE 
typename std::enable_if<std::is_same<typename DerivedV::Scalar,
typename CGAL::Exact_predicates_exact_constructions_kernel::FT>::value, void>::type
igl::cgal::order_facets_around_edges(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedN>& N,
        const Eigen::PlainObjectBase<DerivedE>& E,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const std::vector<std::vector<uE2EType> >& uE2E,
        std::vector<std::vector<uE2oEType> >& uE2oE,
        std::vector<std::vector<uE2CType > >& uE2C ) {

    typedef Eigen::Matrix<typename DerivedN::Scalar, 3, 1> Vector3F;
    typedef Eigen::Matrix<typename DerivedV::Scalar, 3, 1> Vector3E;
    const typename DerivedV::Scalar EPS = 1e-12;
    const size_t num_faces = F.rows();
    const size_t num_undirected_edges = uE.rows();

    auto edge_index_to_face_index = [&](size_t ei) { return ei % num_faces; };
    auto edge_index_to_corner_index = [&](size_t ei) { return ei / num_faces; };

    uE2oE.resize(num_undirected_edges);
    uE2C.resize(num_undirected_edges);

    for(size_t ui = 0;ui<num_undirected_edges;ui++)
    {
        const auto& adj_edges = uE2E[ui];
        const size_t edge_valance = adj_edges.size();
        assert(edge_valance > 0);

        const auto ref_edge = adj_edges[0];
        const auto ref_face = edge_index_to_face_index(ref_edge);
        Vector3F ref_normal = N.row(ref_face);

        const auto ref_corner_o = edge_index_to_corner_index(ref_edge);
        const auto ref_corner_s = (ref_corner_o+1)%3;
        const auto ref_corner_d = (ref_corner_o+2)%3;

        const typename DerivedF::Scalar o = F(ref_face, ref_corner_o);
        const typename DerivedF::Scalar s = F(ref_face, ref_corner_s);
        const typename DerivedF::Scalar d = F(ref_face, ref_corner_d);

        Vector3E exact_edge = V.row(d) - V.row(s);
        exact_edge.array() /= exact_edge.squaredNorm();
        Vector3F edge(
                CGAL::to_double(exact_edge[0]),
                CGAL::to_double(exact_edge[1]),
                CGAL::to_double(exact_edge[2]));
        edge.normalize();

        Eigen::MatrixXd angle_data(edge_valance, 3);
        std::vector<bool> cons(edge_valance);

        for (size_t fei=0; fei<edge_valance; fei++) {
            const auto fe = adj_edges[fei];
            const auto f = edge_index_to_face_index(fe);
            const auto c = edge_index_to_corner_index(fe);
            cons[fei] = (d == F(f, (c+1)%3));
            assert( cons[fei] ||  (d == F(f,(c+2)%3)));
            assert(!cons[fei] || (s == F(f,(c+2)%3)));
            assert(!cons[fei] || (d == F(f,(c+1)%3)));
            Vector3F n = N.row(f);
            angle_data(fei, 0) = ref_normal.cross(n).dot(edge);
            angle_data(fei, 1) = ref_normal.dot(n);
            if (cons[fei]) {
                angle_data(fei, 0) *= -1;
                angle_data(fei, 1) *= -1;
            }
            angle_data(fei, 0) *= -1; // Sort clockwise.
            angle_data(fei, 2) = (cons[fei]?1.:-1.)*(f+1);
        }

        Eigen::VectorXi order;
        igl::sort_angles(angle_data, order);

        auto& ordered_edges = uE2oE[ui];
        auto& consistency = uE2C[ui];

        ordered_edges.resize(edge_valance);
        consistency.resize(edge_valance);
        for (size_t fei=0; fei<edge_valance; fei++) {
            ordered_edges[fei] = adj_edges[order[fei]];
            consistency[fei] = cons[order[fei]];
        }
    }
}

namespace igl {
    namespace cgal {
        namespace order_facets_around_edges_helper {
            template<typename T>
            std::vector<size_t> index_sort(const std::vector<T>& data) {
                const size_t len = data.size();
                std::vector<size_t> order(len);
                for (size_t i=0; i<len; i++) order[i] = i;

                auto comp = [&](size_t i, size_t j) {
                    return data[i] < data[j];
                };
                std::sort(order.begin(), order.end(), comp);
                return order;
            }

            // adj_faces contains signed index starting from +- 1.
            template<
                typename DerivedV,
                typename DerivedF,
                typename DerivedI >
            void order_facets_around_edge(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                size_t s, size_t d, 
                const std::vector<int>& adj_faces,
                Eigen::PlainObjectBase<DerivedI>& order) {

                // Although we only need exact predicates in the algorithm,
                // exact constructions are needed to avoid degeneracies due to
                // casting to double.
                typedef CGAL::Exact_predicates_exact_constructions_kernel K;
                typedef K::Point_3 Point_3;
                typedef K::Plane_3 Plane_3;

                auto get_face_index = [&](int adj_f)->size_t{
                    return abs(adj_f) - 1;
                };

                auto get_opposite_vertex = [&](size_t fid)->size_t {
                    if (F(fid, 0) != s && F(fid, 0) != d) return F(fid, 0);
                    if (F(fid, 1) != s && F(fid, 1) != d) return F(fid, 1);
                    if (F(fid, 2) != s && F(fid, 2) != d) return F(fid, 2);
                    assert(false);
                    return -1;
                };

                // Handle base cases
                if (adj_faces.size() == 0) {
                    order.resize(0);
                    return;
                } else if (adj_faces.size() == 1) {
                    order.resize(1);
                    order[0] = 0;
                    return;
                } else if (adj_faces.size() == 2) {
                    const size_t o1 =
                        get_opposite_vertex(get_face_index(adj_faces[0]));
                    const size_t o2 =
                        get_opposite_vertex(get_face_index(adj_faces[1]));
                    const Point_3 ps(V(s, 0), V(s, 1), V(s, 2));
                    const Point_3 pd(V(d, 0), V(d, 1), V(d, 2));
                    const Point_3 p1(V(o1, 0), V(o1, 1), V(o1, 2));
                    const Point_3 p2(V(o2, 0), V(o2, 1), V(o2, 2));
                    order.resize(2);
                    switch (CGAL::orientation(ps, pd, p1, p2)) {
                        case CGAL::POSITIVE:
                            order[0] = 1;
                            order[1] = 0;
                            break;
                        case CGAL::NEGATIVE:
                            order[0] = 0;
                            order[1] = 1;
                            break;
                        case CGAL::COPLANAR:
                            order[0] = adj_faces[0] < adj_faces[1] ? 0:1;
                            order[1] = adj_faces[0] < adj_faces[1] ? 1:0;
                            break;
                        default:
                            assert(false);
                    }
                    return;
                }

                const size_t num_adj_faces = adj_faces.size();
                const size_t o = get_opposite_vertex(
                        get_face_index(adj_faces[0]));
                const Point_3 p_s(V(s, 0), V(s, 1), V(s, 2));
                const Point_3 p_d(V(d, 0), V(d, 1), V(d, 2));
                const Point_3 p_o(V(o, 0), V(o, 1), V(o, 2));
                const Plane_3 separator(p_s, p_d, p_o);
                assert(!separator.is_degenerate());

                std::vector<Point_3> opposite_vertices;
                for (size_t i=0; i<num_adj_faces; i++) {
                    const size_t o = get_opposite_vertex(
                            get_face_index(adj_faces[i]));
                    opposite_vertices.emplace_back(
                            V(o, 0), V(o, 1), V(o, 2));
                }

                std::vector<int> positive_side;
                std::vector<int> negative_side;
                std::vector<int> tie_positive_oriented;
                std::vector<int> tie_negative_oriented;

                std::vector<size_t> positive_side_index;
                std::vector<size_t> negative_side_index;
                std::vector<size_t> tie_positive_oriented_index;
                std::vector<size_t> tie_negative_oriented_index;

                for (size_t i=0; i<num_adj_faces; i++) {
                    const int f = adj_faces[i];
                    const Point_3& p_a = opposite_vertices[i];
                    auto orientation = separator.oriented_side(p_a);
                    switch (orientation) {
                        case CGAL::ON_POSITIVE_SIDE:
                            positive_side.push_back(f);
                            positive_side_index.push_back(i);
                            break;
                        case CGAL::ON_NEGATIVE_SIDE:
                            negative_side.push_back(f);
                            negative_side_index.push_back(i);
                            break;
                        case CGAL::ON_ORIENTED_BOUNDARY:
                            {
                                const Plane_3 other(p_s, p_d, p_a);
                                const auto target_dir = separator.orthogonal_direction();
                                const auto query_dir = other.orthogonal_direction();
                                if (target_dir == query_dir) {
                                    tie_positive_oriented.push_back(f);
                                    tie_positive_oriented_index.push_back(i);
                                } else if (target_dir == -query_dir) {
                                    tie_negative_oriented.push_back(f);
                                    tie_negative_oriented_index.push_back(i);
                                } else {
                                    assert(false);
                                }
                            }
                            break;
                        default:
                            // Should not be here.
                            assert(false);
                    }
                }

                Eigen::PlainObjectBase<DerivedI> positive_order, negative_order;
                igl::cgal::order_facets_around_edges_helper::order_facets_around_edge(
                        V, F, s, d, positive_side, positive_order);
                igl::cgal::order_facets_around_edges_helper::order_facets_around_edge(
                        V, F, s, d, negative_side, negative_order);
                std::vector<size_t> tie_positive_order =
                    index_sort(tie_positive_oriented);
                std::vector<size_t> tie_negative_order =
                    index_sort(tie_negative_oriented);

                // Copy results into order vector.
                const size_t tie_positive_size = tie_positive_oriented.size();
                const size_t tie_negative_size = tie_negative_oriented.size();
                const size_t positive_size = positive_order.size();
                const size_t negative_size = negative_order.size();

                order.resize(tie_positive_size + positive_size +
                        tie_negative_size + negative_size);

                size_t count=0;
                for (size_t i=0; i<tie_positive_size; i++) {
                    order[count+i] =
                        tie_positive_oriented_index[tie_positive_order[i]];
                }
                count += tie_positive_size;

                for (size_t i=0; i<negative_size; i++) {
                    order[count+i] = negative_side_index[negative_order[i]];
                }
                count += negative_size;

                for (size_t i=0; i<tie_negative_size; i++) {
                    order[count+i] =
                        tie_negative_oriented_index[tie_negative_order[i]];
                }
                count += tie_negative_size;

                for (size_t i=0; i<positive_size; i++) {
                    order[count+i] = positive_side_index[positive_order[i]];
                }
                count += positive_size;
                assert(count == num_adj_faces);

                // Find the correct start point.
                size_t start_idx = 0;
                for (size_t i=0; i<num_adj_faces; i++) {
                    const Point_3& p_a = opposite_vertices[order[i]];
                    const Point_3& p_b =
                        opposite_vertices[order[(i+1)%num_adj_faces]];
                    if (CGAL::orientation(p_s, p_d, p_a, p_b) == CGAL::POSITIVE) {
                        start_idx = (i+1)%num_adj_faces;
                        break;
                    }
                }
                DerivedI circular_order = order;
                for (size_t i=0; i<num_adj_faces; i++) {
                    order[i] = circular_order[(start_idx + i)%num_adj_faces];
                }
            }
        }
    }
}

template<
    typename DerivedV,
    typename DerivedF,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2EType,
    typename uE2oEType,
    typename uE2CType >
IGL_INLINE void igl::cgal::order_facets_around_edges(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        const Eigen::PlainObjectBase<DerivedE>& E,
        const Eigen::PlainObjectBase<DeriveduE>& uE,
        const Eigen::PlainObjectBase<DerivedEMAP>& EMAP,
        const std::vector<std::vector<uE2EType> >& uE2E,
        std::vector<std::vector<uE2oEType> >& uE2oE,
        std::vector<std::vector<uE2CType > >& uE2C ) {

    typedef Eigen::Matrix<typename DerivedV::Scalar, 3, 1> Vector3E;
    const size_t num_faces = F.rows();
    const size_t num_undirected_edges = uE.rows();

    auto edge_index_to_face_index = [&](size_t ei) { return ei % num_faces; };
    auto edge_index_to_corner_index = [&](size_t ei) { return ei / num_faces; };

    uE2oE.resize(num_undirected_edges);
    uE2C.resize(num_undirected_edges);

    for(size_t ui = 0;ui<num_undirected_edges;ui++)
    {
        const auto& adj_edges = uE2E[ui];
        const size_t edge_valance = adj_edges.size();
        assert(edge_valance > 0);

        const auto ref_edge = adj_edges[0];
        const auto ref_face = edge_index_to_face_index(ref_edge);

        const auto ref_corner_o = edge_index_to_corner_index(ref_edge);
        const auto ref_corner_s = (ref_corner_o+1)%3;
        const auto ref_corner_d = (ref_corner_o+2)%3;

        const typename DerivedF::Scalar o = F(ref_face, ref_corner_o);
        const typename DerivedF::Scalar s = F(ref_face, ref_corner_s);
        const typename DerivedF::Scalar d = F(ref_face, ref_corner_d);

        std::vector<bool> cons(edge_valance);
        std::vector<int> adj_faces(edge_valance);
        for (size_t fei=0; fei<edge_valance; fei++) {
            const auto fe = adj_edges[fei];
            const auto f = edge_index_to_face_index(fe);
            const auto c = edge_index_to_corner_index(fe);
            cons[fei] = (d == F(f, (c+1)%3));
            adj_faces[fei] = (f+1) * (cons[fei] ? 1:-1);

            assert( cons[fei] ||  (d == F(f,(c+2)%3)));
            assert(!cons[fei] || (s == F(f,(c+2)%3)));
            assert(!cons[fei] || (d == F(f,(c+1)%3)));
        }

        Eigen::VectorXi order;
        igl::cgal::order_facets_around_edges_helper::order_facets_around_edge(
                V, F, s, d, adj_faces, order);
        assert(order.size() == edge_valance);

        auto& ordered_edges = uE2oE[ui];
        auto& consistency = uE2C[ui];

        ordered_edges.resize(edge_valance);
        consistency.resize(edge_valance);
        for (size_t fei=0; fei<edge_valance; fei++) {
            ordered_edges[fei] = adj_edges[order[fei]];
            consistency[fei] = cons[order[fei]];
        }
    }
}

