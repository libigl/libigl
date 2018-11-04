#include <test_common.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/order_facets_around_edges.h>
#include <igl/unique_edge_map.h>
#include <igl/readDMAT.h>
#include <igl/per_face_normals.h>

namespace {

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

template<typename T>
size_t index_of(const std::vector<T>& array, T val) {
    auto loc = std::find(array.begin(), array.end(), val);
    assert(loc != array.end());
    return loc - array.begin();
}

void assert_consistently_oriented(size_t num_faces,
        const std::vector<int>& expected_face_order,
        const std::vector<int>& e_order) {
    const size_t num_items = expected_face_order.size();
    REQUIRE (e_order.size() == num_items);

    std::vector<int> order(num_items);
    std::transform(e_order.begin(), e_order.end(), order.begin(),
            [=](int val) { return val % num_faces; });

    size_t ref_start = index_of(order, expected_face_order[0]);
    for (size_t i=0; i<num_items; i++) {
        REQUIRE (order[(ref_start + i) % num_items] == expected_face_order[i]);
    }
}

template<typename DerivedV, typename DerivedF>
void assert_order(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        size_t v0, size_t v1,
        std::vector<int> expected_order, const std::string& normal="") {
    Eigen::MatrixXi E, uE, EMAP;
    std::vector<std::vector<int> > uE2E;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);

    std::vector<std::vector<int> > uE2oE;
    std::vector<std::vector<bool> > uE2C;

    if (normal != "") {
        Eigen::MatrixXd N;
        //igl::per_face_normals_stable(V, F, N);
        //igl::per_face_normals(V, F, N);
        test_common::load_matrix(normal, N);
        igl::copyleft::cgal::order_facets_around_edges(
          V, F, N, uE, uE2E, uE2oE, uE2C);
    } else {
        igl::copyleft::cgal::order_facets_around_edges(
          V, F, uE, uE2E, uE2oE, uE2C);
    }

    const size_t num_faces = F.rows();
    const size_t num_uE = uE.rows();
    for (size_t i=0; i<num_uE; i++) {
        const auto& order = uE2oE[i];
        const auto& cons  = uE2C[i];
        const auto ref_edge = uE2E[i][0];
        const auto ref_face = ref_edge % num_faces;
        const auto ref_corner = ref_edge / num_faces;
        const Eigen::Vector2i e{
            F(ref_face, (ref_corner+1)%3),
            F(ref_face, (ref_corner+2)%3) };
        if (order.size() <= 1) continue;
        if (e[0] != v0 && e[0] != v1) continue;
        if (e[1] != v0 && e[1] != v1) continue;
        if (e[0] == v1 && e[1] == v0) {
            std::reverse(expected_order.begin(), expected_order.end());
        }
        assert_consistently_oriented(F.rows(), expected_order, order);
    }
}

} // anonymous namespace

TEST_CASE("copyleft_cgal_order_facets_around_edges: Simple", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(4, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         1.0, 1.0, 0.0;
    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2,
         2, 1, 3;

    assert_order(V, F, 1, 2, {0, 1});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: TripletFaces", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(5, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         1.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    Eigen::MatrixXi F(3, 3);
    F << 0, 1, 2,
         2, 1, 3,
         1, 2, 4;

    assert_order(V, F, 1, 2, {0, 1, 2});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: DuplicatedFaces", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(5, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         1.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    Eigen::MatrixXi F(4, 3);
    F << 0, 1, 2,
         2, 1, 3,
         1, 2, 4,
         4, 1, 2;

    assert_order(V, F, 1, 2, {0, 1, 3, 2});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: MultipleDuplicatedFaces", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(5, 3);
    V << 0.0, 0.0, 0.0,
         1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         1.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    Eigen::MatrixXi F(6, 3);
    F << 0, 1, 2,
         1, 2, 0,
         2, 1, 3,
         1, 3, 2,
         1, 2, 4,
         4, 1, 2;

    assert_order(V, F, 1, 2, {1, 0, 2, 3, 5, 4});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: Debug", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(5, 3);
    V <<
        -44.3205080756887781, 4.22994972382184579e-15, 75,
        -27.933756729740665, -48.382685902179837, 75,
        -55.8675134594812945, -2.81996648254789745e-15, 75,
        -27.933756729740665, -48.382685902179837, 70,
        -31.4903810567666049, -42.2224318643354408, 85;

    Eigen::MatrixXi F(3, 3);
    F << 1, 0, 2,
         2, 3, 1,
         4, 1, 2;

    assert_order(V, F, 1, 2, {0, 2, 1});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: Debug2", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V(5, 3);
    V <<
        -22.160254037844382, 38.3826859021798441, 75,
        -27.9337567297406331, 48.3826859021798654, 75,
        27.9337567297406544, 48.3826859021798512, 75,
        27.9337567297406544, 48.3826859021798512, 70,
        20.8205080756887924, 48.3826859021798512, 85;
    Eigen::MatrixXi F(3, 3);
    F << 1, 0, 2,
         3, 1, 2,
         2, 4, 1;

    assert_order(V, F, 1, 2, {1, 0, 2});
}

TEST_CASE("copyleft_cgal_order_facets_around_edges: NormalSensitivity", "[igl/copyleft/cgal]")
{
    // This example shows that epsilon difference in normal vectors could
    // results in very different ordering of facets.

    Eigen::MatrixXd V;
    test_common::load_matrix("duplicated_faces_V.dmat", V);
    Eigen::MatrixXi F;
    test_common::load_matrix("duplicated_faces_F.dmat", F);

    assert_order(V, F, 223, 224, {2, 0, 3, 1}, "duplicated_faces_N1.dmat");
    assert_order(V, F, 223, 224, {0, 3, 2, 1}, "duplicated_faces_N2.dmat");
}
