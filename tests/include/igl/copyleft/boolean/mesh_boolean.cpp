#include <test_common.h>

#include <vector>

#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/MeshBooleanType.h>
#include <igl/exterior_edges.h>
#include <igl/is_vertex_manifold.h>
#include <igl/unique_edge_map.h>
#include <igl/is_edge_manifold.h>

namespace mesh_boolean_test {

    template<typename DerivedF>
    void assert_no_exterior_edges(
            const Eigen::PlainObjectBase<DerivedF>& F) {
        Eigen::MatrixXi Eb;
        igl::exterior_edges(F, Eb);

        REQUIRE (Eb.rows() == 0);
    }

    template<typename DerivedV, typename DerivedF>
    void assert_is_manifold(
            const Eigen::PlainObjectBase<DerivedV>& V,
            const Eigen::PlainObjectBase<DerivedF>& F) {
        Eigen::MatrixXi B;
        REQUIRE (igl::is_vertex_manifold(F, B));
        REQUIRE (igl::is_edge_manifold(F));
    }

    template<typename DerivedV, typename DerivedF>
    void assert_genus_eq(
            const Eigen::PlainObjectBase<DerivedV>& V,
            const Eigen::PlainObjectBase<DerivedF>& F,
            const int genus) {
        const int num_vertices = V.rows();
        const int num_faces = F.rows();

        Eigen::Matrix<
            typename DerivedF::Scalar,
            Eigen::Dynamic,
            Eigen::Dynamic>
            E, uE, EMAP;
        std::vector<std::vector<size_t> > uE2E;
        igl::unique_edge_map(F, E, uE, EMAP, uE2E);

        const int num_edges = uE.rows();
        const int euler = num_vertices - num_edges + num_faces;
        REQUIRE (2 - 2 * genus == euler);
    }

TEST_CASE("MeshBoolean: TwoCubes", "[igl/copyleft/boolean]")
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    test_common::load_mesh("two-boxes-bad-self-union.ply", V1, F1);

    Eigen::MatrixXd V2(0, 3);
    Eigen::MatrixXi F2(0, 3);

    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;

    igl::copyleft::cgal::mesh_boolean(V1, F1, V2, F2,
            igl::MESH_BOOLEAN_TYPE_UNION,
            Vo, Fo);

    assert_no_exterior_edges(Fo);
    assert_is_manifold(Vo, Fo);
    assert_genus_eq(Vo, Fo, 0);
}

TEST_CASE("MeshBoolean: MinusTest", "[igl/copyleft/boolean]")
    // Many thanks to Eric Yao for submitting this test case.
    Eigen::MatrixXd V1, V2, Vo;
    Eigen::MatrixXi F1, F2, Fo;
    test_common::load_mesh("boolean_minus_test_cube.obj", V1, F1);
    test_common::load_mesh("boolean_minus_test_green.obj", V2, F2);

    igl::copyleft::cgal::mesh_boolean(V1, F1, V2, F2,
            igl::MESH_BOOLEAN_TYPE_MINUS,
            Vo, Fo);

    assert_no_exterior_edges(Fo);
    assert_is_manifold(Vo, Fo);
    assert_genus_eq(Vo, Fo, 1);
}

TEST_CASE("MeshBoolean: IntersectWithSelf", "[igl/copyleft/boolean]")
    Eigen::MatrixXd V1, Vo;
    Eigen::MatrixXi F1, Fo;
    test_common::load_mesh("cube.obj", V1, F1);

    igl::copyleft::cgal::mesh_boolean(V1, F1, V1, F1,
            igl::MESH_BOOLEAN_TYPE_INTERSECT,
            Vo, Fo);

    assert_no_exterior_edges(Fo);
    assert_is_manifold(Vo, Fo);
    assert_genus_eq(Vo, Fo, 0);
}

}
