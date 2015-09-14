#include <test_common.h>

#include <igl/boolean/mesh_boolean.h>
#include <igl/boolean/MeshBooleanType.h>
#include <igl/exterior_edges.h>

TEST(MeshBoolean, TwoCubes) {
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    test_common::load_mesh("two-boxes-bad-self-union.ply", V1, F1);

    Eigen::MatrixXd V2(0, 3);
    Eigen::MatrixXi F2(0, 3);

    Eigen::MatrixXd Vo;
    Eigen::MatrixXi Fo;

    igl::boolean::mesh_boolean(V1, F1, V2, F2,
            igl::boolean::MESH_BOOLEAN_TYPE_UNION,
            Vo, Fo);

    Eigen::MatrixXi Eb;
    igl::exterior_edges(Fo, Eb);

    ASSERT_EQ(0, Eb.rows());
}
