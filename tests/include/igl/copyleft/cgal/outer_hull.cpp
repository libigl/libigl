#include <test_common.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/outer_hull.h>

TEST_CASE("OuterHull: CubeWithFold", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("cube_with_fold.ply", V, F);

    Eigen::MatrixXi G,J,flip;
    // Is this just checking that it doesn't crash?
    igl::copyleft::cgal::outer_hull_legacy(V, F, G, J, flip);
}
