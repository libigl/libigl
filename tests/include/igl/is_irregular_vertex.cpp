#include <test_common.h>
#include <igl/is_irregular_vertex.h>


TEST_CASE("is_irregular_vertex: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Known non-manifold mesh
    test_common::load_mesh("truck.obj", V, F);
    std::vector<bool> vec = igl::is_irregular_vertex(V,F);
    int ret = false;
    for (auto val : vec)
    {
        // true corresponds to 1
        ret += val;
    }
    // some verticies are irregular thus the sum over all verticies should evaluate to a value > 0
    REQUIRE( ret > 0 );

}