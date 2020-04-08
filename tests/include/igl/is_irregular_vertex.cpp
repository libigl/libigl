#include <test_common.h>
#include <igl/is_irregular_vertex.h>


TEST_CASE("is_irregular_vertex: simple", "[igl]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Known "bad" mesh (many boundaries + irregular vertices, non-manifold)
    igl::read_triangle_mesh(test_common::data_path("truck.obj"), V, F);
    std::vector<bool> vec = igl::is_irregular_vertex(V,F);
    // some vertices are irregular thus the sum over all vertices should evaluate to true
    REQUIRE(std::any_of(vec.begin(),vec.end(), [](bool v) { return v; }));

}