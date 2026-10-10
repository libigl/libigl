#include <test_common.h>
#include <igl/embree/reorient_facets_raycast.h>

TEST_CASE(
  "reorient_facets_raycast: cube with every other face flipped",
  "[igl/embree]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // a closed cube, every face pointing out
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  Eigen::MatrixXi G = F;
  for (int f = 0; f < G.rows(); f += 2)
  {
    G.row(f) = F.row(f).reverse().eval();
  }
  // each face alone, by patches, and by parity
  const bool facet_wise[] = {true, false, true};
  const bool use_parity[] = {false, false, true};
  for (int k = 0; k < 3; ++k)
  {
    Eigen::VectorXi I, C;
    igl::embree::reorient_facets_raycast(
      V, G, 100*G.rows(), 10, facet_wise[k], use_parity[k], false, I, C);
    REQUIRE(I.rows() == G.rows());
    Eigen::MatrixXi turned = G;
    for (int f = 0; f < G.rows(); ++f)
    {
      if (I(f)) turned.row(f) = G.row(f).reverse().eval();
    }
    test_common::assert_eq(turned, F);
  }
}
