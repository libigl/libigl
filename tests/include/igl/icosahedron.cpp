#include <test_common.h>
#include <igl/icosahedron.h>
#include <igl/doublearea.h>

TEST_CASE("icosahedron: simple", "[igl]")
{
  Eigen::MatrixXi F(2,3);
  F << 0,1,2,
       1,2,3;
  Eigen::MatrixXd V;
  igl::icosahedron(V,F);
  REQUIRE(V.rows() == 12);
  REQUIRE(V.cols() == 3);
  REQUIRE(F.rows() == 20);
  REQUIRE(F.cols() == 3);
  REQUIRE(((V.rowwise().squaredNorm().array() - 1.0).abs()<1e-15).all());
  Eigen::VectorXd A;
  igl::doublearea(V,F,A);
  const auto std_dev = [](const Eigen::ArrayXd & vec){ return std::sqrt((vec - vec.mean()).square().sum()/(vec.size()-1));};
  REQUIRE(std_dev(A.array()) < 1e-15);
}

