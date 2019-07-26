#include <test_common.h>
#include <igl/path_to_edges.h>

#include <vector>

TEST_CASE("igl_path_to_edges: basic_test", "[igl]")
{
  const Eigen::VectorXi I = (Eigen::VectorXi(6)<<0,1,2,3,4,5).finished();
  const Eigen::MatrixXi Eexpected = (Eigen::MatrixXi(5,2)<<0,1, 1,2, 2,3, 3,4, 4,5).finished();

  Eigen::MatrixXi Eactual;
  igl::path_to_edges(I, Eactual);

  test_common::assert_eq(Eactual, Eexpected);
}

#include <iostream> 

TEST_CASE("igl_path_to_edges: loop_test", "[igl]")
{
  const Eigen::VectorXi I = (Eigen::VectorXi(6)<<0,1,2,3,4,5).finished();
  const Eigen::MatrixXi Eexpected = (Eigen::MatrixXi(6,2)<<0,1, 1,2, 2,3, 3,4, 4,5, 5,0).finished();

  Eigen::MatrixXi Eactual;
  const bool make_loop = true;
  igl::path_to_edges(I, Eactual, make_loop);
  std::cout << Eactual << std::endl;
  test_common::assert_eq(Eactual, Eexpected);
}

TEST_CASE("igl_path_to_edges: vector_basic_test", "[igl]")
{
  const std::vector<int> I{0,1,2,3,4,5};
  const Eigen::MatrixXi Eexpected = (Eigen::MatrixXi(5,2)<<0,1, 1,2, 2,3, 3,4, 4,5).finished();

  Eigen::MatrixXi Eactual;
  igl::path_to_edges(I, Eactual);

  test_common::assert_eq(Eactual, Eexpected);
}


TEST_CASE("igl_path_to_edges: vector_loop_test", "[igl]")
{
  const std::vector<int> I{0,1,2,3,4,5};
  const Eigen::MatrixXi Eexpected = (Eigen::MatrixXi(6,2)<<0,1, 1,2, 2,3, 3,4, 4,5, 5,0).finished();

  Eigen::MatrixXi Eactual;
  const bool make_loop = true;
  igl::path_to_edges(I, Eactual, make_loop);

  test_common::assert_eq(Eactual, Eexpected);
}