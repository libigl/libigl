#include <test_common.h>
#include <igl/triangle_fan.h>
#include <igl/exterior_edges.h>
#include <igl/sortrows.h>
#include <cstdlib>

namespace
{
  // Open "tent": a hexagon of boundary vertices 0..5 around an apex 6.
  void tent(Eigen::MatrixXi & F)
  {
    const int n = 6;
    F.resize(n,3);
    for(int i = 0;i<n;i++)
    {
      F.row(i) << n, i, (i+1)%n;
    }
  }
}

TEST_CASE("triangle_fan: deterministic", "[igl]")
{
  Eigen::MatrixXi F;
  tent(F);
  Eigen::Matrix<int,Eigen::Dynamic,2> E;
  igl::exterior_edges(F,E);
  Eigen::MatrixXi cap0;
  igl::triangle_fan(E,cap0);
  // The output must depend on E only, not on the state of rand().
  for(const unsigned int seed : {1u,2u,3u,4u,5u,6u,7u,8u})
  {
    std::srand(seed);
    Eigen::MatrixXi cap;
    igl::triangle_fan(E,cap);
    test_common::assert_eq(cap,cap0);
  }
}

TEST_CASE("triangle_fan: closes_open_boundary", "[igl]")
{
  Eigen::MatrixXi F;
  tent(F);
  Eigen::Matrix<int,Eigen::Dynamic,2> E;
  igl::exterior_edges(F,E);
  REQUIRE(E.rows() == 6);
  Eigen::MatrixXi cap;
  igl::triangle_fan(E,cap);
  // A simple loop of n edges is closed by n-2 triangles.
  REQUIRE(cap.rows() == 4);
  REQUIRE(cap.cols() == 3);
  // The cap has the same oriented boundary as the patch...
  Eigen::Matrix<int,Eigen::Dynamic,2> capE;
  igl::exterior_edges(cap,capE);
  Eigen::MatrixXi sE,scapE;
  Eigen::VectorXi I;
  igl::sortrows(Eigen::MatrixXi(E),true,sE,I);
  igl::sortrows(Eigen::MatrixXi(capE),true,scapE,I);
  test_common::assert_eq(scapE,sE);
  // ...so the patch together with the reversed cap is closed.
  Eigen::MatrixXi FC(F.rows()+cap.rows(),3);
  FC << F, cap.col(0), cap.col(2), cap.col(1);
  Eigen::Matrix<int,Eigen::Dynamic,2> EC;
  igl::exterior_edges(FC,EC);
  REQUIRE(EC.rows() == 0);
}

TEST_CASE("triangle_fan: empty", "[igl]")
{
  Eigen::Matrix<int,Eigen::Dynamic,2> E(0,2);
  Eigen::MatrixXi cap;
  igl::triangle_fan(E,cap);
  REQUIRE(cap.rows() == 0);
  REQUIRE(cap.cols() == 3);
}
