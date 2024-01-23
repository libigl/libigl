#include <test_common.h>
#include <igl/split_nonmanifold.h>
#include <igl/facet_components.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/triangulated_grid.h>
#include <igl/remove_duplicate_vertices.h>

#include <igl/matlab_format.h>
#include <igl/writePLY.h>
#include <iostream>

void check_same_but_manifold(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & SF,
    const Eigen::VectorXi & I)
{
  // Check that new mesh has no non-manifold edges
  REQUIRE(igl::is_edge_manifold(SF));
  // Check that new mesh has no non-manifold vertices
  REQUIRE(igl::is_vertex_manifold(SF));
  // Check the new mesh references original vertex information
  Eigen::MatrixXi ISF(SF.rows(),SF.cols());
  for(int i = 0;i<ISF.rows();i++)
  {
    for(int j = 0;j<ISF.cols();j++)
    {
      ISF(i,j) = I(SF(i,j));
    }
  }
  test_common::assert_eq(F,ISF);
}

TEST_CASE("split_nonmanifold: edge-fan", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(5,3);
  F<<0,1,3,
     0,3,2,
     0,4,3,
     0,3,5,
     0,3,6;
  Eigen::MatrixXd V(7,3);
  V << 0,0,0,
       1,0,0,
      -1,0,0,
       0,1,0,
       0,0,1,
       0,0,-1,
       1,0,1;
  igl::writePLY("edge-fan.ply",V,F);
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: vertex-boundary", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(5,2);
  V << 0,0,
       1,0,
       0,1,
       2,0,
       1,1;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     1,3,4;
  {
    Eigen::MatrixXd V3(V.rows(),3);
    V3<<V,Eigen::MatrixXd::Zero(V.rows(),1);
    igl::writePLY("vertex-boundary.ply",V3,F);
  }
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  REQUIRE(SVI.rows() == 6);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: edge-disk-flap", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(5,3);
  F<<
    0,1,2,
    0,2,3,
    0,3,4,
    0,4,1,
    0,5,1;
  Eigen::MatrixXd V(6,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    -1,0,0,
    0,-1,0,
    0,0,1;
  igl::writePLY("edge-disk-flap.ply",V,F);
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: edge-disk-tent", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(5,3);
  F<<
    0,1,2,
    0,2,3,
    0,3,1,
    0,4,1,
    1,4,3;
  Eigen::MatrixXd V(5,3);
  V<<
    0,0,0,
    1,0,0,
    -1,1,0,
    0,-1,0,
    0,0,1;
  igl::writePLY("edge-disk-tent.ply",V,F);
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: vertex-kiss", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(6,3);
  F<<
    0,1,3,
    1,2,3,
    2,0,3,
    4,5,3,
    5,6,3,
    6,4,3;
  Eigen::MatrixXd V(7,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1,
    0,0,2,
    1,0,2,
    0,1,2;
  igl::writePLY("vertex-kiss.ply",V,F);
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: non-orientable", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(6,3);
  F<<
    0,2,1,
    2,3,1,
    2,4,3,
    4,5,3,
    4,1,5,
    1,0,5;
  Eigen::MatrixXd V(6,3);
  V<<
     6, 0, 0,
     4, 0, 0,
    -3, 5, 0,
    -2, 4, 0,
    -2,-4, 1,
    -3,-5,-1;
  igl::writePLY("non-orientable.ply",V,F);
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(F,SF,SVI);
  check_same_but_manifold(F,SF,SVI);
}

TEST_CASE("split_nonmanifold: flap", "[igl]")
{
  Eigen::MatrixXd V(12,3);
  V<< 
    1.5,0,0,
    0.75,0,0.433,
    0.75,0,-0.433,
    0,1.5,0,
    0,0.75,0.433,
    0,0.75,-0.433,
    -1.5,0,0,
    -0.75,0,-0.433,
    -0,-1.5,0,
    -0,-0.75,0.433,
    1.5,0,1,
    0,1.5,1;
  Eigen::MatrixXi F(12,3);
  const auto check = [&F,&V]()
  {
    Eigen::MatrixXi SF;
    Eigen::VectorXi SVI;
    igl::split_nonmanifold(F,SF,SVI);
    igl::writePLY("flap.ply",V,F);
    check_same_but_manifold(F,SF,SVI);
    {
      Eigen::VectorXi C;
      const int nc = igl::facet_components(SF,C);
      REQUIRE(nc == 2);
    }
  };
  F<< 
    0,3,1,
    3,4,1,
    2,5,0,
    5,3,0,
    3,6,4,
    5,7,3,
    7,6,3,
    8,0,9,
    0,1,9,
    2,0,8,
    0,3,11,
    0,11,10;
  check();
  F<< 
    0,3,11,
    0,1,9,
    2,5,0,
    7,6,3,
    0,3,1,
    3,6,4,
    5,3,0,
    5,7,3,
    8,0,9,
    3,4,1,
    2,0,8,
    0,11,10;
  check();
}

TEST_CASE("split_nonmanifold: crosses", "[igl]")
{
  for(int p = 1;p<7;p++)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // n = 2^p
    const int n = 1<<p;
    igl::triangulated_grid(n,3,V,F);
    V.array() -= 0.5;
    V.conservativeResize(V.rows(),3);
    V.col(2).setZero();
    Eigen::MatrixXd U(V.rows(),3);
    U << V.col(0),V.col(2),V.col(1);
    Eigen::MatrixXd VV(V.rows()*2,3);
    VV << V,U;
    Eigen::MatrixXi FF(F.rows()*2,3);
    FF << F,F.array()+V.rows();
    {
      Eigen::VectorXi I,J;
      igl::remove_duplicate_vertices(VV,FF,1e-10,V,I,J,F);
    }
    Eigen::MatrixXi SF;
    Eigen::VectorXi SVI;
    igl::split_nonmanifold(F,SF,SVI);
    check_same_but_manifold(F,SF,SVI);
  }
}
