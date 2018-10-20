#include <test_common.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/edge_lengths.h>
#include <igl/triangulated_grid.h>
#include <igl/is_delaunay.h>

TEST(intrinsic_delaunay_triangulation, two_triangles)
{
  const Eigen::MatrixXd V = 
    (Eigen::MatrixXd(4,2)<<
     0,12,
     1,0,
     1,20,
     2,9).finished();
  const Eigen::MatrixXi FN = 
    (Eigen::MatrixXi(2,3)<<
     0,1,2,
     2,1,3).finished();
  Eigen::MatrixXd lN;
  igl::edge_lengths(V,FN,lN);
  Eigen::MatrixXd l;
  Eigen::MatrixXi F;
  igl::intrinsic_delaunay_triangulation( lN, FN, l, F);
  Eigen::MatrixXd lext;
  igl::edge_lengths(V,F,lext);
  test_common::assert_near(l,lext,1e-10);

}

TEST(intrinsic_delaunay_triangulation, skewed_grid)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F_in;
  igl::triangulated_grid(4,4,V,F_in);
  // Skew
  V.col(0) += 1.1*V.col(1);
  Eigen::MatrixXd l_in;
  igl::edge_lengths(V,F_in,l_in);
  Eigen::MatrixXd l;
  Eigen::MatrixXi F;
  igl::intrinsic_delaunay_triangulation( l_in, F_in, l, F);
  Eigen::MatrixXd lext;
  igl::edge_lengths(V,F,lext);
  test_common::assert_near(l,lext,1e-10);
  Eigen::Matrix<bool,Eigen::Dynamic,3> D;
  igl::is_delaunay(V,F,D);
  const Eigen::Matrix<bool,Eigen::Dynamic,3> Dtrue = 
    Eigen::Matrix<bool,Eigen::Dynamic,3>::Constant(F.rows(),F.cols(),true);
  test_common::assert_eq(D,Dtrue);
}
