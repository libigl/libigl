#include <test_common.h>
#include <igl/intrinsic_delaunay_triangulation.h>
#include <igl/edge_lengths.h>
#include <igl/triangulated_grid.h>
#include <igl/is_delaunay.h>
#include <igl/is_intrinsic_delaunay.h>
#include <igl/is_edge_manifold.h>
#include <igl/unique_simplices.h>
#include <igl/matlab_format.h>

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

TEST(intrinsic_delaunay_triangulation, peaks)
{
  Eigen::MatrixXd V2;
  Eigen::MatrixXi F_in;
  igl::triangulated_grid(1000,1000,V2,F_in);
  Eigen::MatrixXd V(V2.rows(),3);
  for(int v=0;v<V.rows();v++)
  {
    const auto x = (V2(v,0)-0.5)*6.0;
    const auto y = (V2(v,1)-0.5)*6.0;
    // peaks.m
    const auto z = 3.*(1.-x)*(1.-x)*std::exp(-(x*x) - (y+1.)*(y+1.)) +
      - 10.*(x/5. - x*x*x - y*y*y*y*y)*std::exp(-x*x-y*y) +
      - 1./3.*std::exp(-(x+1.)*(x+1.) - y*y);
    V(v,0) = x;
    V(v,1) = y;
    V(v,2) = z;
  }
  Eigen::MatrixXd l_in;
  igl::edge_lengths(V,F_in,l_in);
  Eigen::MatrixXd l;
  Eigen::MatrixXi F;
  igl::intrinsic_delaunay_triangulation( l_in, F_in, l, F);
  {
    Eigen::MatrixXi Fu;
    Eigen::VectorXi IA,IC;
    igl::unique_simplices(F,Fu,IA,IC);
    ASSERT_EQ(F.rows(),Fu.rows());
  }
  // Input better have started manifold, otherwise this test doesn't make sense
  assert(igl::is_edge_manifold(F_in));
  ASSERT_TRUE(igl::is_edge_manifold(F));

  //Eigen::MatrixXd lext;
  //igl::edge_lengths(V,F,lext);
  //std::cout<<igl::matlab_format(V,"V")<<std::endl;
  //std::cout<<igl::matlab_format(F_in.array()+1,"F_in")<<std::endl;
  //std::cout<<igl::matlab_format(F.array()+1,"F")<<std::endl;
  //std::cout<<igl::matlab_format(l,"l")<<std::endl;
}

TEST(intrinsic_delaunay_triangulation,unflippable_tet)
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(4,3)<<
    10, 4,7,
     5, 9,0,
     8, 8,8,
     1,10,9).finished();
  const Eigen::MatrixXi F_in = (Eigen::MatrixXi(4,3)<<
     0,1,2,
     0,2,3,
     0,3,1,
     1,3,2).finished();
  const Eigen::Matrix<bool,Eigen::Dynamic,3> Dgt = 
    (Eigen::Matrix<bool,Eigen::Dynamic,3>(4,3)<<
     1,1,1,
     1,0,1,
     1,1,0,
     1,1,1).finished();
  Eigen::Matrix<bool,Eigen::Dynamic,3> D;
  Eigen::MatrixXd l_in;
  igl::edge_lengths(V,F_in,l_in);
  igl::is_intrinsic_delaunay(l_in,F_in,D);
  test_common::assert_eq(D,Dgt);
  Eigen::MatrixXd l;
  Eigen::MatrixXi F;
  igl::intrinsic_delaunay_triangulation( l_in, F_in, l, F);
  // Nothing should have changed: no edges are flippable.
  test_common::assert_eq(F,F_in);
  test_common::assert_eq(l,l_in);
}
