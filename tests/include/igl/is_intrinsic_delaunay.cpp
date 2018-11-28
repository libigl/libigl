#include <test_common.h>
#include <igl/is_intrinsic_delaunay.h>
#include <igl/edge_lengths.h>

TEST_CASE("is_intrinsic_delaunay: two_triangles", "[igl]")
{
  const Eigen::MatrixXd V = 
    (Eigen::MatrixXd(4,2)<<
     0,10,
     1,0,
     1,20,
     2,10).finished();
  const Eigen::MatrixXi FD = 
    (Eigen::MatrixXi(2,3)<<
     0,1,3,
     0,3,2).finished();
  Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> DD,DN;
  Eigen::MatrixXd lD;
  igl::edge_lengths(V,FD,lD);
  igl::is_intrinsic_delaunay(lD,FD,DD);
  for(int f=0;f<DD.rows();f++)
  {
    for(int c=0;c<DD.cols();c++)
    {
      REQUIRE (DD(f,c));
    }
  }
  const Eigen::MatrixXi FN = 
    (Eigen::MatrixXi(2,3)<<
     0,1,2,
     2,1,3).finished();
  Eigen::MatrixXd lN;
  igl::edge_lengths(V,FN,lN);
  igl::is_intrinsic_delaunay(lN,FN,DN);
  REQUIRE (!DN(0,0));
  REQUIRE (!DN(1,2));
}
