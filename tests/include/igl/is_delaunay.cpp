#include <test_common.h>
#include <igl/is_delaunay.h>
#include <igl/matlab_format.h>

TEST(is_delaunay, two_triangles)
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
  igl::is_delaunay(V,FD,DD);
  for(int f=0;f<DD.rows();f++)
  {
    for(int c=0;c<DD.cols();c++)
    {
      ASSERT_TRUE(DD(f,c));
    }
  }
  const Eigen::MatrixXi FN = 
    (Eigen::MatrixXi(2,3)<<
     0,1,2,
     2,1,3).finished();
  igl::is_delaunay(V,FN,DN);
  ASSERT_FALSE(DN(0,0));
  ASSERT_FALSE(DN(1,2));
}
