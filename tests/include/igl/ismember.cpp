#include <test_common.h>
#include <igl/ismember.h>
#include <igl/matlab_format.h>

TEST_CASE("ismember: simple", "[igl]")
{
  Eigen::MatrixXi A(3,4);
  A<<11,12,13,14,21,22,23,24,31,32,33,34;
  Eigen::MatrixXi B(2,3);
  B<<11,13,11,21,22,34;

  Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> IA;
  Eigen::MatrixXi LOCB;
  igl::ismember(A,B,IA,LOCB);
  Eigen::Map<Eigen::VectorXi> vB = 
    Eigen::Map<Eigen::VectorXi>(B.data(),B.rows()*B.cols());
  for(int i = 0;i<A.rows();i++)
  {
    for(int j = 0;j<A.cols();j++)
    {
      // try to find in b
      int bi = 0;
      for(;bi<vB.size();bi++)
      {
        if(A(i,j) == vB(bi))
        {
          break;
        }
      }
      if(IA(i,j))
      {
        REQUIRE (0 <= LOCB(i,j));
        REQUIRE (B.size() > LOCB(i,j));
        REQUIRE (A(i,j) == vB(LOCB(i,j)));
        REQUIRE (bi == LOCB(i,j));
      }else
      {
        REQUIRE (-1 == LOCB(i,j));
        REQUIRE (B.size() == bi);
      }
    }
  }
}

