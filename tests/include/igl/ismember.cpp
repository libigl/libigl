#include <test_common.h>
#include <igl/ismember.h>
#include <igl/matlab_format.h>

TEST(ismember, simple)
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
        ASSERT_GE(LOCB(i,j),0);
        ASSERT_LT(LOCB(i,j),B.size());
        ASSERT_EQ(vB(LOCB(i,j)),A(i,j));
        ASSERT_EQ(LOCB(i,j),bi);
      }else
      {
        ASSERT_EQ(LOCB(i,j),-1);
        ASSERT_EQ(bi,B.size());
      }
    }
  }
}

