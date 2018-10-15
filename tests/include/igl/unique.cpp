#include <test_common.h>
#include <igl/unique_rows.h>
#include <igl/matrix_to_list.h>

TEST(unique,matrix)
{
  Eigen::VectorXi A(12);
  A = (Eigen::VectorXd::Random(A.size(),1).array().abs()*9).cast<int>();
  Eigen::VectorXi C,IA,IC;
  igl::unique_rows(A,C,IA,IC);
  std::vector<bool> inA(A.maxCoeff()+1,false);
  for(int i = 0;i<A.size();i++)
  {
    inA[A(i)] = true;
    ASSERT_EQ(A(i),C(IC(i)));
  }
  std::vector<bool> inC(inA.size(),false);
  // Expect a column vector
  ASSERT_EQ(1,C.cols());
  for(int i = 0;i<C.size();i++)
  {
    // Should be the first time finding this
    ASSERT_FALSE(inC[C(i)]);
    // Mark as found
    inC[C(i)] = true;
    // Should be something also found in A
    ASSERT_TRUE(inA[C(i)]);
    ASSERT_EQ(C(i),A(IA(i)));
  }
  for(int i = 0;i<inC.size();i++)
  {
    ASSERT_EQ(inC[i],inA[i]);
  }
}

TEST(unique_rows,matrix)
{
  Eigen::MatrixXi A(50,4);
  A = (Eigen::MatrixXi::Random(A.rows(),A.cols()).array().abs()*9).cast<int>();
  Eigen::MatrixXi C;
  Eigen::VectorXi IA,IC;
  igl::unique_rows(A,C,IA,IC);
  ASSERT_EQ(A.cols(),C.cols());
  ASSERT_EQ(A.rows(),IC.size());
  ASSERT_EQ(C.rows(),IA.size());
  std::map<std::vector<int>,bool> inA;
  for(int i = 0;i<A.rows();i++)
  {
    Eigen::RowVectorXi Ai = A.row(i);
    std::vector<int> vAi;
    igl::matrix_to_list(Ai,vAi);
    inA[vAi] = true;
    for(int j = 0;j<A.cols();j++)
    {
      ASSERT_EQ(A(i,j),C(IC(i),j));
    }
  }
  std::map<std::vector<int>,bool> inC;
  for(int i = 0;i<C.rows();i++)
  {
    Eigen::RowVectorXi Ci = C.row(i);
    std::vector<int> vCi;
    igl::matrix_to_list(Ci,vCi);
    // Should be the first time finding this
    ASSERT_FALSE(inC[vCi]);
    // Mark as found
    inC[vCi] = true;
    // Should be something also found in A
    ASSERT_TRUE(inA[vCi]);
    for(int j = 0;j<A.cols();j++)
    {
      ASSERT_EQ(C(i,j),A(IA(i),j));
    }
  }
  ASSERT_EQ(inC.size(),inA.size());
  for(const auto pair : inA)
  {
    ASSERT_EQ(inC[pair.first],inA[pair.first]);
  }
  for(const auto pair : inC)
  {
    ASSERT_EQ(inC[pair.first],inA[pair.first]);
  }
}
