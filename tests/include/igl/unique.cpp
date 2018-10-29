#include <test_common.h>
#include <igl/unique_rows.h>
#include <igl/matrix_to_list.h>

TEST_CASE("unique: matrix", "[igl]")
{
  Eigen::VectorXi A(12);
  A = (Eigen::VectorXd::Random(A.size(),1).array().abs()*9).cast<int>();
  Eigen::VectorXi C,IA,IC;
  igl::unique_rows(A,C,IA,IC);
  std::vector<bool> inA(A.maxCoeff()+1,false);
  for(int i = 0;i<A.size();i++)
  {
    inA[A(i)] = true;
    REQUIRE (C(IC(i)) == A(i));
  }
  std::vector<bool> inC(inA.size(),false);
  // Expect a column vector
  REQUIRE (C.cols() == 1);
  for(int i = 0;i<C.size();i++)
  {
    // Should be the first time finding this
    REQUIRE (!inC[C(i)]);
    // Mark as found
    inC[C(i)] = true;
    // Should be something also found in A
    REQUIRE (inA[C(i)]);
    REQUIRE (A(IA(i)) == C(i));
  }
  for(int i = 0;i<inC.size();i++)
  {
    REQUIRE (inA[i] == inC[i]);
  }
}

TEST_CASE("unique_rows: matrix", "[igl]")
{
  Eigen::MatrixXi A(50,4);
  A = (Eigen::MatrixXi::Random(A.rows(),A.cols()).array().abs()*9).cast<int>();
  Eigen::MatrixXi C;
  Eigen::VectorXi IA,IC;
  igl::unique_rows(A,C,IA,IC);
  REQUIRE (C.cols() == A.cols());
  REQUIRE (IC.size() == A.rows());
  REQUIRE (IA.size() == C.rows());
  std::map<std::vector<int>,bool> inA;
  for(int i = 0;i<A.rows();i++)
  {
    Eigen::RowVectorXi Ai = A.row(i);
    std::vector<int> vAi;
    igl::matrix_to_list(Ai,vAi);
    inA[vAi] = true;
    for(int j = 0;j<A.cols();j++)
    {
      REQUIRE (C(IC(i),j) == A(i,j));
    }
  }
  std::map<std::vector<int>,bool> inC;
  for(int i = 0;i<C.rows();i++)
  {
    Eigen::RowVectorXi Ci = C.row(i);
    std::vector<int> vCi;
    igl::matrix_to_list(Ci,vCi);
    // Should be the first time finding this
    REQUIRE (!inC[vCi]);
    // Mark as found
    inC[vCi] = true;
    // Should be something also found in A
    REQUIRE (inA[vCi]);
    for(int j = 0;j<A.cols();j++)
    {
      REQUIRE (A(IA(i),j) == C(i,j));
    }
  }
  REQUIRE (inA.size() == inC.size());
  for(const auto pair : inA)
  {
    REQUIRE (inA[pair.first] == inC[pair.first]);
  }
  for(const auto pair : inC)
  {
    REQUIRE (inA[pair.first] == inC[pair.first]);
  }
}
