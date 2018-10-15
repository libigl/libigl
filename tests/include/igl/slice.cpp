#include <test_common.h>
#include <igl/slice.h>
#include <igl/LinSpaced.h>

TEST(slice, dense_identity)
{
  // https://en.wikipedia.org/wiki/Monkey_testing
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
  Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
  Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
  {
    Eigen::MatrixXd B;
    igl::slice(A,I,J,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::MatrixXd B;
    igl::slice(A,I,1,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::MatrixXd B;
    igl::slice(A,J,2,B);
    test_common::assert_eq(A,B);
  }
}

TEST(slice, sparse_identity)
{
  Eigen::SparseMatrix<double> A = Eigen::MatrixXd::Random(10,9).sparseView();
  Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
  Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
  {
    Eigen::SparseMatrix<double> B;
    igl::slice(A,I,J,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::SparseMatrix<double> B;
    igl::slice(A,I,1,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::SparseMatrix<double> B;
    igl::slice(A,J,2,B);
    test_common::assert_eq(A,B);
  }
}

TEST(slice,density_reverse)
{
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
    Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),A.rows()-1,0);
    Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
    Eigen::MatrixXd B;
    igl::slice(A,I,J,B);
    // reverse rows (i.e., reverse each column vector)
    Eigen::MatrixXd C = A.colwise().reverse().eval();
    test_common::assert_eq(B,C);
  }
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
    Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
    Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),A.cols()-1,0);
    Eigen::MatrixXd B;
    igl::slice(A,I,J,B);
    // reverse cols (i.e., reverse each row vector)
    Eigen::MatrixXd C = A.rowwise().reverse().eval();
    test_common::assert_eq(B,C);
  }
}

TEST(slice,random)
{
  // Test whether unsorted indices are handled correctly by Randomly grow and
  // shrink a matrix by slicing out rows and columns: note that growing will
  // test whether repeated indices are correctly handled
  std::vector<std::pair<int,int> > sizes = {{30,27},{3,4}};
  for(const auto & size : sizes)
  {
    Eigen::MatrixXd A(10,9);
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        A(i,j) = A.rows()*j + i;
      }
    }
    Eigen::VectorXi I = 
      ((Eigen::VectorXd::Random(size.first,1).array()*0.5+0.5)*A.rows()
       ).cast<int>();
    Eigen::VectorXi J = 
      ((Eigen::VectorXd::Random(size.second,1).array()*0.5+0.5)*A.cols()
       ).cast<int>();
    Eigen::MatrixXd B;
    igl::slice(A,I,J,B);
    Eigen::MatrixXd C(I.size(),J.size());
    for(int i = 0;i<I.size();i++)
    {
      for(int j = 0;j<J.size();j++)
      {
        C(i,j) = A.rows()*J(j) + I(i);
      }
    }
    test_common::assert_eq(B,C);
  }
}

