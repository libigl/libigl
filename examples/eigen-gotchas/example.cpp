// Make an effort to verify bugs/gotchas for column and row major
//#define EIGEN_DEFAULT_TO_ROW_MAJOR

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_IM_MAD_AS_HELL_AND_IM_NOT_GOING_TO_TAKE_IT_ANYMORE
#include <Eigen/Sparse>
using namespace Eigen;

#include <cstdio>
#include <iostream>
using namespace std;

#include <igl/print_ijv.h>
using namespace igl;

#if EIGEN_VERSION_AT_LEAST(3,0,92)
#  warning these gotchas have not been verified for your Eigen Version
#else
// Eigen fails to notice at compile time that the inneriterator used to loop
// over the contents of a sparsematrix of type T is a different type
//
// http://www.alecjacobson.com/weblog/?p=2216
void wrong_sparsematrix_inner_iterator_type()
{
  // Fill 10 by 10 matrix with 0.5*(1:10) along the diagonal
  SparseMatrix<double> A(10,10);
  A.reserve(10);
  for(int i = 0;i<10;i++)
  {
    A.insert(i,i) = (double)i/2.0;
  }
  A.finalize();
  cout<<"AIJV=["<<endl;print_ijv(A,1);cout<<endl<<"];"<<endl<<
    "A=sparse(AIJV(:,1),AIJV(:,2),AIJV(:,3),"<<
    A.rows()<<","<<A.cols()<<");"<<endl;

  // Traverse A as *int*
  for(int k = 0;k<A.outerSize();k++)
  {
    // Each entry seems to be cast to int and if it's cast to zero then it's as
    // if it wasn't even there
    for (SparseMatrix<int>::InnerIterator it(A,k); it; ++it)
    {
      printf("A(%d,%d) = %d\n",it.row(),it.col(),it.value());
    }
  }
}

// Eigen can't handle .transpose within expression on right hand side of =
// operator
//
// Temporary solution: Never use Something.transpose() in expression. Always
// first compute Something into a matrix:
//   SparseMatrix S = Something;
// then compute tranpose
//   SparseMatrix ST = Something.transpose();
// then continue
void sparsematrix_transpose_in_rhs_aliasing()
{
  SparseMatrix<double> A(7,2);
  A.reserve(4);
  A.insert(0,0) = -0.5;
  A.insert(2,0) = -0.5;
  A.insert(4,1) = -0.5;
  A.insert(6,1) = -0.5;
  A.finalize();
  cout<<"AIJV=["<<endl;print_ijv(A,1);cout<<endl<<"];"<<endl<<
    "A=sparse(AIJV(:,1),AIJV(:,2),AIJV(:,3),"<<
    A.rows()<<","<<A.cols()<<");"<<endl;

  SparseMatrix<double> B(2,7);
  B.reserve(4);
  B.insert(0,0) = -0.5;
  B.insert(0,2) = -0.5;
  B.insert(1,4) = -0.5;
  B.insert(1,6) = -0.5;
  B.finalize();
  cout<<"BIJV=["<<endl;print_ijv(B,1);cout<<endl<<"];"<<endl<<
    "B=sparse(BIJV(:,1),BIJV(:,2),BIJV(:,3),"<<
    B.rows()<<","<<B.cols()<<");"<<endl;

  SparseMatrix<double> C;

  // Should be empty but isn't
  C = A-B.transpose();
  cout<<"C = A - B.transpose();"<<endl;
  cout<<"CIJV=["<<endl;print_ijv(C,1);cout<<endl<<"];"<<endl<<
    "C=sparse(CIJV(:,1),CIJV(:,2),CIJV(:,3),"<<
    C.rows()<<","<<C.cols()<<");"<<endl;

  // Should be empty but isn't
  C = A-B.transpose().eval();
  cout<<"C = A - B.transpose().eval();"<<endl;
  cout<<"CIJV=["<<endl;print_ijv(C,1);cout<<endl<<"];"<<endl<<
    "C=sparse(CIJV(:,1),CIJV(:,2),CIJV(:,3),"<<
    C.rows()<<","<<C.cols()<<");"<<endl;

  // This works
  SparseMatrix<double> BT = B.transpose();
  C = A-BT;
  cout<<"C = A - BT;"<<endl;
  cout<<"CIJV=["<<endl;print_ijv(C,1);cout<<endl<<"];"<<endl<<
    "C=sparse(CIJV(:,1),CIJV(:,2),CIJV(:,3),"<<
    C.rows()<<","<<C.cols()<<");"<<endl;
}

// Eigen claims the sparseLLT can be constructed using a SparseMatrix but
// at runtime crashes unless it is given exclusively the lower triangle
// 
// Temporary solution, replace:
// SparseLLT<SparseMatrix<T> > A_LLT(A);
// with:
// SparseLLT<SparseMatrix<T> > A_LLT(A.template triangularView<Eigen::Lower>());
//
void sparsellt_needs_triangular_view()
{
  // Succeeds
  SparseMatrix<double> A(2,2);
  A.reserve(4);
  A.insert(0,0) = 1;
  A.insert(0,1) = -0.5;
  A.insert(1,0) = -0.5;
  A.insert(1,1) = 1;
  A.finalize();
  cout<<"AIJV=["<<endl;print_ijv(A,1);cout<<endl<<"];"<<endl<<
    "A=sparse(AIJV(:,1),AIJV(:,2),AIJV(:,3),"<<
    A.rows()<<","<<A.cols()<<");"<<endl;
  SparseLLT<SparseMatrix<double> > A_LLT(A.triangularView<Eigen::Lower>());
  SparseMatrix<double> A_L = A_LLT.matrixL();
  cout<<"A_LIJV=["<<endl;print_ijv(A_L,1);cout<<endl<<"];"<<endl<<
    "A_L=sparse(A_LIJV(:,1),A_LIJV(:,2),A_LIJV(:,3),"<<
    A_L.rows()<<","<<A_L.cols()<<");"<<endl;

  //Crashes
  SparseMatrix<double> B(2,2);
  B.reserve(4);
  B.insert(0,0) = 1;
  B.insert(0,1) = -0.5;
  B.insert(1,0) = -0.5;
  B.insert(1,1) = 1;
  B.finalize();
  cout<<"BIJV=["<<endl;print_ijv(B,1);cout<<endl<<"];"<<endl<<
    "B=sparse(BIJV(:,1),BIJV(:,2),BIJV(:,3),"<<
    B.rows()<<","<<B.cols()<<");"<<endl;
  SparseLLT<SparseMatrix<double> > B_LLT(B);
}

void sparsematrix_nonzeros_after_expression()
{
  SparseMatrix<double> A(2,2);
  A.reserve(4);
  A.insert(0,0) = 1;
  A.insert(0,1) = -0.5;
  A.insert(1,0) = -0.5;
  A.insert(1,1) = 1;
  A.finalize();
  cout<<"AIJV=["<<endl;print_ijv(A,1);cout<<endl<<"];"<<endl<<
    "A=sparse(AIJV(:,1),AIJV(:,2),AIJV(:,3),"<<
    A.rows()<<","<<A.cols()<<");"<<endl;
  // Succeeds
  SparseMatrix<double> AmA = A-A;
  cout<<"(AmA).nonZeros(): "<<AmA.nonZeros()<<endl;
  // Succeeds
  cout<<"(A-A).eval().nonZeros(): "<<(A-A).eval().nonZeros()<<endl;
  // Crashes
  cout<<"(A-A).nonZeros(): "<<(A-A).nonZeros()<<endl;
}

// Eigen's SparseLLT's succeeded() method claims to return whether LLT
// computation was successful. Instead it seems its value is meaningless
//
// Temporary solution: Check for presence NaNs in matrixL()
// bool succeeded = (A_LLT.matrixL()*0).eval().nonZeros() == 0;
void sparsellt_succeeded_is_meaningless()
{
  // Should succeed
  SparseMatrix<double> A(2,2);
  A.reserve(4);
  A.insert(0,0) = 1;
  A.insert(0,1) = -0.5;
  A.insert(1,0) = -0.5;
  A.insert(1,1) = 1;
  A.finalize();
  cout<<"AIJV=["<<endl;print_ijv(A,1);cout<<endl<<"];"<<endl<<
    "A=sparse(AIJV(:,1),AIJV(:,2),AIJV(:,3),"<<
    A.rows()<<","<<A.cols()<<");"<<endl;
  SparseLLT<SparseMatrix<double> > A_LLT(A.triangularView<Eigen::Lower>());
  cout<<"A_LLT.succeeded(): "<<(A_LLT.succeeded()?"TRUE":"FALSE")<<endl;
  SparseMatrix<double> A_L = A_LLT.matrixL();
  // See sparsematrix_nonzeros_after_expression
  cout<<"(A_L*0).eval().nonZeros(): "<<(A_L*0).eval().nonZeros()<<endl;

  cout<<"A_LIJV=["<<endl;print_ijv(A_L,1);cout<<endl<<"];"<<endl<<
    "A_L=sparse(A_LIJV(:,1),A_LIJV(:,2),A_LIJV(:,3),"<<
    A_L.rows()<<","<<A_L.cols()<<");"<<endl;

  // Should not succeed
  SparseMatrix<double> B(2,2);
  B.reserve(4);
  B.insert(0,0) = -1;
  B.insert(0,1) = 0.5;
  B.insert(1,0) = 0.5;
  B.insert(1,1) = -1;
  B.finalize();
  cout<<"BIJV=["<<endl;print_ijv(B,1);cout<<endl<<"];"<<endl<<
    "B=sparse(BIJV(:,1),BIJV(:,2),BIJV(:,3),"<<
    B.rows()<<","<<B.cols()<<");"<<endl;
  SparseLLT<SparseMatrix<double> > B_LLT(B.triangularView<Eigen::Lower>());
  cout<<"B_LLT.succeeded(): "<<(B_LLT.succeeded()?"TRUE":"FALSE")<<endl;
  SparseMatrix<double> B_L = B_LLT.matrixL();
  // See sparsematrix_nonzeros_after_expression
  cout<<"(B_L*0).eval().nonZeros(): "<<(B_L*0).eval().nonZeros()<<endl;
  cout<<"B_LIJV=["<<endl;print_ijv(B_L,1);cout<<endl<<"];"<<endl<<
    "B_L=sparse(B_LIJV(:,1),B_LIJV(:,2),B_LIJV(:,3),"<<
    B_L.rows()<<","<<B_L.cols()<<");"<<endl;
}
#endif

int main(int argc, char * argv[])
{
  //wrong_sparsematrix_inner_iterator_type();
  //sparsematrix_transpose_in_rhs_aliasing();
  //sparsellt_needs_triangular_view();
  //sparsematrix_nonzeros_after_expression();
  //sparsellt_succeeded_is_meaningless();
  return 0;
}
