#include <cstdio>
#include <iostream>
using namespace std;
#include <igl/transpose_blocks.h>
using namespace igl;

int main(int argc,char * argv[])
{
  // hieght of block
  int m = 2;
  // width of block
  int n = 4;
  // number of blocks
  int k = 3;
  // dimension blocks run along
  int dim = 1;

  // Input
  Eigen::MatrixXd A;

  if(dim == 1)
  {
    A.resize(m*k,n);
  }else{
    A.resize(m,n*k);
  }

  // loop over blocks
  for(int b = 0;b<k;b++)
  {
    for(int i = 0;i<m;i++)
    {
      for(int j = 0;j<n;j++)
      {
        if(dim == 1)
        {
          A(b*m+i,j) = 100*b + i*n + j;
        }else// dim == 2
        {
          A(i,b*n+j) = 100*b + i*n + j;
        }
      }
    }
  }
  cout<<"A=["<<endl<<A<<endl<<"];"<<endl;

  Eigen::MatrixXd B;
  transpose_blocks(A,k,dim,B);
  

  cout<<"B=["<<endl<<B<<endl<<"];"<<endl;
  return 0;
}
