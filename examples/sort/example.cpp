// g++ -o main main.cpp -I. -I/usr/local/include/eigen3
#include <Eigen/Core>
#include <iostream>
using namespace std;
#include <igl/sort.h>
using namespace igl;

#include <cstdio>

template <typename T>
void matlab_print(const string name, const T & X)
{
  cout<<name<<"=["<<endl<<X<<endl<<"];"<<endl;
}

int main(int argc, char * argv[])
{
  Eigen::MatrixXd X(2,3);
  X << 3,5,2,1,3,8;
  matlab_print("X",X);

  // sort each row independently
  int dim = 2;
  // sort ascending order
  int ascending = true;
  // Sorted output matrix
  Eigen::MatrixXd Y;
  // sorted indices for sort dimension
  Eigen::MatrixXi IX;
  sort(X,dim,ascending,Y,IX);
  matlab_print<Eigen::MatrixXd>("Y",Y);
  matlab_print<Eigen::MatrixXi>("IX",IX);

  // Verify that IX really does sort X into Y
  int num_outer = (dim == 1 ? X.cols() : X.rows() );
  // get number of rows (or columns)
  int num_inner = (dim == 1 ? X.rows() : X.cols() );
  bool verified = true;
  for(int i = 0;i<num_outer;i++)
  {
    for(int j = 0;j<num_inner;j++)
    {
      if(dim == 1)
      {
        if( Y(j,i) != X(IX(j,i),i))
        {
          printf("Y(%d,%d) = %g != %g = X(%d,%d) = X(IX(%d,%d),%d)\n",
            j,i,Y(j,i),X(IX(j,i),i),IX(j,i),i,j,i,i);
          verified = false;
        }
      }else
      {
        if( Y(i,j) != X(i,IX(i,j)))
        {
          printf("Y(%d,%d) = %g != %g = X(%d,%d) = X(IX(%d,%d),%d)\n",
            i,j,Y(i,j),X(i,IX(i,j)),IX(i,j),i,i,j,i);
          verified = false;
        }
      }
    }
  }
  if(verified)
  {
    printf("Sorting succeeded\n");
  }else
  {
    printf("Sorting failed\n");
  }
  return (verified?0:1);
}
