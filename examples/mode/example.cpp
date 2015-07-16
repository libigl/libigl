// g++ -o main main.cpp -I. -I/usr/local/include/eigen3
#include <igl/mode.h>
#include <Eigen/Core>
#include <iostream>

using namespace std;
using namespace igl;
using namespace Eigen;


template <class T>
void matlab_print(const string name, const T & X)
{
  cout<<name<<"=["<<endl<<X<<endl<<"];"<<endl;
}

int main(int argc, char * argv[])
{
  Eigen::MatrixXd X(3,4);
  X << 
    3,5,4,5,
    1,2,4,2,
    1,1,2,5;
  matlab_print("X",X);

  // Sorted output matrix
  Eigen::Matrix<double,Dynamic,1> M1;
  mode(X,1,M1);
  matlab_print("M1",M1);

  Eigen::Matrix<double,Dynamic,1> M2;
  mode(X,2,M2);
  matlab_print("M2",M2);
}
