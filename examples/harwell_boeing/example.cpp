#define IGL_HEADER_ONLY
#include <igl/harwell_boeing.h>
#include <iostream>

template <typename T>
void print(T & v)
{
  std::cout<<v<<" ";
}

#include <cstdio>
#include <vector>
#include <algorithm>
int main(int argc,char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  SparseMatrix<double> A(5,5);
  A.insert(0,1) = 3;
  A.insert(1,0) = 22;
  A.insert(1,4) = 17;
  A.insert(2,0) = 7;
  A.insert(2,1) = 5;
  A.insert(2,3) = 1;
  A.insert(4,2) = 14;
  A.insert(4,4) = 8;

  vector<double> V;
  vector<int> R,C;

  int nr;
  harwell_boeing(A,nr,V,R,C);
  cout<<"V=[";
  for_each(V.begin(),V.end(),&print<double>);
  cout<<"];"<<endl;
  cout<<"R=[";
  for_each(R.begin(),R.end(),&print<int>);
  cout<<"];"<<endl;
  cout<<"C=[";
  for_each(C.begin(),C.end(),&print<int>);
  cout<<"];"<<endl;

  return 0;
}
