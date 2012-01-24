#ifndef IGL_MODE_H
#define IGL_MODE_H
#include <Eigen/Dense>
namespace igl
{
  // Takes mode of coefficients in a matrix along a given dension
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   X  m by n original matrix
  //   d  dension along which to take mode, m or n
  // Outputs:
  //   M  vector containing mode along dension d, if d==1 then this will be a
  //     n-long vector if d==2 then this will be a m-long vector
  template <typename T>
  inline void mode(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
    const int d, 
    Eigen::Matrix<T,Eigen::Dynamic,1> & M);
}

// Implementation 
#include <vector>

template <typename T>
inline void igl::mode(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
  const int d, 
  Eigen::Matrix<T,Eigen::Dynamic,1> & M)
{
  assert(d==1 || d==2);
  using namespace std;
  int m = X.rows();
  int n = X.cols();
  M.resize((d==1)?n:m,1);
  for(int i = 0;i<((d==2)?m:n);i++)
  {
    vector<int> counts(((d==2)?n:m),0);
    for(int j = 0;j<((d==2)?n:m);j++)
    {
      T v = (d==2)?X(i,j):X(j,i);
      for(int k = 0;k<((d==2)?n:m);k++)
      {
        T u = (d==2)?X(i,k):X(k,i);
        if(v == u)
        {
          counts[k]++;
        }
      }
    }
    assert(counts.size() > 0);
    int max_count = -1;
    int max_count_j = -1;
    int j =0;
    for(vector<int>::iterator it = counts.begin();it<counts.end();it++)
    {
      if(max_count < *it)
      {
        max_count = *it;
        max_count_j = j;
      }
      j++;
    }
    M(i,0) = (d==2)?X(i,max_count_j):X(max_count_j,i);
  }
}

#endif
