#include "median.h"
#include "matrix_to_list.h"

#include <vector>
#include <algorithm>

IGL_INLINE bool igl::median(const Eigen::VectorXd & V, double & m)
{
  using namespace std;
  using namespace igl;
  if(V.size() == 0)
  {
    return false;
  }
  vector<double> vV;
  matrix_to_list(V,vV);
  // http://stackoverflow.com/a/1719155/148668
  size_t n = vV.size()/2;
  nth_element(vV.begin(),vV.begin()+n,vV.end());
  m = vV[n];
  return true;
}
