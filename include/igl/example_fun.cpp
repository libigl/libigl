#include "igl/example_fun.h"
#include <iostream>

template <typename Printable>
IGL_INLINE bool igl::example_fun(const Printable & input)
{
  using namespace std;
  cout<<"example_fun: "<<input<<endl;
  return true;
}

#ifndef IGL_HEADER_ONLY
template bool igl::example_fun<double>(const double& input);
template bool igl::example_fun<int>(const int& input);
#endif
