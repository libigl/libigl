#include "example_fun.h"
#include <iostream>

template <typename Printable>
IGL_INLINE bool igl::example_fun(const Printable & input)
{
  using namespace std;
  cout<<"example_fun: "<<input<<endl;
  return true;
}

#ifndef IGL_HEADER_ONLY
// List all useful instanciations here
namespace igl_explicit_instanciations
{
  bool (*example_fun_A)(const double &) = &igl::example_fun<double>;
  bool (*example_fun_B)(const int &) = &igl::example_fun<int>;
}
#endif
