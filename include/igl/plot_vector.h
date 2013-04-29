#ifndef IGL_PLOT_VECTOR_H
#define IGL_PLOT_VECTOR_H
#include "igl_inline.h"

#include <vector>
namespace igl 
{
  // Not clear what these are supposed to be doing. Currently they print
  // vectors to standard error...
  template <typename T>
  IGL_INLINE void plot_vector( std::vector<T>& v);
  template <typename T>
  IGL_INLINE void plot_vector( std::vector< std::vector<T> >& v);
  template <typename T>
  IGL_INLINE void plot_vector(std::vector< std::vector< std::vector<T> > >& v);
}

#ifdef IGL_HEADER_ONLY
#  include "plot_vector.cpp"
#endif

#endif
