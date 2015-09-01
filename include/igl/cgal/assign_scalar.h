#ifndef IGL_CGAL_ASSIGN_SCALAR_H
#define IGL_CGAL_ASSIGN_SCALAR_H
#include "../igl_inline.h"
#include "CGAL_includes.hpp"
namespace igl
{
  namespace cgal
  {
    // Inputs:
    //   cgal  cgal scalar
    // Outputs:
    //   d  output scalar
    IGL_INLINE void assign_scalar(
      const typename CGAL::Epeck::FT & cgal,
      CGAL::Epeck::FT & d);
    IGL_INLINE void assign_scalar(
      const typename CGAL::Epeck::FT & cgal,
      double & d);
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "assign_scalar.cpp"
#endif
#endif
