#include "assign_scalar.h"

IGL_INLINE void igl::cgal::assign_scalar(
  const typename CGAL::Epeck::FT & cgal,
  CGAL::Epeck::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::cgal::assign_scalar(
  const typename CGAL::Epeck::FT & cgal,
  double & d)
{
  d = CGAL::to_double(cgal);
}
