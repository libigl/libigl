#include "point_segment_squared_distance.h"

template < typename Kernel>
IGL_INLINE void igl::copyleft::cgal::point_segment_squared_distance(
  const CGAL::Point_3<Kernel> & P1,
  const CGAL::Segment_3<Kernel> & S2,
  CGAL::Point_3<Kernel> & P2,
  typename Kernel::FT & d)
{
  if(S2.is_degenerate())
  {
    P2 = S2.source();
    d = (P1-P2).squared_length();
    return;
  }
  // http://stackoverflow.com/a/1501725/148668
  const auto sqr_len = S2.squared_length();
  assert(sqr_len != 0);
  const auto & V = S1.source();
  const auto & W = S2.target();
  const auto t = (P1-V).dot(W-V)/sqr_len;
  if(t<0)
  {
    P2 = V;
  }else if(t>1)
  {
    P2 = W;
  }else
  {
    P2 = V  + t*(W-V);
  }
  d = (P1-P2).squared_length();
}

