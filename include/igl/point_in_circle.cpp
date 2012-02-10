#include "point_in_circle.h"

IGL_INLINE bool igl::point_in_circle(
  const double qx, 
  const double qy,
  const double cx, 
  const double cy,
  const double r)
{
  return (qx-cx)*(qx-cx) + (qy-cy)*(qy-cy) - r*r < 0;
}
