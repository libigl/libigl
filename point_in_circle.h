#ifndef IGL_POINT_IN_CIRCLE
#define IGL_POINT_IN_CIRCLE

namespace igl
{
  // Determine if 2d point is in a circle
  // Inputs:
  //   qx  x-coordinate of query point
  //   qy  y-coordinate of query point
  //   cx  x-coordinate of circle center
  //   cy  y-coordinate of circle center
  //   r  radius of circle
  // Returns true if query point is in circle, false otherwise
  inline bool point_in_circle(
    const double qx, 
    const double qy,
    const double cx, 
    const double cy,
    const double r);
}

// Implementation

inline bool igl::point_in_circle(
  const double qx, 
  const double qy,
  const double cx, 
  const double cy,
  const double r)
{
  return (qx-cx)*(qx-cx) + (qy-cy)*(qy-cy) - r*r < 0;
}
#endif
