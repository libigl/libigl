#ifndef IGL_LINE_SEGMENT_IN_RECTANGLE_H
#define IGL_LINE_SEGMENT_IN_RECTANGLE_H
#include <Eigen/Core>
namespace igl
{
  // Determine whether a line segment overlaps with a rectangle.
  //
  // Inputs:
  //   s  source point of line segment
  //   d  dest point of line segment
  //   A  first corner of rectangle
  //   B  opposite corner of rectangle
  // Returns true if line segment is at all inside rectangle
  bool line_segment_in_rectangle(
    const Eigen::Vector2d & s,
    const Eigen::Vector2d & d,
    const Eigen::Vector2d & A,
    const Eigen::Vector2d & B);
}

#ifdef IGL_HEADER_ONLY
#  include "line_segment_in_rectangle.cpp"
#endif

#endif
