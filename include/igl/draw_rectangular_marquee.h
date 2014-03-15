#ifndef IGL_DRAW_RECTANGULAR_MARQUEE_H
#define IGL_DRAW_RECTANGULAR_MARQUEE_H
namespace igl
{
  // Draw a rectangular marquee (selection box) in screen space. This sets up
  // screen space using current viewport.
  //
  // Inputs:
  //   from_x  x coordinate of from point
  //   from_y  y coordinate of from point
  //   to_x  x coordinate of to point
  //   to_y  y coordinate of to point
  void draw_rectangular_marquee(
    const int from_x,
    const int from_y,
    const int to_x,
    const int to_y);
}

#ifdef IGL_HEADER_ONLY
#  include "draw_rectangular_marquee.cpp"
#endif

#endif
