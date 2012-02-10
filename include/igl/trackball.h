#ifndef IGL_TRACKBALL_H
#define IGL_TRACKBALL_H
#include "igl_inline.h"

namespace igl
{
  // Applies a trackball drag to a given rotation
  // Inputs:
  //   w  width of the trackball context
  //   h  height of the trackball context
  //   speed_factor  controls how fast the trackball feels, 1 is normal
  //   down_quat  rotation at mouse down, i.e. the rotation we're applying the
  //     trackball motion to (as quaternion)
  //   down_mouse_x  x position of mouse down
  //   down_mouse_y  y position of mouse down
  //   mouse_x  current x position of mouse
  //   mouse_y  current y position of mouse
  // Outputs:
  //   quat  the resulting rotation (as quaternion)
  template <typename Q_type>
  IGL_INLINE void trackball(
    const int w,
    const int h,
    const Q_type speed_factor,
    const Q_type * down_quat,
    const int down_mouse_x,
    const int down_mouse_y,
    const int mouse_x,
    const int mouse_y,
    Q_type * quat);
}

#ifdef IGL_HEADER_ONLY
#  include "trackball.cpp"
#endif

#endif
