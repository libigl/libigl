#ifndef IGL_CAMERA_H
#define IGL_CAMERA_H
#include "igl_inline.h"

// Simple Camera class
namespace igl
{
  class Camera
  {
    // Public Fields
    public:
      // Rotation as quaternion (0,0,0,1) is identity
      double rotation[4];
      // Translation of origin
      double pan[3];
      // Zoom scale
      double zoom;
      // Perspective angle
      double angle;
    // Public functions
    public:
      Camera();
      // Copy constructor
      // 
      // Inputs:
      //   that  other Camera to be copied
      Camera(const Camera & that);
      ~Camera(){}
  };
};

#ifdef IGL_HEADER_ONLY
#  include "Camera.cpp"
#endif

#endif
