#include "Camera.h"
#include <igl/canonical_quaternions.h>
#include <algorithm>

igl::Camera::Camera():
  zoom(1.0),
  angle(45)
{
  using namespace igl;
  using namespace std;
  // Defaults
  // canonical (X,Y) view
  copy(XY_PLANE_QUAT_D,XY_PLANE_QUAT_D+4,rotation);
  pan[0] = 0.0;
  pan[1] = 0.0;
  pan[2] = 0.0;
}

igl::Camera::Camera(const Camera & that):
  zoom(that.zoom),
  angle(that.angle)
{
  pan[0] = that.pan[0];
  pan[1] = that.pan[1];
  pan[2] = that.pan[2];
  for(int i = 0; i<4; i++)
  {
    rotation[i] = that.rotation[i];
  }
}
