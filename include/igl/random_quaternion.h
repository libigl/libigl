#ifndef IGL_RANDOM_QUATERNION_H
#define IGL_RANDOM_QUATERNION_H
#include "igl_inline.h"
#include <Eigen/Geometry>
namespace igl
{
  // Return a random quaternion via uniform sampling of the 4-sphere
  template <typename Scalar>
  IGL_INLINE Eigen::Quaternion<Scalar> random_quaternion();
}
#ifndef IGL_STATIC_LIBRARY
#include "random_quaternion.cpp"
#endif
#endif
