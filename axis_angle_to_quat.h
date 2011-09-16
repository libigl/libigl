#include <EPS.h>
namespace igl
{
  // Convert axis angle representation of a rotation to a quaternion
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   axis  3d vector
  //   angle  scalar
  // Outputs:
  //   quaternion
  inline void axis_angle_to_quat(
    const double *axis, 
    const double angle,
    double *out);
}

// Implementation

// http://www.antisphere.com/Wiki/tools:anttweakbar
inline void igl::axis_angle_to_quat(
  const double *axis, 
  const double angle,
  double *out)
{
    double n = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    if( fabs(n)>igl::DOUBLE_EPS )
    {
        double f = 0.5*angle;
        out[3] = cos(f);
        f = sin(f)/sqrt(n);
        out[0] = axis[0]*f;
        out[1] = axis[1]*f;
        out[2] = axis[2]*f;
    }
    else
    {
        out[3] = 1.0;
        out[0] = out[1] = out[2] = 0.0;
    }
}
