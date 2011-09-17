#ifndef IGL_QUAT_MULT_H
#define IGL_QUAT_MULT_H

namespace igl
{
  // Computes out = q1 * q2 with quaternion multiplication
  // A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
  // such that q = x*i + y*j + z*k + w
  // Inputs:
  //   q1  left quaternion
  //   q2  right quaternion
  // Outputs:
  //   out  result of multiplication
  inline void quat_mult(
    const double *q1, 
    const double *q2,
    double *out);
};

// Implementation
// http://www.antisphere.com/Wiki/tools:anttweakbar
inline void igl::quat_mult(
  const double *q1, 
  const double *q2,
  double *out)
{
    out[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    out[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
    out[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
    out[3] = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
}

#endif
