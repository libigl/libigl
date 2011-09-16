#ifndef IGL_CROSS_H
#define IGL_CROSS_H
namespace igl
{
  // Computes out = cross(a,b)
  // Inputs:
  //   a  left 3d vector
  //   b  right 3d vector
  // Outputs:
  //   out  result 3d vector
  inline void cross(
    const double *a, 
    const double *b,
    double *out);
}

// Implementation
// http://www.antisphere.com/Wiki/tools:anttweakbar
inline void igl::cross(
  const double *a, 
  const double *b,
  double *out)
{
  out[0] = a[1]*b[2]-a[2]*b[1];
  out[1] = a[2]*b[0]-a[0]*b[2];
  out[2] = a[0]*b[1]-a[1]*b[0];
}
#endif
