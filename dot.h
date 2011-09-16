#ifndef IGL_DOT_H
#define IGL_DOT_H
namespace igl
{
  // Computes out = dot(a,b)
  // Inputs:
  //   a  left 3d vector
  //   b  right 3d vector
  // Returns scalar dot product
  inline double dot(
    const double *a, 
    const double *b);
}

// Implementation
// http://www.antisphere.com/Wiki/tools:anttweakbar
inline double igl::dot(
  const double *a, 
  const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
#endif
