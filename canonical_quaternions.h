#ifndef IGL_CANONICAL_QUATERNIONS_H
#define IGL_CANONICAL_QUATERNIONS_H
// Define some canonical quaternions for floats and doubles
// A Quaternion, q, is defined here as an arrays of four scalars (x,y,z,w),
// such that q = x*i + y*j + z*k + w
namespace igl
{
#  define SQRT_2_OVER_2 0.707106781f

  // Float versions
  // Identity
  const float IDENTITY_QUAT_F[4] = {0,0,0,1};
  // The following match the Matlab canonical views
  // X point right, Y pointing up and Z point out
  const float XY_PLANE_QUAT_F[4] = {0,0,0,1};
  // X points right, Y points *in* and Z points up
  const float XZ_PLANE_QUAT_F[4] = {-SQRT_2_OVER_2,0,0,SQRT_2_OVER_2};
  // X points out, Y points right, and Z points up
  const float YZ_PLANE_QUAT_F[4] = {-0.5,-0.5,-0.5,0.5};
  const float CANONICAL_VIEW_QUAT_F[][4] = 
    {
      {             0,             0,             0,             1},
      {             0,             0, SQRT_2_OVER_2, SQRT_2_OVER_2},
      {             0,             0,             1,             0},
      {             0,             0, SQRT_2_OVER_2,-SQRT_2_OVER_2},
  
      {             0,            -1,             0,             0},
      {-SQRT_2_OVER_2, SQRT_2_OVER_2,             0,             0},
      {            -1,             0,             0,             0},
      {-SQRT_2_OVER_2,-SQRT_2_OVER_2,             0,             0},
  
      {          -0.5,          -0.5,          -0.5,           0.5},
      {             0,-SQRT_2_OVER_2,             0, SQRT_2_OVER_2},
      {           0.5,          -0.5,           0.5,           0.5},
      { SQRT_2_OVER_2,             0, SQRT_2_OVER_2,             0},
  
      { SQRT_2_OVER_2,             0,-SQRT_2_OVER_2,             0},
      {           0.5,           0.5,          -0.5,           0.5},
      {             0, SQRT_2_OVER_2,             0, SQRT_2_OVER_2},
      {          -0.5,           0.5,           0.5,           0.5},
  
      {             0, SQRT_2_OVER_2, SQRT_2_OVER_2,             0},
      {          -0.5,           0.5,           0.5,          -0.5},
      {-SQRT_2_OVER_2,             0,             0,-SQRT_2_OVER_2},
      {          -0.5,          -0.5,          -0.5,          -0.5},
  
      {-SQRT_2_OVER_2,             0,             0, SQRT_2_OVER_2},
      {          -0.5,          -0.5,           0.5,           0.5},
      {             0,-SQRT_2_OVER_2, SQRT_2_OVER_2,             0},
      {           0.5,          -0.5,           0.5,          -0.5}
    };
#  undef SQRT_2_OVER_2

#  define SQRT_2_OVER_2 0.707106781186548f
  // Double versions
  // Identity
  const double IDENTITY_QUAT_D[4] = {0,0,0,1};
  // The following match the Matlab canonical views
  // X point right, Y pointing up and Z point out
  const double XY_PLANE_QUAT_D[4] = {0,0,0,1};
  // X points right, Y points *in* and Z points up
  const double XZ_PLANE_QUAT_D[4] = {-SQRT_2_OVER_2,0,0,SQRT_2_OVER_2};
  // X points out, Y points right, and Z points up
  const double YZ_PLANE_QUAT_D[4] = {-0.5,-0.5,-0.5,0.5};
  const double CANONICAL_VIEW_QUAT_D[][4] = 
    {
      {             0,             0,             0,             1},
      {             0,             0, SQRT_2_OVER_2, SQRT_2_OVER_2},
      {             0,             0,             1,             0},
      {             0,             0, SQRT_2_OVER_2,-SQRT_2_OVER_2},
  
      {             0,            -1,             0,             0},
      {-SQRT_2_OVER_2, SQRT_2_OVER_2,             0,             0},
      {            -1,             0,             0,             0},
      {-SQRT_2_OVER_2,-SQRT_2_OVER_2,             0,             0},
  
      {          -0.5,          -0.5,          -0.5,           0.5},
      {             0,-SQRT_2_OVER_2,             0, SQRT_2_OVER_2},
      {           0.5,          -0.5,           0.5,           0.5},
      { SQRT_2_OVER_2,             0, SQRT_2_OVER_2,             0},
  
      { SQRT_2_OVER_2,             0,-SQRT_2_OVER_2,             0},
      {           0.5,           0.5,          -0.5,           0.5},
      {             0, SQRT_2_OVER_2,             0, SQRT_2_OVER_2},
      {          -0.5,           0.5,           0.5,           0.5},
  
      {             0, SQRT_2_OVER_2, SQRT_2_OVER_2,             0},
      {          -0.5,           0.5,           0.5,          -0.5},
      {-SQRT_2_OVER_2,             0,             0,-SQRT_2_OVER_2},
      {          -0.5,          -0.5,          -0.5,          -0.5},
  
      {-SQRT_2_OVER_2,             0,             0, SQRT_2_OVER_2},
      {          -0.5,          -0.5,           0.5,           0.5},
      {             0,-SQRT_2_OVER_2, SQRT_2_OVER_2,             0},
      {           0.5,          -0.5,           0.5,          -0.5}
    };

  const size_t NUM_CANONICAL_VIEW_QUAT = 24;

  // NOTE: I want to rather be able to return a Q_type[][] but C++ is not
  // making it easy. So instead I've written a per-element accessor
  
  // Return element [i][j] of the corresponding CANONICAL_VIEW_QUAT_* of the
  // given templated type
  // Inputs:
  //   i  index of quaternion
  //   j  index of coordinate in quaternion i
  // Returns values of CANONICAL_VIEW_QUAT_*[i][j]
  template <typename Q_type> 
  inline const Q_type CANONICAL_VIEW_QUAT(size_t i, size_t j);
  // Template specializations for float and double
  template <> 
  inline const float CANONICAL_VIEW_QUAT<float>(size_t i, size_t j);
  template <> 
  inline const double CANONICAL_VIEW_QUAT<double>(size_t i, size_t j);

#  undef SQRT_2_OVER_2
}

// Implementation

template <> inline const float igl::CANONICAL_VIEW_QUAT<float>(size_t i, size_t j)
{
  return igl::CANONICAL_VIEW_QUAT_F[i][j];
}
template <> inline const double igl::CANONICAL_VIEW_QUAT<double>(size_t i, size_t j)
{
  return igl::CANONICAL_VIEW_QUAT_D[i][j];
}

#endif
