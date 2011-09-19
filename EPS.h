#ifndef EPS_H
#define EPS_H
// Define a standard value for double epsilon
namespace igl
{
  const double DOUBLE_EPS    = 1.0e-14;
  const double DOUBLE_EPS_SQ = 1.0e-28;
  const float FLOAT_EPS    = 1.0e-7;
  const float FLOAT_EPS_SQ = 1.0e-14;
  // Function returning EPS for corresponding type
  template <typename S_type> inline S_type EPS();
  // Template specializations for float and double
  template <> inline float EPS<float>();
  template <> inline double EPS<double>();
}

// Implementation
template <> inline float igl::EPS()
{
  return igl::FLOAT_EPS;
}

template <> inline double igl::EPS()
{
  return igl::DOUBLE_EPS;
}
#endif
