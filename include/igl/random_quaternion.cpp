#include "random_quaternion.h"

template <typename Scalar>
IGL_INLINE Eigen::Quaternion<Scalar> igl::random_quaternion()
{
  // http://mathproofs.blogspot.com/2005/05/uniformly-distributed-random-unit.html
  const auto & unit_rand = []()->Scalar
  {
    return ((Scalar)rand() / (Scalar)RAND_MAX);
  };
  const Scalar t0 = 2.*M_PI*unit_rand();
  const Scalar t1 = acos(1.-2.*unit_rand());
  const Scalar t2 = 0.5*(M_PI*unit_rand() + acos(unit_rand()));
  return Eigen::Quaternion<Scalar>(
    1.*sin(t0)*sin(t1)*sin(t2),
    1.*cos(t0)*sin(t1)*sin(t2),
    1.*cos(t1)*sin(t2),
    1.*cos(t2));
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template Eigen::Quaternion<double, 0> igl::random_quaternion<double>();
#endif
