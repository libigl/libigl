#include "ray_box_intersect.h"
template <
  typename Derivedsource,
  typename Deriveddir,
  typename Scalar>
IGL_INLINE bool igl::ray_box_intersect(
  const Eigen::PlainObjectBase<Derivedsource> & origin,
  const Eigen::PlainObjectBase<Deriveddir> & dir,
  const Eigen::AlignedBox<Scalar,3> & box,
  const Scalar & t0,
  const Scalar & t1)
{
  using namespace Eigen;
  // This should be precomputed and provided as input
  typedef Matrix<Scalar,1,3>  RowVector3S;
  const RowVector3S inv_dir( 1./dir(0),1./dir(1),1./dir(2));
  const std::vector<bool> sign = { inv_dir(0)<0, inv_dir(1)<0, inv_dir(2)<0};
  // http://people.csail.mit.edu/amy/papers/box-jgt.pdf
  // "An Efficient and Robust Ray–Box Intersection Algorithm"
  Scalar tmin, tmax, tymin, tymax, tzmin, tzmax;
  std::vector<RowVector3S> bounds = {box.min(),box.max()};
  tmin = ( bounds[sign[0]](0)   - origin(0)) * inv_dir(0);
  tmax = ( bounds[1-sign[0]](0) - origin(0)) * inv_dir(0);
  tymin = (bounds[sign[1]](1)   - origin(1)) * inv_dir(1);
  tymax = (bounds[1-sign[1]](1) - origin(1)) * inv_dir(1);
  if ( (tmin > tymax) || (tymin > tmax) )
  {
    return false;
  }
  if (tymin > tmin)
  {
    tmin = tymin;
  }
  if (tymax < tmax)
  {
    tmax = tymax;
  }
  tzmin = (bounds[sign[2]](2) - origin(2))   * inv_dir(2);
  tzmax = (bounds[1-sign[2]](2) - origin(2)) * inv_dir(2);
  if ( (tmin > tzmax) || (tzmin > tmax) )
  {
    return false;
  }
  if (tzmin > tmin)
  {
    tmin = tzmin;
  }
  if (tzmax < tmax)
  {
    tmax = tzmax;
  }
  return ( (tmin < t1) && (tmax > t0) );
}
