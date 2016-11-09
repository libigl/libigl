#include "row_to_point.h"

template <
  typename Kernel,
  typename DerivedV>
IGL_INLINE CGAL::Point_2<Kernel> igl::copyleft::cgal::row_to_point(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const typename DerivedV::Index & i)
{
  return CGAL::Point_2<Kernel>(V(i,0),V(i,1));
}
