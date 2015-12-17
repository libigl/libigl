#ifndef IGL_RAY_MESH_INTERSECT_H
#define IGL_RAY_MESH_INTERSECT_H
#include "igl_inline.h"
#include "Hit.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{
  // Shoot a ray against a mesh (V,F) and collect all hits.
  //
  // Inputs:
  //   source  3-vector origin of ray
  //   dir  3-vector direction of ray
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of mesh face indices into V
  // Outputs:
  //    hits  **sorted** list of hits
  // Returns true if there were any hits (hits.size() > 0)
  //
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV, 
    typename DerivedF> 
  IGL_INLINE bool ray_mesh_intersect(
    const Eigen::PlainObjectBase<Derivedsource> & source,
    const Eigen::PlainObjectBase<Deriveddir> & dir,
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    std::vector<igl::Hit> & hits);
  // Outputs:
  //   hit  first hit, set only if it exists
  // Returns true if there was a hit
  template <
    typename Derivedsource,
    typename Deriveddir,
    typename DerivedV, 
    typename DerivedF> 
  IGL_INLINE bool ray_mesh_intersect(
    const Eigen::PlainObjectBase<Derivedsource> & source,
    const Eigen::PlainObjectBase<Deriveddir> & dir,
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    igl::Hit & hit);
}
#ifndef IGL_STATIC_LIBRARY
#  include "ray_mesh_intersect.cpp"
#endif
#endif
