#include "ray_mesh_intersect.h"

extern "C"
{
#include "raytri.c"
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV, 
  typename DerivedF> 
IGL_INLINE bool igl::ray_mesh_intersect(
  const Eigen::PlainObjectBase<Derivedsource> & s,
  const Eigen::PlainObjectBase<Deriveddir> & dir,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  std::vector<igl::Hit> & hits)
{
  using namespace Eigen;
  using namespace std;
  // Should be but can't be const 
  Vector3d s_d = s.template cast<double>();
  Vector3d dir_d = dir.template cast<double>();
  hits.clear();
  // loop over all triangles
  for(int f = 0;f<F.rows();f++)
  {
    // Should be but can't be const 
    RowVector3d v0 = V.row(F(f,0)).template cast<double>();
    RowVector3d v1 = V.row(F(f,1)).template cast<double>();
    RowVector3d v2 = V.row(F(f,2)).template cast<double>();
    // shoot ray, record hit
    double t,u,v;
    if(intersect_triangle1(
      s_d.data(), dir_d.data(), v0.data(), v1.data(), v2.data(), &t, &u, &v) &&
      t>0)
    {
      hits.push_back({(int)f,(int)-1,(float)u,(float)v,(float)t});
    }
  }
  // Sort hits based on distance
  std::sort(
    hits.begin(),
    hits.end(),
    [](const Hit & a, const Hit & b)->bool{ return a.t < b.t;});
  return hits.size() > 0;
}

template <
  typename Derivedsource,
  typename Deriveddir,
  typename DerivedV, 
  typename DerivedF> 
IGL_INLINE bool igl::ray_mesh_intersect(
  const Eigen::PlainObjectBase<Derivedsource> & source,
  const Eigen::PlainObjectBase<Deriveddir> & dir,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  igl::Hit & hit)
{
  std::vector<igl::Hit> hits;
  ray_mesh_intersect(source,dir,V,F,hits);
  if(hits.size() > 0)
  {
    hit = hits.front();
    return true;
  }else
  {
    return false;
  }
}
