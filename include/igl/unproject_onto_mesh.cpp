// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject_onto_mesh.h"
#include "unproject.h"
#include "unproject_ray.h"
#include "ray_mesh_intersect.h"
#include <vector>

template < typename DerivedV, typename DerivedF, typename Derivedbc>
IGL_INLINE bool igl::unproject_onto_mesh(
  const Eigen::Vector2f& pos,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  int & fid,
  Eigen::PlainObjectBase<Derivedbc> & bc)
{
  using namespace std;
  using namespace Eigen;
  const auto & shoot_ray = [&V,&F](
    const Eigen::Vector3f& s,
    const Eigen::Vector3f& dir,
    igl::Hit & hit)->bool
  {
    std::vector<igl::Hit> hits;
    if(!ray_mesh_intersect(s,dir,V,F,hits))
    {
      return false;
    }
    hit = hits[0];
    return true;
  };
  return unproject_onto_mesh(pos,model,proj,viewport,shoot_ray,fid,bc);
}

template <typename Derivedbc>
IGL_INLINE bool igl::unproject_onto_mesh(
  const Eigen::Vector2f& pos,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const std::function<
    bool(
      const Eigen::Vector3f&,
      const Eigen::Vector3f&,
      igl::Hit &)
      > & shoot_ray,
  int & fid,
  Eigen::PlainObjectBase<Derivedbc> & bc)
{
  using namespace std;
  using namespace Eigen;
  Vector3f s,dir;
  unproject_ray(pos,model,proj,viewport,s,dir);
  Hit hit;
  if(!shoot_ray(s,dir,hit))
  {
    return false;
  }
  bc.resize(3);
  bc << 1.0-hit.u-hit.v, hit.u, hit.v;
  fid = hit.id;
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif

