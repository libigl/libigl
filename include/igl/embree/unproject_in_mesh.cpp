// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject_in_mesh.h"
#include "EmbreeIntersector.h"
#include "../unproject.h"
#include <vector>

template <typename Derivedobj>
IGL_INLINE int igl::embree::unproject_in_mesh(
  const Eigen::Vector2f& pos,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj,
  std::vector<igl::embree::Hit > & hits)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  // Source and direction on screen
  Vector3f win_s(pos(0),pos(1),0);
  Vector3f win_d(pos(0),pos(1),1);
  // Source, destination and direction in world
  Vector3f s,d,dir;
  s = igl::unproject(win_s,model,proj,viewport);
  d = igl::unproject(win_d,model,proj,viewport);
  dir = d-s;

  // Shoot ray, collect all hits (could just collect first two)
  int num_rays_shot;
  hits.clear();
  ei.intersectRay(s,dir,hits,num_rays_shot);
  switch(hits.size())
  {
    case 0:
      break;
    case 1:
    {
      obj = (s + dir*hits[0].t).cast<typename Derivedobj::Scalar>();
      break;
    }
    case 2:
    default:
    {
      obj = 0.5*((s + dir*hits[0].t) + (s + dir*hits[1].t)).cast<typename Derivedobj::Scalar>();
      break;
    }
  }
  return hits.size();
}

template <typename Derivedobj>
IGL_INLINE int igl::embree::unproject_in_mesh(
    const Eigen::Vector2f& pos,
    const Eigen::Matrix4f& model,
    const Eigen::Matrix4f& proj,
    const Eigen::Vector4f& viewport,
    const EmbreeIntersector & ei,
    Eigen::PlainObjectBase<Derivedobj> & obj)
{
  std::vector<igl::embree::Hit> hits;
  return unproject_in_mesh(pos,model,proj,viewport,ei,obj,hits);
}


#ifdef IGL_STATIC_LIBRARY
template int igl::embree::unproject_in_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, igl::embree::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::vector<igl::embree::Hit, std::allocator<igl::embree::Hit> >&);
template int igl::embree::unproject_in_mesh<Eigen::Matrix<double, 1, 3, 1, 1, 3> >(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, igl::embree::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&, std::vector<igl::embree::Hit, std::allocator<igl::embree::Hit> >&);
template int igl::embree::unproject_in_mesh<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, igl::embree::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
#endif
