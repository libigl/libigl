// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unproject_in_mesh.h"
#include "EmbreeIntersector.h"
#include <igl/unproject.h>
#include <vector>

#ifndef IGL_OPENGL_4
template <
  typename Derivedobj>
IGL_INLINE int igl::unproject_in_mesh(
  const double x,
  const double y,
  const igl::EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  std::vector<igl::Hit> hits;
  return igl::unproject_in_mesh(x,y,ei,obj,hits);
}

template <
  typename Derivedobj>
IGL_INLINE int igl::unproject_in_mesh(
  const double x,
  const double y,
  const igl::EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj,
  std::vector<igl::Hit > & hits)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  // Source and direction on screen
  Vector3f win_s = Vector3f(x,y,0);
  Vector3f win_d(x,y,1);
  // Source, destination and direction in world
  Vector3f s,d,dir;
  unproject(win_s,s);
  unproject(win_d,d);
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
#endif

template <typename Derivedobj>
IGL_INLINE int igl::unproject_in_mesh(
  const Eigen::Vector2f& pos,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj,
  std::vector<igl::Hit > & hits)
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


#ifdef IGL_STATIC_LIBRARY
#  ifndef IGL_OPENLGL_4
template int igl::unproject_in_mesh<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(double, double, igl::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&, std::vector<igl::Hit, std::allocator<igl::Hit> >&);
template int igl::unproject_in_mesh<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(double, double, igl::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
#  endif
template int igl::unproject_in_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 4, 0, 4, 4> const&, Eigen::Matrix<float, 4, 1, 0, 4, 1> const&, igl::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, std::vector<igl::Hit, std::allocator<igl::Hit> >&);
#endif
