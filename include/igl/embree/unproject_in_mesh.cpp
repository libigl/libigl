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

template <
  typename Derivedobj>
int igl::unproject_in_mesh(
  const int x,
  const int y,
  const igl::EmbreeIntersector & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  std::vector<igl::Hit> hits;
  return igl::unproject_in_mesh(x,y,ei,obj,hits);
}

template <
  typename Derivedobj>
int igl::unproject_in_mesh(
  const int x,
  const int y,
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

#ifndef IGL_HEADER_ONLY
template int igl::unproject_in_mesh<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(int, int, igl::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&, std::vector<igl::Hit, std::allocator<igl::Hit> >&);
template int igl::unproject_in_mesh<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(int, int, igl::EmbreeIntersector const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
#endif
