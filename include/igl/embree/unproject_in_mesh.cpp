#include "unproject_in_mesh.h"
#include "EmbreeIntersector.h"
#include <igl/unproject.h>
#include <vector>

template <
  typename Scalar,
  typename Index,
  typename Derivedobj>
int igl::unproject_in_mesh(
  const int x,
  const int y,
  const igl::EmbreeIntersector<Scalar,Index> & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj)
{
  std::vector<igl::Hit> hits;
  return igl::unproject_in_mesh(x,y,ei,obj,hits);
}

template <
  typename Scalar,
  typename Index,
  typename Derivedobj>
int igl::unproject_in_mesh(
  const int x,
  const int y,
  const igl::EmbreeIntersector<Scalar,Index> & ei,
  Eigen::PlainObjectBase<Derivedobj> & obj,
  std::vector<igl::Hit > & hits)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  // Source and direction on screen
  Vector3d win_s = Vector3d(x,y,0);
  Vector3d win_d(x,y,1);
  // Source, destination and direction in world
  Vector3d s,d,dir;
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
      obj = s + dir*hits[0].t;
      break;
    }
    case 2:
    default:
    {
      obj = 0.5*((s + dir*hits[0].t) + (s + dir*hits[1].t));
      break;
    }
  }
  return hits.size();
}

#ifndef IGL_HEADER_ONLY
#endif
