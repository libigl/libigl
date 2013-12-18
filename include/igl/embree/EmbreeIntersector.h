// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
// igl function interface for Embree2.0
//
// Necessary changes to switch from previous Embree versions:
// * Use igl:Hit instead of embree:Hit (where id0 -> id)
// * Embree2.0 finds now also back face intersections

#ifndef IGL_EMBREE_INTERSECTOR_H
#define IGL_EMBREE_INTERSECTOR_H

#include "Hit.h"
#include <Eigen/Core>
#include "Embree_convenience.h"
#include <vector>

namespace igl
{
  class EmbreeIntersector
  {
  public:
    // Initialize embree engine. This will be called on instance `init()`
    // calls. If already inited then this function does nothing: it is harmless
    // to call more than once.
    static inline void global_init();
  private:
    // Deinitialize the embree engine. This should probably never be called by
    // the user. Hence it's private. Do you really want to do this?
    static inline void global_deinit();
  public:
    typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> PointMatrixType;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>  FaceMatrixType;
    typedef Eigen::Matrix<float,1,3> RowVector3;
  public:
    inline EmbreeIntersector();
  private:
    // Copying and assignment are not allowed.
    inline EmbreeIntersector(const EmbreeIntersector & that);
    inline EmbreeIntersector & operator=(const EmbreeIntersector &);
  public:
    virtual inline ~EmbreeIntersector();
      
    // Initialize with a given mesh.
    //
    // Inputs:
    //   V  #V by 3 list of vertex positions
    //   F  #F by 3 list of Oriented triangles
    // Side effects:
    //   The first time this is ever called the embree engine is initialized.
    inline void init(
      const PointMatrixType & V,
      const FaceMatrixType & F,
      const char* structure = "default",
      const char* builder = "default",
      const char* traverser = "default");
    // Deinitialize embree datasctructures for current mesh.  Also called on
    // destruction: no need to call if you just want to init() once and
    // destroy.
    inline void deinit();
  
    // Given a ray find the first hit
    // 
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    // Output:
    //   hit        information about hit
    // Returns true if and only if there was a hit
    inline bool intersectRay(
      const RowVector3& origin, 
      const RowVector3& direction,
      Hit& hit,
      float tnear = 0,
      float tfar = embree::inf) const;

    // Given a ray find the all hits in order
    // 
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    // Output:
    //   hit        information about hit
    //   num_rays   number of rays shot (at least one)
    // Returns true if and only if there was a hit
    inline bool intersectRay(
      const RowVector3& origin,
      const RowVector3& direction,
      std::vector<Hit > &hits,
      int& num_rays,
      float tnear = 0,
      float tfar = embree::inf) const;

    // Given a ray find the first hit
    // 
    // Inputs:
    //   a    3d first end point of segment
    //   ab   3d vector from a to other endpoint b
    // Output:
    //   hit  information about hit
    // Returns true if and only if there was a hit
    inline bool intersectSegment(const RowVector3& a, const RowVector3& ab, Hit &hit) const;
    
  private:
    embree::RTCGeometry* mesh;
    embree::RTCTriangle* triangles;
    embree::RTCVertex *vertices;
    embree::Intersector1 *intersector;
    bool initialized;
  };
}

// Implementation
#include <igl/EPS.h>
// This unfortunately cannot be a static field of EmbreeIntersector because it
// would depend on the template and then we might end up with initializing
// embree twice. If only there was a way to ask embree if it's already
// initialized...
namespace igl
{
  // Keeps track of whether the **Global** Embree intersector has been
  // initialized. This should never been done at the global scope.
  static bool EmbreeIntersector_inited = false;
}

template <typename RowVector3>
inline embree::Vector3f toVector3f(const RowVector3 &p) { return embree::Vector3f((float)p[0], (float)p[1], (float)p[2]); }

inline void igl::EmbreeIntersector::global_init()
{
  if(!EmbreeIntersector_inited)
  {
    embree::rtcInit();
#ifdef IGL_VERBOSE
    embree::rtcSetVerbose(3);
#endif
    embree::rtcStartThreads();
    EmbreeIntersector_inited = true;
  }
}

inline void igl::EmbreeIntersector::global_deinit()
{
  EmbreeIntersector_inited = false;
  embree::rtcStopThreads();
  embree::rtcExit();
  embree::rtcFreeMemory();
}

inline igl::EmbreeIntersector::EmbreeIntersector()
  :
  mesh(NULL),
  triangles(NULL),
  vertices(NULL),
  intersector(NULL),
  initialized(false)
{
}

inline igl::EmbreeIntersector::EmbreeIntersector(
  const EmbreeIntersector & /*that*/)
  :// To make -Weffc++ happy
  mesh(NULL),
  triangles(NULL),
  vertices(NULL),
  intersector(NULL),
  initialized(false)
{
  assert(false && "Copying EmbreeIntersector is not allowed");
}

inline igl::EmbreeIntersector & igl::EmbreeIntersector::operator=(
  const EmbreeIntersector & /*that*/)
{
  assert(false && "Assigning an EmbreeIntersector is not allowed");
  return *this;
}


inline void igl::EmbreeIntersector::init(
  const PointMatrixType & V,
  const FaceMatrixType & F,
  const char* structure,
  const char* builder,
  const char* traverser)
{
  
  if (initialized)
    deinit();
  
  using namespace std;
  global_init();

  if(V.size() == 0 || F.size() == 0)
  {
    return;
  }
  
  mesh = embree::rtcNewTriangleMesh(F.rows(),V.rows(),structure);

  // fill vertex buffer
  vertices = embree::rtcMapPositionBuffer(mesh);
  for(int i=0;i<(int)V.rows();i++)
  {
    vertices[i] = embree::RTCVertex((float)V(i,0),(float)V(i,1),(float)V(i,2));
  }
  embree::rtcUnmapPositionBuffer(mesh);

  // fill triangle buffer
  triangles = embree::rtcMapTriangleBuffer(mesh);
  for(int i=0;i<(int)F.rows();i++)
  {
    triangles[i] = embree::RTCTriangle((int)F(i,0),(int)F(i,1),(int)F(i,2),i);
  }
  embree::rtcUnmapTriangleBuffer(mesh);

  embree::rtcBuildAccel(mesh,builder);
  embree::rtcCleanupGeometry(mesh);
  
  intersector = embree::rtcQueryIntersector1(mesh,traverser);
  
  initialized = true;
}

igl::EmbreeIntersector
::~EmbreeIntersector()
{
  if (initialized)
    deinit();
}

void igl::EmbreeIntersector::deinit()
{
  embree::rtcDeleteIntersector1(intersector);
  embree::rtcDeleteGeometry(mesh);
}

inline bool igl::EmbreeIntersector::intersectRay(
  const RowVector3& origin,
  const RowVector3& direction,
  Hit& hit,
  float tnear,
  float tfar) const
{
  embree::Ray ray(toVector3f(origin), toVector3f(direction), tnear, tfar);
  intersector->intersect(ray);
  
  if(ray)
  {
    hit.id = ray.id0;
    hit.u = ray.u;
    hit.v = ray.v;
    hit.t = ray.tfar;
    return true;
  }

  return false;
}

inline bool 
igl::EmbreeIntersector
::intersectRay(
  const RowVector3& origin, 
  const RowVector3& direction,
  std::vector<Hit > &hits,
  int& num_rays,
  float tnear,
  float tfar) const
{
  using namespace std;
  num_rays = 0;
  hits.clear();
  int last_id0 = -1;
  double self_hits = 0;
  // This epsilon is directly correleated to the number of missed hits, smaller
  // means more accurate and slower
  //const double eps = DOUBLE_EPS;
  const double eps = FLOAT_EPS;
  double min_t = tnear;
  bool large_hits_warned = false;
  embree::Ray ray(toVector3f(origin),toVector3f(direction));

  while(true)
  {
    ray.tnear = min_t;
    ray.tfar = tfar;
    ray.id0 = -1;
    num_rays++;
    intersector->intersect(ray);
    if(ray)
    {
      // Hit self again, progressively advance
      if(ray.id0 == last_id0 || ray.tfar <= min_t)
      {
        // push min_t a bit more
        //double t_push = pow(2.0,self_hits-4)*(hit.t<eps?eps:hit.t);
        double t_push = pow(2.0,self_hits)*eps;
#ifdef IGL_VERBOSE
        cout<<"  t_push: "<<t_push<<endl;
#endif
        //o = o+t_push*d;
        min_t += t_push;
        self_hits++;
      }
      else
      {
        Hit hit;
        hit.id = ray.id0;
        hit.u = ray.u;
        hit.v = ray.v;
        hit.t = ray.tfar;
        hits.push_back(hit);
#ifdef IGL_VERBOSE
        cout<<"  t: "<<hit.t<<endl;
#endif
        // Instead of moving origin, just change min_t. That way calculations
        // all use exactly same origin values
        min_t = ray.tfar;

        // reset t_scale
        self_hits = 0;
      }
      last_id0 = ray.id0;
    }
    else
      break; // no more hits
    
    if(hits.size()>1000 && !large_hits_warned)
    {
      cerr<<"Warning: Large number of hits..."<<endl;
      cerr<<"[ ";
      for(vector<Hit>::iterator hit = hits.begin(); hit != hits.end();hit++)
      {
        cerr<<(hit->id+1)<<" ";
      }
      
      cerr.precision(std::numeric_limits< double >::digits10);
      cerr<<"[ ";
      
      for(vector<Hit>::iterator hit = hits.begin(); hit != hits.end(); hit++)
      {
        cerr<<(hit->t)<<endl;;
      }

      cerr<<"]"<<endl;
      large_hits_warned = true;

      return hits.empty();
    }
  }

  return hits.empty();
}

inline bool 
igl::EmbreeIntersector
::intersectSegment(const RowVector3& a, const RowVector3& ab, Hit &hit) const
{
  embree::Ray ray(toVector3f(a), toVector3f(ab), embree::zero, embree::one);
  intersector->intersect(ray);

  if(ray)
  {
    hit.id = ray.id0;
    hit.u = ray.u;
    hit.v = ray.v;
    hit.t = ray.tfar;
    return true;
  }

  return false;
}

#endif //EMBREE_INTERSECTOR_H
