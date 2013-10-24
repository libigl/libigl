#ifndef IGL_EMBREE_INTERSECTOR_H
#define IGL_EMBREE_INTERSECTOR_H

#include "include/embree.h"
#include "include/intersector1.h"
#include "common/ray.h"
#include <vector>

namespace igl
{
  struct Hit
  {
    int id; // primitive id
    float u,v; // barycentric coordinates
    float t; // distance = direction*t to intersection
  };

  template <
  typename Scalar,
  typename Index>
  class EmbreeIntersector
  {
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> PointMatrixType;
    typedef Eigen::Matrix<Index,Eigen::Dynamic,Eigen::Dynamic>  FaceMatrixType;
    typedef Eigen::Matrix<Scalar,1,3> RowVector3;
  public:
    // V  #V by 3 list of vertex positions
    // F  #F by 3 list of Oriented triangles
    EmbreeIntersector();
    EmbreeIntersector(
      const PointMatrixType & V,
      const FaceMatrixType & F,
      const char* structure = "default",
      const char* builder = "default",
      const char* traverser = "default");
    virtual ~EmbreeIntersector();
  
    // Given a ray find the first hit
    // 
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    // Output:
    //   hit        information about hit
    // Returns true if and only if there was a hit
    bool intersectRay(
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
    bool intersectRay(
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
    bool intersectSegment(const RowVector3& a, const RowVector3& ab, Hit &hit) const;
    
  private:
    embree::RTCGeometry* mesh;
    embree::RTCTriangle* triangles;
    embree::RTCVertex *vertices;
    embree::Intersector1 *intersector;
  };
}

// Implementation
#include <igl/EPS.h>

template <typename RowVector3>
inline embree::Vector3f toVector3f(const RowVector3 &p) { return embree::Vector3f((float)p[0], (float)p[1], (float)p[2]); }

template <
typename Scalar,
typename Index>
igl::EmbreeIntersector < Scalar, Index>
::EmbreeIntersector()
{
  static bool inited = false;
  if(!inited)
  {
    embree::rtcInit();
    embree::rtcStartThreads();
    inited = true;
  }
}

template <
typename Scalar,
typename Index>
igl::EmbreeIntersector < Scalar, Index>
::EmbreeIntersector(const PointMatrixType & V,
                    const FaceMatrixType & F,
                    const char* structure,
                    const char* builder,
                    const char* traverser)
{
  static bool inited = false;
  if(!inited)
  {
    embree::rtcInit();
#ifdef VERBOSE
    embree::rtcSetVerbose(3);
#endif
    embree::rtcStartThreads();
    inited = true;
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
}

template <
typename Scalar,
typename Index>
igl::EmbreeIntersector < Scalar, Index>
::~EmbreeIntersector()
{
  embree::rtcDeleteIntersector1(intersector);
  embree::rtcDeleteGeometry(mesh);
//  embree::rtcStopThreads();
//  embree::rtcExit();
//  embree::rtcFreeMemory();
}

template <
typename Scalar,
typename Index>
bool 
igl::EmbreeIntersector< Scalar, Index>
::intersectRay(
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

template <
typename Scalar,
typename Index>
bool 
igl::EmbreeIntersector < Scalar, Index>
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
#ifdef VERBOSE
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
#ifdef VERBOSE
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

template <
typename Scalar,
typename Index>
bool 
igl::EmbreeIntersector < Scalar, Index>
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
