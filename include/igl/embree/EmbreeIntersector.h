#ifndef IGL_EMBREE_INTERSECTOR_H
#define IGL_EMBREE_INTERSECTOR_H

#undef interface
#undef near
#undef far

#include "common/intersector.h"
#include "common/accel.h"
#include <vector>
//#include "types.h"

namespace igl
{
  template <
  typename PointMatrixType,
  typename FaceMatrixType,
  typename RowVector3>
  class EmbreeIntersector
  {
  public:
    // V  #V by 3 list of vertex positions
    // F  #F by 3 list of Oriented triangles
    //
    // Note: this will only find front-facing hits. To consider all hits then
    // pass [F;fliplr(F)]
    EmbreeIntersector(const PointMatrixType & V, const FaceMatrixType & F);
    virtual ~EmbreeIntersector();
  
    // Given a ray find the first *front-facing* hit
    // 
    // Inputs:
    //   origin  3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    // Output:
    //   hit  embree information about hit
    // Returns true if and only if there was a hit
    bool intersectRay(
      const RowVector3& origin, 
      const RowVector3& direction,
      embree::Hit &hit) const;
    // Given a ray find the all *front-facing* hits in order
    // 
    // Inputs:
    //   origin  3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    // Output:
    //   hit  embree information about hit
    // Returns true if and only if there was a hit
    bool intersectRay(
      const RowVector3& origin, 
      const RowVector3& direction, 
      std::vector<embree::Hit > &hits) const;
  
    // Given a ray find the first *front-facing* hit
    // 
    // Inputs:
    //   a  3d first end point of segment
    //   ab  3d vector from a to other endpoint b
    // Output:
    //   hit  embree information about hit
    // Returns true if and only if there was a hit
    bool intersectSegment(const RowVector3& a, const RowVector3& ab, embree::Hit &hit) const;
    
  private:
    embree::BuildTriangle *triangles;
    embree::BuildVertex *vertices;
    embree::Ref<embree::Accel> _accel;
    embree::Ref<embree::Intersector> _intersector;
  };
}
#ifdef IGL_HEADER_ONLY
#  include "EmbreeIntersector.cpp"
#endif

#endif //EMBREE_INTERSECTOR_H
