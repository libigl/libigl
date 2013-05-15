#ifndef EMBREE_INTERSECTOR_H
#define EMBREE_INTERSECTOR_H

#undef interface
#undef near
#undef far

#include "common/intersector.h"
#include "common/accel.h"
//#include "types.h"

template <
typename PointMatrixType,
typename FaceMatrixType,
typename RowVector3>
class EmbreeIntersector
{
public:
  EmbreeIntersector(const PointMatrixType & V, const FaceMatrixType & F);
  virtual ~EmbreeIntersector();
  
  bool intersectRay(const RowVector3& origin, const RowVector3& direction, embree::Hit &hit) const;
  bool intersectSegment(const RowVector3& a, const RowVector3& ab, embree::Hit &hit) const;
  
private:
  embree::BuildTriangle *triangles;
  embree::BuildVertex *vertices;
  embree::Ref<embree::Accel> _accel;
  embree::Ref<embree::Intersector> _intersector;
};

#endif //EMBREE_INTERSECTOR_H
