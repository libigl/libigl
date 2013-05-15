#include "EmbreeIntersector.h"

template <typename RowVector3>
inline embree::Vec3f toVec3f(const RowVector3 &p) { return embree::Vec3f((float)p[0], (float)p[1], (float)p[2]); }

template <
typename PointMatrixType,
typename FaceMatrixType,
typename RowVector3>
EmbreeIntersector < PointMatrixType, FaceMatrixType, RowVector3>
::EmbreeIntersector(const PointMatrixType & V, const FaceMatrixType & F)
{
  static bool inited = false;
  if(!inited)
  {
	//embree::TaskScheduler::start();//init();
    inited = true;
  }
  
  size_t numVertices = 0;
  size_t numTriangles = 0;
  
  triangles = (embree::BuildTriangle*) embree::rtcMalloc(sizeof(embree::BuildTriangle) * F.rows());
  vertices = (embree::BuildVertex*) embree::rtcMalloc(sizeof(embree::BuildVertex) * V.rows());

  for(int i = 0; i < (int)V.rows(); ++i)
  {
    vertices[numVertices++] = embree::BuildVertex((float)V(i,0),(float)V(i,1),(float)V(i,2));
  }
  
  for(int i = 0; i < (int)F.rows(); ++i)
  {
    triangles[numTriangles++] = embree::BuildTriangle((int)F(i,0),(int)F(i,1),(int)F(i,2),i);
  }
  
  _accel = embree::rtcCreateAccel("default", "default", triangles, numTriangles, vertices, numVertices);
  _intersector = _accel->queryInterface<embree::Intersector>();
}

template <
typename PointMatrixType,
typename FaceMatrixType,
typename RowVector3>
EmbreeIntersector < PointMatrixType, FaceMatrixType, RowVector3>
::~EmbreeIntersector()
{
	embree::rtcFreeMemory();
}

template <
typename PointMatrixType,
typename FaceMatrixType,
typename RowVector3>
bool 
EmbreeIntersector < PointMatrixType, FaceMatrixType, RowVector3>
::intersectRay(const RowVector3& origin, const RowVector3& direction, embree::Hit &hit) const
{
	embree::Ray ray(toVec3f(origin), toVec3f(direction), 1e-4f);
	_intersector->intersect(ray, hit);
	return hit ; 
}

template <
typename PointMatrixType,
typename FaceMatrixType,
typename RowVector3>
bool 
EmbreeIntersector < PointMatrixType, FaceMatrixType, RowVector3>
::intersectSegment(const RowVector3& a, const RowVector3& ab, embree::Hit &hit) const
{
	embree::Ray ray(toVec3f(a), toVec3f(ab), embree::zero, embree::one);
	_intersector->intersect(ray, hit);
	return hit ; 
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
#include <Eigen/Core>
template class EmbreeIntersector<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 3, 1, 0, 3, 1> >;
#endif
