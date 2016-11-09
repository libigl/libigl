#ifndef IGL_COPYLEFT_CGAL_POINT_TRIANGLE_SQUARED_DISTANCE_H
#define IGL_COPYLEFT_CGAL_POINT_TRIANGLE_SQUARED_DISTANCE_H
#include <igl/igl_inline.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Point_3.h>
namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Given a point P1 and triangle T2 find the points on each of closest
      // approach and the squared distance thereof.
      // 
      // Inputs:
      //   P1  point
      //   T2  triangle
      // Outputs:
      //   P2  point on T2 closest to P1
      //   d  distance betwee P1 and T2
      template < typename Kernel>
      IGL_INLINE void point_triangle_squared_distance(
        const CGAL::Point_3<Kernel> & P1,
        const CGAL::Triangle_3<Kernel> & T2,
        CGAL::Point_3<Kernel> & P2,
        typename Kernel::FT & d
        );

    }
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "point_triangle_squared_distance.cpp"
#endif

#endif


