#include "test_common.h"
#include <igl/tri_tri_intersect.h>

TEST_CASE("tri_tri_intersection_test_3d intersect", "[igl]")
{
  // try with oblique interecion plane
  for(double shift=-1000;shift<=1000;shift+=100.0)
  {
    Eigen::RowVector3d p1(0,0,0),    q1(1,0,0),r1(0,1,0);
    Eigen::RowVector3d p2(shift,0,1),q2(1,1,0),r2(-shift,0,-1);

    // should intersect along the vector (0,0,0) -> (0.5,0.5,0)
    Eigen::RowVector3d s,t;
    bool coplanar;
    REQUIRE( igl::tri_tri_intersection_test_3d(p1,q1,r1, p2,q2,r2, coplanar, s, t) );
    REQUIRE( !coplanar);
    
    if(s.norm()<1e-5)
    {
      Eigen::RowVector3d t_correct(0.5,0.5,0);
      Eigen::RowVector3d s_correct(0,0,0);
      test_common::assert_near( t, t_correct, 1e-10);
      test_common::assert_near( s, s_correct, 1e-10);
    } else {
      Eigen::RowVector3d s_correct(0.5,0.5,0);
      Eigen::RowVector3d t_correct(0,0,0);
      test_common::assert_near( t, t_correct, 1e-10);
      test_common::assert_near( s, s_correct, 1e-10);
    }
  }
}

TEST_CASE("tri_tri_intersection_test_3d not intersect", "[igl]")
{
  // positive test that two triangles intersect
  Eigen::RowVector3d p1(0,0,0),q1(1,0,0),r1(0,1,0);
  Eigen::RowVector3d p2(0,0,1),q2(1,0,1),r2(0,1,1);

  // should intersect along the vector (0,0,0) -> (sqrt(2),sqrt(2),0)
  Eigen::RowVector3d s,t;
  bool coplanar;
  REQUIRE( !igl::tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2, coplanar, s, t) );
}


TEST_CASE("tri_tri_intersection_test_3d coplanar", "[igl]")
{
  // positive test that two triangles intersect
  Eigen::RowVector3d p1(0,0,0),q1(1,0,0),r1(0,1,0);
  Eigen::RowVector3d p2(0,0,0),q2(0.5,0,0),r2(0,0.5,0);

  // should intersect along the vector (0,0,0) -> (sqrt(2),sqrt(2),0)
  Eigen::RowVector3d s,t;
  bool coplanar;
  REQUIRE( igl::tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2, coplanar, s, t) );
  REQUIRE(coplanar);
}
