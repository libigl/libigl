#include "test_common.h"
#include <igl/predicates/triangle_triangle_intersect.h>

TEST_CASE("triangle_triangle_intersect intersect", "[igl/predicates]")
{
  // try with oblique interecion plane
  for(double shift=-1000;shift<=1000;shift+=100.0)
  {
    Eigen::RowVector3d p1(0,0,0),    q1(1,0,0),r1(0,1,0);
    Eigen::RowVector3d p2(shift,0,1),q2(1,1,0),r2(-shift,0,-1);

    bool coplanar;
    REQUIRE( igl::predicates::triangle_triangle_intersect(p1,q1,r1, p2,q2,r2,coplanar));
  }
}

TEST_CASE("triangle_triangle_intersect not intersect", "[igl/predicates]")
{
  // positive test that two triangles intersect
  Eigen::RowVector3d p1(0,0,0),q1(1,0,0),r1(0,1,0);
  Eigen::RowVector3d p2(0,0,1),q2(1,0,1),r2(0,1,1);

  bool coplanar;
  REQUIRE( !igl::predicates::triangle_triangle_intersect(p1,q1,r1,p2,q2,r2,coplanar) );
}


TEST_CASE("triangle_triangle_intersect coplanar", "[igl/predicates]")
{
  // positive test that two triangles intersect
  Eigen::RowVector3d p1(0,0,0),q1(1,0,0),r1(0,1,0);
  Eigen::RowVector3d p2(0,0,0),q2(0.5,0,0),r2(0,0.5,0);

  // should intersect along the vector (0,0,0) -> (sqrt(2),sqrt(2),0)
  bool coplanar;
  REQUIRE( igl::predicates::triangle_triangle_intersect(p1,q1,r1,p2,q2,r2,coplanar) );
  REQUIRE( coplanar );
}

TEST_CASE("triangle_triangle_intersect: non-intersecting", "[igl]")
{
  Eigen::Matrix<double,Eigen::Dynamic,3> V(6,3);
  V<< 
    0.39234799146652199,0.38443601131439198,0.44744500517845198,
    0.38924700021743752,0.385637506842613,0.45762349665164948,
    0.38700349628925301,0.38789276033639897,0.45634675025939975,
    0.39079749584197976,0.38503675907850249,0.45253425091505073,
    0.38769650459289529,0.38623825460672351,0.46271274238824822,
    0.39279749989509577,0.37299175560474374,0.45553924888372399;
  Eigen::MatrixXi F(2,3);
  F<< 
    0,3,2,
    1,5,4;

  bool coplanar;
  REQUIRE(!igl::predicates::triangle_triangle_intersect(
        V.row(F(0,0)).eval(),
        V.row(F(0,1)).eval(),
        V.row(F(0,2)).eval(),
        V.row(F(1,0)).eval(),
        V.row(F(1,1)).eval(),
        V.row(F(1,2)).eval(),
        coplanar));
}


TEST_CASE("triangle_triangle_intersect: coplanar", "[igl]")
{
  Eigen::Matrix<double,Eigen::Dynamic,3> V(6,3);
  V<< 
  0.30947001278400399,0.80785250663757346,0.47595100104808802,
  0.299046009778976,0.79801350831985496,0.47843150794506101,
  0.32418551295995701,0.79794725775718722,0.48203899711370451,
  0.30425801128148999,0.80293300747871421,0.47719125449657451,
  0.31897351145744302,0.79302775859832797,0.48327925056219101,
  0.31376150995492902,0.78810825943946872,0.4845195040106775;
  Eigen::MatrixXi F(2,3);
  F<< 
    0,4,2,
    1,5,3;

  bool coplanar;
  REQUIRE(!igl::predicates::triangle_triangle_intersect(
        V.row(F(0,0)).eval(),
        V.row(F(0,1)).eval(),
        V.row(F(0,2)).eval(),
        V.row(F(1,0)).eval(),
        V.row(F(1,1)).eval(),
        V.row(F(1,2)).eval(),coplanar));
}
