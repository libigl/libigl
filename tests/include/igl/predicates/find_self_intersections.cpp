#include "test_common.h"
#include <igl/predicates/find_self_intersections.h>
#include <igl/upsample.h>
#include <igl/triangle_triangle_intersect.h>
#include <igl/combine.h>
#include <igl/sortrows.h>
#include <igl/unique.h>
#include <igl/matlab_format.h>
#include <igl/PI.h>
#include <iostream>
#include <iomanip>


TEST_CASE("find_self_intersections: cube", "[igl/predicates]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE( !igl::predicates::find_self_intersections(V,F,true,IF,CP) );
}

TEST_CASE("find_self_intersections: cube-triangle", "[igl/predicates]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  auto dims=(V.colwise().maxCoeff()-V.colwise().minCoeff())*2;
  auto cz=(V.col(2).minCoeff()+V.col(2).maxCoeff())/2;

  // add a triangle intersecting cube 
  Eigen::MatrixXd Vp(3,3);
  Vp << V.col(0).minCoeff()-dims(0),   V.col(1).minCoeff()-dims(1),   cz,
        V.col(0).minCoeff()-dims(0),   V.col(1).maxCoeff()+dims(1)*3, cz,
        V.col(0).maxCoeff()+dims(0)*3, V.col(1).maxCoeff()-dims(1),   cz;
  Eigen::MatrixXi Fp(1,3);
  Fp << 0,1,2;
  
  Eigen::MatrixXd V_;
  Eigen::MatrixXi F_;
  igl::combine<Eigen::MatrixXd,Eigen::MatrixXi>({V,Vp},{F,Fp}, V_,F_);

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE( igl::predicates::find_self_intersections(V_,F_,false,IF,CP) );
  {
    Eigen::VectorXi I;
    igl::unique(IF,I);
    // should have 9 triangles that are intersecting each other
    REQUIRE( I.size()==9 );
    // plane that we added intersects other triangles
    REQUIRE((IF.col(1).array()==(F_.rows()-1)).all() );
  }
}


TEST_CASE("find_self_intersections: rose", "[igl/predicates]")
{
  Eigen::MatrixXd V(10,3);
  Eigen::MatrixXi F(9,3);
  for(int i=0;i<9;i++)
  {
    const double theta_i = 4.0*igl::PI*double(i)/9.0;
    V.row(i) << std::cos(theta_i), std::sin(theta_i), 1;
    F.row(i)<<9,i,(i+1)%9;
  }
  V.row(9) << 0,0,0;

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE( igl::predicates::find_self_intersections(V,F,false,IF,CP));
  Eigen::MatrixXi IF_gt(9,2);
  IF_gt<<
    0,4,
    0,5,
    1,5,
    1,6,
    2,6,
    2,7,
    3,7,
    3,8,
    4,8;
  {
    Eigen::MatrixXi sIF;
    Eigen::VectorXi _;
    igl::sortrows(Eigen::MatrixXi(IF),true,sIF,_);
    test_common::assert_eq(IF_gt,sIF);
  }
}

TEST_CASE("find_self_intersections: shared-edge", "[igl/predicates]")
{
  Eigen::MatrixXd V(4,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    -1,0,0;
  Eigen::MatrixXi F(2,3);
  F <<
    0,1,2,
    0,2,3;
  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE(!igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );
  V.row(3) << 2.0,0,0;
  REQUIRE( igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==1 );
  REQUIRE( CP.rows()==1 );
  REQUIRE( CP(0) );
}

TEST_CASE("find_self_intersections: coplanar", "[igl/predicates]")
{
  Eigen::MatrixXd V(6,3);
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

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE( !igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );
}

TEST_CASE("find_self_intersections: non-intersecting", "[igl/predicates]")
{
  Eigen::MatrixXd V(6,3);
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

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE( !igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );
}

TEST_CASE("find_self_intersections: upsampled-knight", "[igl/predicates]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"),V,F);

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;

  REQUIRE( !igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );

  igl::upsample(Eigen::MatrixXd(V),Eigen::MatrixXi(F),V,F);


  REQUIRE( !igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );

  igl::upsample(Eigen::MatrixXd(V),Eigen::MatrixXi(F),V,F);

  REQUIRE( !igl::predicates::find_self_intersections(V,F,false,IF,CP));
  REQUIRE( IF.rows()==0 );
  REQUIRE( CP.rows()==0 );

}

TEST_CASE("find_self_intersections: extract", "[igl/predicates]")
{
  Eigen::MatrixXd V(4+3,3);
  V << 0,0,0,
       1,0,0,
       1,1,0,
       0,1,0,
       1,0,1,
       0,1,1,
       0.5,0.5,-1;
  Eigen::MatrixXi F(3,3);
  F << 0,1,2,
       0,2,3,
    4,5,6;
  Eigen::MatrixXi IF,EE; Eigen::MatrixXd EV;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  bool found = igl::predicates::find_self_intersections(V,F,false,IF,CP);

  igl::triangle_triangle_intersect(V,F,IF,EV,EE);
  REQUIRE( (CP==false).all() );
  REQUIRE( found );
  REQUIRE( IF.rows() == 2);
  Eigen::MatrixXi IF_gt(2,2);
  IF_gt << 
    0,2,
    1,2;
  test_common::assert_near_rows(IF,IF_gt,0);

  Eigen::MatrixXd EV_gt(3,3);
  EV_gt<<
    0.5 ,0.5,0,
    0.25,0.75,0,
    0.75,0.25,0;
  test_common::assert_near_rows(EV,EV_gt,0);
  REQUIRE( EE.rows() == 2);
  
}
