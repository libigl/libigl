#include "test_common.h"
#include <igl/fast_find_self_intersections.h>
#include <igl/combine.h>
#include <igl/sortrows.h>
#include <igl/unique.h>
#include <igl/matlab_format.h>
#include <igl/PI.h>
#include <iostream>
#include <iomanip>


TEST_CASE("fast_find_self_intersections: cube", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  REQUIRE( !igl::fast_find_self_intersections(V,F,true,true,IF,EV,EE,EI) );

  REQUIRE(IF.rows()==0);
  REQUIRE(EV.rows()==0);
  REQUIRE(EE.rows()==0);
  REQUIRE(EI.size()==0);
}

TEST_CASE("fast_find_self_intersections: cube-triangle", "[igl]")
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


  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  REQUIRE( igl::fast_find_self_intersections(V_,F_,false,false,IF,EV,EE,EI) );
  {
    Eigen::VectorXi I;
    igl::unique(IF,I);
    // should have 9 triangles that are intersecting each other
    REQUIRE( I.size()==9 );
    // plane that we added intersects other triangles
    REQUIRE((IF.col(1).array()==(F_.rows()-1)).all() );
  }
  // and 8 edges
  REQUIRE( EE.rows()==8 );
  REQUIRE( EV.rows()==2*8 );
  REQUIRE( EI.size()==8 );
}


TEST_CASE("fast_find_self_intersections: rose", "[igl]")
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

  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  REQUIRE( igl::fast_find_self_intersections(V,F,false,false,IF,EV,EE,EI) );
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

TEST_CASE("fast_find_self_intersections: shared-edge", "[igl]")
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
  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  REQUIRE(!igl::fast_find_self_intersections(V,F,false,false,IF,EV,EE,EI));
  V.row(3) << 2.0,0,0;
  REQUIRE( igl::fast_find_self_intersections(V,F,false,false,IF,EV,EE,EI));
  REQUIRE( IF.rows()==1 );
  REQUIRE( EE.rows()==0 );
  REQUIRE( EV.rows()==0 );
  REQUIRE( EI.rows()==0 );
}
