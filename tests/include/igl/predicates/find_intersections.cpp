#include "test_common.h"
#include <igl/predicates/find_intersections.h>
#include <igl/combine.h>
#include <igl/triangle_triangle_intersect.h>
#include <igl/unique.h>
#include <igl/sortrows.h>
#include <igl/matlab_format.h>
#include <iostream>



TEST_CASE("find_intersections: cube-triangle", "[igl/predicates]")
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
  

  Eigen::MatrixXi IF;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  igl::predicates::find_intersections(V,F,Vp,Fp,false,IF,CP);
  {
    Eigen::VectorXi I;
    igl::unique(IF,I);
    // should have 8 triangles that are intersecting plane (so 9 unique)
    REQUIRE( I.size()== 9 );
  }

  // all triangles from the cube intersect the same plane
  REQUIRE( (IF.col(1).array()==0).all() );
}

TEST_CASE("find_intersections: extract", "[igl/predicates]")
{
  Eigen::MatrixXd V1(4,3);
  V1 << 0,0,0,
       1,0,0,
       1,1,0,
       0,1,0;
  Eigen::MatrixXi F1(2,3);
  F1 << 0,1,2,
       0,2,3;
  Eigen::MatrixXd V2(3,3);
  V2 << 1,0,1,
       0,1,1,
       0.5,0.5,-1;
  Eigen::MatrixXi F2(1,3);
  F2 << 0,1,2;

  Eigen::MatrixXi IF,EE; Eigen::MatrixXd EV;
  Eigen::Array<bool,Eigen::Dynamic,1> CP;
  bool found = igl::predicates::find_intersections(V1,F1,V2,F2,false,IF,CP);
  REQUIRE( (CP==false).all() );
  REQUIRE( found );
  REQUIRE( IF.rows() == 2);
  Eigen::MatrixXi IF_gt(2,2);
  IF_gt << 
    0,0,
    1,0;
  test_common::assert_near_rows(IF,IF_gt,0);

  igl::triangle_triangle_intersect(V1,F1,V2,F2,IF,EV,EE);
  Eigen::MatrixXd EV_gt(3,3);
  EV_gt<<
    0.5 ,0.5,0,
    0.25,0.75,0,
    0.75,0.25,0;
  test_common::assert_near_rows(EV,EV_gt,0);
  REQUIRE( EE.rows() == 2);
}
