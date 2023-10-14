#include "test_common.h"
#include <igl/fast_find_intersections.h>
#include <igl/fast_find_self_intersections.h>
#include <igl/combine.h>
#include <igl/unique.h>
#include <igl/matlab_format.h>
#include <iostream>



TEST_CASE("fast_find_intersections: cube-triangle", "[igl]")
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
  

  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  igl::fast_find_intersections(V,F,Vp,Fp,false,false,IF,EV,EE,EI);
  {
    Eigen::VectorXi I;
    igl::unique(IF,I);
    // should have 8 triangles that are intersecting plane (so 9 unique)
    REQUIRE( I.size()== 9 );
  }

  // all triangles from the cube intersect the same plane
  REQUIRE( (IF.col(1).array()==0).all() );
  
  // and 8 line edges
  REQUIRE( EE.rows()==8 );
  REQUIRE( EV.rows()==2*8 );
  REQUIRE( EI.size()==8 );
  // TODO: check if the edges are at the right place
}

TEST_CASE("fast_find_intersections: coplanar", "[igl]")
{
}
