#include <test_common.h>
#include <igl/fast_find_mesh_selfintersect.h>
#include <igl/fast_find_mesh_intersect.h>
#include <igl/combine.h>

TEST_CASE("fast_find_mesh_selfintersect: negative", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi I;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  REQUIRE (! igl::fast_find_mesh_selfintersect(V,F,I) );

  REQUIRE ( I.sum()==0);
}

TEST_CASE("fast_find_mesh_selfintersect: positive", "[igl]")
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

  Eigen::VectorXi I;
  Eigen::MatrixXd edges;

  REQUIRE ( igl::fast_find_mesh_selfintersect(V_,F_,I,edges) );
  // should have 9 triangles that are intersecting each other
  REQUIRE ( I.sum()==9);
  // and 16 line edges
  REQUIRE ( edges.rows()==16);
  // plane that we added intersects other triangles
  REQUIRE ( I(F_.rows()-1)!=0 );
  // TODO: check if the edges are at the right place (?) 
}

TEST_CASE("fast_find_mesh_intersect:", "[igl]")
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
  
  Eigen::MatrixXi I;
  Eigen::MatrixXd edges;

  igl::fast_find_mesh_intersect(V,F,Vp,Fp,I,edges);

  // should have 8 triangles that are intersecting plane
  REQUIRE ( I.rows()==8);

  // all triangles from the cube intersect the same plane
  REQUIRE( (I.col(1).array()==0).all());
  
  // and 16 line edges
  REQUIRE ( edges.rows()==16);
  // TODO: check if the edges are at the right place

}
