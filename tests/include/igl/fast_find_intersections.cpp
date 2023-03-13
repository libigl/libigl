#include <test_common.h>
#include <igl/tri_tri_intersect.h>
#include <igl/fast_find_intersections.h>
#include <igl/fast_find_self_intersections.h>
#include <igl/combine.h>

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


TEST_CASE("fast_find_self_intersections: negative", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi I;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  REQUIRE (! igl::fast_find_self_intersections(V,F,I) );

  REQUIRE ( I.sum()==0);
}

TEST_CASE("fast_find_self_intersections: positive", "[igl]")
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

  REQUIRE ( igl::fast_find_self_intersections(V_,F_,I,edges) );
  // should have 9 triangles that are intersecting each other
  REQUIRE ( I.sum()==9);
  // and 16 line edges
  REQUIRE ( edges.rows()==16);
  // plane that we added intersects other triangles
  REQUIRE ( I(F_.rows()-1)!=0 );
  // TODO: check if the edges are at the right place (?) 
}

TEST_CASE("fast_find_intersections:", "[igl]")
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

  igl::fast_find_intersections(V,F,Vp,Fp,I,edges);

  // should have 8 triangles that are intersecting plane
  REQUIRE ( I.rows()==8);

  // all triangles from the cube intersect the same plane
  REQUIRE( (I.col(1).array()==0).all());
  
  // and 16 line edges
  REQUIRE ( edges.rows()==16);
  // TODO: check if the edges are at the right place
}
