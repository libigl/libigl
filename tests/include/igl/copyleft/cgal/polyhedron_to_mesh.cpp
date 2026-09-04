#include <test_common.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/join_coplanar_neighboring_facets.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

TEST_CASE(
  "igl_copyleft_cgal_polyhedron_to_mesh: cube",
  "[igl/copyleft/cgal/]")
{
  Eigen::MatrixXd V(8,3);
  V << 0,0,0,
       1,0,0,
       1,1,0,
       0,1,0,
       0,0,1,
       1,0,1,
       1,1,1,
       0,1,1;
  Eigen::MatrixXi F(12,3);
  F <<
    0,1,2,
    0,2,3,
    4,7,5,
    5,7,6,
    0,4,1,
    1,4,5,
    1,5,2,
    2,5,6,
    2,6,3,
    3,6,7,
    3,7,0,
    0,7,4;
  CGAL::Polyhedron_3<
    CGAL::Simple_cartesian<double>, 
    CGAL::Polyhedron_items_with_id_3, 
    CGAL::HalfedgeDS_default, 
    std::allocator<int> > 
    poly;
  REQUIRE( igl::copyleft::cgal::mesh_to_polyhedron(V,F,poly));
  Eigen::MatrixXd V_out;
  Eigen::MatrixXi F_out;
  igl::copyleft::cgal::polyhedron_to_mesh(poly,V_out,F_out);
  REQUIRE( V_out.rows() == 8 );
  REQUIRE( F_out.rows() == 12 );
  test_common::assert_eq(V,V_out);
  test_common::assert_eq(F,F_out);
  igl::copyleft::cgal::join_coplanar_neighboring_facets(poly);
  REQUIRE( poly.size_of_facets() == 6 );
  Eigen::VectorXi I,C;
  igl::copyleft::cgal::polyhedron_to_mesh(poly,V_out,I,C);
  REQUIRE( V_out.rows() == 8 );
  REQUIRE( C.size() == 6+1 );
  auto sizes = (C.tail(C.size()-1) - C.head(C.size()-1)).eval();
  REQUIRE((sizes.array() == 4).all());
}
