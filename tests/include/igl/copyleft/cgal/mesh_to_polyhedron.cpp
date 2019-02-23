#include <test_common.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

TEST_CASE(
  "igl_copyleft_cgal_mesh_to_polyhedron: positive",
  "[igl/copyleft/cgal/]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh(param, V, F);
    CGAL::Polyhedron_3<
      CGAL::Simple_cartesian<double>, 
      CGAL::Polyhedron_items_with_id_3, 
      CGAL::HalfedgeDS_default, 
      std::allocator<int> > 
      poly;
    REQUIRE ( igl::copyleft::cgal::mesh_to_polyhedron(V,F,poly) );
  };

  test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}

TEST_CASE(
  "igl_copyleft_cgal_mesh_to_polyhedron: negative",
  "[igl/copyleft/cgal/]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  test_common::load_mesh("truck.obj", V, F);
  CGAL::Polyhedron_3<
    CGAL::Simple_cartesian<double>, 
    CGAL::Polyhedron_items_with_id_3, 
    CGAL::HalfedgeDS_default, 
    std::allocator<int> > 
    poly;
  REQUIRE (! igl::copyleft::cgal::mesh_to_polyhedron(V,F,poly) );
}
