#include <test_common.h>
#include <Eigen/Dense>

#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>
#include <igl/copyleft/cgal/assign.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

TEST_CASE("remesh_self_intersections: CubeWithFold", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path("cube_with_fold.ply"), V, F);

    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
    typedef Eigen::Matrix<K::FT, Eigen::Dynamic, Eigen::Dynamic> MatrixXe;

    MatrixXe VV;
    Eigen::MatrixXi FF, IF;
    Eigen::VectorXi J, IM;
    igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
    igl::copyleft::cgal::remesh_self_intersections(V, F, param, VV, FF, IF, J, IM);
}

TEST_CASE(
  "remesh_self_intersections: exact",
  "[igl/copyleft/cgal/]")
{
  const auto test_case = [](const std::string & param)
  {
    using namespace igl::copyleft::cgal;
    // Load and cast to exact
    typedef CGAL::Epeck::FT EScalar;
    typedef Eigen::Matrix<EScalar,Eigen::Dynamic,3> MatrixXE;
    Eigen::MatrixXi F;
    MatrixXE Ve;
    {
      Eigen::MatrixXd Vd;
      igl::read_triangle_mesh(test_common::data_path(param),Vd,F);
      igl::copyleft::cgal::assign(Vd,Ve);
    }
    MatrixXE Vsu;
    Eigen::MatrixXi Fsu;
    // resolve intersections
    {
      RemeshSelfIntersectionsParam params(false,false,true);
      Eigen::MatrixXi IF;
      Eigen::VectorXi Jsu,IM;
      remesh_self_intersections(Ve,F,params,Vsu,Fsu,IF,Jsu,IM);
    }
    // detect intersections (there should be none)
    {
      RemeshSelfIntersectionsParam params(true,false,false);
      MatrixXE _VV;
      Eigen::MatrixXi _FF, IF;
      Eigen::VectorXi _J,_IM;
      remesh_self_intersections(Vsu,Fsu,params,_VV,_FF,IF,_J,_IM);
      REQUIRE(IF.rows() == 0);
    }
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}

