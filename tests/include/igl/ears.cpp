#include <test_common.h>
#include <igl/ears.h>
#include <igl/triangulated_grid.h>
#include <igl/is_boundary_edge.h>

void test_F(const Eigen::MatrixXi & F)
{
  Eigen::VectorXi ear,ear_opp;
  igl::ears(F,ear,ear_opp);

  Eigen::Array<bool,Eigen::Dynamic,1> B;
  Eigen::MatrixXi E;
  Eigen::VectorXi EMAP;
  igl::is_boundary_edge(F,B,E,EMAP);

  for(int e = 0;e<ear.size();e++)
  {
    int ue1 = EMAP(ear(e) + F.rows()*((ear_opp(e)+1)%3));
    int ue2 = EMAP(ear(e) + F.rows()*((ear_opp(e)+2)%3));
    REQUIRE(B(ue1));
    REQUIRE(B(ue2));
  }
}

TEST_CASE("ears: grid", "[igl]" )
{
  Eigen::MatrixXi F;
  igl::triangulated_grid(3,3,F);
  test_F(F);
}

TEST_CASE("ears: two-boundary", "[igl]" )
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Load example mesh: GetParam() will be name of mesh file
    igl::read_triangle_mesh(test_common::data_path(param), V, F);
    test_F(F);
  };

  test_common::run_test_cases(
      {
      "bunny.off",
      "elephant.off",
      "hemisphere.obj"
      }, test_case);

}
