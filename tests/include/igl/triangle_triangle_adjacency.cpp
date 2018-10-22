
#include <test_common.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Geometry>


TEST_CASE("triangle_triangle_adjacency: dot", "[igl]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F,TT,TTi;
    // Load example mesh: GetParam() will be name of mesh file
    test_common::load_mesh(param, V, F);
    igl::triangle_triangle_adjacency(F,TT,TTi);
    REQUIRE (TT.rows() == F.rows());
    REQUIRE (TTi.rows() == F.rows());
    REQUIRE (TT.cols() == F.cols());
    REQUIRE (TTi.cols() == F.cols());
    for(int f = 0;f<F.rows();f++)
    {
      for(int c = 0;c<3;c++)
      {
        if(TT(f,c) >= 0)
        {
          REQUIRE (F.rows() > TT(f,c));
          REQUIRE (0 <= TTi(f,c));
          REQUIRE (3 > TTi(f,c));
          REQUIRE (f == TT(TT(f,c),TTi(f,c)));
        }
      }
    }
    // REQUIRE (b == a);
    // REQUIRE (a==b);
    // REQUIRE(a == Approx(b).margin(1e-15))
    // REQUIRE (1e-12 > a);
  };

  test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}
