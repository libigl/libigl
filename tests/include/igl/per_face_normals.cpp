
#include <test_common.h>
#include <igl/per_face_normals.h>
#include <Eigen/Geometry>

TEST_CASE("per_face_normals: dot", "[igl]")
{
  const auto test_case = [](const std::string &param)
  {
	  Eigen::MatrixXd V,N;
	  Eigen::MatrixXi F;
	  // Load example mesh: GetParam() will be name of mesh file
	  test_common::load_mesh(param, V, F);
	  igl::per_face_normals(V,F,N);
	  REQUIRE (N.rows() == F.rows());
	  for(int f = 0;f<N.rows();f++)
	  {
	    for(int c = 0;c<3;c++)
	    {
	      // Every half-edge dot the normal should be 0
	      REQUIRE(std::abs((V.row(F(f,c))-V.row(F(f,(c+1)%3))).dot(N.row(f))) < 1e-12);
	    }
	  }
	  // REQUIRE (b == a);
	  // REQUIRE (a==b);
	  // ASSERT_NEAR(a,b,1e-15)
	  // REQUIRE (1e-12 > a);
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}
