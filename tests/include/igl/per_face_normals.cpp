
#include <test_common.h>
#include <igl/per_face_normals.h>
#include <Eigen/Geometry>

class per_face_normals : public ::testing::TestWithParam<std::string> {};

TEST_P(per_face_normals, dot)
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  // Load example mesh: GetParam() will be name of mesh file
  test_common::load_mesh(GetParam(), V, F);
  igl::per_face_normals(V,F,N);
  ASSERT_EQ(F.rows(),N.rows());
  for(int f = 0;f<N.rows();f++)
  {
    for(int c = 0;c<3;c++)
    {
      // Every half-edge dot the normal should be 0
      ASSERT_LT(
        std::abs((V.row(F(f,c))-V.row(F(f,(c+1)%3))).dot(N.row(f))),
        1e-12);
    }
  }
  // ASSERT_EQ(a,b);
  // ASSERT_TRUE(a==b);
  // ASSERT_NEAR(a,b,1e-15)
  // ASSERT_LT(a,1e-12);
}

INSTANTIATE_TEST_CASE_P
(
  all_meshes,
  per_face_normals,
  ::testing::ValuesIn(test_common::all_meshes()),
  test_common::string_test_name
);
