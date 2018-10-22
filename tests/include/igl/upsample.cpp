
#include <test_common.h>
#include <igl/upsample.h>


TEST_CASE("upsample: single_triangle", "[igl]")
{
  Eigen::MatrixXi NF_groundtruth(4,3);
  NF_groundtruth << 0,3,5 ,1,4,3 ,3,4,5 ,4,2,5;
  Eigen::MatrixXd NV_groundtruth(6,2);
  NV_groundtruth <<0,0 ,1,0 ,0,1 ,0.5,0 ,0.5,0.5 ,0,0.5;
  Eigen::MatrixXd S_groundtruth(6,3);
  S_groundtruth<<1,0,0 ,0,1,0 ,0,0,1 ,0.5,0.5,0 ,0,0.5,0.5 ,0.5,0,0.5;

  Eigen::MatrixXi F(1,3);
  F<<0,1,2;
  Eigen::MatrixXd V(3,2);
  V<<0,0,1,0,0,1;
  Eigen::MatrixXi NF;
  Eigen::MatrixXd NV;
  Eigen::SparseMatrix<double> S;
  igl::upsample(V.rows(),F,S,NF);
  test_common::assert_eq(NF_groundtruth,NF);
  test_common::assert_eq(S_groundtruth,Eigen::MatrixXd(S));
  igl::upsample(V,F,NV,NF);
  test_common::assert_eq(NF_groundtruth,NF);
  test_common::assert_eq(NV_groundtruth,NV);
}

TEST_CASE("upsample: V_comes_first_F_ordering", "[igl]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V,NV;
    Eigen::MatrixXi F,NF;
    // Load example mesh: GetParam() will be name of mesh file
    test_common::load_mesh(param, V, F);
    igl::upsample(V,F,NV,NF);
    REQUIRE (V.rows() <= NV.rows());
    REQUIRE (4*F.rows() == NF.rows());
    // V should be first part of V
    test_common::assert_eq(V,NV.topLeftCorner(V.rows(),V.cols()));
    // Expect a particular order
    for(int f = 0;f<F.rows();f++)
    {
      REQUIRE (NF((f*4)+0,0) == F(f,0));
      REQUIRE (NF((f*4)+1,0) == F(f,1));
      REQUIRE (NF((f*4)+3,1) == F(f,2));
    }
  };

  test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}
