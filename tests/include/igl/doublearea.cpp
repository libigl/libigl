#include <test_common.h>
#include <igl/doublearea.h>


TEST_CASE("doublearea: VF_vs_ABC", "[igl]")
{
    auto test_case = [](const std::string &param)
    {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      test_common::load_mesh(param, V, F);

      // Check that computing double area with (V,F) is the same as computing
      // double area with (V1,V2,V2)
      Eigen::VectorXd A1,A2;
      igl::doublearea(V,F,A1);
      Eigen::MatrixXd A(F.rows(),3);
      Eigen::MatrixXd B(F.rows(),3);
      Eigen::MatrixXd C(F.rows(),3);
      for(int f = 0;f<F.rows();f++)
      {
        A.row(f) = V.row(F(f,0));
        B.row(f) = V.row(F(f,1));
        C.row(f) = V.row(F(f,2));
      }
      igl::doublearea(A,B,C,A2);
      REQUIRE (A2.size() == A1.size());
      for(int a = 0;a<A1.size();a++)
      {
        REQUIRE (A2(a) == Approx (A1(a)).margin(1e-15));
      }
    };

    test_common::run_test_cases(test_common::all_meshes(), test_case);
}
