#include <test_common.h>
#include <igl/edge_flaps.h>

TEST_CASE("edge_flaps: verify", "[igl]" "[slow]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path(param), V, F);

    Eigen::MatrixXi efE,efEF,efEI;
    Eigen::VectorXi efEMAP;
    igl::edge_flaps(F,efE,efEMAP,efEF,efEI);
    REQUIRE (efEF.rows() == efE.rows());
    REQUIRE (2 == efE.cols());
    REQUIRE (efEF.cols() == efE.cols());
    // for each edge, make sure edge appears in face
    for(int e = 0;e<efE.rows();e++)
    {
      for(int fe = 0;fe<2;fe++)
      {
        const int f = efEF(e,fe);
        // index of corner
        const int c = efEI(e,fe);
        REQUIRE (f<F.rows());
        // only check if not on boundary
        if(f >= 0)
        {
          // Either efE(e,[1 2]) = [i,j] appears after vertex c of face f
          // Or  efE(e,[2 1]) = [j,i] appears after vertex c of face f
          CHECK((
            ((efE(e,0) == F(f,(c+1)%3)) && (efE(e,1) == F(f,(c+2)%3))) ||
            ((efE(e,1) == F(f,(c+1)%3)) && (efE(e,0) == F(f,(c+2)%3)))));
        }
      }
    }
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}
