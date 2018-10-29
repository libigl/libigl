#include <test_common.h>
#include <igl/intrinsic_delaunay_cotmatrix.h>
#include <igl/EPS.h>
#include <igl/triangulated_grid.h>
#include <igl/is_border_vertex.h>

class intrinsic_delaunay_cotmatrix : public ::testing::TestWithParam<std::string> {};

TEST(intrinsic_delaunay_cotmatrix,skewed_grid)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  const int s = 7;
  igl::triangulated_grid(s,s,V,F);
  // Skew against diagonal direction
  V.col(0) -= 1.1*V.col(1);
  Eigen::SparseMatrix<double> L;
  Eigen::MatrixXi F_intrinsic;
  Eigen::MatrixXd l_intrinsic;
  igl::intrinsic_delaunay_cotmatrix(V,F,L,l_intrinsic,F_intrinsic);
  Eigen::VectorXi LI,LJ;
  Eigen::VectorXd LV;
  igl::find(L,LI,LJ,LV);
  const auto is_boundary_edge = [](const int i, const int j, const int s)->bool
  {
    const auto is_boundary_vertex = [](const int i, const int s)->bool
    {
      return (i<s) || (i>=s*s-s) || (i%s == 0) || (i%s == s-1);
    };
    return is_boundary_vertex(i,s) && is_boundary_vertex(j,s);
  };
  // Off diagonals should be all non-positive
  for(int k = 0;k<LI.size();k++)
  {
    if(LI(k) != LJ(k) && !is_boundary_edge(LI(k),LJ(k),s))
    {
      ASSERT_GT(LV(k),-igl::EPS<double>());
    }
  }
}

TEST_P(intrinsic_delaunay_cotmatrix,manifold_meshes)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  test_common::load_mesh(GetParam(), V, F);
  Eigen::SparseMatrix<double> L;
  Eigen::MatrixXi F_intrinsic;
  Eigen::MatrixXd l_intrinsic;
  igl::intrinsic_delaunay_cotmatrix(V,F,L,l_intrinsic,F_intrinsic);
  Eigen::VectorXi LI,LJ;
  Eigen::VectorXd LV;
  igl::find(L,LI,LJ,LV);

  const std::vector<bool> is_boundary_vertex = igl::is_border_vertex(F);
  // Off diagonals should be all non-positive
  for(int k = 0;k<LI.size();k++)
  {
    if(LI(k) != LJ(k) && 
      !(is_boundary_vertex[LI(k)] && is_boundary_vertex[LJ(k)]))
    {
      ASSERT_GT(LV(k),-igl::EPS<double>());
    }
  }
}

INSTANTIATE_TEST_CASE_P
(
 manifold_meshes,
 intrinsic_delaunay_cotmatrix,
 ::testing::ValuesIn(test_common::manifold_meshes()),
 test_common::string_test_name
);
