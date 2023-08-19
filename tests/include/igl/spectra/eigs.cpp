#include <test_common.h>
#include <igl/spectra/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/triangulated_grid.h>

#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("eigs: grid", "[igl/spectra]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::triangulated_grid(10,10,V,F);
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  Eigen::VectorXd S_matlab(5,1);
  S_matlab<<
    -0.00000000000000600,
    -9.76979543268284445,
    -9.76979543268286221,
    -19.53959086536570311,
    -37.90080021472548566;
  const bool success = igl::spectra::eigs(L,M,S_matlab.size(),igl::EIGS_TYPE_SM,U,S);
  REQUIRE(success);


  test_common::assert_near(S,S_matlab,1e-4);
}

