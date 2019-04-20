#include <test_common.h>
#include <igl/grad_intrinsic.h>
#include <igl/triangulated_grid.h>
#include <igl/edge_lengths.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/EPS.h>

TEST_CASE("grad_intrinsic: laplace_grid", "[igl]")
{
  Eigen::MatrixXd V2;
  Eigen::MatrixXi F;
  igl::triangulated_grid(3,3,V2,F);
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(V2.rows(),3);
  V.topLeftCorner(V2.rows(),2) = V2;
  V.col(2) = V.col(1).array() * V.col(0).array() + V.col(1).array();
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  Eigen::MatrixXd l;
  igl::edge_lengths(V,F,l);
  Eigen::SparseMatrix<double> G;
  igl::grad_intrinsic(l,F,G);
  Eigen::VectorXd dblA;
  igl::doublearea(l,0,dblA);
  Eigen::SparseMatrix<double> GTAG =
    G.transpose() * (dblA.colwise().replicate(2).asDiagonal()) * G;
  test_common::assert_near(
    Eigen::MatrixXd(L),Eigen::MatrixXd(-0.5*GTAG),igl::EPS<double>());
}
