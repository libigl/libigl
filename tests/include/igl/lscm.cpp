#include <test_common.h>
#include <igl/grad.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/lscm.h>
#include <igl/EPS.h>

double compute_lscm_energy(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXd& W_flat)
{
  const int nV = V.rows();
  const int nF = F.rows();

  assert(V.cols() == 3);
  assert(F.cols() == 3);

  // Compute gradient
  Eigen::SparseMatrix<double> G;
  igl::grad(V, F, G);

  Eigen::VectorXd nablaU_flat = G * W_flat.head(nV);
  Eigen::VectorXd nablaV_flat = G * W_flat.tail(nV);

  Eigen::MatrixXd nablaU(nF, 3);
  Eigen::MatrixXd nablaV(nF, 3);
  nablaU << nablaU_flat.segment(0, nF), nablaU_flat.segment(nF, nF), nablaU_flat.segment(2 * nF, nF);
  nablaV << nablaV_flat.segment(0, nF), nablaV_flat.segment(nF, nF), nablaV_flat.segment(2 * nF, nF);

  // Compute face normal
  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, N);

  // Compute area
  Eigen::VectorXd dblA;
  igl::doublearea(V, F, dblA);

  // Sum up
  double sum = 0;

  for (int i = 0; i < nF; ++i) {
    Eigen::Vector3d N_i = N.row(i);
    Eigen::Vector3d nablaU_i = nablaU.row(i);
    Eigen::Vector3d nablaV_i = nablaV.row(i);

    // Rotate nablaU_i about N_i by 90 degrees ccw
    nablaU_i = N_i.cross(nablaU_i);

    sum += 0.5 * dblA[i] * (nablaU_i - nablaV_i).squaredNorm();
  }

  return 0.5 * sum;
}

TEST_CASE("lscm: lscm_energy_check", "[igl]")
{
  // Load a mesh in OFF format
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("camelhead.off"), V, F);

  // Fix two points on the boundary
  Eigen::VectorXi bnd, b(2, 1);
  igl::boundary_loop(F, bnd);
  b(0) = bnd(0);
  b(1) = bnd(bnd.size() / 2);
  Eigen::MatrixXd bc(2, 2);
  bc << 0, 0, 1, 0;

  // LSCM parametrization
  Eigen::MatrixXd V_uv;
  Eigen::SparseMatrix<double> Q;
  igl::lscm(V, F, b, bc, V_uv, Q);

  Eigen::VectorXd W_flat(2 * V.rows());
  W_flat << V_uv.col(0), V_uv.col(1);

  // Consistency check
  double e1 = compute_lscm_energy(V, F, W_flat);
  double e2 = 0.5 * W_flat.transpose() * Q * W_flat;

  REQUIRE(std::abs(e1 - e2) < igl::EPS<double>());
}
