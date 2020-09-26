#include <test_common.h>
#include <igl/cat.h>
#include <igl/triangle/scaf.h>
#include <igl/slice.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/flip_avoiding_line_search.h>

// Assert that building the SCAF equation system via scaf_system() and
// solving it manually yields the same result as calling scaf_solve() directly.
TEST_CASE("scaf_system: Test scaf_system() vs scaf_solve()", "[igl]")
{
  // Read cube mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  // Cut to disk
  F = F.topRows(F.rows() - 1).eval();

  // Compute Tutte's embedding
  Eigen::MatrixXd uv_init;
  {
      Eigen::MatrixXd bnd_uv;
      std::vector<int> boundary;
      igl::boundary_loop(F, boundary);
      Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(boundary.data(), boundary.size());
      igl::map_vertices_to_circle(V, bnd, bnd_uv);
      igl::harmonic(F, bnd, bnd_uv ,1, uv_init);
      uv_init.conservativeResize(V.rows(), 2);
  }

  // Init empty constraint vectors
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;

  // Run one scaf iteration as reference
  igl::triangle::SCAFData s_ref;
  {
      igl::triangle::scaf_precompute(V, F, uv_init, s_ref, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);
      igl::triangle::scaf_solve(s_ref, 1);
  }

  // Obtain SCAF linear system perform iteration manually
  igl::triangle::SCAFData s_test;
  {
      // Set up system
      igl::triangle::scaf_precompute(V, F, uv_init, s_test, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);
      Eigen::SparseMatrix<double> L;
      Eigen::VectorXd rhs;
      igl::triangle::scaf_system(s_test, L, rhs);

      // Solve
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
      auto x = solver.compute(L).solve(rhs);

      // Write result to uv_target
      Eigen::MatrixXd uv_target = s_test.w_uv;
      const int n_v_var = s_test.w_uv.rows() - s_test.frame_ids.size();
      for (int i = 0; i < n_v_var; ++i)
      {
          int row = i;
          if (i >= s_test.v_num)
              row += s_test.frame_ids.size();

          uv_target(row, 0) = x[i];
          uv_target(row, 1) = x[i + n_v_var];
      }

      // Line search
      auto E = [&s_test](Eigen::MatrixXd &uv) { return igl::triangle::scaf::compute_energy(s_test, uv, true); };
      Eigen::MatrixXi w_T;
      igl::cat(1, s_test.m_T, s_test.s_T, w_T);
      igl::flip_avoiding_line_search(w_T, s_test.w_uv, uv_target, E, -1);
  }

  // Compare both results
  REQUIRE(s_test.w_uv == s_ref.w_uv);
}
