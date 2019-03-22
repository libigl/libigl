#include <igl/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOBJ.h>
#include <igl/Timer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>

#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
igl::Timer timer;
igl::SCAFData scaf_data;

bool show_uv = false;
float uv_scale = 0.2;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
    show_uv = false;
  else if (key == '2')
    show_uv = true;

  if (key == ' ')
  {
    timer.start();
    igl::scaf_solve(scaf_data, 1);
    std::cout << "time = " << timer.getElapsedTime() << std::endl;
  }

  const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
  if (show_uv)
  {
    viewer.data().clear();
    viewer.data().set_mesh(V_uv,F);
    viewer.data().set_uv(V_uv);
    viewer.core().align_camera_center(V_uv,F);
  }
  else
  {
    viewer.data().set_mesh(V,F);
    viewer.data().set_uv(V_uv);
    viewer.core().align_camera_center(V,F);
  }

  viewer.data().compute_normals();

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  // Load a mesh in OFF format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/camel_b.obj", V, F);

  Eigen::MatrixXd bnd_uv, uv_init;

  Eigen::VectorXd M;
  igl::doublearea(V, F, M);
  std::vector<std::vector<int>> all_bnds;
  igl::boundary_loop(F, all_bnds);

  // Heuristic primary boundary choice: longest
  auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int> &a, const std::vector<int> &b) { return a.size()<b.size(); });

  Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());

  igl::map_vertices_to_circle(V, bnd, bnd_uv);
  bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
  if (all_bnds.size() == 1)
  {
    if (bnd.rows() == V.rows()) // case: all vertex on boundary
    {
      uv_init.resize(V.rows(), 2);
      for (int i = 0; i < bnd.rows(); i++)
        uv_init.row(bnd(i)) = bnd_uv.row(i);
    }
    else
    {
      igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
      if (igl::flipped_triangles(uv_init, F).size() != 0)
        igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
    }
  }
  else
  {
    // if there is a hole, fill it and erase additional vertices.
    all_bnds.erase(primary_bnd);
    Eigen::MatrixXi F_filled;
    igl::topological_hole_fill(F, bnd, all_bnds, F_filled);
    igl::harmonic(F_filled, bnd, bnd_uv ,1, uv_init);
    uv_init = uv_init.topRows(V.rows());
  }

  Eigen::VectorXi b; Eigen::MatrixXd bc;
  igl::scaf_precompute(V, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
  viewer.data().set_uv(V_uv);
  viewer.callback_key_down = &key_down;

  // Enable wireframe
  viewer.data().show_lines = true;

  // Draw checkerboard texture
  viewer.data().show_texture = true;


  std::cerr << "Press space for running an iteration." << std::endl;
  std::cerr << "Press 1 for Mesh 2 for UV" << std::endl;

  // Launch the viewer
  viewer.launch();
}
