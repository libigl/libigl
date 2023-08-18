#include "param_2d_demo_iter.h"
#include "check_mesh_for_issues.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/flipped_triangles.h>
#include <igl/MappingEnergyType.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/slim.h>
#include <igl/Timer.h>

void param_2d_demo_iter(
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    double & uv_scale_param,
    bool & first_iter,
    igl::SLIMData& sData,
    igl::Timer & timer,
    igl::opengl::glfw::Viewer& viewer)
{
  using namespace std;
  if (first_iter) {
    timer.start();
    igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/face.obj", V, F);
    check_mesh_for_issues(V,F);
    cout << "\tMesh is valid!" << endl;

    Eigen::MatrixXd uv_init;
    Eigen::VectorXi bnd; Eigen::MatrixXd bnd_uv;
    igl::boundary_loop(F,bnd);
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    igl::harmonic(V,F,bnd,bnd_uv,1,uv_init);
    if (igl::flipped_triangles(uv_init,F).size() != 0) {
      igl::harmonic(F,bnd,bnd_uv,1,uv_init); // use uniform laplacian
    }

    cout << "initialized parametrization" << endl;

    sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
    Eigen::VectorXi b; Eigen::MatrixXd bc;
    slim_precompute(V,F,uv_init,sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b,bc,0);

    uv_scale_param = 15 * (1./sqrt(sData.mesh_area));
    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V,F);
    viewer.data().set_uv(sData.V_o*uv_scale_param);
    viewer.data().compute_normals();
    viewer.data().show_texture = true;

    first_iter = false;
  } else {
    timer.start();
    slim_solve(sData,1); // 1 iter
    viewer.data().set_uv(sData.V_o*uv_scale_param);
  }
  cout << "time = " << timer.getElapsedTime() << endl;
  cout << "energy = " << sData.energy << endl;
}

