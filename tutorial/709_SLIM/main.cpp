#include "check_mesh_for_issues.h"
#include "get_soft_constraint_for_circle.h"
#include "get_cube_corner_constraints.h"
#include "param_2d_demo_iter.h"


#include <igl/slim.h>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/Timer.h>

#include <igl/MappingEnergyType.h>
#include <igl/serialize.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/barycenter.h>

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

void soft_const_demo_iter(igl::opengl::glfw::Viewer& viewer);
void deform_3d_demo_iter(igl::opengl::glfw::Viewer& viewer);
void display_3d_mesh(igl::opengl::glfw::Viewer& viewer);

Eigen::MatrixXd V;
Eigen::MatrixXi F;
bool first_iter = true;
igl::SLIMData sData;
igl::Timer timer;

double uv_scale_param;

enum DEMO_TYPE {
  PARAM_2D,
  SOFT_CONST,
  DEFORM_3D
};
DEMO_TYPE demo_type;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
  if (key == ' ') {
    switch (demo_type) {
      case PARAM_2D: {
        param_2d_demo_iter(V,F,uv_scale_param,first_iter,sData,timer,viewer);
        break;
      }
      case SOFT_CONST: {
        soft_const_demo_iter(viewer);
        break;
      }
      case DEFORM_3D: {
        deform_3d_demo_iter(viewer);
        break;
      }
      default:
        break;
    }
  }

  return false;
}

void soft_const_demo_iter(igl::opengl::glfw::Viewer& viewer) {
  using namespace std;
  if (first_iter) {

    igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/circle.obj", V, F);

    check_mesh_for_issues(V,F);
    cout << "\tMesh is valid!" << endl;
    Eigen::MatrixXd V_0 = V.block(0,0,V.rows(),2);

    Eigen::VectorXi b; Eigen::MatrixXd bc;
    get_soft_constraint_for_circle(V_0,F,b,bc);
    double soft_const_p = 1e5;
    slim_precompute(V,F,V_0,sData,igl::MappingEnergyType::SYMMETRIC_DIRICHLET,b,bc,soft_const_p);

    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V,F);
    viewer.data().compute_normals();
    viewer.data().show_lines = true;

    first_iter = false;

  } else {
    slim_solve(sData,1); // 1 iter
    viewer.data().set_mesh(sData.V_o, F);
  }
}

void deform_3d_demo_iter(igl::opengl::glfw::Viewer& viewer) {
  using namespace std;
  if (first_iter) {
    timer.start();
    igl::readOBJ(TUTORIAL_SHARED_PATH "/cube_40k.obj", V, F);

    Eigen::MatrixXd V_0 = V;
    Eigen::VectorXi b; Eigen::MatrixXd bc;
    get_cube_corner_constraints(V_0,F,b,bc);

    double soft_const_p = 1e5;
    sData.exp_factor = 5.0;
    slim_precompute(V,F,V_0,sData,igl::MappingEnergyType::EXP_CONFORMAL,b,bc,soft_const_p);
    //cout << "precomputed" << endl;

    first_iter = false;
    display_3d_mesh(viewer);

  } else {
    timer.start();
    slim_solve(sData,1); // 1 iter
    display_3d_mesh(viewer);
  }
  cout << "time = " << timer.getElapsedTime() << endl;
  cout << "energy = " << sData.energy << endl;
}

void display_3d_mesh(igl::opengl::glfw::Viewer& viewer) {
  using namespace std;
  using namespace Eigen;
  MatrixXd V_temp; MatrixXi F_temp;
  Eigen::MatrixXd Barycenters;

  igl::barycenter(sData.V,sData.F,Barycenters);
  //cout << "Barycenters.rows() = " << Barycenters.rows() << endl;
  //double t = double((key - '1')+1) / 9.0;
  double view_depth = 10.;
  double t = view_depth/9.;

  VectorXd v = Barycenters.col(2).array() - Barycenters.col(2).minCoeff();
  v /= v.col(0).maxCoeff();

  vector<int> s;

  for (unsigned i=0; i<v.size();++i)
    if (v(i) < t)
      s.push_back(i);

  V_temp.resize(s.size()*4,3);
  F_temp.resize(s.size()*4,3);

  for (unsigned i=0; i<s.size();++i){
    V_temp.row(i*4+0) = sData.V_o.row(sData.F(s[i],0));
    V_temp.row(i*4+1) = sData.V_o.row(sData.F(s[i],1));
    V_temp.row(i*4+2) = sData.V_o.row(sData.F(s[i],2));
    V_temp.row(i*4+3) = sData.V_o.row(sData.F(s[i],3));
    F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
    F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
    F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
    F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
  }
  viewer.data().set_mesh(V_temp,F_temp);
  viewer.core().align_camera_center(V_temp,F_temp);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = true;
}

int main(int argc, char *argv[]) {
  using namespace std;
  using namespace Eigen;

  cerr << "Press space for running an iteration." << std::endl;
  cerr << "Syntax: " << argv[0] << " demo_number (1 to 3)" << std::endl;
  cerr << "1. 2D unconstrained parametrization" << std::endl;
  cerr << "2. 2D deformation with positional constraints" << std::endl;
  cerr << "3. 3D mesh deformation with positional constraints" << std::endl;

  demo_type = PARAM_2D;

   if (argc == 2) {
     switch (std::atoi(argv[1])) {
       case 1: {
         demo_type = PARAM_2D;
         break;
       } case 2: {
         demo_type = SOFT_CONST;
         break;
       } case 3: {
         demo_type = DEFORM_3D;
         break;
       }
       default: {
         cerr << "Wrong demo number - Please choose one between 1-3" << std:: endl;
         exit(1);
       }
     }
   }


  // Launch the viewer
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;

  // Disable wireframe
  viewer.data().show_lines = false;

  // Draw checkerboard texture
  viewer.data().show_texture = false;

  // First iteration
  key_down(viewer, ' ', 0);

  viewer.launch();

  return 0;
}



