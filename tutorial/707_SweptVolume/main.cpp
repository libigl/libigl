#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/material_colors.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/swept_volume.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/PI.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F,SF;
  Eigen::MatrixXd V,SV,VT;
  bool show_swept_volume = false;
  // Define a rigid motion
  const auto & transform = [](const double t)->Eigen::Affine3d
  {
    Eigen::Affine3d T = Eigen::Affine3d::Identity();
    T.rotate(Eigen::AngleAxisd(t*2.*igl::PI,Eigen::Vector3d(0,1,0)));
    T.translate(Eigen::Vector3d(0,0.125*cos(2.*igl::PI*t),0));
    return T;
  };
  // Read in inputs as double precision floating point meshes
  read_triangle_mesh(
      TUTORIAL_SHARED_PATH "/bunny.off",V,F);
  cout<<R"(Usage:
[space]  Toggle between transforming original mesh and swept volume
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  viewer.core.is_animating = !show_swept_volume;
  const int grid_size = 50;
  const int time_steps = 200;
  const double isolevel = 0.1;
  std::cerr<<"Computing swept volume...";
  igl::copyleft::swept_volume(
    V,F,transform,time_steps,grid_size,isolevel,SV,SF);
  std::cerr<<" finished."<<std::endl;

  viewer.callback_pre_draw =
    [&](igl::opengl::glfw::Viewer & viewer)->bool
    {
      if(!show_swept_volume)
      {
        Eigen::Affine3d T = transform(0.25*igl::get_seconds());
        VT = V*T.matrix().block(0,0,3,3).transpose();
        Eigen::RowVector3d trans = T.matrix().block(0,3,3,1).transpose();
        VT = ( VT.rowwise() + trans).eval();
        viewer.data().set_vertices(VT);
        viewer.data().compute_normals();
      }
      return false;
    };
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case ' ':
          show_swept_volume = !show_swept_volume;
          viewer.data().clear();
          if(show_swept_volume)
          {
            viewer.data().set_mesh(SV,SF);
            Eigen::Vector3d ambient = Eigen::Vector3d(SILVER_AMBIENT[0], SILVER_AMBIENT[1], SILVER_AMBIENT[2]);
            Eigen::Vector3d diffuse = Eigen::Vector3d(SILVER_DIFFUSE[0], SILVER_DIFFUSE[1], SILVER_DIFFUSE[2]);
            Eigen::Vector3d specular = Eigen::Vector3d(SILVER_SPECULAR[0], SILVER_SPECULAR[1], SILVER_SPECULAR[2]);
            viewer.data().uniform_colors(ambient,diffuse,specular);
          }
          else
          {
            viewer.data().set_mesh(V,F);
          }
          viewer.core.is_animating = !show_swept_volume;
          viewer.data().set_face_based(true);
          break;
      }
      return true;
    };
  viewer.launch();
}
