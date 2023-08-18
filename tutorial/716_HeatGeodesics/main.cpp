#include "isolines_colormap.h"
#include "update.h"
#include <igl/read_triangle_mesh.h>
#include <igl/heat_geodesics.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

int main(int argc, char *argv[])
{
  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  igl::read_triangle_mesh( argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off",V,F);
  double t = std::pow(igl::avg_edge_length(V,F),2);

  // Precomputation
  igl::HeatGeodesicsData<double> data;
  const auto precompute = [&]()
  {
    if(!igl::heat_geodesics_precompute(V,F,t,data))
    {
      std::cerr<<"Error: heat_geodesics_precompute failed."<<std::endl;
      exit(EXIT_FAILURE);
    };
  };
  precompute();

  igl::opengl::glfw::Viewer viewer;
  bool down_on_mesh = false;
  const auto update = [&]()->bool
  {
    const double x = viewer.current_mouse_x;
    const double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Eigen::VectorXd D;
    if(::update(
      V,F,t,x,y,
      viewer.core().view,viewer.core().proj,viewer.core().viewport,
      data,
      D))
    {
      viewer.data().set_data(D);
      return true;
    }
    return false;
  };
  viewer.callback_mouse_down =
    [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    if(update())
    {
      down_on_mesh = true;
      return true;
    }
    return false;
  };
  viewer.callback_mouse_move =
    [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
    {
      if(down_on_mesh)
      {
        update();
        return true;
      }
      return false;
    };
  viewer.callback_mouse_up =
    [&down_on_mesh](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    down_on_mesh = false;
    return false;
  };
  std::cout<<R"(Usage:
  [click]  Click on shape to pick new geodesic distance source
  ,/.      Decrease/increase t by factor of 10.0
  D,d      Toggle using intrinsic Delaunay discrete differential operators

)";

  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer& /*viewer*/, unsigned int key, int mod)->bool
  {
    switch(key)
    {
    default:
      return false;
    case 'D':
    case 'd':
      data.use_intrinsic_delaunay = !data.use_intrinsic_delaunay;
      std::cout<<(data.use_intrinsic_delaunay?"":"not ")<<
        "using intrinsic delaunay..."<<std::endl;
      precompute();
      update();
      break;
    case '.':
    case ',':
      t *= (key=='.'?10.0:0.1);
      precompute();
      update();
      std::cout<<"t: "<<t<<std::endl;
      break;
    }
    return true;
  };

  // Show mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(Eigen::VectorXd::Zero(V.rows()));
  viewer.data().set_colormap(isolines_colormap());
  viewer.data().show_lines = false;
  viewer.launch();

}
