#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/tetrahedralized_grid.h>
#include <igl/marching_tets.h>
#include <igl/signed_distance.h>
#include <igl/writeOBJ.h>
#include <igl/get_seconds.h>
#include <Eigen/Core>

#include "tutorial_shared_path.h"

#include <igl/marching_cubes.h>

int main(int argc, char * argv[])
{
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };

  // Load a surface mesh which is a cube
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);

  Eigen::RowVector3i side;
  Eigen::MatrixXd TV;
  tictoc();
  igl::voxel_grid(V,0,100,1,TV,side);
  printf("igl::voxel_grid               %g secs\n",tictoc());
  Eigen::MatrixXi TT5,TT6;
  tictoc();
  igl::tetrahedralized_grid(TV,side,igl::TETRAHEDRALIZED_GRID_TYPE_5,TT5);
  printf("igl::tetrahedralized_grid     %g secs\n",tictoc());
  igl::tetrahedralized_grid(TV,side,igl::TETRAHEDRALIZED_GRID_TYPE_6_ROTATIONAL,TT6);
  printf("igl::tetrahedralized_grid     %g secs\n",tictoc());

  tictoc();
  Eigen::VectorXd S;
  {
    Eigen::VectorXi I;
    Eigen::MatrixXd C,N;
    igl::signed_distance(
      TV,V,F,
      igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
      std::numeric_limits<double>::min(),
      std::numeric_limits<double>::max(),
      S,I,C,N);
  }
  printf("igl::signed_distance          %g secs\n",tictoc());
  
  std::vector<Eigen::MatrixXd> SV(3);
  std::vector<Eigen::MatrixXi> SF(3);
  igl::marching_tets(TV,TT5,S,0,SV[0],SF[0]);
  printf("igl::marching_tets5           %g secs\n",tictoc());
  igl::marching_tets(TV,TT6,S,0,SV[1],SF[1]);
  printf("igl::marching_tets6           %g secs\n",tictoc());
  tictoc();
  igl::marching_cubes(S,TV,side(0),side(1),side(2),0,SV[2],SF[2]);
  printf("igl::marching_cubes           %g secs\n",tictoc());


  //igl::writeOBJ("tmc.obj",SV,SF);

  // Draw the mesh stored in (SV, SF)
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
  vr.data().is_visible = false;
  vr.append_mesh();
  int sel = 0;
  const auto update = [&]()
  {
    vr.data().clear();
    vr.data().set_mesh(SV[sel],SF[sel]);
  };
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ',': 
      case '.': 
        sel = (sel+SV.size()+(key=='.'?1:-1))%SV.size();
        update();
        return true;
      case ' ': 
        vr.data_list[0].is_visible = !vr.data_list[0].is_visible;
        vr.data_list[1].is_visible = !vr.data_list[1].is_visible;
        vr.selected_data_index = (vr.selected_data_index+1)%2;
        return true;
    }
    return false;
  };
  update();
  vr.launch();
}
