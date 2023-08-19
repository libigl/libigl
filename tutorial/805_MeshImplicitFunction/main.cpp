#include "contours.h"
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char * argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXd VV;
  Eigen::MatrixXi FF;
  Eigen::MatrixXd NN;
  Eigen::MatrixXd mcV;
  Eigen::MatrixXi mcF;
  Eigen::MatrixXd mcN;
  contours(V,E,VV,FF,NN,mcV,mcF,mcN);

  igl::opengl::glfw::Viewer vr;
  bool show_edges = true;
  bool use_dc = true;
  const auto update = [&]()
  {
    const bool was_face_based = vr.data().face_based ;
    vr.data().clear();
    if(use_dc)
    {
      vr.data().set_mesh(VV,FF);
      vr.data().show_lines = false;
      vr.data().set_normals(NN);
      if(show_edges)
      {
        vr.data().clear_edges();
        vr.data().set_edges(V,E,Eigen::RowVector3d(0,0,0));
      }
    }else
    {
      vr.data().set_mesh(mcV,mcF);
      vr.data().set_normals(mcN);
      vr.data().show_lines = show_edges;
    }
    vr.data().face_based = was_face_based;
  };
  update();
  vr.data().face_based = true;
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ' ': use_dc=!use_dc; update();return true;
      case 'L': case 'l': show_edges=!show_edges; update();return true;
    }
    return false;
  };
  std::cout<<R"(
[space]  Toggle between dual contouring and marching cubes
)";
  vr.launch();
}

