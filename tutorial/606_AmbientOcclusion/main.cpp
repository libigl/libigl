#include <igl/avg_edge_length.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOFF.h>
#include <igl/embree/ambient_occlusion.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

#include "tutorial_shared_path.h"

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::VectorXd AO;

// It allows to change the degree of the field when a number is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace Eigen;
  using namespace std;
  const RowVector3d color(0.9,0.85,0.9);
  switch(key)
  {
    case '1':
      // Show the mesh without the ambient occlusion factor
      viewer.data().set_colors(color);
      break;
    case '2':
    {
      // Show the mesh with the ambient occlusion factor
      MatrixXd C = color.replicate(V.rows(),1);
      for (unsigned i=0; i<C.rows();++i)
        C.row(i) *= AO(i);//std::min<double>(AO(i)+0.2,1);
      viewer.data().set_colors(C);
      break;
    }
    case '.':
      viewer.core.lighting_factor += 0.1;
      break;
    case ',':
      viewer.core.lighting_factor -= 0.1;
      break;
    default: break;
  }
  viewer.core.lighting_factor = 
    std::min(std::max(viewer.core.lighting_factor,0.f),1.f);

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  using namespace Eigen;
  cout<<
    "Press 1 to turn off Ambient Occlusion"<<endl<<
    "Press 2 to turn on Ambient Occlusion"<<endl<<
    "Press . to turn up lighting"<<endl<<
    "Press , to turn down lighting"<<endl;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V, F);

  MatrixXd N;
  igl::per_vertex_normals(V,F,N);

  // Compute ambient occlusion factor using embree
  igl::embree::ambient_occlusion(V,F,V,N,500,AO);
  AO = 1.0 - AO.array();

  // Show mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  key_down(viewer,'2',0);
  viewer.data().show_lines = false;
  viewer.core.lighting_factor = 0.0f;
  viewer.launch();
}
