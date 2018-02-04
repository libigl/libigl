#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "tutorial_shared_path.h"
#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

    // Draw the scene in the buffers
    viewer.core.draw_buffer(
      viewer.data(),false,R,G,B,A);

    // Save it to a PNG
    igl::png::writePNG(R,G,B,A,"out.png");
  }

  if (key == '2')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;

    // Read the PNG
    igl::png::readPNG("out.png",R,G,B,A);

    // Replace the mesh with a triangulated square
    Eigen::MatrixXd V(4,3);
    V <<
      -0.5,-0.5,0,
       0.5,-0.5,0,
       0.5, 0.5,0,
      -0.5, 0.5,0;
    Eigen::MatrixXi F(2,3);
    F <<
      0,1,2,
      2,3,0;
    Eigen::MatrixXd UV(4,2);
    UV <<
      0,0,
      1,0,
      1,1,
      0,1;

    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_uv(UV);
    viewer.core.align_camera_center(V);
    viewer.data().show_texture = true;

    // Use the image as a texture
    viewer.data().set_texture(R,G,B);

  }


  return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  std::cerr << "Press 1 to render the scene and save it in a png." << std::endl;
  std::cerr << "Press 2 to load the saved png and use it as a texture." << std::endl;

  // Plot the mesh and register the callback
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}
