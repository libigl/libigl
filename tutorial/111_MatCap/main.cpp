#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/stb/read_image.h>

int main(int argc, char *argv[])
{
  igl::opengl::glfw::Viewer v;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
  igl::stb::read_image(argc>2?argv[2]: TUTORIAL_SHARED_PATH "/jade.png",R,G,B,A);
  v.data().set_mesh(V,F);
  v.data().set_texture(R,G,B,A);
  v.data().use_matcap = true;
  v.data().show_lines = false;
  v.launch();
}
