#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/list_to_matrix.h>
#include <igl/matlab_format.h>
#include <igl/AABB.h>
#include <igl/screen_space_selection.h>

#include <igl/opengl/glfw/imgui/SelectionPlugin.h>

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);

  // Plot the mesh
  igl::opengl::glfw::Viewer vr;
  igl::opengl::glfw::imgui::SelectionPlugin plugin;
  Eigen::VectorXd W = Eigen::VectorXd::Zero(V.rows());
  Eigen::Array<double,Eigen::Dynamic,1> and_visible = 
    Eigen::Array<double,Eigen::Dynamic,1>::Zero(V.rows());
  const Eigen::MatrixXd CM = (Eigen::MatrixXd(2,3)<<
      0.3,0.3,0.5,                 
      255.0/255.0,228.0/255.0,58.0/255.0).finished();
  bool only_visible = false;
  const auto update = [&]()
  {
    const bool was_face_based = vr.data().face_based;
    Eigen::VectorXd S = W;
    if(only_visible) { S.array() *= and_visible; }
    vr.data().set_data(S,0,1,igl::COLOR_MAP_TYPE_PLASMA,2);
    vr.data().face_based = was_face_based;
    vr.data().set_colormap(CM);
  };
  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V,F);
  plugin.callback = [&]()
  {
    screen_space_selection(V,F,tree,vr.core().view,vr.core().proj,vr.core().viewport,plugin.L,W,and_visible);
    update();
  };
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ' ': only_visible = !only_visible; update(); return true;
      case 'D': case 'd': W.setZero(); update(); return true;
    }
    return false;
  };
  std::cout<<R"(
Usage:
  [space]  Toggle whether to take visibility into account
  D,d      Clear selection
)";
  vr.plugins.push_back(&plugin);
  vr.data().set_mesh(V,F);
  vr.data().set_face_based(true);
  vr.core().background_color.head(3) = CM.row(0).head(3).cast<float>();
  vr.data().line_color.head(3) = (CM.row(0).head(3)*0.5).cast<float>();
  vr.data().show_lines = F.rows() < 20000;
  update();
  vr.launch();
}
