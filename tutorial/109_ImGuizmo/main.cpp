#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuizmoPlugin.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imguizmo/ImGuizmo.h>

int main(int argc, char *argv[])
{
  // Load a mesh from file
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/cow.off",V,F);
  // Set up viewer
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
  // Custom menu
  igl::opengl::glfw::imgui::ImGuizmoPlugin plugin;
  vr.plugins.push_back(&plugin);
  // Initialize ImGuizmo at mesh centroid
  plugin.T.block(0,3,3,1) = 
    0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff()).transpose().cast<float>();
  // Update can be applied relative to this remembered initial transform
  const Eigen::Matrix4f T0 = plugin.T;
  // Attach callback to apply imguizmo's transform to mesh
  plugin.callback = [&](const Eigen::Matrix4f & T)
  {
    const Eigen::Matrix4d TT = (T*T0.inverse()).cast<double>().transpose();
    vr.data().set_vertices(
      (V.rowwise().homogeneous()*TT).rowwise().hnormalized());
    vr.data().compute_normals();
  };
  // Maya-style keyboard shortcuts for operation
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ' ': plugin.visible = !plugin.visible; return true;
      case 'W': case 'w': plugin.operation = ImGuizmo::TRANSLATE; return true;
      case 'E': case 'e': plugin.operation = ImGuizmo::ROTATE;    return true;
      case 'R': case 'r': plugin.operation = ImGuizmo::SCALE;     return true;
    }
    return false;
  };
  std::cout<<R"(
W,w  Switch to translate operation
E,e  Switch to rotate operation
R,r  Switch to scale operation
)";
  vr.launch();
}
