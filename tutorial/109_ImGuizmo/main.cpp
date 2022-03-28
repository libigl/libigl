#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <GLFW/glfw3.h>

int main(int argc, char *argv[])
{
  // Load a mesh from file
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/cow.off",V,F);
  // Set up viewer
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
 
  igl::opengl::glfw::imgui::ImGuiPlugin imgui_plugin;
  vr.plugins.push_back(&imgui_plugin);

  // Add a 3D gizmo plugin
  igl::opengl::glfw::imgui::ImGuizmoWidget gizmo;
  imgui_plugin.widgets.push_back(&gizmo);
  // Initialize ImGuizmo at mesh centroid
  gizmo.T.block(0,3,3,1) =
    0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff()).transpose().cast<float>();
  // Update can be applied relative to this remembered initial transform
  const Eigen::Matrix4f T0 = gizmo.T;
  // Attach callback to apply imguizmo's transform to mesh
  gizmo.callback = [&](const Eigen::Matrix4f & T)
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
      case ' ': gizmo.visible = !gizmo.visible; return true;
      case 'W': case 'w': gizmo.operation = ImGuizmo::TRANSLATE; return true;
      case 'E': case 'e': gizmo.operation = ImGuizmo::ROTATE;    return true;
      case 'R': case 'r': gizmo.operation = ImGuizmo::SCALE;     return true;
    }
    return false;
  };

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  imgui_plugin.widgets.push_back(&menu);

  std::cout<<R"(
W,w  Switch to translate operation
E,e  Switch to rotate operation
R,r  Switch to scale operation
)";
  vr.launch();
}
