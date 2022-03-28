#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>



int main(int argc, char *argv[])
{
  Eigen::MatrixXd V, V_uv;
  Eigen::MatrixXi F;
  // Load a mesh in OFF format
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/camelhead.off", V, F);

  // Find the open boundary
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);

  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  // Harmonic parametrization for the internal vertices
  igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

  // Scale UV to make the texture more clear
  V_uv *= 5;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_uv(V_uv);
  // Attach callback to allow toggling between 3D and 2D view
  viewer.callback_key_pressed = 
    [&V,&V_uv,&F](igl::opengl::glfw::Viewer& viewer, unsigned int key, int /*mod*/)
  {
    if(key == '3' || key == '2')
    {
      // Plot the 3D mesh or 2D UV coordinates
      viewer.data().set_vertices(key=='3'?V:V_uv);
      viewer.data().compute_normals();
      viewer.core().align_camera_center(key=='3'?V:V_uv,F);
      // key press was used
      return true;
    }
    // key press not used
    return false;
  };

  // Disable wireframe
  viewer.data().show_lines = false;

  // Draw default checkerboard texture
  viewer.data().show_texture = true;

  std::cout<<R"(
3  Show 3D model
2  Show 2D parametrization
)";

  // Launch the viewer
  viewer.launch();
}
