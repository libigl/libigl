#include <igl/readOFF.h>
#include <igl/combine.h>

#include <igl/opengl/glfw/Viewer.h>
//#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <igl/fast_find_self_intersections.h>

Eigen::MatrixXd V1,V2,V;
Eigen::MatrixXi F1,F2,F;

int main(int argc, char *argv[])
{

  // Load two meshes 
  igl::readOFF(TUTORIAL_SHARED_PATH "/planexy.off", V1, F1);
  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off",     V2, F2);

  // Combine into one mesh (will produce self-intersections)
  igl::combine<Eigen::MatrixXd,Eigen::MatrixXi>({V1,V2},{F1,F2}, V,F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  Eigen::VectorXi I;
  Eigen::MatrixXd edges;

  if(igl::fast_find_self_intersections(V,F,I,edges))
  {
    std::cout<<"Found "<<I.sum()<<" self intersections"<<std::endl;

    // plot edge vertices
    //viewer.data().add_points(edges, Eigen::RowVector3d(1,0,0));

    // Plot the edges of the self intersects
    for (unsigned i=0;i<edges.rows(); i+=2)
    {
      viewer.data().add_edges
      (
        edges.row(i),
        edges.row(i+1),
        Eigen::RowVector3d(1,0,0)
      );
    }
    std::cout<<std::endl;
  }
  viewer.data().set_data(I.cast<double>());
  viewer.data().double_sided=true;

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  //plugin.widgets.push_back(&menu);
  menu.callback_draw_viewer_window = [](){};

  // Launch the viewer
  viewer.launch();
}
