#include <igl/readOFF.h>
#include <igl/combine.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <igl/fast_find_intersections.h>

Eigen::MatrixXd V1,V2;
Eigen::MatrixXi F1,F2;

igl::AABB<Eigen::MatrixXd,3> tree;

double min_z,max_z;
double slice_z;


void update_visualization(igl::opengl::glfw::Viewer & viewer)
{
  Eigen::MatrixXi I;
  Eigen::MatrixXd edges;

  //shifted intersection object
  Eigen::MatrixXd V2_(V2.rows(),V2.cols());
  V2_<< V2.col(0), V2.col(1), V2.col(2).array()+slice_z;

  igl::fast_find_intersections(tree, V1,F1, V2_,F2, I,edges);
  Eigen::MatrixXi edges_link=Eigen::MatrixXi::NullaryExpr(edges.rows()/2,2, [](int i,int j) { return i*2+j;});
 
  // Plot the edges of the intersects
  viewer.data().set_edges ( edges, edges_link, Eigen::RowVector3d(1,0,0));
  
  // show faces which are intersected
  Eigen::VectorXd face_data=Eigen::VectorXd::Zero(F1.rows());

  for(int i=0; i<I.rows(); ++i)
    face_data(I(i,0)) = 1.0;
  
  viewer.data().set_data(face_data);
}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)
{
  switch(key)
  {
    default:
      return false;
    case '.':
      slice_z = std::min(slice_z+0.01,max_z);
      break;
    case ',':
      slice_z = std::max(slice_z-0.01,min_z);
      break;
  }
  update_visualization(viewer);
  return true;
}


int main(int argc, char *argv[])
{
  // Load two meshes 
  igl::readOFF(TUTORIAL_SHARED_PATH "/cow.off",     V1, F1);
  igl::readOFF(TUTORIAL_SHARED_PATH "/planexy.off", V2, F2);

  // initialize AABB tree with first mesh (it doesn't move)
  tree.init(V1,F1);

  std::cout<<"Usage:"<<std::endl;
  std::cout<<"'.'/','  push back/pull forward slicing plane."<<std::endl;
  std::cout<<std::endl;

  min_z=V1.col(2).minCoeff();
  max_z=V1.col(2).maxCoeff();
  slice_z=(max_z+min_z)/2;

  //center slicing object
  V2.col(2).array() -= V2.col(2).array().mean();

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V1, F1);
  

  update_visualization(viewer);

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  //plugin.widgets.push_back(&menu);
  menu.callback_draw_viewer_window = [](){};
  viewer.callback_key_down = &key_down;

  // show the cow closer
  viewer.core().camera_zoom = 2.0;

  // Launch the viewer
  viewer.launch();
}
