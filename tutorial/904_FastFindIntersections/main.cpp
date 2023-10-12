#include <igl/read_triangle_mesh.h>
#include <igl/AABB.h>
#include <igl/unique.h>

#include <igl/opengl/glfw/Viewer.h>

#include <igl/fast_find_intersections.h>

Eigen::MatrixXd V1,V2;
Eigen::MatrixXi F1,F2;

igl::AABB<Eigen::MatrixXd,3> tree;

double min_z,max_z;
double slice_z;


void update_visualization(igl::opengl::glfw::Viewer & viewer)
{
  //shifted intersection object
  Eigen::MatrixXd V2_(V2.rows(),V2.cols());
  V2_<< V2.col(0), V2.col(1), V2.col(2).array()+slice_z;

  Eigen::MatrixXi IF,EE;
  Eigen::MatrixXd EV;
  Eigen::VectorXi EI;
  igl::fast_find_intersections(tree, V1,F1, V2_,F2,false,false,IF,EV,EE,EI);
 
  // Plot the edges of the intersects
  viewer.data().set_edges( EV,EE, Eigen::RowVector3d(1,0,0));
  
  // show faces which are intersected
  Eigen::VectorXi I;
  igl::unique(IF,I);
  Eigen::VectorXd D = Eigen::MatrixXd::Zero(F1.rows(),1);
  D(I).setConstant(1.0);
  viewer.data().set_data(D,0,1,igl::COLOR_MAP_TYPE_PARULA);
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
  igl::read_triangle_mesh(
    argc<=1?TUTORIAL_SHARED_PATH "/cow.off"    :argv[1],V1,F1);
  igl::read_triangle_mesh(
    argc<=2?TUTORIAL_SHARED_PATH "/planexy.off":argv[2],V2,F2);

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
  viewer.data().set_face_based(true);
  

  update_visualization(viewer);

  viewer.callback_key_down = &key_down;

  // show the cow closer
  viewer.core().camera_zoom = 2.0;

  // Launch the viewer
  viewer.launch();
}
