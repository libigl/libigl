#include <igl/read_triangle_mesh.h>
#include <igl/unique.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/fast_find_self_intersections.h>

int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc<=1?TUTORIAL_SHARED_PATH "/cow.off":argv[1], V, F);


  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  Eigen::VectorXi EI;
  Eigen::MatrixXd EV;
  Eigen::MatrixXi IF,EE;

  if(igl::fast_find_self_intersections(V,F,false,false,IF,EV,EE,EI))
  {
    std::cout<<"Found "<<IF.rows()<<" self intersecting pairs"<<std::endl;
    // plot edge vertices
    viewer.data().set_edges(EV,EE, Eigen::RowVector3d(1,0,0));
  }
  Eigen::VectorXi I;
  igl::unique(IF,I);
  Eigen::VectorXd D = Eigen::MatrixXd::Zero(F.rows(),1);
  D(I).setConstant(1.0);
  viewer.data().set_data(D,0,1,igl::COLOR_MAP_TYPE_PARULA);
  viewer.data().set_face_based(true);
  viewer.data().double_sided=true;

  // Launch the viewer
  viewer.launch();
}
