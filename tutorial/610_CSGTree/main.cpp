#include "get_mesh.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <Eigen/Core>


int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  cout<<R"(
[,]  Toggle between boolean sub-tree operations
)";

  MatrixXi FA,FB,FC,FD,FE;
  MatrixXd VA,VB,VC,VD,VE;
  // Read in inputs as double precision floating point meshes
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/cube.obj"     ,VA,FA);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/sphere.obj"   ,VB,FB);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/xcylinder.obj",VC,FC);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/ycylinder.obj",VD,FD);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/zcylinder.obj",VE,FE);
  igl::opengl::glfw::Viewer viewer;

  int num_views = 5+4;
  int view_id = num_views-1;
  const auto & update = [&]()
  {
    viewer.data().clear();
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd I;
    get_mesh(VA,FA,VB,FB,VC,FC,VD,FD,VE,FE,view_id,V,F,I);
    viewer.data().set_mesh(V,F);
    MatrixXd C;
    jet(I,1,5,C);
    viewer.data().set_colors(C);
  };
  update();

  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)->bool
    {
      switch(key)
      {
        case ']':
          view_id = (view_id+1)%num_views;
          break;
        case '[':
          view_id = (view_id+num_views-1)%num_views;
          break;
        default:
          return false;
      }
      update();
      return true;
    };
  viewer.launch();
}
