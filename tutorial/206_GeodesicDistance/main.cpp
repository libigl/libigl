#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/exact_geodesic.h>
#include <igl/colormap.h>
#include <igl/unproject_onto_mesh.h>

#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

void plotMeshDistance(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& d, const double strip_size )
{
    // Rescale the function depending on the strip size
    Eigen::VectorXd f = (d/strip_size);

    // The function should be 1 on each integer coordinate
    f = (f*M_PI).array().sin().abs();

    // Compute per-vertex colors
    Eigen::MatrixXd C;
    igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,f,false,C);

    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OFF format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);


  igl::opengl::glfw::Viewer viewer;
  // Plot a distance when a vertex is picked
  viewer.callback_mouse_down =
  [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
      int fid;
      Eigen::Vector3f bc;
      // Cast a ray in the view direction starting from the mouse position
      double x = viewer.current_mouse_x;
      double y = viewer.core.viewport(3) - viewer.current_mouse_y;
      if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view,
                                  viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
      {
          int max;
          bc.maxCoeff(&max);
          int vid = F(fid,max);
          Eigen::VectorXi VS,FS,VT,FT;
          // The selected vertex is the source
          VS.resize(1);
          VS << vid;
          // All vertices are the targets
          VT.setLinSpaced(V.rows(),0,V.rows()-1);
          Eigen::VectorXd d;
          igl::exact_geodesic(V,F,VS,FS,VT,FT,d);

          plotMeshDistance(viewer,V,F,d,0.05);
      }
      return false;
  };
  viewer.data().set_mesh(V,F);

  cout << "Press [space] to smooth." << endl;;
  cout << "Press [r] to reset." << endl;;
  return viewer.launch();
}
