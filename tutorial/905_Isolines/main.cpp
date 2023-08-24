#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/material_colors.h>
#include <igl/isolines.h>

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
      argc>1?argv[1]:TUTORIAL_SHARED_PATH "/fertility.off",V,F);
  // Use y-coordinate as scalar field
  Eigen::VectorXd f = V.col(1);

  Eigen::MatrixXd iV;
  Eigen::MatrixXi iE;
  {
    // How many isolines in the range (min,max)?
    const int n = argc>2?atoi(argv[2]):128;
    // This is actually unnecessary since isolines will not output degenerate
    // segments.
    //Eigen::VectorXd vals = Eigen::VectorXd::LinSpaced(n+2,f.minCoeff(),f.maxCoeff()).segment(1,n);
    // Instead just use all n+2 and if the min-,max-value isolines are
    // non-degneerate then we'll compute them, too.
    Eigen::VectorXd vals = Eigen::VectorXd::LinSpaced(n+2,f.minCoeff(),f.maxCoeff());
    {
      Eigen::VectorXi I;
      igl::isolines(V,F,f,vals,iV,iE,I);
    }
  }

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().label_size = 10;
  viewer.data().set_face_based(true);
  viewer.data().show_faces = true;
  viewer.data().show_lines = false;
  viewer.data().uniform_colors(
    Eigen::Vector3d(0.94*viewer.core().background_color.head<3>().cast<double>()),
    Eigen::Vector3d(0.05*viewer.core().background_color.head<3>().cast<double>()),
    Eigen::Vector3d(0.01*viewer.core().background_color.head<3>().cast<double>()));

  viewer.core().lighting_factor = 0.5;
  viewer.data().set_edges(iV, iE, 
      Eigen::RowVector3d(igl::GOLD_DIFFUSE[0], igl::GOLD_DIFFUSE[1], igl::GOLD_DIFFUSE[2]));
  viewer.data().line_width = 1;
  
  viewer.launch();
}

