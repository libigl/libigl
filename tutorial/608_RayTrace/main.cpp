#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/read_triangle_mesh.h>
#include <igl/png/writePNG.h>
#include <igl/PI.h>
#include <Eigen/Geometry>
// embree
#include <igl/embree/EmbreeRenderer.h>

#include "tutorial_shared_path.h"

#include <iostream>

int main(int argc, char *argv[])
{
  int width=640;
  int height=480;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  const char *mesh_file=argc > 1 ? argv[1] : TUTORIAL_SHARED_PATH "/fertility.off";
  const char *png_file=argc > 2 ? argv[2]: "fertility_curvature.png";

  // Load mesh 
  igl::read_triangle_mesh(mesh_file, V,F);
  
  Eigen::VectorXd K;
  // Compute integral of Gaussian curvature
  igl::gaussian_curvature(V,F,K);
  // Compute mass matrix
  Eigen::SparseMatrix<double> M,Minv;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  igl::invert_diag(M,Minv);
  // Divide by area to get integral average
  K = (Minv*K).eval();

  // embree object
  igl::embree::EmbreeRenderer er;
  er.set_mesh(V,F,true);

  //er.set_uniform_color(Eigen::RowVector3d(0.3,0.3,0.3));
  er.set_data(K,igl::COLOR_MAP_TYPE_JET);

  Eigen::Matrix3d rot_matrix;

  // specify rotation
  rot_matrix =  Eigen::AngleAxisd( 10*igl::PI/180.0, Eigen::Vector3d::UnitX())
              * Eigen::AngleAxisd(  5*igl::PI/180.0, Eigen::Vector3d::UnitY())
              * Eigen::AngleAxisd(  4*igl::PI/180.0, Eigen::Vector3d::UnitZ());
  er.set_rot(rot_matrix);

  er.set_zoom(1.5);
  er.set_orthographic(false);
     
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(width, height);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(width, height);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(width, height);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(width, height);

  // render view using embree
  er.render_buffer(R,G,B,A);

  std::cout<<"Rendered scene saved to "<<png_file<<std::endl;

  // save to PNG file
  igl::png::writePNG(R,G,B,A,png_file);
  return 0;
}
