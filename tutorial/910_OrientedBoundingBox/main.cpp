#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/oriented_bounding_box.h>
#include <igl/copyleft/cgal/convex_hull.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/oriented_bounding_box.h>
#include <igl/moments.h>
#include <igl/super_fibonacci.h>
#include <igl/matlab_format.h>
#include <igl/write_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/parallel_for.h>
#include <Eigen/Core>
#include <iostream>
#include <limits>


int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;

  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  // Read in inputs as double precision floating point meshes
  igl::read_triangle_mesh(
      argc>1?argv[1]: TUTORIAL_SHARED_PATH "/hand.mesh",V,F);

  V.array() += 1;

  tictoc();
  {
    Eigen::RowVector3d min_corner = V.colwise().minCoeff();
    Eigen::RowVector3d max_corner = V.colwise().maxCoeff();
  };
  double t_aabb = tictoc();

  // PCA
  tictoc();
  Eigen::RowVector3d mean = V.colwise().mean();
  Eigen::MatrixXd V_centered = V.rowwise() - mean;
  Eigen::Matrix3d cov = (V_centered.transpose() * V_centered) / (V.rows() - 1);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
  Eigen::Matrix3d PR = eig.eigenvectors();
  double t_pca = tictoc();


  Eigen::Matrix3d R;
  tictoc();
  igl::copyleft::cgal::oriented_bounding_box(V, R);
  double t_cgal = tictoc();

  tictoc();
  Eigen::MatrixXd W;
  Eigen::MatrixXi G;
  igl::copyleft::cgal::convex_hull(V,W,G);
  Eigen::Matrix3d BR;
  igl::oriented_bounding_box(W, 50000, igl::ORIENTED_BOUNDING_BOX_MINIMIZE_SURFACE_AREA, BR);
  std::cout<<igl::matlab_format(BR, "saBR")<<std::endl;
  igl::oriented_bounding_box(W, 50000, igl::ORIENTED_BOUNDING_BOX_MINIMIZE_VOLUME, BR);
  std::cout<<igl::matlab_format(BR, "vBR")<<std::endl;
  double t_best = tictoc();

  // Cube mesh
  Eigen::MatrixXd CV(8,3);
  CV<< 
    0,0,0,
    0,1,0,
    1,0,0,
    1,1,0,
    1,0,1,
    1,1,1,
    0,0,1,
    0,1,1;
  Eigen::MatrixXi CF(12,3);
  CF<< 
    1,2,0,
    1,3,2,
    3,4,2,
    3,5,4,
    0,4,6,
    0,2,4,
    7,3,1,
    7,5,3,
    7,0,6,
    7,1,0,
    5,6,4,
    5,7,6;
  Eigen::MatrixXi CE(12,2);
  CE<< 
    0,1,
    0,2,
    0,6,
    1,3,
    1,7,
    2,3,
    2,4,
    3,5,
    4,5,
    4,6,
    5,7,
    6,7;


  const auto transform_cube = [](const Eigen::MatrixXd & V, const Eigen::MatrixXd & CV)
  {
    Eigen::RowVector3d min_corner = V.colwise().minCoeff();
    Eigen::RowVector3d max_corner = V.colwise().maxCoeff();
    Eigen::RowVector3d diff = max_corner - min_corner;
    Eigen::MatrixXd AV = CV;
    AV.array().rowwise() *= diff.array();
    AV.rowwise() += min_corner;
    return AV;
  };
  Eigen::MatrixXd AV = transform_cube(V, CV);

  Eigen::MatrixXd OV = transform_cube((V*R).eval(), CV);
  OV = (OV*R.transpose()).eval();

  Eigen::MatrixXd PV = transform_cube((V*PR).eval(), CV);
  PV = (PV*PR.transpose()).eval();

  Eigen::MatrixXd BV = transform_cube((V*BR).eval(), CV);
  BV = (BV*BR.transpose()).eval();

  const auto volume = [](const Eigen::MatrixXd & V, const Eigen::MatrixXi & F)
  {
    Eigen::Vector3d m1;
    Eigen::Matrix3d m2;
    double vol = 0;
    igl::moments(V,F,vol,m1,m2);
    return vol;
  };
  double Avol = volume(AV,CF);
  double Pvol = volume(PV,CF);
  double Ovol = volume(OV,CF);
  double Bvol = volume(BV,CF);
  printf("method: %20s %20s  %20s %20s\n", "AABB", "PCA", "CGAL", "super_fibonacci");
  printf("volume: %20g %20g  %20g %20g\n", Avol, Pvol, Ovol, Bvol);
  printf("time:   %20g %20g  %20g %20g\n", t_aabb, t_pca, t_cgal, t_best);


  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);

  viewer.append_mesh();
  viewer.data().set_mesh(AV,CF);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(Eigen::RowVector4d(0.43064,0.77072,0.51013,0.5));
  viewer.data().set_edges(AV,CE,Eigen::RowVector3d(0,0,0));
  viewer.data().show_lines = false;

  viewer.append_mesh();
  viewer.data().set_mesh(PV,CF);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(Eigen::RowVector4d(0.9374,0.54164,0.66868,0.5));
  viewer.data().set_edges(PV,CE,Eigen::RowVector3d(0,0,0));
  viewer.data().show_lines = false;

  viewer.append_mesh();
  viewer.data().set_mesh(OV,CF);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(Eigen::RowVector4d(0.91191,0.60119,0.32995,0.5));
  viewer.data().show_lines = false;

  viewer.append_mesh();
  viewer.data().set_mesh(BV,CF);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(Eigen::RowVector4d(0.34481,0.72126,0.96535,0.5));
  viewer.data().show_lines = false;
  viewer.data().set_edges(OV,CE,Eigen::RowVector3d(0,0,0));

  // set background to white
  viewer.core().background_color = Eigen::Vector4f(1,1,1,1);

  viewer.data().show_lines = false;
  viewer.launch();
}

