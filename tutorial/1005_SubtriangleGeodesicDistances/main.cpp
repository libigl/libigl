#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/isolines_map.h>
#include <igl/parula.h>
#include <igl/upsample.h>
#include <igl/get_seconds.h>
#include <igl/exact_geodesic.h>
#include <igl/blue_noise.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/triangle/remesh_at_points.h>
#include <igl/opengl/glfw/Viewer.h>
#include <unordered_set>
#include <cassert>

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;

  Eigen::MatrixXd V; Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh(
      argc>1 ? argv[1] : TUTORIAL_SHARED_PATH "/decimated-knight.off", V, F))
    { std::cerr << "Failed to load mesh.\n"; return 1; }

      
  const int n_desired = argc>2?atoi(argv[2]):100;
  // Heuristic to  determine radius from desired number 
  const double r = [&V,&F](const int n)
  {
    Eigen::VectorXd A;
    igl::doublearea(V,F,A);
    return sqrt(((A.sum()*0.5/(n*0.6162910373))/igl::PI));
  }(n_desired);

  Eigen::MatrixXd B;
  Eigen::VectorXi I;
  Eigen::MatrixXd P;
  igl::blue_noise(V,F,r,B,I,P);
  const double tol = r*0.5;
  // Aggressively snap points to vertices and edges
  for(int i = 0;i<B.rows();i++)
  {
    for(int j = 0;j<B.cols();j++)
    {
      if(B(i,j) < tol) { B(i,j) = 0; }
    }
    const double sum = B.row(i).sum();
    B.row(i) /= sum;
  }
  // non-zeros in each row of B
  const auto nnz = (B.array()>0).rowwise().count();
  printf("counts: %d %d %d\n", (nnz.array()==1).count(), (nnz.array()==2).count(), (nnz.array()==3).count());

  Eigen::MatrixXd VV;
  Eigen::MatrixXi FF;
  Eigen::VectorXi J,K;
  igl::triangle::remesh_at_points(V,F,B,I,VV,FF,J,K);
  igl::write_triangle_mesh("remesh_at_points_output.ply",VV,FF);

  // Upsample for vis
  igl::upsample(VV,FF,4);

  Eigen::VectorXd D;
  igl::exact_geodesic(VV,FF,
      K,Eigen::VectorXi(),
      Eigen::VectorXi(Eigen::VectorXi::LinSpaced(VV.rows(),0,VV.rows()-1)),
      Eigen::VectorXi(),
      D);


  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(VV, FF);
  vr.data().set_face_based(false);
  vr.data().show_lines = false;
  vr.data().set_data(D);
  Eigen::MatrixXd CM;
  igl::parula(Eigen::VectorXd::LinSpaced(11,0,1).eval(),false,CM);
  igl::isolines_map(Eigen::MatrixXd(CM),CM);
  vr.data().set_colormap(CM);

  vr.data().point_size = 8;
  vr.data().add_points(P, Eigen::RowVector3d(1,0,0));
  vr.data().add_points(VV(K,Eigen::placeholders::all), Eigen::RowVector3d(0,1,0));
  vr.launch();
}

