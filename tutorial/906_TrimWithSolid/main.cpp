#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/sort_triangles.h>
#include <igl/material_colors.h>
#include <igl/copyleft/cgal/trim_with_solid.h>

int main(int argc, char *argv[])
{
  using namespace igl;
  IGL_TICTOC_LAMBDA;
  Eigen::MatrixXd VA, VB;
  Eigen::MatrixXi FA, FB;
  // Load a hot mess of a mesh
  igl::read_triangle_mesh(argc>1?argv[1]:TUTORIAL_SHARED_PATH "/truck.obj", VA, FA);
  // Load a solid mesh
  igl::read_triangle_mesh(argc>2?argv[2]:TUTORIAL_SHARED_PATH "/bunny.off", VB, FB);
  if(argc<2)
  {
    // resize bunny
    VB.rowwise() -= VB.colwise().mean();
    VB /= (VB.colwise().maxCoeff()-VB.colwise().minCoeff()).maxCoeff();
    VB *= 1.15;
    VB.rowwise() += VA.colwise().mean();
  }

  Eigen::VectorXi J;
  Eigen::MatrixXd VC;
  Eigen::MatrixXi FC;
  Eigen::Array<bool, Eigen::Dynamic, 1> D;
  using namespace igl::copyleft::cgal;
  tictoc();
  // More patches, less intersection handling
  igl::copyleft::cgal::trim_with_solid(VA, FA, VB, FB, CHECK_EACH_PATCH, VC, FC, D, J);
  printf("CHECK_EACH_PATCH: %g secs, |FC| = %d, |D| = %d\n", tictoc(),FC.rows(),D.count());
  // More intersection handling, fewer patches
  igl::copyleft::cgal::trim_with_solid(VA, FA, VB, FB, RESOLVE_BOTH_AND_RESTORE_THEN_CHECK_EACH_PATCH, VC, FC, D, J);
  printf("RESOLVE_BOTH_...: %g secs, |FC| = %d, |D| = %d\n", tictoc(),FC.rows(),D.count());

  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(VC, FC);
  // Turn on double sided lighting
  vr.data().double_sided = true;
  vr.data().set_face_based(true);
  vr.data().set_data(D.cast<double>());
  Eigen::MatrixXd CM = (Eigen::MatrixXd(2,3)<< 
      1,1,1,
      GOLD_DIFFUSE[0],GOLD_DIFFUSE[1],GOLD_DIFFUSE[2]
      ).finished();
  vr.data().set_colormap(CM);

  vr.append_mesh();
  vr.data().set_mesh(VB, FB);
  vr.data().show_lines = false;
  // Make a semi-transparent orange matcap
  {
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> 
      R(256,256), G(256,256), B(256,256), A(256,256);
    for(int i = 0;i<R.rows();i++)
    {
      for(int j = 0;j<R.cols();j++)
      {
        R(i,j) = 255;
        G(i,j) = 110;
        B(i,j) = 20;
        // distance to middle of image
        const double x = (i+0.5-R.rows()/2.0)/(R.rows()/2.0);
        const double y = (j+0.5-R.cols()/2.0)/(R.cols()/2.0);
        const double r = std::min(1.0,sqrt(x*x+y*y));
        A(i,j) = 255*(1-sqrt(1-r*r));
      }
    }
    vr.data().set_texture(R,G,B,A);
  }
  vr.data().use_matcap = true;
  // On mouse up resort for better transparency
  vr.callback_mouse_up = 
    [&](igl::opengl::glfw::Viewer &, int button, int mod)
  {
    Eigen::VectorXi _;
    igl::sort_triangles(
      VB,Eigen::MatrixXi(FB), vr.core().view, vr.core().proj,FB,_);
    vr.data_list[1].set_mesh(VB, FB);
    return false;
  };
  // set first mesh as selected
  vr.selected_data_index = 0;

  vr.launch_init();
  vr.core().draw(vr.data(),true);
  vr.callback_mouse_up(vr,0,0);
  vr.launch_rendering(true);
  vr.launch_shut();
}
