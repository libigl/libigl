#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <igl/readDMAT.h>
#include <igl/cubic.h>
#include <igl/cycodebase/box_cubic.h>

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  Eigen::MatrixXd P;
  Eigen::MatrixXi C;
  igl::readDMAT(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/libigl-cursive-P.dmat",P);
  igl::readDMAT(argc>1?argv[2]: TUTORIAL_SHARED_PATH "/libigl-cursive-C.dmat",C);

  igl::opengl::glfw::Viewer vr;
  // Cheezy conversion of (P,C) to polyline using 50 samples per curve
  const int ns = 50;
  Eigen::MatrixXd V(ns*C.rows(),P.cols());
  Eigen::MatrixXi E((ns-1)*C.rows(),2);
  for(int c = 0;c<C.rows();c++)
  {
    for(int i = 0;i<ns;i++)
    {
      const double t = double(i)/(ns-1);
      // Evaluate cubic Bezier at parameter t
      Eigen::RowVectorXd Vct;
      igl::cubic(P(C.row(c),Eigen::all),t,Vct);
      V.row(c*ns + i) = Vct;
      if(i>0)
      {
        E.row(c*(ns-1) + (i-1))<< c*ns + (i-1), c*ns + i;
      }
    }
  }

  Eigen::MatrixXd B1,B2;
  igl::cycodebase::box_cubic(P,C,B1,B2);
  Eigen::MatrixXd BV(C.rows()*4,2);
  BV<<
    B1.col(0), B1.col(1),
    B1.col(0), B2.col(1),
    B2.col(0), B2.col(1),
    B2.col(0), B1.col(1);
  Eigen::MatrixXi BE(C.rows()*4,2);
  for(int c = 0;c<C.rows();c++)
  {
    BE.row(c*4 + 0)<< 0*C.rows() + c, 1*C.rows() + c;
    BE.row(c*4 + 1)<< 1*C.rows() + c, 2*C.rows() + c;
    BE.row(c*4 + 2)<< 2*C.rows() + c, 3*C.rows() + c;
    BE.row(c*4 + 3)<< 3*C.rows() + c, 0*C.rows() + c;
  }



  vr.data().set_mesh(V,E(Eigen::all,{0,1,1}).eval());
  vr.data().set_points(P,Eigen::RowVector3d(1,0.7,0.2));
  Eigen::MatrixXi CE(C.rows()*2,2);
  CE<< C(Eigen::all,{0,1}),
    C(Eigen::all,{2,3});
  vr.data().set_edges(P,CE,Eigen::RowVector3d(1,0.7,0.2));
  vr.data().point_size = 10;
  //vr.data().set_points(BV,Eigen::RowVector3d(0,1,0));
  //vr.data().set_edges(BV,BE,Eigen::RowVector3d(0,1,1));
  vr.append_mesh();
  vr.data().set_mesh(BV,BE(Eigen::all,{0,1,1}).eval());
  vr.data().line_color = Eigen::RowVector4f(1,1,1,1);


  vr.launch();
  return 0;
}
