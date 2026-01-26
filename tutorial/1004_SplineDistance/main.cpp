#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <igl/find.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/triangulated_grid.h>
#include <igl/cubic.h>
#include <igl/cycodebase/box_cubic.h>
#include <igl/cycodebase/point_spline_squared_distance.h>
#include <igl/cycodebase/spline_eytzinger_aabb.h>
#include <igl/predicates/spline_winding_number.h>

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  Eigen::MatrixXd P;
  Eigen::MatrixXi C;
  igl::readDMAT(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/libigl-cursive-P.dmat",P);
  igl::readDMAT(argc>1?argv[2]: TUTORIAL_SHARED_PATH "/libigl-cursive-C.dmat",C);

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


  Eigen::MatrixXd GV;
  Eigen::MatrixXi GF;
  {
    Eigen::RowVector2d min_corner = P.colwise().minCoeff();
    Eigen::RowVector2d max_corner = P.colwise().maxCoeff();
    // samples on x-axis
    const int sx = 512;
    const int sy = sx * (max_corner(1)-min_corner(1)) / (max_corner(0)-min_corner(0));
    Eigen::RowVector2d diag = (max_corner - min_corner).eval();
    const double bbd = diag.norm();
    diag /= diag.norm();
    min_corner -= 0.1*bbd*diag;
    max_corner += 0.1*bbd*diag;

    igl::triangulated_grid(sx,sy,GV,GF);
    // Scale and translate grid to fit bounding box
    GV.col(0) = (GV.col(0).array() * (max_corner(0)-min_corner(0))) + min_corner(0);
    GV.col(1) = (GV.col(1).array() * (max_corner(1)-min_corner(1))) + min_corner(1);
  }


  Eigen::VectorXd sqrD,S;
  Eigen::VectorXi I;
  Eigen::MatrixXd K;
  Eigen::VectorXd W;

  Eigen::Matrix<double,Eigen::Dynamic,2,Eigen::RowMajor> B1,B2;
  Eigen::VectorXi leaf;
  tictoc();
  igl::cycodebase::spline_eytzinger_aabb(P, C, B1, B2,leaf);
  printf("AABB build time: %g secs\n",tictoc());
  igl::writeDMAT("P.dmat", Eigen::MatrixXd(P));
  igl::writeDMAT("C.dmat", Eigen::MatrixXi(C));
  igl::writeDMAT("B1.dmat",Eigen::MatrixXd(B1));
  igl::writeDMAT("B2.dmat",Eigen::MatrixXd(B2));

  {
    tictoc();
    igl::cycodebase::point_spline_squared_distance(
        GV, P, C, B1,B2,leaf, sqrD, I, S, K);
    printf("sqrD: %d points in %g secs\n",GV.rows(),tictoc());
    tictoc();
    igl::predicates::spline_winding_number(P,C,B1,B2,leaf,GV,W);
    printf("Wind: %d points in %g secs\n",GV.rows(),tictoc());
  }
  Eigen::VectorXd unsigned_D = sqrD.array().sqrt();
  Eigen::VectorXd signed_D = unsigned_D.array() * (0.5 - W.array().abs()) * 2.0;

  // Just leaves
  //Eigen::MatrixXd B1,B2;
  //igl::cycodebase::box_cubic(P,C,B1,B2);

  const auto nonempty = igl::find((leaf.array() != -2).eval());
  const int nb = nonempty.size();
  Eigen::MatrixXd BV(nb*4,2);
  BV<<
    B1(nonempty,0), B1(nonempty,1),
    B1(nonempty,0), B2(nonempty,1),
    B2(nonempty,0), B2(nonempty,1),
    B2(nonempty,0), B1(nonempty,1);
  Eigen::MatrixXi BE(nb*4,2);
  for(int c = 0;c<nb;c++)
  {
    BE.row(c*4 + 0)<< 0*nb + c, 1*nb + c;
    BE.row(c*4 + 1)<< 1*nb + c, 2*nb + c;
    BE.row(c*4 + 2)<< 2*nb + c, 3*nb + c;
    BE.row(c*4 + 3)<< 3*nb + c, 0*nb + c;
  }

  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(GV,GF);
  vr.data().show_lines = false;
  const int g_index = vr.selected_data_index;
  vr.append_mesh();
  vr.data().set_mesh(V,E(Eigen::all,{0,1,1}).eval());
  vr.data().set_points(P,Eigen::RowVector3d(1,0.7,0.2));
  Eigen::MatrixXi CE(C.rows()*2,2);
  CE<< C(Eigen::all,{0,1}),
    C(Eigen::all,{2,3});
  vr.data().set_edges(P,CE,Eigen::RowVector3d(1,0.7,0.2));
  vr.data().point_size = 10;
  vr.append_mesh();
  vr.data().set_edges(BV,BE,Eigen::RowVector3d(0,0,0));

  enum Quantity
  {
    SIGNED_DISTANCE,
    UNSIGNED_DISTANCE,
    WINDING_NUMBER
  } quantity = SIGNED_DISTANCE;

  const auto update = [&]()
  {
    vr.core().lighting_factor = 0;
    const double dmax = unsigned_D.maxCoeff();
    const double dmin = -dmax;
    switch(quantity)
    {
      case SIGNED_DISTANCE:
        vr.data_list[g_index].set_data(signed_D,dmin,dmax,igl::COLOR_MAP_TYPE_ZOE,28);
        break;
      case UNSIGNED_DISTANCE:
        vr.data_list[g_index].set_data(unsigned_D,dmin,dmax,igl::COLOR_MAP_TYPE_ZOE,28);
        break;
      case WINDING_NUMBER:
        vr.data_list[g_index].set_data(W);
        break;
    }
  };
  update();

  vr.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)->bool
  {
    switch(key)
    {
      default:
        return false;
      case 'W':
      case 'w':
        quantity = WINDING_NUMBER;
        break;
      case 'U':
      case 'u':
        quantity = UNSIGNED_DISTANCE;
        break;
      case 'S':
      case 's':
        quantity = SIGNED_DISTANCE;
        break;
    }
    update();
    return true;
  };
  std::cout<<R"(
  S,s     Signed distance
  W,w     Winding number
  U,u     Unsigned distance
  N,n     Toggle naive / accelerated
  )"<<std::endl;
  vr.launch();
  return 0;
}
