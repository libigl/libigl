#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <igl/readDMAT.h>
#include <igl/triangulated_grid.h>
#include <igl/box_simplices.h>
#include <igl/eytzinger_aabb.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/eytzinger_aabb_winding_number_tree.h>
#include <igl/eytzinger_aabb_winding_number.h>
double sdCapsule(
    const Eigen::RowVector2d & p, 
    const Eigen::RowVector2d & a, 
    const Eigen::RowVector2d & b, 
    const double & r )
{
  using Scalar = double;
  Eigen::RowVector2d pa = p - a, ba = b - a;
  Scalar h = std::max(Scalar(0), std::min(Scalar(1), pa.dot(ba)/ba.dot(ba)));
  return (pa - ba*h).norm() - r;
}

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  igl::readDMAT(argc>1?argv[1]: TUTORIAL_SHARED_PATH "/libigl-cursive-V.dmat",V);
  igl::readDMAT(argc>1?argv[2]: TUTORIAL_SHARED_PATH "/libigl-cursive-E.dmat",E);
  Eigen::MatrixXi F(E.rows(),3);
  F<<E,E.col(1);

  Eigen::MatrixXd GV;
  Eigen::MatrixXi GF;
  {
    Eigen::RowVector2d min_corner = V.colwise().minCoeff();
    Eigen::RowVector2d max_corner = V.colwise().maxCoeff();
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

  tictoc();
  Eigen::Matrix<double,Eigen::Dynamic,2,Eigen::RowMajor> EB1,EB2;
  igl::box_simplices(V,E,EB1,EB2);
  Eigen::Matrix<double,Eigen::Dynamic,2,Eigen::RowMajor> B1,B2;
  Eigen::VectorXi leaf;
  igl::eytzinger_aabb(EB1,EB2,B1,B2,leaf);
  printf("AABB build: %f s\n",tictoc());

  Eigen::VectorXi I,C;
  igl::eytzinger_aabb_winding_number_tree(E,leaf,I,C);
  printf("Winding number tree build: %f s\n",tictoc());


  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(GV,GF);

  enum Quantity
  {
    SIGNED_DISTANCE,
    UNSIGNED_DISTANCE,
    WINDING_NUMBER
  } quantity = SIGNED_DISTANCE;
  bool naive = false;

  const auto update = [&]()
  {
    tictoc();
    Eigen::VectorXd D(GV.rows());
    for(int i=0;i<GV.rows();i++)
    {
      Eigen::RowVector2d p = GV.row(i);
      double di,wi;
      if(naive)
      {
        di = std::numeric_limits<double>::infinity();
        wi = 0;
        for(int e = 0;e<E.rows();e++)
        {
          if(quantity == UNSIGNED_DISTANCE || quantity == SIGNED_DISTANCE)
          {
            double d = sdCapsule(
              p,
              V.row(E(e,0)),
              V.row(E(e,1)),
              0.0);
            if(d < di){ di = d; }
          }
          if(quantity == WINDING_NUMBER || quantity == SIGNED_DISTANCE)
          {
            double angle = std::atan2(
              (V(E(e,0),0)-p(0))*(V(E(e,1),1)-p(1)) - (V(E(e,0),1)-p(1))*(V(E(e,1),0)-p(0)),
              (V(E(e,0),0)-p(0))*(V(E(e,1),0)-p(0)) + (V(E(e,0),1)-p(1))*(V(E(e,1),1)-p(1))
            );
            wi += angle;
          }
        }
      }else
      {
        const std::function<double(const int)> primitive_p = [&](const int e)-> double
        {
          double d = sdCapsule(
            p,
            V.row(E(e,0)),
            V.row(E(e,1)),
            0.0);
          return d;
        };
        if(quantity == UNSIGNED_DISTANCE || quantity == SIGNED_DISTANCE)
        {
          igl::eytzinger_aabb_sdf<false>(p,primitive_p,B1,B2,leaf,di);
        }
        if(quantity == WINDING_NUMBER || quantity == SIGNED_DISTANCE)
        {
          igl::eytzinger_aabb_winding_number(p,V,E,B1,B2,leaf,I,C,wi);
        }
      }
      switch(quantity)
      {
        case SIGNED_DISTANCE:
          D(i) = di * (std::abs(wi) > 0.5?-1:1);
          break;
        case UNSIGNED_DISTANCE:
          D(i) = di;
          break;
        case WINDING_NUMBER:
          D(i) = wi;
          break;
      }
    }
    printf("Compute %s %s: %f s\n",
      naive?"naive":"aabb",
      quantity==SIGNED_DISTANCE?"signed distance":
      quantity==UNSIGNED_DISTANCE?"unsigned distance":"winding number",
      tictoc());
    viewer.data().set_data(D);
  };

  update();
  viewer.data().show_lines = false;
  viewer.data().line_width = 2;
  viewer.core().lighting_factor = 0;
  viewer.data().set_edges(V,E,Eigen::RowVector3d(0,0,0));
  // key
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)->bool
  {
    switch(key)
    {
      default:
        return false;
      case 'N':
      case 'n':
        naive = !naive;
        break;
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
 


  viewer.launch();
  return 0;
}
