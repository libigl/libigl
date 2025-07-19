#include <igl/read_triangle_mesh.h>
#include <igl/eytzinger_aabb.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/unique.h>
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>
#include <limits>

#include <Eigen/Geometry>
template <typename T>
inline T sign(const T & x)
{
  return (x > T(0)) - (x < T(0));
}

double dot2(const Eigen::RowVector3d & v)
{
  return v.dot(v);
}

double clamp(const double & x, const double & a, const double & b)
{
  return std::max(a, std::min(b, x));
}

// https://iquilezles.org/articles/triangledistance/
template <typename T1, typename T2, typename T3, typename Tp>
typename Tp::Scalar udTriangle(
    const T1 & v1, const T2 & v2, const T3 & v3, const Tp & p)
{
  using vec3 = Eigen::Matrix<typename Tp::Scalar, 3, 1>;
  // prepare data    
  vec3 v21 = v2 - v1; vec3 p1 = p - v1;
  vec3 v32 = v3 - v2; vec3 p2 = p - v2;
  vec3 v13 = v1 - v3; vec3 p3 = p - v3;
  vec3 nor =  v21.cross( v13 );

    return sqrt( // inside/outside test    
                 (sign(v21.cross(nor).dot(p1)) + 
                  sign(v32.cross(nor).dot(p2)) + 
                  sign(v13.cross(nor).dot(p3))<2.0) 
                  ?
                  // 3 edges    
                  std::min( std::min( 
                  dot2(v21*clamp(v21.dot(p1)/dot2(v21),0.0,1.0)-p1), 
                  dot2(v32*clamp(v32.dot(p2)/dot2(v32),0.0,1.0)-p2) ), 
                  dot2(v13*clamp(v13.dot(p3)/dot2(v13),0.0,1.0)-p3) )
                  :
                  // 1 face    
                  nor.dot(p1)*nor.dot(p1)/dot2(nor) );
}

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  // Read mesh
  Eigen::Matrix<double,Eigen::Dynamic,3> V;
  Eigen::Matrix<int,Eigen::Dynamic,3> F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }

  V.rowwise() -= 0.5* (V.colwise().maxCoeff() + V.colwise().minCoeff());
  V /= V.array().abs().maxCoeff();
  V *= 0.8;


  tictoc();
  Eigen::Matrix<double,Eigen::Dynamic,3> PB1,PB2;
  PB1.setConstant(F.rows(),V.cols(),std::numeric_limits<double>::infinity());
  PB2.setConstant(F.rows(),V.cols(),-std::numeric_limits<double>::infinity());
  for(int f = 0; f < F.rows(); f++) 
  {
    for(int c = 0; c < F.cols(); c++) 
    {
      for(int d = 0; d < V.cols(); d++) 
      {
        PB1(f,d) = std::min(PB1(f,d),V(F(f,c),d));
        PB2(f,d) = std::max(PB2(f,d),V(F(f,c),d));
      }
    }
  }
  printf("%-20s: %g secs\n","PB",tictoc());

  const int ns = 64;
  Eigen::RowVector3i res(ns,ns,ns);
  Eigen::Matrix<double,Eigen::Dynamic,3> GV;
  igl::grid(res,GV);
  GV *= 2;
  GV.array() -= 1;

  Eigen::Matrix<double,Eigen::Dynamic,3> B1,B2;
  Eigen::VectorXi leaf;
  igl::eytzinger_aabb(PB1,PB2,B1,B2,leaf);
  printf("%-20s: %g secs\n","eytzinger_aabb",tictoc());
  {
    Eigen::VectorXi U;
    igl::unique(leaf,U);
    printf("%d → %d\n",PB1.rows(),U.size()-2);
  }

    //Eigen::RowVector3d p(0,0,0);
    //double f ;
    //const std::function<double(const int i)> primitive = [&](const int i)
    //{
    //  return udTriangle( V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), p);
    //};
    //igl::eytzinger_aabb_sdf(p,primitive,B1,B2,leaf,f);
    //printf("sdf(%g,%g,%g) = %g\n",p(0),p(1),p(2),f);

  tictoc();
  Eigen::VectorXd S;
  const std::function<double(const Eigen::RowVector3d &,const int i)> primitive = 
    [&](const Eigen::RowVector3d & p, const int i)
  {
    return udTriangle( V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), p);
  };
  igl::eytzinger_aabb_sdf(GV,primitive,B1,B2,leaf,S);
  printf("%-20s: %g secs\n","eytzinger_aabb_sdf",tictoc());

  Eigen::Matrix<double,Eigen::Dynamic,3> mV;
  Eigen::Matrix<int,Eigen::Dynamic,3> mF;
  const double isovalue = (2.0/(res(0)-1))*2;
  igl::marching_cubes(S,GV,res(0),res(1),res(2),isovalue,mV,mF);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = false;
  viewer.append_mesh();
  viewer.data().set_mesh(mV,mF);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = false;
  viewer.data().set_colors(Eigen::RowVector3d(0.8,0.2,0.2));

  viewer.launch();
  return 0;
}
