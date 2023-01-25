#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <algorithm>

int main(int argc, char *argv[])
{
  igl::opengl::glfw::Viewer v;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);

  // Create a floor
  Eigen::MatrixXd fV, fU;
  Eigen::MatrixXi fF;
  {
    igl::triangulated_grid(2,2,fU,fF);
    fV = fU;
    fV.array() -= 0.5;
    fV.array() *= 2 * 2 * 
      (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
    fV = (fV * (Eigen::Matrix<double,2,3>()<<1,0,0,0,0,-1).finished() ).eval();
    fV.col(1).array() += V.col(1).minCoeff();
  }

  const int s = 16;
  const int f = 100;
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> X(s*f,s*f);
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(s*f,s*f);
  for(int i = 0;i<s*f;i++)
  {
    const double x = double(i)/double(s*f-1)*2-1;
    for(int j = 0;j<s*f;j++)
    {
      const int u = i/f;
      const int v = j/f;
      const double y = double(j)/double(s*f-1)*2-1;
      const double r = std::min(std::max( (1.0 - sqrt(x*x+y*y))*3.0 ,0.0),1.0);
      //const double a = 3*r*r - 2*r*r*r;
      const auto smooth_step = [](const double w)
      {
        return ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) ;
      };
      double a = smooth_step(r);
      X(i,j) = u%2 == v%2 ? 245 : 235;
      A(i,j) = a * 255;
    }
  }

  v.data().set_mesh(fV,fF);
  v.data().set_uv(fU);
  v.data().uniform_colors(Eigen::Vector3d(0.1,0.1,0.1),Eigen::Vector3d(1,1,1),Eigen::Vector3d(0,0,0));
  v.data().set_texture(X,X,X,A);
  v.data().show_texture = true;
  v.data().show_lines = false;

  v.append_mesh();
  v.data().set_mesh(V,F);
  v.data().show_lines = false;

  v.launch();
}
