#include <igl/copyleft/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  MatrixXi F;
  MatrixXd V;
  // Read in inputs as double precision floating point meshes
  read_triangle_mesh(
      TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  // number of vertices on the largest side
  const int s = 50;
  const RowVector3d Vmin = V.colwise().minCoeff();
  const RowVector3d Vmax = V.colwise().maxCoeff();
  const double h = (Vmax-Vmin).maxCoeff()/(double)s;
  const RowVector3i res = (s*((Vmax-Vmin)/(Vmax-Vmin).maxCoeff())).cast<int>();
  // create grid
  cout<<"Creating grid..."<<endl;
  MatrixXd GV(res(0)*res(1)*res(2),3);
  for(int zi = 0;zi<res(2);zi++)
  {
    const auto lerp = [&](const int di, const int d)->double
      {return Vmin(d)+(double)di/(double)(res(d)-1)*(Vmax(d)-Vmin(d));};
    const double z = lerp(zi,2);
    for(int yi = 0;yi<res(1);yi++)
    {
      const double y = lerp(yi,1);
      for(int xi = 0;xi<res(0);xi++)
      {
        const double x = lerp(xi,0);
        GV.row(xi+res(0)*(yi + res(1)*zi)) = RowVector3d(x,y,z);
      }
    }
  }
  // compute values
  cout<<"Computing distances..."<<endl;
  VectorXd S,B;
  {
    VectorXi I;
    MatrixXd C,N;
    signed_distance(GV,V,F,SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);
    // Convert distances to binary inside-outside data --> aliasing artifacts
    B = S;
    for_each(B.data(),B.data()+B.size(),[](double& b){b=(b>0?1:(b<0?-1:0));});
  }
  cout<<"Marching cubes..."<<endl;
  MatrixXd SV,BV;
  MatrixXi SF,BF;
  igl::copyleft::marching_cubes(S,GV,res(0),res(1),res(2),SV,SF);
  igl::copyleft::marching_cubes(B,GV,res(0),res(1),res(2),BV,BF);

  cout<<R"(Usage:
'1'  Show original mesh.
'2'  Show marching cubes contour of signed distance.
'3'  Show marching cubes contour of indicator function.
)";
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(SV,SF);
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
          viewer.data().clear();
          viewer.data().set_mesh(V,F);
          break;
        case '2':
          viewer.data().clear();
          viewer.data().set_mesh(SV,SF);
          break;
        case '3':
          viewer.data().clear();
          viewer.data().set_mesh(BV,BF);
          break;
      }
      viewer.data().set_face_based(true);
      return true;
    };
  viewer.launch();
}
