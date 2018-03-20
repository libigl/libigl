#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V,U;
Eigen::MatrixXi F;
int c=0;
double bbd = 1;
bool twod = 0;
int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  VectorXd D;
  if(!read_triangle_mesh(
     argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off",V,F))
  {
    cout<<"failed to load mesh"<<endl;
  }
  twod = V.col(2).minCoeff()==V.col(2).maxCoeff();
  bbd = (V.colwise().maxCoeff()-V.colwise().minCoeff()).norm();
  SparseMatrix<double> L,M;
  cotmatrix(V,F,L);
  L = (-L).eval();
  massmatrix(V,F,MASSMATRIX_TYPE_DEFAULT,M);
  const size_t k = 5;
  if(!eigs(L,M,k+1,EIGS_TYPE_SM,U,D))
  {
    cout<<"failed."<<endl;
  }
  // Normalize
  U = ((U.array()-U.minCoeff())/(U.maxCoeff()-U.minCoeff())).eval();

  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer & viewer,unsigned char key,int)->bool
  {
    switch(key)
    {
      default:
        return false;
      case ' ':
      {
        U = U.rightCols(k).eval();
        // Rescale eigen vectors for visualization
        VectorXd Z =
          bbd*0.5*U.col(c);
        Eigen::MatrixXd C;
        igl::parula(U.col(c).eval(),false,C);
        c = (c+1)%U.cols();
        if(twod)
        {
          V.col(2) = Z;
        }
        viewer.data().set_mesh(V,F);
        viewer.data().compute_normals();
        viewer.data().set_colors(C);
        return true;
      }
    }
  };
  viewer.callback_key_down(viewer,' ',0);
  viewer.data().show_lines = false;
  std::cout<<
R"(
  [space] Cycle through eigen modes
)";
  viewer.launch();
}
