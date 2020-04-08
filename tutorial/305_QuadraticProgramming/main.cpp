#include <igl/active_set.h>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>
#include "tutorial_shared_path.h"
  
Eigen::VectorXi b;
Eigen::VectorXd B,bc,lx,ux,Beq,Bieq,Z;
Eigen::SparseMatrix<double> Q,Aeq,Aieq;

void solve(igl::opengl::glfw::Viewer &viewer)
{
  using namespace std;
  igl::active_set_params as;
  as.max_iter = 8;
  igl::active_set(Q,B,b,bc,Aeq,Beq,Aieq,Bieq,lx,ux,as,Z);
  viewer.data().set_data(Z);
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)
{
  switch(key)
  {
    case '.':
      Beq(0) *= 2.0;
      solve(viewer);
      return true;
    case ',':
      Beq(0) /= 2.0;
      solve(viewer);
      return true;
    case ' ':
      solve(viewer);
      return true;
    default:
      return false;
  }
}


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off",V,F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().show_lines = false;
  viewer.callback_key_down = &key_down;

  // One fixed point
  b.resize(1,1);
  // point on belly.
  b<<2556;
  bc.resize(1,1);
  bc<<1;

  // Construct Laplacian and mass matrix
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  //M = (M/M.diagonal().maxCoeff()).eval();
  igl::invert_diag(M,Minv);
  // Bi-Laplacian
  Q = L.transpose() * (Minv * L);
  // Zero linear term
  B = VectorXd::Zero(V.rows(),1);

  // Lower and upper bound
  lx = VectorXd::Zero(V.rows(),1);
  ux = VectorXd::Ones(V.rows(),1);

  // Equality constraint constrain solution to sum to 1
  Beq.resize(1,1);
  Beq(0) = 0.08;
  Aeq = M.diagonal().sparseView().transpose();
  // (Empty inequality constraints)
  solve(viewer);
  cout<<
    "Press '.' to increase scale and resolve."<<endl<<
    "Press ',' to decrease scale and resolve."<<endl;

  viewer.launch();
}
