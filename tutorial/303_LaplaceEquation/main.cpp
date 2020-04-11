#include <igl/boundary_facets.h>
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/jet.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/readOFF.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/unique.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off",V,F);
  // Find boundary edges
  MatrixXi E;
  igl::boundary_facets(F,E);
  // Find boundary vertices
  VectorXi b,IA,IC;
  igl::unique(E,b,IA,IC);
  // List of all vertex indices
  VectorXi all,in;
  igl::colon<int>(0,V.rows()-1,all);
  // List of interior indices
  igl::setdiff(all,b,in,IA);

  // Construct and slice up Laplacian
  SparseMatrix<double> L,L_in_in,L_in_b;
  igl::cotmatrix(V,F,L);
  igl::slice(L,in,in,L_in_in);
  igl::slice(L,in,b,L_in_b);

  // Dirichlet boundary conditions from z-coordinate
  VectorXd bc;
  VectorXd Z = V.col(2);
  igl::slice(Z,b,bc);

  // Solve PDE
  SimplicialLLT<SparseMatrix<double > > solver(-L_in_in);
  VectorXd Z_in = solver.solve(L_in_b*bc);
  // slice into solution
  igl::slice_into(Z_in,in,Z);

  // Alternative, short hand
  igl::min_quad_with_fixed_data<double> mqwf;
  // Linear term is 0
  VectorXd B = VectorXd::Zero(V.rows(),1);
  // Empty constraints
  VectorXd Beq;
  SparseMatrix<double> Aeq;
  // Our cotmatrix is _negative_ definite, so flip sign
  igl::min_quad_with_fixed_precompute((-L).eval(),b,Aeq,true,mqwf);
  igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().show_lines = false;
  viewer.data().set_data(Z);
  viewer.launch();
}
