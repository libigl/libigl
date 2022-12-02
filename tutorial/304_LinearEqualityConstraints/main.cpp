#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <iostream>

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off",V,F);

  // Two fixed points
  VectorXi b(2,1);
  // Left hand, left foot
  b<<4331,5957;
  VectorXd bc(2,1);
  bc<<1,-1;

  // Construct Laplacian and mass matrix
  SparseMatrix<double> L,M,Minv,Q;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);
  // Bi-Laplacian
  Q = L * (Minv * L);
  // Zero linear term
  VectorXd B = VectorXd::Zero(V.rows(),1);

  VectorXd Z,Z_const;
  {
    // Alternative, short hand
    igl::min_quad_with_fixed_data<double> mqwf;
    // Empty constraints
    VectorXd Beq;
    SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(Q,b,Aeq,true,mqwf);
    igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);
  }

  {
    igl::min_quad_with_fixed_data<double> mqwf;
    // Constraint forcing difference of two points to be 0
    SparseMatrix<double> Aeq(1,V.rows());
    // Right hand, right foot
    Aeq.insert(0,6074) = 1;
    Aeq.insert(0,6523) = -1;
    Aeq.makeCompressed();
    VectorXd Beq(1,1);
    Beq(0) = 0;
    igl::min_quad_with_fixed_precompute(Q,b,Aeq,true,mqwf);
    igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z_const);
  }

  // Use same color axes
  const double min_z = std::min(Z.minCoeff(),Z_const.minCoeff());
  const double max_z = std::max(Z.maxCoeff(),Z_const.maxCoeff());

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().show_lines = false;
  viewer.data().set_data(Z,min_z,max_z);

  viewer.callback_key_down =
    [&Z,&Z_const,&min_z,&max_z](igl::opengl::glfw::Viewer& viewer,unsigned char key,int mod)->bool
    {
      if(key == ' ')
      {
        static bool toggle = true;
        viewer.data().set_data(toggle?Z_const:Z,min_z,max_z);
        toggle = !toggle;
        return true;
      }else
      {
        return false;
      }
    };
  cout<<
    "Press [space] to toggle between unconstrained and constrained."<<endl;
  viewer.launch();
}
