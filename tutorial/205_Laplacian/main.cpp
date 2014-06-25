#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/repdiag.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
#include <igl/viewer/Viewer.h>

#include <iostream>

Eigen::MatrixXd V,U;
Eigen::MatrixXi F;
Eigen::SparseMatrix<double> L;
igl::Viewer viewer;

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OFF format
  igl::readOFF("../shared/cow.off", V, F);

  // Compute Laplace-Beltrami operator: #V by #V
  igl::cotmatrix(V,F,L);

  // Alternative construction of same Laplacian
  SparseMatrix<double> G,K;
  // Gradient/Divergence
  igl::grad(V,F,G);
  // Diagonal per-triangle "mass matrix"
  VectorXd dblA;
  igl::doublearea(V,F,dblA);
  // Place areas along diagonal #dim times
  const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
  // Laplacian K built as discrete divergence of gradient or equivalently
  // discrete Dirichelet energy Hessian
  K = -G.transpose() * T * G;
  cout<<"|K-L|: "<<(K-L).norm()<<endl;

  const auto &key_down = [](igl::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case 'r':
      case 'R':
        U = V;
        break;
      case ' ':
      {
        // Recompute just mass matrix on each step
        SparseMatrix<double> M;
        igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_VORONOI,M);
        // Solve (M-delta*L) U = M*U
        const auto & S = (M - 0.001*L);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
        assert(solver.info() == Eigen::Success);
        U = solver.solve(M*U).eval();
        // Normalize to unit surface area (important for numerics)
        U.array() /= sqrt(M.diagonal().array().sum());
        break;
      }
      default:
        return false;
    }
    // Send new positions, update normals, recenter
    viewer.set_vertices(U);
    viewer.compute_normals();
    viewer.align_camera_center();
    return true;
  };


  // Use original normals as pseudo-colors
  MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  MatrixXd C = N.rowwise().normalized().array()*0.5+0.5;

  // Initialize smoothing with base mesh
  U = V;
  viewer.set_mesh(U, F);
  viewer.set_colors(C);
  viewer.callback_key_down = key_down;

  cout<<"Press [space] to smooth."<<endl;;
  cout<<"Press [r] to reset."<<endl;;
  return viewer.launch();
}
