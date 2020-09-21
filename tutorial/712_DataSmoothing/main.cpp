#include <igl/read_triangle_mesh.h>
#include <igl/hessian_energy.h>
#include <igl/curved_hessian_energy.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/isolines_map.h>
#include <igl/parula.h>
#include <igl/vertex_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/heat_geodesics.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <iostream>
#include <set>
#include <limits>
#include <stdlib.h>

#include "tutorial_shared_path.h"



int main(int argc, char * argv[])
{
  typedef Eigen::SparseMatrix<double> SparseMat;
  srand(57);
  
  //Read our mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh(
                              argc>1?argv[1]: TUTORIAL_SHARED_PATH "/beetle.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }
  
  //Constructing an exact function to smooth
  igl::HeatGeodesicsData<double> hgData;
  igl::heat_geodesics_precompute(V, F, hgData);
  Eigen::VectorXd heatDist;
  Eigen::VectorXi gamma(1); gamma << 1947; //1631;
  igl::heat_geodesics_solve(hgData, gamma, heatDist);
  Eigen::VectorXd zexact =
  0.1*(heatDist.array() + (-heatDist.maxCoeff())).pow(2)
  + 3*V.block(0,1,V.rows(),1).array().cos();
  
  //Make the exact function noisy
  const double s = 0.1*(zexact.maxCoeff() - zexact.minCoeff());
  Eigen::VectorXd znoisy = zexact + s*Eigen::VectorXd::Random(zexact.size());
  
  //Constructing the squared Laplacian and squared Hessian energy
  SparseMat L, M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
  Eigen::SimplicialLDLT<SparseMat> solver(M);
  SparseMat MinvL = solver.solve(L);
  SparseMat QL = L.transpose()*MinvL;
  SparseMat QH;
  igl::hessian_energy(V, F, QH);
  SparseMat QcH;
  igl::curved_hessian_energy(V, F, QcH);
  
  //Solve to find Laplacian-smoothed Hessian-smoothed, and
  // curved-Hessian-smoothed solutions
  const double al = 3e-7;
  Eigen::SimplicialLDLT<SparseMat> lapSolver(al*QL + (1.-al)*M);
  Eigen::VectorXd zl = lapSolver.solve(al*M*znoisy);
  const double ah = 2e-7;
  Eigen::SimplicialLDLT<SparseMat> hessSolver(ah*QH + (1.-ah)*M);
  Eigen::VectorXd zh = hessSolver.solve(ah*M*znoisy);
  const double ach = 3e-7;
  Eigen::SimplicialLDLT<SparseMat> curvedHessSolver(al*QcH + (1.-ach)*M);
  Eigen::VectorXd zch = curvedHessSolver.solve(ach*M*znoisy);
  
  //Viewer that shows all functions: zexact, znoisy, zl, zh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.callback_key_down =
  [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
  {
    //Graduate result to show isolines, then compute color matrix
    const Eigen::VectorXd* z;
    switch(key) {
      case '1':
        z = &zexact;
        break;
      case '2':
        z = &znoisy;
        break;
      case '3':
        z = &zl;
        break;
      case '4':
        z = &zh;
        break;
      case '5':
        z = &zch;
        break;
      default:
        return false;
    }
    viewer.data().set_data(*z);
    return true;
  };
  std::cout << R"(Smoothing a noisy function.
Usage:
1  Show original function
2  Show noisy function
3  Biharmonic smoothing (zero Neumann boundary)
4  Biharmonic smoothing (natural planar Hessian boundary)
5  Biharmonic smoothing (natural curved Hessian boundary)

)";
  Eigen::MatrixXd CM;
  igl::parula(Eigen::VectorXd::LinSpaced(21,0,1).eval(),false,CM);
  igl::isolines_map(Eigen::MatrixXd(CM),CM);
  viewer.data().set_colormap(CM);
  viewer.data().set_data(znoisy);
  viewer.launch();
  
  
  //Constructing a step function to smooth
  Eigen::VectorXd zstep = Eigen::VectorXd::Zero(V.rows());
  for(int i=0; i<V.rows(); ++i) {
    zstep(i) = V(i,2)<-0.25 ? 1. : (V(i,2)>0.31 ? 2. : 0);
  }
  
  //Smooth that function
  const double sl = 2e-5;
  Eigen::SimplicialLDLT<SparseMat> stepLapSolver(sl*QL + (1.-sl)*M);
  Eigen::VectorXd stepzl = stepLapSolver.solve(al*M*zstep);
  const double sh = 6e-6;
  Eigen::SimplicialLDLT<SparseMat> stepHessSolver(sh*QH + (1.-sh)*M);
  Eigen::VectorXd stepzh = stepHessSolver.solve(ah*M*zstep);
  const double sch = 2e-5;
  Eigen::SimplicialLDLT<SparseMat> stepCurvedHessSolver(sl*QcH + (1.-sch)*M);
  Eigen::VectorXd stepzch = stepCurvedHessSolver.solve(ach*M*zstep);
  
  //Display functions
  igl::opengl::glfw::Viewer viewer2;
  viewer2.data().set_mesh(V,F);
  viewer2.data().show_lines = false;
  viewer2.callback_key_down =
  [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
  {
    //Graduate result to show isolines, then compute color matrix
    const Eigen::VectorXd* z;
    switch(key) {
      case '1':
        z = &zstep;
        break;
      case '2':
        z = &stepzl;
        break;
      case '3':
        z = &stepzh;
        break;
      case '4':
        z = &stepzch;
        break;
      default:
        return false;
    }
    viewer.data().set_data(*z);
    return true;
  };
  std::cout << R"(Smoothing a step function.
Usage:
1  Show step function
2  Biharmonic smoothing (zero Neumann boundary)
3  Biharmonic smoothing (natural planar Hessian boundary)
4  Biharmonic smoothing (natural curved Hessian boundary)

)";
  
  viewer2.data().set_colormap(CM);
  viewer2.data().set_data(zstep);
  viewer2.launch();
  
  return 0;
}
