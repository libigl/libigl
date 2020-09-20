#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/orient_halfedges.h>
#include <igl/cr_vector_laplacian.h>
#include <igl/cr_vector_mass.h>
#include <igl/edge_midpoints.h>
#include <igl/edge_vectors.h>
#include <igl/average_from_edges_onto_vertices.h>
#include <igl/PI.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/Geometry>

#include <iostream>
#include <set>
#include <limits>
#include <stdlib.h>

#include "tutorial_shared_path.h"



int main(int argc, char * argv[])
{
  typedef Eigen::SparseMatrix<double> SparseMat;
  
  //Constants used for smoothing
  const double howMuchToSmoothBy = 1e-1;
  const int howManySmoothingInterations = 50;
  
  //Read our mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/elephant.obj",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }
  
  //Orient edges for plotting
  Eigen::MatrixXi E, oE;
  igl::orient_halfedges(F, E, oE);
  
  //Compute edge midpoints & edge vectors
  Eigen::MatrixXd edgeMps, parVec, perpVec;
  igl::edge_midpoints(V, F, E, oE, edgeMps);
  igl::edge_vectors(V, F, E, oE, parVec, perpVec);
  
  //Constructing a function to add noise to
  const auto zraw_function = [] (const Eigen::Vector3d& x) {
    return Eigen::Vector3d(0.2*x(1) + cos(2*x(1)+0.2),
                           0.5*x(0) + 0.15,
                           0.3*cos(0.2+igl::PI*x(2)));
  };
  
  Eigen::VectorXd zraw(2*edgeMps.rows());
  for(int i=0; i<edgeMps.rows(); ++i) {
    const Eigen::Vector3d f = zraw_function(edgeMps.row(i));
    zraw(i) = f.dot(parVec.row(i));
    zraw(i+edgeMps.rows()) = f.dot(perpVec.row(i));
  }
  
  //Add noise
  srand(71);
  const double l = 15;
  Eigen::VectorXd znoisy = zraw + l*Eigen::VectorXd::Random(zraw.size());
  
  //Denoise function using the vector Dirichlet energy
  Eigen::VectorXd zsmoothed = znoisy;
  for(int i=0; i<howManySmoothingInterations; ++i) {
    //Compute Laplacian and mass matrix
    SparseMat L, M;
    igl::cr_vector_mass(V, F, E, oE, M);
    igl::cr_vector_laplacian(V, F, E, oE, L);
    
    //Implicit step
    Eigen::SimplicialLDLT<SparseMat> rhsSolver(M + howMuchToSmoothBy*L);
    zsmoothed = rhsSolver.solve(M*zsmoothed);
  }
  
  //Convert vector fields for plotting
  const auto cr_result_to_vecs_and_colors = [&]
  (const Eigen::VectorXd& z, Eigen::MatrixXd& vecs, Eigen::MatrixXd& colors) {
    vecs.resize(edgeMps.rows(), 3);
    for(int i=0; i<edgeMps.rows(); ++i) {
      vecs.row(i) = z(i)*parVec.row(i)
      + z(i+edgeMps.rows())*perpVec.row(i);
    }
    igl::average_from_edges_onto_vertices
    (F, E, oE, vecs.rowwise().norm(), colors);
  };
  Eigen::MatrixXd noisyvecs, noisycolors, smoothedvecs, smoothedcolors,
  rawvecs, rawcolors;
  cr_result_to_vecs_and_colors(znoisy, noisyvecs, noisycolors);
  cr_result_to_vecs_and_colors(zsmoothed, smoothedvecs, smoothedcolors);
  cr_result_to_vecs_and_colors(zraw, rawvecs, rawcolors);
  
  
  //Viewer that shows noisy and denoised functions
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.callback_key_down =
  [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
  {
    const Eigen::MatrixXd *vecs, *colors;
    switch(key) {
      case '1':
        vecs = &rawvecs;
        colors = &rawcolors;
        break;
      case '2':
        vecs = &noisyvecs;
        colors = &noisycolors;
        break;
      case '3':
        vecs = &smoothedvecs;
        colors = &smoothedcolors;
        break;
      default:
        return false;
    }
    viewer.data().set_data(*colors);
    viewer.data().clear_edges();
    const double s = 0.08; //How much to scale vectors during plotting
    viewer.data().add_edges(edgeMps, edgeMps + s*(*vecs),
                            Eigen::RowVector3d(0.1, 0.1, 0.1));
    return true;
  };
  std::cout << R"(Usage:
1  Show raw function
2  Show noisy function
3  Show smoothed function

)";
  Eigen::MatrixXd CM;
  igl::parula(Eigen::VectorXd::LinSpaced
              (500,znoisy.minCoeff(), znoisy.maxCoeff()).eval(), true, CM);
  viewer.data().set_colormap(CM);
  viewer.callback_key_down(viewer, '1', 0);
  viewer.launch();
  
  return 0;
}
