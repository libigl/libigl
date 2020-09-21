#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/remove_unreferenced.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/orient_halfedges.h>
#include <igl/cr_vector_laplacian.h>
#include <igl/cr_vector_mass.h>
#include <igl/crouzeix_raviart_cotmatrix.h>
#include <igl/crouzeix_raviart_massmatrix.h>
#include <igl/edge_midpoints.h>
#include <igl/edge_vectors.h>
#include <igl/average_from_edges_onto_vertices.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/heat_geodesics.h>

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
  typedef Eigen::Matrix<double, 1, 1> Vector1d;
  typedef Eigen::Matrix<int, 1, 1> Vector1i;
  
  //Constants used for smoothing
  const double howMuchToSmoothBy = 1e-1;
  const int howManySmoothingInterations = 50;
  
  //Read our mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/cheburashka.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }
  
  //Compute vector Laplacian and mass matrix
  Eigen::MatrixXi E, oE;//Compute Laplacian and mass matrix
  SparseMat vecL, vecM;
  igl::cr_vector_mass(V, F, E, oE, vecM);
  igl::cr_vector_laplacian(V, F, E, oE, vecL);
  const int m = vecL.rows()/2; //The number of edges in the mesh
  
  //Convert the E / oE matrix format to list of edges / EMAP format required
  // by the functions constructing scalar Crouzeix-Raviart functions
  Eigen::MatrixXi Elist(m,2), EMAP(3*F.rows(),1);
  for(int i=0; i<F.rows(); ++i) {
    for(int j=0; j<3; ++j) {
      const int e = E(i,j);
      EMAP(i+j*F.rows()) = e;
      if(oE(i,j)>0) {
        Elist.row(e) << F(i, (j+1)%3), F(i, (j+2)%3);
      }
    }
  }
  SparseMat scalarL, scalarM;
  igl::crouzeix_raviart_massmatrix(V, F, Elist, EMAP, scalarM);
  igl::crouzeix_raviart_cotmatrix(V, F, Elist, EMAP, scalarL);
  
  //Compute edge midpoints & edge vectors
  Eigen::MatrixXd edgeMps, parVec, perpVec;
  igl::edge_midpoints(V, F, E, oE, edgeMps);
  igl::edge_vectors(V, F, E, oE, parVec, perpVec);
  
  //Perform the vector heat method
  const int initialIndex = 14319;
  const double initialPara=0.95, initialPerp=0.08;
  const double t = 0.01;
  
  SparseMat Aeq;
  Eigen::VectorXd Beq;
  Eigen::VectorXi known = Eigen::Vector2i(initialIndex, initialIndex+m);
  Eigen::VectorXd knownVals = Eigen::Vector2d(initialPara, initialPerp);
  Eigen::VectorXd Y0 = Eigen::VectorXd::Zero(2*m), Yt;
  Y0(initialIndex) = initialPara; Y0(initialIndex+m) = initialPerp;
  igl::min_quad_with_fixed
  (SparseMat(vecM+t*vecL), Eigen::VectorXd(-vecM*Y0), known, knownVals,
   Aeq, Beq, false, Yt);
  
  Eigen::VectorXd u0 = Eigen::VectorXd::Zero(m), ut;
  u0(initialIndex) = sqrt(initialPara*initialPara + initialPerp*initialPerp);
  Eigen::VectorXi knownScal = Vector1i(initialIndex);
  Eigen::VectorXd knownScalVals = Vector1d(u0(initialIndex));
  igl::min_quad_with_fixed
  (SparseMat(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*u0), knownScal,
   knownScalVals, Aeq, Beq, false, ut);
  
  Eigen::VectorXd phi0 = Eigen::VectorXd::Zero(m), phit;
  phi0(initialIndex) = 1;
  Eigen::VectorXd knownScalValsPhi = Vector1d(1);
  igl::min_quad_with_fixed
  (SparseMat(scalarM+t*scalarL), Eigen::VectorXd(-scalarM*phi0), knownScal,
   knownScalValsPhi, Aeq, Beq, false, phit);
  
  Eigen::ArrayXd Xtfactor = ut.array() /
  (phit.array() * (Yt.array().segment(0,m)*Yt.array().segment(0,m)
                   + Yt.array().segment(m,m)*Yt.array().segment(m,m)).sqrt());
  Eigen::VectorXd Xt(2*m);
  Xt.segment(0,m) = Xtfactor * Yt.segment(0,m).array();
  Xt.segment(m,m) = Xtfactor * Yt.segment(m,m).array();
  
  
  //Compute scalar heat colors
  igl::HeatGeodesicsData<double> hgData;
  igl::heat_geodesics_precompute(V, F, hgData);
  Eigen::VectorXd heatColor;
  Eigen::VectorXi gamma = Elist.row(initialIndex);
  igl::heat_geodesics_solve(hgData, gamma, heatColor);
  
  
  //Convert vector field for plotting
  Eigen::MatrixXd vecs(m, 3);
  for(int i=0; i<edgeMps.rows(); ++i) {
    vecs.row(i) = Xt(i)*parVec.row(i) + Xt(i+edgeMps.rows())*perpVec.row(i);
  }
  
  
  //Viewer that shows parallel transported vector
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().show_lines = false;
  viewer.data().set_data(heatColor.maxCoeff()-heatColor.array(), //invert colormap
                         igl::COLOR_MAP_TYPE_VIRIDIS);
  const double s = 0.012; //How much to scale vectors during plotting
  Eigen::MatrixXd vecColors(m, 3);
  for(int i=0; i<m; ++i) {
    vecColors.row(i) << 0.1, 0.1, 0.1;
  }
  vecColors.row(initialIndex) << 0.9, 0.1, 0.1;
  viewer.data().add_edges(edgeMps, edgeMps + s*vecs, vecColors);
  
  std::cout << R"(The red vector is parallel transported to every point on the surface.
The surface is shaded by geodesic distance from the red vector.
)"
  << std::endl;
  viewer.launch();
  
  return 0;
}
