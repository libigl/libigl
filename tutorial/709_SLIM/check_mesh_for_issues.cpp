#include "check_mesh_for_issues.h"
#include <iostream>

#include <igl/doublearea.h>
#include <igl/vertex_components.h>
#include <igl/euler_characteristic.h>
#include <igl/is_edge_manifold.h>
#include <igl/adjacency_matrix.h>
#include <Eigen/Sparse>

void check_mesh_for_issues(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {

  using namespace std;
  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F,A);

  Eigen::MatrixXi C, Ci;
  igl::vertex_components(A, C, Ci);

  int connected_components = Ci.rows();
  if (connected_components!=1) {
    cout << "Error! Input has multiple connected components" << endl; exit(1);
  }
  int euler_char = igl::euler_characteristic(F);
  if (euler_char!=1) 
  {
    cout << 
      "Error! Input does not have a disk topology, it's euler char is " << 
      euler_char << endl; 
    exit(1);
  }
  bool is_edge_manifold = igl::is_edge_manifold(F);
  if (!is_edge_manifold) {
    cout << "Error! Input is not an edge manifold" << endl; exit(1);
  }

  Eigen::VectorXd areas; igl::doublearea(V,F,areas);
  const double eps = 1e-14;
  for (int i = 0; i < areas.rows(); i++) {
    if (areas(i) < eps) {
      cout << "Error! Input has zero area faces" << endl; exit(1);
    }
  }
}
