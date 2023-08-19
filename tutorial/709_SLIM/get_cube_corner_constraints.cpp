
#include "get_cube_corner_constraints.h"
#include <igl/PI.h>
#include <igl/cat.h>
#include <set>

void get_cube_corner_constraints(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc) {
  using namespace std;
  double min_x,max_x,min_y,max_y,min_z,max_z;
  min_x = V.col(0).minCoeff(); max_x = V.col(0).maxCoeff();
  min_y = V.col(1).minCoeff(); max_y = V.col(1).maxCoeff();
  min_z = V.col(2).minCoeff(); max_z = V.col(2).maxCoeff();


  // get all cube corners
  b.resize(8,1); bc.resize(8,3);
  int x;
  for (int i = 0; i < V.rows(); i++) {
    if (V.row(i) == Eigen::RowVector3d(min_x,min_y,min_z)) b(0) = i;
    if (V.row(i) == Eigen::RowVector3d(min_x,min_y,max_z)) b(1) = i;
    if (V.row(i) == Eigen::RowVector3d(min_x,max_y,min_z)) b(2) = i;
    if (V.row(i) == Eigen::RowVector3d(min_x,max_y,max_z)) b(3) = i;
    if (V.row(i) == Eigen::RowVector3d(max_x,min_y,min_z)) b(4) = i;
    if (V.row(i) == Eigen::RowVector3d(max_x,max_y,min_z)) b(5) = i;
    if (V.row(i) == Eigen::RowVector3d(max_x,min_y,max_z)) b(6) = i;
    if (V.row(i) == Eigen::RowVector3d(max_x,max_y,max_z)) b(7) = i;
  }

  // get all cube edges
  std::set<int> cube_edge1; Eigen::VectorXi cube_edge1_vec;
  for (int i = 0; i < V.rows(); i++) {
    if ((V(i,0) == min_x && V(i,1) == min_y)) {
      cube_edge1.insert(i);
    }
  }
  Eigen::VectorXi edge1;
  int_set_to_eigen_vector(cube_edge1, edge1);

  std::set<int> cube_edge2; Eigen::VectorXi edge2;
  for (int i = 0; i < V.rows(); i++) {
    if ((V(i,0) == max_x && V(i,1) == max_y)) {
      cube_edge2.insert(i);
    }
  }
  int_set_to_eigen_vector(cube_edge2, edge2);
  b = igl::cat(1,edge1,edge2);

  std::set<int> cube_edge3; Eigen::VectorXi edge3;
  for (int i = 0; i < V.rows(); i++) {
    if ((V(i,0) == max_x && V(i,1) == min_y)) {
      cube_edge3.insert(i);
    }
  }
  int_set_to_eigen_vector(cube_edge3, edge3);
  b = igl::cat(1,b,edge3);

  std::set<int> cube_edge4; Eigen::VectorXi edge4;
  for (int i = 0; i < V.rows(); i++) {
    if ((V(i,0) == min_x && V(i,1) == max_y)) {
      cube_edge4.insert(i);
    }
  }
  int_set_to_eigen_vector(cube_edge4, edge4);
  b = igl::cat(1,b,edge4);

  bc.resize(b.rows(),3);
  Eigen::Matrix3d m; m = Eigen::AngleAxisd(0.3 * igl::PI, Eigen::Vector3d(1./sqrt(2.),1./sqrt(2.),0.)/*Eigen::Vector3d::UnitX()*/);
  int i = 0;
  for (; i < cube_edge1.size(); i++) {
    Eigen::RowVector3d edge_rot_center(min_x,min_y,(min_z+max_z)/2.);
    bc.row(i) = (V.row(b(i)) - edge_rot_center) * m + edge_rot_center;
  }
  for (; i < cube_edge1.size() + cube_edge2.size(); i++) {
    Eigen::RowVector3d edge_rot_center(max_x,max_y,(min_z+max_z)/2.);
    bc.row(i) = (V.row(b(i)) - edge_rot_center) * m.transpose() + edge_rot_center;
  }
  for (; i < cube_edge1.size() + cube_edge2.size() + cube_edge3.size(); i++) {
    bc.row(i) = 0.75*V.row(b(i));
  }
  for (; i < b.rows(); i++) {
    bc.row(i) = 0.75*V.row(b(i));
  }
}

void int_set_to_eigen_vector(const std::set<int>& int_set, Eigen::VectorXi& vec) {
  vec.resize(int_set.size()); int idx = 0;
  for(auto f : int_set) {
      vec(idx) = f; idx++;
    }
}
