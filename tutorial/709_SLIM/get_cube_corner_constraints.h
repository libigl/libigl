
#include <Eigen/Core>
#include <set>
void get_cube_corner_constraints(Eigen::MatrixXd& V_o, Eigen::MatrixXi& F, Eigen::VectorXi& b, Eigen::MatrixXd& bc);
void int_set_to_eigen_vector(const std::set<int>& int_set, Eigen::VectorXi& vec);
