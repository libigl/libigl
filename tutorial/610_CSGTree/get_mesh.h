#include <Eigen/Core>

void get_mesh(
    const Eigen::MatrixXd &VA,
    const Eigen::MatrixXi &FA,
    const Eigen::MatrixXd &VB,
    const Eigen::MatrixXi &FB,
    const Eigen::MatrixXd &VC,
    const Eigen::MatrixXi &FC,
    const Eigen::MatrixXd &VD,
    const Eigen::MatrixXi &FD,
    const Eigen::MatrixXd &VE,
    const Eigen::MatrixXi &FE,
    const int view_id,
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F,
    Eigen::VectorXd &I);

