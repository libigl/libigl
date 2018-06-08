#ifndef IGL_FAST_WINDING_NUMBER
#define IGL_FAST_WINDING_NUMBER
#include <Eigen/Core>
namespace igl
{
    void fast_winding_number_precompute(const Eigen::MatrixXd & P,
                                        const Eigen::MatrixXd & N,
                                        const Eigen::MatrixXd & A,
                                        const std::vector<std::vector<int> > & point_indices,
                                        const std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                                        Eigen::MatrixXd & CM,
                                        Eigen::VectorXd & R,
                                        Eigen::MatrixXd & EC
                    );
    
    void fast_winding_number(const Eigen::MatrixXd & P,
                             const Eigen::MatrixXd & N,
                             const Eigen::MatrixXd & A,
                             const std::vector<std::vector<int> > & point_indices,
                             const std::vector<Eigen::Matrix<int,8,1>,               Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                             const Eigen::MatrixXd & CM,
                             const Eigen::VectorXd & R,
                             const Eigen::MatrixXd & EC,
                             const Eigen::MatrixXd & Q,
                             Eigen::VectorXd & WN
                             );
}
#ifndef IGL_STATIC_LIBRARY
#  include "fast_winding_number.cpp"
#endif

#endif

