#ifndef IGL_BUILD_OCTREE
#define IGL_BUILD_OCTREE
#include <Eigen/Core>
//namespace igl
//{
    void build_octree(const Eigen::MatrixXd & P,
                    std::vector<std::vector<int> > & point_indices,
                    std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                    std::vector<Eigen::RowVector3d, Eigen::aligned_allocator<Eigen::RowVector3d>> & centers,
                    std::vector<double> & widths
                    );
//}
#ifndef IGL_STATIC_LIBRARY
#  include "build_octree.cpp"
#endif

#endif

