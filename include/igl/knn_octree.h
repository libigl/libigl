#ifndef IGL_KNN_OCTREE
#define IGL_KNN_OCTREE
#include <Eigen/Core>
namespace igl
{
void knn_octree(const Eigen::MatrixXd & P,
                const int & k,
                const std::vector<std::vector<int> > & point_indices,
                const std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                const std::vector<Eigen::RowVector3d, Eigen::aligned_allocator<Eigen::RowVector3d>> & centers,
                const std::vector<double> & widths,
                Eigen::MatrixXi & I
                );
}
#ifndef IGL_STATIC_LIBRARY
#  include "knn_octree.cpp"
#endif
#endif

