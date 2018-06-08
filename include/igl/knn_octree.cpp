#include "point_areas_and_normals.h"
#include <cmath>
#include <igl/parallel_for.h>

double distance_to_width_one_cube(Eigen::RowVector3d point){
    return std::sqrt(std::pow(std::max(std::abs(point(0))-1,0.0),2)
                     + std::pow(std::max(std::abs(point(1))-1,0.0),2)
                     + std::pow(std::max(std::abs(point(2))-1,0.0),2));
}

double distance_to_cube(Eigen::RowVector3d point, Eigen::RowVector3d cube_center, double cube_width){
    Eigen::RowVector3d transformed_point = (point-cube_center)/cube_width;
    return cube_width*distance_to_width_one_cube(transformed_point);
}

//    template <typename DerivedP, typename DerivedI, typename DerivedO>
//    IGL_INLINE void point_areas_and_normals(
//                               const Eigen::MatrixBase<DerivedP>& P,
//                               const Eigen::MatrixBase<DerivedI>& I,
//                               const Eigen::MatrixBase<DerivedO>& O,
//                               Eigen::MatrixBase<DerivedA> & A,
//                               Eigen::MatrixBase<DerivedN> & N);

namespace igl {
    void knn_octree(const Eigen::MatrixXd & P,
                    const int & k,
                    const std::vector<std::vector<int> > & point_indices,
                    const std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                    const std::vector<Eigen::RowVector3d, Eigen::aligned_allocator<Eigen::RowVector3d>> & centers,
                    const std::vector<double> & widths,
                    Eigen::MatrixXi & I
                    )
    {
        const int n = P.rows();
        const int real_k = std::min(n,k);
        I.resize(n,real_k);
        
      igl::parallel_for(n,[&](int i)
      {
          int points_found = 0;
          Eigen::RowVector3d point_of_interest = P.row(i);
          
          //To make my priority queue take both points and octree cells, I use the indices 0 to n-1
          //for the n points, and the indices n to n+m-1 for the m octree cells
          
          // Using lambda to compare elements.
          auto cmp = [&point_of_interest, &P, &centers, &widths, &n](int left, int right) {
              double leftdistance, rightdistance;
              if(left < n){ //left is a point index
                  leftdistance = (P.row(left) - point_of_interest).norm();
              } else { //left is an octree cell
                  leftdistance = distance_to_cube(point_of_interest,centers.at(left-n),widths.at(left-n));
              }
              
              if(right < n){ //left is a point index
                  rightdistance = (P.row(right) - point_of_interest).norm();
              } else { //left is an octree cell
                  rightdistance = distance_to_cube(point_of_interest,centers.at(right-n),widths.at(right-n));
              }
              return leftdistance >= rightdistance;
          };
          std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);
          queue.push(n); //This is the 0th octree cell (ie the root)
          while(points_found < real_k){
              int curr_cell_or_point = queue.top();
              queue.pop();
              if(curr_cell_or_point < n){ //current index is for is a point
                  I(i,points_found) = curr_cell_or_point;
                  points_found++;
              } else {
                  int curr_cell = curr_cell_or_point - n;
                  if(children.at(curr_cell)(0) == -1){ //In the case of a leaf
                      if(point_indices.at(curr_cell).size() > 0){ //Assumption: Leaves either have one point, or none
                          queue.push(point_indices.at(curr_cell).at(0)); //push the point (pardon the pun)
                      }
                  } else { //Not a leaf
                      for(int j = 0; j < 8; j++){
                          queue.push(children.at(curr_cell)(j)+n); //+n to adjust for the octree cells
                      }
                  }
              }
          }
//          points_found = 0;
//          std::list<int> l = {n}; //This is the 0th octree cell (ie the root)
//          while(points_found < real_k){
//              int curr_cell_or_point = l.front();
//              l.pop_front();
//              if(curr_cell_or_point < n){ //current index is for is a point
//                  I(i,points_found) = curr_cell_or_point;
//                  points_found++;
//              } else {
//                  int curr_cell = curr_cell_or_point - n;
//                  if(children.at(curr_cell)(0) == -1){ //In the case of a leaf
//                      if(point_indices.at(curr_cell).size() > 0){ //Assumption: Leaves either have one point, or none
//                          l.emplace_back(point_indices.at(curr_cell).at(0)); //push the point (pardon the pun)
//                      }
//                  } else { //Not a leaf
//                      for(int j = 0; j < 8; j++){
//                          l.emplace_back(children.at(curr_cell)(j)+n); //+n to adjust for the octree cells
//                      }
//                  }
//                  l.sort(cmp);
//                  int new_size = std::min<int>(l.size(),real_k-points_found);
//                  l.resize(new_size);
//              }
//          }
//
      },1000);
    }
}




