#ifndef IGL_IN_ELEMENT_H
#define IGL_IN_ELEMENT_H

#include <igl/igl_inline.h>
#include "InElementAABB.h"
#include <Eigen/Core>

namespace igl
{
  // Determine whether each point in a list of points is in the elements of a
  // mesh.
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions. 
  //   Ele  #Ele by dim+1 list of mesh indices into #V. 
  //   Q  #Q by dim list of query point positions
  //   aabb  axis-aligned bounding box tree object (see InElementAABB.h)
  // Outputs:
  //   I  #Q list of indices into Ele of first containing element (-1 means no
  //     containing element)
  IGL_INLINE void in_element(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele,
    const Eigen::MatrixXd & Q,
    const InElementAABB & aabb,
    Eigen::VectorXi & I);
  // Outputs:
  //   I  #Q by #Ele sparse matrix revealing whether each element contains each
  //     point: I(q,e) means point q is in element e
  IGL_INLINE void in_element(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & Ele,
    const Eigen::MatrixXd & Q,
    const InElementAABB & aabb,
    Eigen::SparseMatrix<double> & I);
  //
  // Example:
  //   InElementAABB aabb;
  //   aabb.init(V,Ele);
};

#ifndef IGL_STATIC_LIBRARY
#include "in_element.cpp"
#endif

#endif
