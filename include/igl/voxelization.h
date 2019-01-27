#ifndef IGL_VOXELIZATION_H
#define IGL_VOXELIZATION_H
#include "igl_inline.h"
#include <Eigen/Dense>

namespace igl
{
  // This function provides voxelization from triangle meshes.

  // It uses igl::winding_number for inside-outside segmentation,
  // so it's robust for all kinds of meshes, and requires oriented
  // triangle information.
  //
  // Inputs:
  //   V        #V by 3 matrix of vertex coordinates
  //   F        #F by 3 matrix of indices of simplex corners into V
  //   maxn maximum length of voxels in all dimensions
  // Outputs:
  //   X        number of voxels in the X dimension
  //   Y        number of voxels in the Y dimension
  //   Z        number of voxels in the Z dimension 
  //   min_cord the left-bottom coordinate of all vertices
  //   voxel    X*Y*Z vector of voxelization (0 or 1)
	template <typename DerivedV, typename DerivedF, typename DerivedMC>
	IGL_INLINE void voxelization(
		const Eigen::MatrixBase<DerivedV> &V,
		const Eigen::MatrixBase<DerivedF> &F,
		const int maxn,
		int &X,
		int &Y,
		int &Z,
		Eigen::PlainObjectBase<DerivedMC> &min_cord,
		double &resolution,
		std::vector<int> &voxel);
}

#ifndef IGL_STATIC_LIBRARY
#  include "voxelization.cpp"
#endif

#endif
