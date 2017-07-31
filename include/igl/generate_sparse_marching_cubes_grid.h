#pragma once

#include <Eigen/Core>

namespace igl {

    // generate_sparse_marching_cubes_grid( p0, scalarFunc, eps, CV, CS, CI )
    //
    // Given a point, p0, on an isosurface, construct a shell of epsilon sized cubes surrounding the surface.
    // These cubes can be used as the input to marching cubes.
    //
    // Input:
    //  p0  A 3D point on the isosurface surface defined by scalarFunc(x) = 0
    //  scalarFunc  A scalar function from R^3 to R -- points which map to 0 lie
    //      on the surface, points which are negative lie inside the surface,
    //      and points which are positive lie outside the surface
    //  eps  The edge length of the cubes surrounding the surface
    //
    // Output:
    //   CS  #cube-vertices by 3 list of scalar values at the cube vertices
    //   CV  #cube-vertices by 3 list of cube vertex positions
    //   CI  #number of cubes by 8 list of indexes into CS and CV. Each row represents a cube
    //
    template <typename DerivedS, typename DerivedV, typename DerivedI>
    void generate_sparse_marching_cubes_grid(const Eigen::RowVector3d& p0, std::function<double(const Eigen::RowVector3d&)> scalarFunc, double eps,
                             Eigen::PlainObjectBase<DerivedS>& CS, Eigen::PlainObjectBase<DerivedV>& CV, Eigen::PlainObjectBase<DerivedI>& CI);

}
#ifndef IGL_STATIC_LIBRARY
#    include "generate_sparse_marching_cubes_grid.cpp"
#endif
