#ifndef IGL_SWEPT_VOLUME_H
#define IGL_SWEPT_VOLUME_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
namespace igl
{
  /// Compute the surface of the swept volume of a solid object with surface
  /// (V,F) mesh under going rigid motion.
  /// 
  /// @param[in] V  #V by 3 list of mesh positions in reference pose
  /// @param[in] F  #F by 3 list of mesh indices into V
  /// @param[in] transform  function handle so that transform(t) returns the rigid
  ///     transformation at time t∈[0,1]
  /// @param[in] steps  number of time steps: steps=3 --> t∈{0,0.5,1}
  /// @param[in] grid_res  number of grid cells on the longest side containing the
  ///     motion (isolevel+1 cells will also be added on each side as padding)
  /// @param[in] isolevel  distance level to be contoured as swept volume
  /// @param[out] SV  #SV by 3 list of mesh positions of the swept surface
  /// @param[out] SF  #SF by 3 list of mesh faces into SV
  IGL_INLINE void swept_volume(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::function<Eigen::Affine3d(const double t)> & transform,
    const size_t steps,
    const size_t grid_res,
    const size_t isolevel,
    Eigen::MatrixXd & SV,
    Eigen::MatrixXi & SF);
  /// \brief Pass stepped transforms directly. Passed to the function above with
  /// piecewise constant interpolation. So that, if steps=3, then [0,0.25) →
  /// transforms[0], [0.25,0.75) → transforms[1], [0.75,1] → transforms[2]. In
  /// other words,
  ///
  /// transform(t) = transforms[std::min( (size_t)(t*(steps-1) + 0.5), steps-1)]
  ///
  /// @param[in] transforms  #steps list of rigid transformations at each time
  /// step
  IGL_INLINE void swept_volume(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::vector<
      Eigen::Affine3d,Eigen::aligned_allocator<Eigen::Affine3d> > & transforms,
    const size_t grid_res,
    const size_t isolevel,
    Eigen::MatrixXd & SV,
    Eigen::MatrixXi & SF);
  
}

#ifndef IGL_STATIC_LIBRARY
#  include "swept_volume.cpp"
#endif

#endif
