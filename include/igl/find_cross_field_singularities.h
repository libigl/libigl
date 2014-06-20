#ifndef IGL_FIND_CROSS_FIELD_SINGULARITIES_H
#define IGL_FIND_CROSS_FIELD_SINGULARITIES_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  //todo
  // Creates a quad mesh from a triangular mesh and a set of two directions
  // per face, using the algorithm described in the paper
  // "Mixed-Integer Quadrangulation" by D. Bommes, H. Zimmer, L. Kobbelt
  // ACM SIGGRAPH 2009, Article No. 77 (http://dl.acm.org/citation.cfm?id=1531383)
  
  // Inputs:
  //   Vin        #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F          #F by 4 eigen Matrix of face (quad) indices
  //   maxIter    maximum numbers of iterations
  //   threshold  minimum allowed threshold for non-planarity
  // Output:
  //   Vout       #V by 3 eigen Matrix of planar mesh vertex 3D positions
  //
  
  
  template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedO>
  IGL_INLINE void find_cross_field_singularities(const Eigen::PlainObjectBase<DerivedV> &V,
                                                 const Eigen::PlainObjectBase<DerivedF> &F,
                                                 const Eigen::PlainObjectBase<DerivedM> &Handle_MMatch,
                                                 Eigen::PlainObjectBase<DerivedO> &isSingularity,
                                                 Eigen::PlainObjectBase<DerivedO> &singularityIndex);

  
  // TODO: this returns singularity index modulo 4. It may need to be modified to cover indices
  template <typename DerivedV, typename DerivedF, typename DerivedO>
  IGL_INLINE void find_cross_field_singularities(const Eigen::PlainObjectBase<DerivedV> &V,
                                                 const Eigen::PlainObjectBase<DerivedF> &F,
                                                 const Eigen::PlainObjectBase<DerivedV> &PD1,
                                                 const Eigen::PlainObjectBase<DerivedV> &PD2,
                                                 Eigen::PlainObjectBase<DerivedO> &isSingularity,
                                                 Eigen::PlainObjectBase<DerivedO> &singularityIndex);
}
#ifdef IGL_HEADER_ONLY
#include "find_cross_field_singularities.cpp"
#endif

#endif
