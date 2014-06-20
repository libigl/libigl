#ifndef IGL_MIXED_INTEGER_QUADRANGULATE_H
#define IGL_MIXED_INTEGER_QUADRANGULATE_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

namespace igl
{
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
  
  template <typename DerivedV, typename DerivedF, typename DerivedU>
  IGL_INLINE void mixed_integer_quadrangulate(const Eigen::PlainObjectBase<DerivedV> &V,
                                              const Eigen::PlainObjectBase<DerivedF> &F,
                                              const Eigen::PlainObjectBase<DerivedV> &PD1,
                                              const Eigen::PlainObjectBase<DerivedV> &PD2,
                                              Eigen::PlainObjectBase<DerivedU> &UV,
                                              Eigen::PlainObjectBase<DerivedF> &FUV,
                                              double GradientSize = 30.0,
                                              double Stiffness = 5.0,
                                              bool DirectRound = false,
                                              int iter = 5,
                                              int localIter = 5, bool DoRound = true,
                                              std::vector<int> roundVertices = std::vector<int>(),
                                              std::vector<std::vector<int> > hardFeatures = std::vector<std::vector<int> >());

  template <typename DerivedV, typename DerivedF, typename DerivedU>
  IGL_INLINE void mixed_integer_quadrangulate(const Eigen::PlainObjectBase<DerivedV> &V,
                                              const Eigen::PlainObjectBase<DerivedF> &F,
                                              const Eigen::PlainObjectBase<DerivedV> &PD1_combed,
                                              const Eigen::PlainObjectBase<DerivedV> &PD2_combed,
                                              const Eigen::PlainObjectBase<DerivedV> &BIS1_combed,
                                              const Eigen::PlainObjectBase<DerivedV> &BIS2_combed,
                                              const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_MMatch,
                                              const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_Singular,
                                              const Eigen::Matrix<int, Eigen::Dynamic, 1> &Handle_SingularDegree,
                                              const Eigen::Matrix<int, Eigen::Dynamic, 3> &Handle_Seams,
                                              Eigen::PlainObjectBase<DerivedU> &UV,
                                              Eigen::PlainObjectBase<DerivedF> &FUV,
                                              double GradientSize = 30.0,
                                              double Stiffness = 5.0,
                                              bool DirectRound = false,
                                              int iter = 5,
                                              int localIter = 5, bool DoRound = true,
                                              std::vector<int> roundVertices = std::vector<int>(),
                                              std::vector<std::vector<int> > hardFeatures = std::vector<std::vector<int> >());
}
#ifdef IGL_HEADER_ONLY
#include "mixed_integer_quadrangulate.cpp"
#endif

#endif
