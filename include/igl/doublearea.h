#ifndef IGL_DOUBLEAREA_H
#define IGL_DOUBLEAREA_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // DOUBLEAREA computes twice the area for each input triangle
  //
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   DeriveddblA  derived type of eigen matrix for dblA (e.g. derived from
  //     MatrixXd)
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by simplex_size list of mesh faces (must be triangles)
  // Outputs:
  //   dblA  #F list of triangle double areas
  template <typename DerivedV, typename DerivedF, typename DeriveddblA>
  IGL_INLINE void doublearea( 
    const Eigen::PlainObjectBase<DerivedV> & V, 
    const Eigen::PlainObjectBase<DerivedF> & F, 
    Eigen::PlainObjectBase<DeriveddblA> & dblA);
  // Same as above but use instrinsic edge lengths rather than (V,F) mesh
  // Templates:
  //   DerivedV  derived type of eigen matrix for V (e.g. derived from
  //     MatrixXd)
  //   DerivedF  derived type of eigen matrix for F (e.g. derived from
  //     MatrixXi)
  //   DeriveddblA  derived type of eigen matrix for dblA (e.g. derived from
  //     MatrixXd)
  // Inputs:
  //   l  #F by dim list of edge lengths using 
  //     for triangles, columns correspond to edges 23,31,12
  // Outputs:
  //   dblA  #F list of triangle double areas
  template <typename Derivedl, typename DeriveddblA>
  IGL_INLINE void doublearea( 
    const Eigen::PlainObjectBase<Derivedl> & l, 
    Eigen::PlainObjectBase<DeriveddblA> & dblA);
}

#ifdef IGL_HEADER_ONLY
#  include "doublearea.h"
#endif

#endif
