
#ifndef IGL_HALFEDGE_FLIP_H
#define IGL_HALFEDGE_FLIP_H

#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
    /*
     * Get information of halfedge flips of triangle meshes
     *
     *  Templates:
     *      DerivedF derived from face indices matrix type: i.e. MatrixXi
     *      DerivedHE derived from halfedge flip matrix type: i.e. MatrixXi
     *
     *  Inputs:
     *      F #F by 3 list of mesh faces (must be triangles)
     *
     *  Outputs:
     *      HE #F by 3 list of id of flipped halfedge, HE(i, j) = 3 * u + v stands
     *          for flip relation between halfedge (3 * i + j) and (3 * u + v).
     *          The id (3 * i + j) of a halfedge depends on the tirangle id (i) 
     *          it belongs to and corresponding vertex rank (j \in {0, 1, 2})
     */
    template <typename DerivedF, typename DerivedHE>
    IGL_INLINE void halfedge_flip(
            const Eigen::MatrixBase<DerivedF>& F,
            Eigen::PlainObjectBase<DerivedHE>& HE);
}

#ifndef IGL_STATIC_LIBRARY
#include "halfedge_flip.cpp"
#endif


#endif
