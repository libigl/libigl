#ifndef IS_INSIDE
#define IS_INSIDE

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
    namespace cgal {

        // Determine if mesh (V1, F1) is inside of mesh (V2, F2).
        //
        // Precondition:
        // Both mesh must represent closed, self-intersection free,
        // non-degenerated surfaces that
        // are the boundary of 3D volumes. In addition, (V1, F1) must not
        // intersect with (V2, F2).
        //
        // Inputs:
        //   V1  #V1 by 3 list of vertex position of mesh 1
        //   F1  #F1 by 3 list of triangles indices into V1
        //   I1  #I1 list of indices into F1, faces to consider.
        //   V2  #V2 by 3 list of vertex position of mesh 2
        //   F2  #F2 by 3 list of triangles indices into V2
        //   I2  #I2 list of indices into F2, faces to consider.
        //
        // Outputs:
        //   return true iff (V1, F1) is entirely inside of (V2, F2).
        template<typename DerivedV, typename DerivedF, typename DerivedI>
            IGL_INLINE bool is_inside(
                    const Eigen::PlainObjectBase<DerivedV>& V1,
                    const Eigen::PlainObjectBase<DerivedF>& F1,
                    const Eigen::PlainObjectBase<DerivedI>& I1,
                    const Eigen::PlainObjectBase<DerivedV>& V2,
                    const Eigen::PlainObjectBase<DerivedF>& F2,
                    const Eigen::PlainObjectBase<DerivedI>& I2);

        template<typename DerivedV, typename DerivedF>
            IGL_INLINE bool is_inside(
                    const Eigen::PlainObjectBase<DerivedV>& V1,
                    const Eigen::PlainObjectBase<DerivedF>& F1,
                    const Eigen::PlainObjectBase<DerivedV>& V2,
                    const Eigen::PlainObjectBase<DerivedF>& F2);

        // Determine if queries points are inside of mesh (V, F).
        //
        // Precondition:
        // The input mesh must be a closed, self-intersection free,
        // non-degenerated surface.  Queries points must be either inside or
        // outside of the mesh.
        //
        // Inputs:
        //   V  #V by 3 array of vertex positions.
        //   F  #F by 3 array of triangles.
        //   I  #I list of triangle indices to consider.
        //   P  #P by 3 array of query points.
        //
        // Outputs:
        //   inside  #P list of booleans that is true iff the corresponding
        //           query point is inside of the mesh.
        template<typename DerivedV, typename DerivedF, typename DerivedI,
            typename DerivedP, typename DerivedB>
            IGL_INLINE void is_inside_batch(
                    const Eigen::PlainObjectBase<DerivedV>& V,
                    const Eigen::PlainObjectBase<DerivedF>& F,
                    const Eigen::PlainObjectBase<DerivedI>& I,
                    const Eigen::PlainObjectBase<DerivedP>& P,
                    Eigen::PlainObjectBase<DerivedB>& inside);

        template<typename DerivedV, typename DerivedF, typename DerivedP,
            typename DerivedB>
            IGL_INLINE void is_inside_batch(
                    const Eigen::PlainObjectBase<DerivedV>& V,
                    const Eigen::PlainObjectBase<DerivedF>& F,
                    const Eigen::PlainObjectBase<DerivedP>& P,
                    Eigen::PlainObjectBase<DerivedB>& inside);
    }
}

#ifndef IGL_STATIC_LIBRARY
#include "is_inside.cpp"
#endif
#endif
