#ifndef CLOSEST_FACET_H
#define CLOSEST_FACET_H

#include "../igl_inline.h"
#include <Eigen/Core>

namespace igl {
    namespace cgal {

        // Determine the closest facet for each of the input points.
        //
        // Inputs:
        //   V  #V by 3 array of vertices.
        //   F  #F by 3 array of faces.
        //   I  #I list of triangle indices to consider.
        //   P  #P by 3 array of query points.
        //
        // Outputs:
        //   R  #P list of closest facet indices.
        //   S  #P list of bools indicating on which side of the closest facet
        //      each query point lies.
        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedI,
            typename DerivedP,
            typename DerivedR,
            typename DerivedS >
        IGL_INLINE void closest_facet(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedI>& I,
                const Eigen::PlainObjectBase<DerivedP>& P,
                Eigen::PlainObjectBase<DerivedR>& R,
                Eigen::PlainObjectBase<DerivedS>& S);

        template<
            typename DerivedV,
            typename DerivedF,
            typename DerivedP,
            typename DerivedR,
            typename DerivedS >
        IGL_INLINE void closest_facet(
                const Eigen::PlainObjectBase<DerivedV>& V,
                const Eigen::PlainObjectBase<DerivedF>& F,
                const Eigen::PlainObjectBase<DerivedP>& P,
                Eigen::PlainObjectBase<DerivedR>& R,
                Eigen::PlainObjectBase<DerivedS>& S);
    }
}

#ifndef IGL_STATIC_LIBRARY
#include "closest_facet.cpp"
#endif
#endif
