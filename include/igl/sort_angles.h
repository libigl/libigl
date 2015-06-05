#ifndef SORT_ANGLES_H
#define SORT_ANGLES_H

#include "igl_inline.h"

namespace igl {
    // Sort angles in ascending order in a numerically robust way.
    //
    // Instead of computing angles using atan2(y, x), sort directly on (y, x).
    //
    // Inputs:
    //   M: m by n matrix of scalars. (n >= 2).  Assuming the first column of M
    //      contains values for y, and the second column is x.  Using the rest
    //      of the columns as tie-breaker.
    //   R: an array of m indices.  M.row(R[i]) contains the i-th smallest
    //      angle.
    template<typename DerivedM, typename DerivedR>
    IGL_INLINE void sort_angles(
            const Eigen::PlainObjectBase<DerivedM>& M,
            Eigen::PlainObjectBase<DerivedR>& R);
}

#ifndef IGL_STATIC_LIBRARY
#include "sort_angles.cpp"
#endif

#endif
