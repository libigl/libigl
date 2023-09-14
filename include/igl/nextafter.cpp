// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "nextafter.h"
#include <cmath>
#include <limits>


// see "Robust BVH Ray Traversal" by Thiago Ize, section 3:
// for why we need this
template <typename Derived>
IGL_INLINE void igl::nextafter(
    Eigen::MatrixBase<Derived>& inout,
    int it
    )
{
    typedef typename Derived::Scalar Scalar;

    for (int i = 0; i < inout.rows(); ++i) {
        for (int j = 0; j < inout.cols(); ++j) {
            for (int k = 0; k < it; ++k) {
                inout(i, j) = std::nextafter(inout(i, j), std::numeric_limits<Scalar>::infinity());
            }
        }
    }
}
