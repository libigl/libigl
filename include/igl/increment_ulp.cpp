// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "increment_ulp.h"
#include <cmath>
#include <limits>


// see "Robust BVH Ray Traversal" by Thiago Ize, section 3:
// for why we need this
template <typename Derived>
IGL_INLINE void igl::increment_ulp(
    Eigen::MatrixBase<Derived>& inout,
    int it
    )
{
    typedef typename Derived::Scalar Scalar;

    inout = inout.unaryExpr([&it](Scalar v){
              for (int k = 0; k < it; ++k) {
                v = std::nextafter(v, std::signbit(v) ? -std::numeric_limits<Scalar>::infinity(): std::numeric_limits<Scalar>::infinity());
              }
              return v;
            });
}
