// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Qingnan Zhou <qnzhou@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/predicates/predicates.h>
#include <predicates.h>

namespace igl {
namespace predicates {

using REAL = double;

template<typename Vector2D>
IGL_INLINE Orientation orient2d(
        const Eigen::MatrixBase<Vector2D>& pa,
        const Eigen::MatrixBase<Vector2D>& pb,
        const Eigen::MatrixBase<Vector2D>& pc) {

    using Point = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>;
    const Point& a = pa.template cast<REAL>();
    const Point& b = pb.template cast<REAL>();
    const Point& c = pc.template cast<REAL>();

    auto r = ::orient2d(
            const_cast<REAL*>(a.data()),
            const_cast<REAL*>(b.data()),
            const_cast<REAL*>(c.data()));

    if (r > 0) return Orientation::POSITIVE;
    else if (r < 0) return Orientation::NEGATIVE;
    else return Orientation::COLLINEAR;
}

template<typename Vector3D>
IGL_INLINE Orientation orient3d(
        const Eigen::MatrixBase<Vector3D>& pa,
        const Eigen::MatrixBase<Vector3D>& pb,
        const Eigen::MatrixBase<Vector3D>& pc,
        const Eigen::MatrixBase<Vector3D>& pd) {

    using Point = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>;
    const Point& a = pa.template cast<REAL>();
    const Point& b = pb.template cast<REAL>();
    const Point& c = pc.template cast<REAL>();
    const Point& d = pd.template cast<REAL>();

    auto r = ::orient3d(
            const_cast<REAL*>(a.data()),
            const_cast<REAL*>(b.data()),
            const_cast<REAL*>(c.data()),
            const_cast<REAL*>(d.data()));

    if (r > 0) return Orientation::POSITIVE;
    else if (r < 0) return Orientation::NEGATIVE;
    else return Orientation::COPLANAR;
}

template<typename Vector2D>
IGL_INLINE Orientation incircle(
        const Eigen::MatrixBase<Vector2D>& pa,
        const Eigen::MatrixBase<Vector2D>& pb,
        const Eigen::MatrixBase<Vector2D>& pc,
        const Eigen::MatrixBase<Vector2D>& pd) {

    using Point = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>;
    const Point& a = pa.template cast<REAL>();
    const Point& b = pb.template cast<REAL>();
    const Point& c = pc.template cast<REAL>();
    const Point& d = pd.template cast<REAL>();

    auto r = ::incircle(
            const_cast<REAL*>(a.data()),
            const_cast<REAL*>(b.data()),
            const_cast<REAL*>(c.data()),
            const_cast<REAL*>(d.data()));

    if (r > 0) return Orientation::INSIDE;
    else if (r < 0) return Orientation::OUTSIDE;
    else return Orientation::COCIRCULAR;
}

template<typename Vector3D>
IGL_INLINE Orientation insphere(
        const Eigen::MatrixBase<Vector3D>& pa,
        const Eigen::MatrixBase<Vector3D>& pb,
        const Eigen::MatrixBase<Vector3D>& pc,
        const Eigen::MatrixBase<Vector3D>& pd,
        const Eigen::MatrixBase<Vector3D>& pe) {

    using Point = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>;
    const Point& a = pa.template cast<REAL>();
    const Point& b = pb.template cast<REAL>();
    const Point& c = pc.template cast<REAL>();
    const Point& d = pd.template cast<REAL>();
    const Point& e = pe.template cast<REAL>();

    auto r = ::insphere(
            const_cast<REAL*>(a.data()),
            const_cast<REAL*>(b.data()),
            const_cast<REAL*>(c.data()),
            const_cast<REAL*>(d.data()),
            const_cast<REAL*>(e.data()));

    if (r > 0) return Orientation::INSIDE;
    else if (r < 0) return Orientation::OUTSIDE;
    else return Orientation::COSPHERICAL;
}

}
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template
igl::predicates::Orientation
igl::predicates::orient2d<Eigen::Matrix<double, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient2d<Eigen::Matrix<float, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient2d<Eigen::Matrix<double, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient2d<Eigen::Matrix<float, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&
);


template
igl::predicates::Orientation
igl::predicates::orient3d<Eigen::Matrix<double, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient3d<Eigen::Matrix<float, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient3d<Eigen::Matrix<double, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::orient3d<Eigen::Matrix<float, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&
);


template
igl::predicates::Orientation
igl::predicates::incircle<Eigen::Matrix<double, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::incircle<Eigen::Matrix<float, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::incircle<Eigen::Matrix<double, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::incircle<Eigen::Matrix<float, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&
);


template
igl::predicates::Orientation
igl::predicates::insphere<Eigen::Matrix<double, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::insphere<Eigen::Matrix<float, -1, -1, 1, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 1, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::insphere<Eigen::Matrix<double, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&
);
template
igl::predicates::Orientation
igl::predicates::insphere<Eigen::Matrix<float, -1, -1, 0, -1, -1>>
(
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&,
const Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1>>&
);

#endif
