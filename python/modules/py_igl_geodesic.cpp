// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "../python_shared.h"

#include <igl/geodesic/exact_geodesic.h>

void python_export_igl_geodesic(py::module &me) {

  py::module m = me.def_submodule(
    "geodesic", "Wrappers for libigl functions that use geodesic");

  #include "../py_igl/geodesic/py_exact_geodesic.cpp"

}
