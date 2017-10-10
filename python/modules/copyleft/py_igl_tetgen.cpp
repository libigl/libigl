// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/tetgen/tetrahedralize.h>


void python_export_igl_tetgen(py::module &me) {

  py::module m = me.def_submodule(
    "tetgen", "Wrappers for libigl functions that use tetgen");

  #include "../../py_igl/copyleft/tetgen/py_tetrahedralize.cpp"

}
