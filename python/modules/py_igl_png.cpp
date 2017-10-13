// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "../python_shared.h"

#include <igl/png/readPNG.h>
#include <igl/png/writePNG.h>


void python_export_igl_png(py::module &me) {

  py::module m = me.def_submodule(
    "png", "Wrappers for libigl functions that use png");

  #include "../py_igl/png/py_readPNG.cpp"
  #include "../py_igl/png/py_writePNG.cpp"

}
