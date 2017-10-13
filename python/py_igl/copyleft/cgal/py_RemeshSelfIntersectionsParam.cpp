// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::class_<igl::copyleft::cgal::RemeshSelfIntersectionsParam > RemeshSelfIntersectionsParam(m, "RemeshSelfIntersectionsParam");

RemeshSelfIntersectionsParam
.def("__init__", [](igl::copyleft::cgal::RemeshSelfIntersectionsParam &m)
{
    new (&m) igl::copyleft::cgal::RemeshSelfIntersectionsParam();
    m.detect_only = false;
    m.first_only = false;
    m.stitch_all = false;
})
.def_readwrite("detect_only", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::detect_only)
.def_readwrite("first_only", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::first_only)
.def_readwrite("stitch_all", &igl::copyleft::cgal::RemeshSelfIntersectionsParam::stitch_all)
;
