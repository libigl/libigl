// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::enum_<igl::MeshBooleanType>(m, "MeshBooleanType")
    .value("MESH_BOOLEAN_TYPE_UNION", igl::MESH_BOOLEAN_TYPE_UNION)
    .value("MESH_BOOLEAN_TYPE_INTERSECT", igl::MESH_BOOLEAN_TYPE_INTERSECT)
    .value("MESH_BOOLEAN_TYPE_MINUS", igl::MESH_BOOLEAN_TYPE_MINUS)
    .value("MESH_BOOLEAN_TYPE_XOR", igl::MESH_BOOLEAN_TYPE_XOR)
    .value("MESH_BOOLEAN_TYPE_RESOLVE", igl::MESH_BOOLEAN_TYPE_RESOLVE)
    .value("NUM_MESH_BOOLEAN_TYPES", igl::NUM_MESH_BOOLEAN_TYPES)
    .export_values();


