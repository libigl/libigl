// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::enum_<igl::ARAPEnergyType>(m, "ARAPEnergyType")
    .value("ARAP_ENERGY_TYPE_SPOKES", igl::ARAP_ENERGY_TYPE_SPOKES)
    .value("ARAP_ENERGY_TYPE_SPOKES_AND_RIMS", igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS)
    .value("ARAP_ENERGY_TYPE_ELEMENTS", igl::ARAP_ENERGY_TYPE_ELEMENTS)
    .value("ARAP_ENERGY_TYPE_DEFAULT", igl::ARAP_ENERGY_TYPE_DEFAULT)
    .value("NUM_ARAP_ENERGY_TYPES", igl::NUM_ARAP_ENERGY_TYPES)
    .export_values();
