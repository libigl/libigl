// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

py::enum_<igl::MappingEnergyType>(m, "MappingEnergyType")
    .value("MAPPING_ENERGY_TYPE_ARAP", igl::MappingEnergyType::ARAP)
    .value("MAPPING_ENERGY_TYPE_LOG_ARAP", igl::MappingEnergyType::LOG_ARAP)
    .value("MAPPING_ENERGY_TYPE_SYMMETRIC_DIRICHLET", igl::MappingEnergyType::SYMMETRIC_DIRICHLET)
    .value("MAPPING_ENERGY_TYPE_CONFORMAL", igl::MappingEnergyType::CONFORMAL)
    .value("MAPPING_ENERGY_TYPE_EXP_CONFORMAL", igl::MappingEnergyType::EXP_CONFORMAL)
    .value("MAPPING_ENERGY_TYPE_EXP_SYMMETRIC_DIRICHLET", igl::MappingEnergyType::EXP_SYMMETRIC_DIRICHLET)
    .export_values();