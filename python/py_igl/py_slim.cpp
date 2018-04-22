// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

py::enum_<igl::SLIMData::SLIM_ENERGY>(m, "SLIMEnergyType")
    .value("SLIM_ENERGY_TYPE_ARAP", igl::SLIMData::ARAP)
    .value("SLIM_ENERGY_TYPE_LOG_ARAP", igl::SLIMData::LOG_ARAP)
    .value("SLIM_ENERGY_TYPE_SYMMETRIC_DIRICHLET", igl::SLIMData::SYMMETRIC_DIRICHLET)
    .value("SLIM_ENERGY_TYPE_CONFORMAL", igl::SLIMData::CONFORMAL)
    .value("SLIM_ENERGY_TYPE_EXP_CONFORMAL", igl::SLIMData::EXP_CONFORMAL)
    .value("SLIM_ENERGY_TYPE_EXP_SYMMETRIC_DIRICHLET", igl::SLIMData::EXP_SYMMETRIC_DIRICHLET)
    .export_values();

py::class_<igl::SLIMData>(m, "SLIMData");

m.def("slim_precompute", &igl::slim_precompute, __doc_igl_slim_precompute, py::arg("V"), py::arg("F"), py::arg("V_init"), 
py::arg("data"), py::arg("slim_energy"), py::arg("b"), py::arg("bc"), py::arg("soft_p"));

m.def("slim_solve", &igl::slim_solve, __doc_igl_slim_solve, py::arg("data"), py::arg("iter_num"));