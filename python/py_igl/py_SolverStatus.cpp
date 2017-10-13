// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::enum_<igl::SolverStatus>(m, "SolverStatus")
    .value("SOLVER_STATUS_CONVERGED", igl::SOLVER_STATUS_CONVERGED)
    .value("SOLVER_STATUS_MAX_ITER", igl::SOLVER_STATUS_MAX_ITER)
    .value("SOLVER_STATUS_ERROR", igl::SOLVER_STATUS_ERROR)
    .value("NUM_SOLVER_STATUSES", igl::NUM_SOLVER_STATUSES)
    .export_values();
