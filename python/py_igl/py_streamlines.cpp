// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::class_<igl::StreamlineData> StreamlineData(m, "StreamlineData");
StreamlineData
.def(py::init<>())
.def_readwrite("TT", &igl::StreamlineData::TT)
.def_readwrite("E", &igl::StreamlineData::E)
.def_readwrite("F2E", &igl::StreamlineData::F2E)
.def_readwrite("E2F", &igl::StreamlineData::E2F)
.def_readwrite("field", &igl::StreamlineData::field)
.def_readwrite("match_ab", &igl::StreamlineData::match_ab)
.def_readwrite("match_ba", &igl::StreamlineData::match_ba)
.def_readwrite("nsample", &igl::StreamlineData::nsample)
.def_readwrite("degree", &igl::StreamlineData::degree)
;

py::class_<igl::StreamlineState> StreamlineState(m, "StreamlineState");
StreamlineState
.def(py::init<>())
.def_readwrite("start_point", &igl::StreamlineState::start_point)
.def_readwrite("end_point", &igl::StreamlineState::end_point)
.def_readwrite("current_face", &igl::StreamlineState::current_face)
.def_readwrite("current_direction", &igl::StreamlineState::current_direction)
.def("copy", [](const igl::StreamlineState &m) { return igl::StreamlineState(m); })
;

m.def("streamlines_init", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& temp_field,
  const bool treat_as_symmetric,
  igl::StreamlineData &data,
  igl::StreamlineState &state,
  double percentage
)
{
  return igl::streamlines_init(V, F, temp_field, treat_as_symmetric, data, state, percentage);

},__doc_igl_streamlines_init,
py::arg("V"), py::arg("F"), py::arg("temp_field"), py::arg("treat_as_symmetric"),
py::arg("data"), py::arg("state"), py::arg("percentage")=0.3);

m.def("streamlines_next", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::StreamlineData &data,
  igl::StreamlineState &state
)
{
  return igl::streamlines_next(V, F, data, state);

},__doc_igl_streamlines_next,
py::arg("V"), py::arg("F"), py::arg("data"), py::arg("state"));
