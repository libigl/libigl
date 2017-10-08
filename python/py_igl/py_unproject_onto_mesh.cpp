// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// COMPLETE BINDINGS ========================

m.def("unproject_onto_mesh", []
(
  const Eigen::MatrixXd & pos,
  const Eigen::MatrixXd & model,
  const Eigen::MatrixXd & proj,
  const Eigen::MatrixXd & viewport,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& fid, // TODO: Can we replace this with integer object reference?
  Eigen::MatrixXd& bc
)
{
  assert_is_Vector2("pos", pos);
  Eigen::Vector2f posv;
  if (pos.size() != 0)
    posv = Eigen::Vector2f(pos.cast<float>());
  assert_is_Matrix4("model", model);
  Eigen::Matrix4f modelm;
  if (model.size() != 0)
    modelm = model.cast<float>();
  assert_is_Matrix4("proj", proj);
  Eigen::Matrix4f projm;
  if (proj.size() != 0)
    projm = proj.cast<float>();
  assert_is_Vector4("viewport", viewport);
  Eigen::Vector4f viewportv;
  if (viewport.size() != 0)
    viewportv = Eigen::Vector4f(viewport.cast<float>());

  Eigen::VectorXd bcv;
  int fidi;
  bool ret = igl::unproject_onto_mesh(posv, modelm, projm, viewportv, V, F, fidi, bcv);
  fid(0, 0) = fidi;
  bc = bcv;
  return ret;
}, __doc_igl_unproject_onto_mesh,
py::arg("pos"), py::arg("model"), py::arg("proj"), py::arg("viewport"), py::arg("V"), py::arg("F"), py::arg("fid"), py::arg("bc"));

// INCOMPLETE BINDINGS ========================

//m.def("unproject_onto_mesh", []
//(
//  Eigen::Vector2f & pos,
//  Eigen::Matrix4f & model,
//  Eigen::Matrix4f & proj,
//  Eigen::Vector4f & viewport,
//  std::function<bool (const Eigen::Vector3f &, const Eigen::Vector3f &, igl::Hit &)> & shoot_ray,
//  int & fid,
//  Eigen::MatrixXd& bc
//)
//{
//  return igl::unproject_onto_mesh(pos, model, proj, viewport, shoot_ray, fid, bc);
//}, __doc_igl_unproject_onto_mesh,
//py::arg("pos"), py::arg("model"), py::arg("proj"), py::arg("viewport"), py::arg("shoot_ray"), py::arg("fid"), py::arg("bc"));

