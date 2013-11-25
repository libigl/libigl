// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "unique_simplices.h"
#include "sort.h"
#include "unique.h"

IGL_INLINE void igl::unique_simplices(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & FF)
{
  using namespace Eigen;
  using namespace igl;
  // Sort each face
  MatrixXi sortF, unusedI;
  igl::sort(F,2,1,sortF,unusedI);
  // Find unique faces
  VectorXi IA,IC;
  MatrixXi C;
  igl::unique_rows(sortF,C,IA,IC);
  FF.resize(IA.size(),F.cols());
  // Copy into output
  for(int i = 0;i<IA.rows();i++)
  {
    FF.row(i) = F.row(IA(i));
  }
}
