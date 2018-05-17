// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MAPPING_ENERGY_WITH_JACOBIANS_H
#define IGL_MAPPING_ENERGY_WITH_JACOBIANS_H

#include "igl_inline.h"
#include <Eigen/Dense>
#include "MappingEnergyType.h"

namespace igl
{
   IGL_INLINE double mapping_energy_with_jacobians(const Eigen::MatrixXd &Ji, 
                                                  const Eigen::VectorXd &areas, 
                                                  igl::MappingEnergyType slim_energy, 
                                                  double exp_factor);
  
}
#ifndef IGL_STATIC_LIBRARY
#  include "mapping_energy_with_jacobians.cpp"
#endif

#endif