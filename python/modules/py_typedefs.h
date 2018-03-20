// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

typedef std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;
PYBIND11_MAKE_OPAQUE(RotationList)

//typedef std::vector<Eigen::Vector3d> TranslationList;
//PYBIND11_MAKE_OPAQUE(TranslationList);

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>)

