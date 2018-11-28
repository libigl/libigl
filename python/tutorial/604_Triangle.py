#!/usr/bin/env python
#
# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import check_dependencies

dependencies = ["triangle", "glfw"]
check_dependencies(dependencies)


# Input polygon
V = igl.eigen.MatrixXd([[-1, -1], [1, -1], [1, 1], [-1, 1], [-2, -2], [2, -2], [2, 2], [-2, 2]])
E = igl.eigen.MatrixXd([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6,7], [7,4]]).castint()
H = igl.eigen.MatrixXd([[0, 0]])

# Triangulated Interior
V2 = igl.eigen.MatrixXd()
F2 = igl.eigen.MatrixXi()

igl.triangle.triangulate(V, E, H, "a0.005q", V2, F2)

# Plot the mesh
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(V2, F2)
viewer.launch()
