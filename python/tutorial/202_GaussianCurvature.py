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

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["glfw"]
check_dependencies(dependencies)


# Load mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "bumpy.off", V, F)

# Compute Gaussian curvature
K = igl.eigen.MatrixXd()
igl.gaussian_curvature(V, F, K)

# Compute pseudocolor
C = igl.eigen.MatrixXd()
igl.jet(K, True, C)

# Plot the mesh with pseudocolors
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(V, F)
viewer.data().set_colors(C)
viewer.launch()
