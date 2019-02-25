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


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "cheburashka.off", V, F)

# Read scalar function values from a file, U: #V by 1
U = igl.eigen.MatrixXd()
igl.readDMAT(TUTORIAL_SHARED_PATH + "cheburashka-scalar.dmat", U)
U = U.col(0)

# Compute gradient operator: #F*3 by #V
G = igl.eigen.SparseMatrixd()
igl.grad(V, F, G)

# Compute gradient of U
GU = (G * U).MapMatrix(F.rows(), 3)

# Compute gradient magnitude
GU_mag = GU.rowwiseNorm()

viewer = igl.glfw.Viewer()
viewer.data().set_mesh(V, F)

# Compute pseudocolor for original function
C = igl.eigen.MatrixXd()

igl.jet(U, True, C)

# Or for gradient magnitude
# igl.jet(GU_mag,True,C)

viewer.data().set_colors(C)

# Average edge length divided by average gradient (for scaling)
max_size = igl.avg_edge_length(V, F) / GU_mag.mean()

# Draw a black segment in direction of gradient at face barycenters
BC = igl.eigen.MatrixXd()
igl.barycenter(V, F, BC)

black = igl.eigen.MatrixXd([[0.0, 0.0, 0.0]])
viewer.data().add_edges(BC, BC + max_size * GU, black)

# Hide wireframe
viewer.data().show_lines = False

viewer.launch()
