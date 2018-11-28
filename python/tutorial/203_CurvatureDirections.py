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
igl.read_triangle_mesh(TUTORIAL_SHARED_PATH + "fertility.off", V, F)

# Alternative discrete mean curvature
HN = igl.eigen.MatrixXd()
L = igl.eigen.SparseMatrixd()
M = igl.eigen.SparseMatrixd()
Minv = igl.eigen.SparseMatrixd()

igl.cotmatrix(V, F, L)
igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_VORONOI, M)

igl.invert_diag(M, Minv)

# Laplace-Beltrami of position
HN = -Minv * (L * V)

# Extract magnitude as mean curvature
H = HN.rowwiseNorm()

# Compute curvature directions via quadric fitting
PD1 = igl.eigen.MatrixXd()
PD2 = igl.eigen.MatrixXd()

PV1 = igl.eigen.MatrixXd()
PV2 = igl.eigen.MatrixXd()

igl.principal_curvature(V, F, PD1, PD2, PV1, PV2)

# Mean curvature
H = 0.5 * (PV1 + PV2)

viewer = igl.glfw.Viewer()
viewer.data().set_mesh(V, F)

# Compute pseudocolor
C = igl.eigen.MatrixXd()
igl.parula(H, True, C)

viewer.data().set_colors(C)

# Average edge length for sizing
avg = igl.avg_edge_length(V, F)

# Draw a blue segment parallel to the minimal curvature direction
red = igl.eigen.MatrixXd([[0.8, 0.2, 0.2]])
blue = igl.eigen.MatrixXd([[0.2, 0.2, 0.8]])

viewer.data().add_edges(V + PD1 * avg, V - PD1 * avg, blue)

# Draw a red segment parallel to the maximal curvature direction
viewer.data().add_edges(V + PD2 * avg, V - PD2 * avg, red)

# Hide wireframe
viewer.data().show_lines = False

viewer.launch()
