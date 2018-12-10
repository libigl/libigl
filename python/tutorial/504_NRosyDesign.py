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
from math import atan2, pi, cos, sin

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["comiso", "glfw"]
check_dependencies(dependencies)


# Mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Constrained faces id
b = igl.eigen.MatrixXi()

# Constrained faces representative vector
bc = igl.eigen.MatrixXd()

# Degree of the N-RoSy field
N = 4


# Converts a representative vector per face in the full set of vectors that describe
# an N-RoSy field
def representative_to_nrosy(V, F, R, N, Y):
    B1 = igl.eigen.MatrixXd()
    B2 = igl.eigen.MatrixXd()
    B3 = igl.eigen.MatrixXd()

    igl.local_basis(V, F, B1, B2, B3)

    Y.resize(F.rows() * N, 3)

    for i in range(0, F.rows()):
        x = R.row(i) * B1.row(i).transpose()
        y = R.row(i) * B2.row(i).transpose()
        angle = atan2(y[0], x[0])

        for j in range(0, N):
            anglej = angle + 2 * pi * j / float(N)
            xj = cos(anglej)
            yj = sin(anglej)
            Y.setRow(i * N + j, xj * B1.row(i) + yj * B2.row(i))


# Plots the mesh with an N-RoSy field and its singularities on top
# The constrained faces (b) are colored in red.
def plot_mesh_nrosy(viewer, V, F, N, PD1, S, b):
    # Clear the mesh
    viewer.data().clear()
    viewer.data().set_mesh(V, F)

    # Expand the representative vectors in the full vector set and plot them as lines
    avg = igl.avg_edge_length(V, F)
    Y = igl.eigen.MatrixXd()
    representative_to_nrosy(V, F, PD1, N, Y)

    B = igl.eigen.MatrixXd()
    igl.barycenter(V, F, B)

    Be = igl.eigen.MatrixXd(B.rows() * N, 3)
    for i in range(0, B.rows()):
        for j in range(0, N):
            Be.setRow(i * N + j, B.row(i))

    viewer.data().add_edges(Be, Be + Y * (avg / 2), igl.eigen.MatrixXd([[0, 0, 1]]))

    # Plot the singularities as colored dots (red for negative, blue for positive)
    for i in range(0, S.size()):
        if S[i] < -0.001:
            viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[1, 0, 0]]))
        elif S[i] > 0.001:
            viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[0, 1, 0]]));

    # Highlight in red the constrained faces
    C = igl.eigen.MatrixXd.Constant(F.rows(), 3, 1)
    for i in range(0, b.size()):
        C.setRow(b[i], igl.eigen.MatrixXd([[1, 0, 0]]))
    viewer.data().set_colors(C)


# It allows to change the degree of the field when a number is pressed
def key_down(viewer, key, modifier):
    global N
    if ord('1') <= key <= ord('9'):
        N = key - ord('0')

    R = igl.eigen.MatrixXd()
    S = igl.eigen.MatrixXd()

    igl.comiso.nrosy(V, F, b, bc, igl.eigen.MatrixXi(), igl.eigen.MatrixXd(), igl.eigen.MatrixXd(), N, 0.5, R, S)
    plot_mesh_nrosy(viewer, V, F, N, R, S, b)

    return False


# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "bumpy.off", V, F)

# Threshold faces with high anisotropy
b = igl.eigen.MatrixXd([[0]]).castint()
bc = igl.eigen.MatrixXd([[1, 1, 1]])

viewer = igl.glfw.Viewer()

# Interpolate the field and plot
key_down(viewer, ord('4'), 0)

# Plot the mesh
viewer.data().set_mesh(V, F)
viewer.callback_key_down = key_down

# Disable wireframe
viewer.data().show_lines = False

# Launch the viewer
viewer.launch()
