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
V_uv = igl.eigen.MatrixXd()
initial_guess = igl.eigen.MatrixXd()

show_uv = False


def key_down(viewer, key, modifier):
    global show_uv, V_uv
    if key == ord('1'):
        show_uv = False
    elif key == ord('2'):
        show_uv = True
    elif key == ord('q'):
        V_uv = initial_guess

    if show_uv:
        viewer.data().set_mesh(V_uv, F)
        viewer.core.align_camera_center(V_uv, F)
    else:
        viewer.data().set_mesh(V, F)
        viewer.core.align_camera_center(V, F)

    viewer.data().compute_normals()
    return False


# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "camelhead.off", V, F)

# Compute the initial solution for ARAP (harmonic parametrization)
bnd = igl.eigen.MatrixXi()
igl.boundary_loop(F, bnd)
bnd_uv = igl.eigen.MatrixXd()
igl.map_vertices_to_circle(V, bnd, bnd_uv)

igl.harmonic(V, F, bnd, bnd_uv, 1, initial_guess)

# Add dynamic regularization to avoid to specify boundary conditions
arap_data = igl.ARAPData()
arap_data.with_dynamics = True
b = igl.eigen.MatrixXi.Zero(0, 0)
bc = igl.eigen.MatrixXd.Zero(0, 0)

# Initialize ARAP
arap_data.max_iter = 100

# 2 means that we're going to *solve* in 2d
igl.arap_precomputation(V, F, 2, b, arap_data)

# Solve arap using the harmonic map as initial guess
V_uv = igl.eigen.MatrixXd(initial_guess)  # important, make a copy of it!

igl.arap_solve(bc, arap_data, V_uv)

# Scale UV to make the texture more clear
V_uv *= 20

# Plot the mesh
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(V, F)
viewer.data().set_uv(V_uv)
viewer.callback_key_down = key_down

# Disable wireframe
viewer.data().show_lines = False

# Draw checkerboard texture
viewer.data().show_texture = True

# Launch the viewer
viewer.launch()
