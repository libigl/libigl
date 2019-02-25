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

dependencies = ["embree", "glfw"]
check_dependencies(dependencies)


# Mesh + AO values + Normals
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
AO = igl.eigen.MatrixXd()
N = igl.eigen.MatrixXd()


viewer = igl.glfw.Viewer()


def key_down(viewer, key, modifier):

    color = igl.eigen.MatrixXd([[0.9, 0.85, 0.9]])

    if key == ord('1'):
        # Show the mesh without the ambient occlusion factor
        viewer.data().set_colors(color)
    elif key == ord('2'):
        # Show the mesh with the ambient occlusion factor
        C = color.replicate(V.rows(), 1)
        for i in range(C.rows()):
            C.setRow(i, C.row(i) * AO[i, 0])
        viewer.data().set_colors(C)
    elif key == ord('.'):
        viewer.core.lighting_factor += 0.1
    elif key == ord(','):
        viewer.core.lighting_factor -= 0.1
    else:
        return False

    viewer.core.lighting_factor = min(max(viewer.core.lighting_factor, 0.0), 1.0)
    return True


print("Press 1 to turn off Ambient Occlusion\nPress 2 to turn on Ambient Occlusion\nPress . to turn up lighting\nPress , to turn down lighting")

# Load a surface mesh
igl.readOFF(TUTORIAL_SHARED_PATH + "fertility.off", V, F)

# Calculate vertex normals
igl.per_vertex_normals(V, F, N)

# Compute ambient occlusion factor using embree
igl.embree.ambient_occlusion(V, F, V, N, 500, AO)
AO = 1.0 - AO

# Plot the generated mesh
viewer.data().set_mesh(V, F)
key_down(viewer, ord('2'), 0)
viewer.callback_key_down = key_down
viewer.data().show_lines = False
viewer.core.lighting_factor = 0.0
viewer.launch()
