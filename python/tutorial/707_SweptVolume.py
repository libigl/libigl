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
from math import pi, cos

sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies, print_usage

dependencies = ["copyleft", "glfw"]
check_dependencies(dependencies)


def key_down(viewer, key, modifier):
    global show_swept_volume, SV, SF, V, F
    if key == ord(' '):
        show_swept_volume = not show_swept_volume
        viewer.data().clear()

        if show_swept_volume:
            viewer.data().set_mesh(SV, SF)
            viewer.data().uniform_colors(igl.eigen.MatrixXd([0.2, 0.2, 0.2]), igl.eigen.MatrixXd([1.0, 1.0, 1.0]), igl.eigen.MatrixXd([1.0, 1.0, 1.0])) # TODO replace with constants from cpp
        else:
            viewer.data().set_mesh(V, F)

        viewer.core.is_animating = not show_swept_volume
        viewer.data().set_face_based(True)

    return True


def pre_draw(viewer):
    global show_swept_volume, V
    if not show_swept_volume:
        T = transform(0.25 * igl.get_seconds())
        VT = V * T.matrix().block(0, 0, 3, 3).transpose()
        trans = T.matrix().block(0, 3, 3, 1).transpose()
        Vtrans = igl.eigen.MatrixXd(VT.rows(), VT.cols())
        Vtrans.rowwiseSet(trans)
        VT += Vtrans
        viewer.data().set_vertices(VT)
        viewer.data().compute_normals()
    return False


# Define a rigid motion
def transform(t):
    T = igl.eigen.Affine3d()
    T.setIdentity()
    T.rotate(t * 2 * pi, igl.eigen.MatrixXd([0, 1, 0]))
    T.translate(igl.eigen.MatrixXd([0, 0.125 * cos(2 * pi * t), 0]))
    return T


if __name__ == "__main__":
    keys = {"space": "toggle between transforming original mesh and swept volume"}
    print_usage(keys)

    V = igl.eigen.MatrixXd()
    SV = igl.eigen.MatrixXd()
    VT = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()
    SF = igl.eigen.MatrixXi()
    show_swept_volume = False
    grid_size = 50
    time_steps = 200
    isolevel = 1

    igl.read_triangle_mesh(TUTORIAL_SHARED_PATH + "bunny.off", V, F)

    print("Computing swept volume...")
    igl.copyleft.swept_volume(V, F, transform, time_steps, grid_size, isolevel, SV, SF)
    print("...finished.")

    # Plot the generated mesh
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(V, F)
    viewer.data().set_face_based(True)
    viewer.core.is_animating = not show_swept_volume
    viewer.callback_pre_draw = pre_draw
    viewer.callback_key_down = key_down
    viewer.launch()
