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
import numpy as np

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl
from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["glfw"]
check_dependencies(dependencies)

# Input mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

data = igl.StreamlineData()
state = igl.StreamlineState()

treat_as_symmetric = True

# animation params
anim_t = 0
anim_t_dir = 1


def representative_to_nrosy(V, F, R, N, Y):
    B1 = igl.eigen.MatrixXd()
    B2 = igl.eigen.MatrixXd()
    B3 = igl.eigen.MatrixXd()

    igl.local_basis(V, F, B1, B2, B3)

    Y.resize(F.rows(), 3 * N)
    for i in range(0, F.rows()):
        x = R.row(i) * B1.row(i).transpose()
        y = R.row(i) * B2.row(i).transpose()
        angle = np.arctan2(y, x)

        for j in range(0, N):
            anglej = angle + np.pi * float(j) / float(N)
            xj = float(np.cos(anglej))
            yj = float(np.sin(anglej))
            Y.setBlock(i, j * 3, 1, 3, xj * B1.row(i) + yj * B2.row(i))


def pre_draw(viewer):
    if not viewer.core.is_animating:
        return False

    global anim_t
    global start_point
    global end_point

    igl.streamlines_next(V, F, data, state)

    value = (anim_t % 100) / 100.0

    if value > 0.5:
        value = 1 - value
    value /= 0.5
    r, g, b = igl.parula(value)
    viewer.data().add_edges(state.start_point, state.end_point, igl.eigen.MatrixXd([[r, g, b]]))

    anim_t += anim_t_dir

    return False


def key_down(viewer, key, modifier):
    if key == ord(' '):
        viewer.core.is_animating = not viewer.core.is_animating
        return True

    return False


def main():
    # Load a mesh in OFF format
    igl.readOFF(TUTORIAL_SHARED_PATH + "bumpy.off", V, F)

    # Create a Vector Field
    temp_field = igl.eigen.MatrixXd()
    b = igl.eigen.MatrixXi([[0]])
    bc = igl.eigen.MatrixXd([[1, 1, 1]])
    S = igl.eigen.MatrixXd()  # unused

    degree = 3
    igl.comiso.nrosy(V, F, b, bc, igl.eigen.MatrixXi(), igl.eigen.MatrixXd(), igl.eigen.MatrixXd(), 1, 0.5, temp_field, S)
    temp_field2 = igl.eigen.MatrixXd()
    representative_to_nrosy(V, F, temp_field, degree, temp_field2)

    # Initialize tracer
    igl.streamlines_init(V, F, temp_field2, treat_as_symmetric, data, state)

    # Setup viewer
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(V, F)
    viewer.callback_pre_draw = pre_draw
    viewer.callback_key_down = key_down

    viewer.core.show_lines = False

    viewer.core.is_animating = False
    viewer.core.animation_max_fps = 30.0

    # Paint mesh grayish
    C = igl.eigen.MatrixXd()
    C.setConstant(viewer.data().V.rows(), 3, .9)
    viewer.data().set_colors(C)

    # Draw vector field on sample points
    state0 = state.copy()

    igl.streamlines_next(V, F, data, state0)
    v = state0.end_point - state0.start_point
    v = v.rowwiseNormalized()

    viewer.data().add_edges(state0.start_point,
                          state0.start_point + 0.059 * v,
                          igl.eigen.MatrixXd([[1.0, 1.0, 1.0]]))

    print("Press [space] to toggle animation")
    viewer.launch()

if __name__ == "__main__":
    main()
