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
from math import sin, cos, pi

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["glfw"]
check_dependencies(dependencies)


sea_green = igl.eigen.MatrixXd([[70. / 255., 252. / 255., 167. / 255.]])

V = igl.eigen.MatrixXd()
U = igl.eigen.MatrixXd()

F = igl.eigen.MatrixXi()

S = igl.eigen.MatrixXd()
b = igl.eigen.MatrixXi()

mid = igl.eigen.MatrixXd()

anim_t = 0.0
anim_t_dir = 0.03
arap_data = igl.ARAPData()


def pre_draw(viewer):
    global anim_t

    bc = igl.eigen.MatrixXd(b.size(), V.cols())
    for i in range(0, b.size()):
        bc.setRow(i, V.row(b[i]))
        if S[b[i]] == 0:
            r = mid[0] * 0.25
            bc[i, 0] += r * sin(0.5 * anim_t * 2. * pi)
            bc[i, 1] = bc[i, 1] - r + r * cos(pi + 0.5 * anim_t * 2. * pi)
        elif S[b[i]] == 1:
            r = mid[1] * 0.15
            bc[i, 1] = bc[i, 1] + r + r * cos(pi + 0.15 * anim_t * 2. * pi)
            bc[i, 2] -= r * sin(0.15 * anim_t * 2. * pi)
        elif S[b[i]] == 2:
            r = mid[1] * 0.15
            bc[i, 2] = bc[i, 2] + r + r * cos(pi + 0.35 * anim_t * 2. * pi)
            bc[i, 0] += r * sin(0.35 * anim_t * 2. * pi)

    igl.arap_solve(bc, arap_data, U)
    viewer.data().set_vertices(U)
    viewer.data().compute_normals()

    if viewer.core.is_animating:
        anim_t += anim_t_dir

    return False


def key_down(viewer, key, mods):
    if key == ord(' '):
        viewer.core.is_animating = not viewer.core.is_animating
        return True
    return False


igl.readOFF(TUTORIAL_SHARED_PATH + "decimated-knight.off", V, F)
U = igl.eigen.MatrixXd(V)
igl.readDMAT(TUTORIAL_SHARED_PATH + "decimated-knight-selection.dmat", S)

# Vertices in selection

b = igl.eigen.MatrixXd([[t[0] for t in [(i, S[i]) for i in range(0, V.rows())] if t[1] >= 0]]).transpose().castint()

# Centroid
mid = 0.5 * (V.colwiseMaxCoeff() + V.colwiseMinCoeff())

# Precomputation
arap_data.max_iter = 100
igl.arap_precomputation(V, F, V.cols(), b, arap_data)

# Set color based on selection
C = igl.eigen.MatrixXd(F.rows(), 3)
purple = igl.eigen.MatrixXd([[80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0]])
gold = igl.eigen.MatrixXd([[255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0]])

for f in range(0, F.rows()):
    if S[F[f, 0]] >= 0 and S[F[f, 1]] >= 0 and S[F[f, 2]] >= 0:
        C.setRow(f, purple)
    else:
        C.setRow(f, gold)

# Plot the mesh with pseudocolors
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(U, F)
viewer.data().set_colors(C)
viewer.callback_pre_draw = pre_draw
viewer.callback_key_down = key_down
viewer.core.is_animating = True
viewer.core.animation_max_fps = 30.
print("Press [space] to toggle animation")
viewer.launch()
