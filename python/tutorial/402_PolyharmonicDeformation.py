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


z_max = 1.0
z_dir = -0.03
k = 2
resolve = True

V = igl.eigen.MatrixXd()
U = igl.eigen.MatrixXd()

Z = igl.eigen.MatrixXd()

F = igl.eigen.MatrixXi()
b = igl.eigen.MatrixXi()

bc = igl.eigen.MatrixXd()


def pre_draw(viewer):
    global z_max, z_dir, k, resolve, V, U, Z, F, b, bc

    if resolve:
        igl.harmonic(V, F, b, bc, k, Z)
        resolve = False

    U.setCol(2, z_max * Z)
    viewer.data().set_vertices(U)
    viewer.data().compute_normals()

    if viewer.core.is_animating:
        z_max += z_dir
        z_dir *= (-1.0 if z_max >= 1.0 or z_max <= 0.0 else 1.0)

    return False


def key_down(viewer, key, mods):
    global z_max, z_dir, k, resolve, V, U, Z, F, b, bc

    if key == ord(' '):
        viewer.core.is_animating = not viewer.core.is_animating
    elif key == ord('.'):
        k = k + 1
        k = (4 if k > 4 else k)
        resolve = True
    elif key == ord(','):
        k = k - 1
        k = (1 if k < 1 else k)
        resolve = True
    return True


igl.readOBJ(TUTORIAL_SHARED_PATH + "bump-domain.obj", V, F)
U = igl.eigen.MatrixXd(V)

# Find boundary vertices outside annulus

Vrn = V.rowwiseNorm()
is_outer = [Vrn[i] - 1.00 > -1e-15 for i in range(0, V.rows())]
is_inner = [Vrn[i] - 0.15 < 1e-15 for i in range(0, V.rows())]
in_b = [is_outer[i] or is_inner[i] for i in range(0, len(is_outer))]

b = igl.eigen.MatrixXd([[i for i in range(0, V.rows()) if (in_b[i])]]).transpose().castint()

bc.resize(b.size(), 1)

for bi in range(0, b.size()):
    bc[bi] = (0.0 if is_outer[b[bi]] else 1.0)

# Pseudo-color based on selection
C = igl.eigen.MatrixXd(F.rows(), 3)
purple = igl.eigen.MatrixXd([[80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0]])
gold = igl.eigen.MatrixXd([[255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0]])

for f in range(0, F.rows()):
    if in_b[F[f, 0]] and in_b[F[f, 1]] and in_b[F[f, 2]]:
        C.setRow(f, purple)
    else:
        C.setRow(f, gold)

# Plot the mesh with pseudocolors
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(U, F)
viewer.data().show_lines = False
viewer.data().set_colors(C)
viewer.core.trackball_angle = igl.eigen.Quaterniond(0.81,-0.58,-0.03,-0.03)
viewer.core.trackball_angle.normalize()
viewer.callback_pre_draw = pre_draw
viewer.callback_key_down = key_down
viewer.core.is_animating = True
viewer.core.animation_max_fps = 30.0
print("Press [space] to toggle animation.")
print("Press '.' to increase k.")
print("Press ',' to decrease k.")
viewer.launch()
