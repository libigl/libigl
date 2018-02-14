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

bc_frac = 1.0
bc_dir = -0.03
deformation_field = False

V = igl.eigen.MatrixXd()
U = igl.eigen.MatrixXd()
V_bc = igl.eigen.MatrixXd()
U_bc = igl.eigen.MatrixXd()

# Z = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
b = igl.eigen.MatrixXi()


def pre_draw(viewer):
    global bc_frac, bc_dir, deformation_field, V, U, V_bc, U_bc, F, b
    # Determine boundary conditions
    if (viewer.core.is_animating):
        bc_frac += bc_dir
        bc_dir *= (-1.0 if bc_frac >= 1.0 or bc_frac <= 0.0 else 1.0)

    U_bc_anim = V_bc + bc_frac * (U_bc - V_bc)

    if (deformation_field):
        D = igl.eigen.MatrixXd()
        D_bc = U_bc_anim - V_bc
        igl.harmonic(V, F, b, D_bc, 2, D)
        U = V + D
    else:
        igl.harmonic(V, F, b, U_bc_anim, 2, U)

    viewer.data().set_vertices(U)
    viewer.data().compute_normals()
    return False


def key_down(viewer, key, mods):
    global bc_frac, bc_dir, deformation_field, V, U, V_bc, U_bc, F, b

    if key == ord(' '):
        viewer.core.is_animating = not viewer.core.is_animating
        return True
    if key == ord('D') or key == ord('d'):
        deformation_field = not deformation_field
        return True
    return False


igl.readOBJ(TUTORIAL_SHARED_PATH + "decimated-max.obj", V, F)
U = igl.eigen.MatrixXd(V)

# S(i) = j: j<0 (vertex i not in handle), j >= 0 (vertex i in handle j)
S = igl.eigen.MatrixXd()
igl.readDMAT(TUTORIAL_SHARED_PATH + "decimated-max-selection.dmat", S)

S = S.castint()

b = igl.eigen.MatrixXd([[t[0] for t in [(i, S[i]) for i in range(0, V.rows())] if t[1] >= 0]]).transpose().castint()

# Boundary conditions directly on deformed positions
U_bc.resize(b.rows(), V.cols())
V_bc.resize(b.rows(), V.cols())

for bi in range(0, b.rows()):
    V_bc.setRow(bi, V.row(b[bi]))

    if S[b[bi]] == 0:
        # Don't move handle 0
        U_bc.setRow(bi, V.row(b[bi]))
    elif S[b[bi]] == 1:
        # Move handle 1 down
        U_bc.setRow(bi, V.row(b[bi]) + igl.eigen.MatrixXd([[0, -50, 0]]))
    else:
        # Move other handles forward
        U_bc.setRow(bi, V.row(b[bi]) + igl.eigen.MatrixXd([[0, 0, -25]]))

# Pseudo-color based on selection
C = igl.eigen.MatrixXd(F.rows(), 3)
purple = igl.eigen.MatrixXd([[80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0]])
gold = igl.eigen.MatrixXd([[255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0]])

for f in range(0, F.rows()):
    if (S[F[f, 0]]) >= 0 and S[F[f, 1]] >= 0 and S[F[f, 2]] >= 0:
        C.setRow(f, purple)
    else:
        C.setRow(f, gold)

# Plot the mesh with pseudocolors
viewer = igl.glfw.Viewer()
viewer.data().set_mesh(U, F)
viewer.data().show_lines = False
viewer.data().set_colors(C)
# viewer.core.trackball_angle = igl.eigen.Quaterniond(sqrt(2.0),0,sqrt(2.0),0)
# viewer.core.trackball_angle.normalize()

viewer.callback_pre_draw = pre_draw
viewer.callback_key_down = key_down

viewer.core.animation_max_fps = 30.0
print("Press [space] to toggle deformation.")
print("Press 'd' to toggle between biharmonic surface or displacements.")
viewer.launch()
