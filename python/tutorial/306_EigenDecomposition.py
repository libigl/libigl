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
U = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

c = 0
bbd = 1.0
twod = False

if not igl.read_triangle_mesh(TUTORIAL_SHARED_PATH + "beetle.off", V, F):
    print("failed to load mesh")

twod = V.col(2).minCoeff() == V.col(2).maxCoeff()
bbd = (V.colwiseMaxCoeff() - V.colwiseMinCoeff()).norm()

L = igl.eigen.SparseMatrixd()
M = igl.eigen.SparseMatrixd()

igl.cotmatrix(V, F, L)
L = -L
igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_DEFAULT, M)
k = 5

D = igl.eigen.MatrixXd()
if not igl.eigs(L, M, k + 1, igl.EIGS_TYPE_SM, U, D):
    print("Eigs failed.")

U = (U - U.minCoeff()) / (U.maxCoeff() - U.minCoeff())

viewer = igl.glfw.Viewer()


def key_down(viewer, key, mod):
    global U, c

    if key == ord(' '):
        U = U.rightCols(k)

        # Rescale eigen vectors for visualization
        Z = bbd * 0.5 * U.col(c)
        C = igl.eigen.MatrixXd()
        igl.parula(U.col(c), False, C)
        c = (c + 1) % U.cols()

        if twod:
            V.setcol(2, Z)

        viewer.data().set_mesh(V, F)
        viewer.data().compute_normals()
        viewer.data().set_colors(C)
        return True
    return False


viewer.callback_key_down = key_down
viewer.callback_key_down(viewer, ord(' '), 0)
viewer.data().show_lines = False
viewer.launch()
