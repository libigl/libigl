# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import sys, os
import math

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["glfw"]
check_dependencies(dependencies)


V = igl.eigen.MatrixXd()
U = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

L = igl.eigen.SparseMatrixd()
viewer = igl.glfw.Viewer()

# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "cow.off", V, F)

# Compute Laplace-Beltrami operator: #V by #V
igl.cotmatrix(V, F, L)

# Alternative construction of same Laplacian
G = igl.eigen.SparseMatrixd()
K = igl.eigen.SparseMatrixd()

# Gradient/Divergence
igl.grad(V, F, G)

# Diagonal per-triangle "mass matrix"
dblA = igl.eigen.MatrixXd()
igl.doublearea(V, F, dblA)

# Place areas along diagonal #dim times

T = (dblA.replicate(3, 1) * 0.5).asDiagonal() * 1

# Laplacian K built as discrete divergence of gradient or equivalently
# discrete Dirichelet energy Hessian

temp = -G.transpose()
K = -G.transpose() * T * G
print("|K-L|: ", (K - L).norm())


def key_pressed(viewer, key, modifier):
    global V, U, F, L

    if key == ord('r') or key == ord('R'):
        U = V
        print("RESET")

    elif key == ord(' '):

        # Recompute just mass matrix on each step
        M = igl.eigen.SparseMatrixd()

        igl.massmatrix(U, F, igl.MASSMATRIX_TYPE_BARYCENTRIC, M)

        # Solve (M-delta*L) U = M*U
        S = (M - 0.001 * L)

        solver = igl.eigen.SimplicialLLTsparse(S)

        U = solver.solve(M * U)

        # Compute centroid and subtract (also important for numerics)
        dblA = igl.eigen.MatrixXd()
        igl.doublearea(U, F, dblA)

        print(dblA.sum())

        area = 0.5 * dblA.sum()
        BC = igl.eigen.MatrixXd()
        igl.barycenter(U, F, BC)
        centroid = igl.eigen.MatrixXd([[0.0, 0.0, 0.0]])

        for i in range(0, BC.rows()):
            centroid += 0.5 * dblA[i, 0] / area * BC.row(i)

        U -= centroid.replicate(U.rows(), 1)

        # Normalize to unit surface area (important for numerics)
        U = U / math.sqrt(area)
    else:
        return False

    # Send new positions, update normals, recenter
    viewer.data().set_vertices(U)
    viewer.data().compute_normals()
    viewer.core.align_camera_center(U, F)
    return True


# Use original normals as pseudo-colors
N = igl.eigen.MatrixXd()
igl.per_vertex_normals(V, F, N)
C = N.rowwiseNormalized() * 0.5 + 0.5

# Initialize smoothing with base mesh
U = V
viewer.data().set_mesh(U, F)
viewer.data().set_colors(C)
viewer.callback_key_pressed = key_pressed

print("Press [space] to smooth.")
print("Press [r] to reset.")

viewer.launch()
