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

from shared import TUTORIAL_SHARED_PATH, check_dependencies, print_usage

dependencies = ["glfw"]
check_dependencies(dependencies)


def append_mesh(C_vis, F_vis, V_vis, V, F, color):
    F_vis.conservativeResize(F_vis.rows() + F.rows(), 3)
    F_vis.setBottomRows(F.rows(), F + V_vis.rows())
    V_vis.conservativeResize(V_vis.rows() + V.rows(), 3)
    V_vis.setBottomRows(V.rows(), V)
    C_vis.conservativeResize(C_vis.rows() + F.rows(), 3)
    colorM = igl.eigen.MatrixXd(F.rows(), C_vis.cols())
    colorM.rowwiseSet(color)
    C_vis.setBottomRows(F.rows(), colorM)


def update(viewer):
    global V, F, T, W, slice_z, overlay
    plane = igl.eigen.MatrixXd([0, 0, 1, -((1 - slice_z) * V.col(2).minCoeff() + slice_z * V.col(2).maxCoeff())])
    V_vis = igl.eigen.MatrixXd()
    F_vis = igl.eigen.MatrixXi()
    J = igl.eigen.MatrixXi()
    bary = igl.eigen.SparseMatrixd()
    igl.marching_tets(V, T, plane, V_vis, F_vis, J, bary)
    W_vis = igl.eigen.MatrixXd()
    igl.slice(W, J, W_vis)
    C_vis = igl.eigen.MatrixXd()
    igl.parula(W_vis, False, C_vis)

    if overlay == 1:  # OVERLAY_INPUT
        append_mesh(C_vis, F_vis, V_vis, V, F, igl.eigen.MatrixXd([[1., 0.894, 0.227]]))
    elif overlay == 2:  # OVERLAY_OUTPUT
        append_mesh(C_vis, F_vis, V_vis, V, F, igl.eigen.MatrixXd([[0.8, 0.8, 0.8]]))

    viewer.data().clear()
    viewer.data().set_mesh(V_vis, F_vis)
    viewer.data().set_colors(C_vis)
    viewer.data().set_face_based(True)


def key_down(viewer, key, modifier):
    global overlay, slice_z

    if key == ord(' '):
        overlay = (overlay + 1) % 3
    elif key == ord('.'):
        slice_z = min(slice_z + 0.01, 0.99)
    elif key == ord(','):
        slice_z = max(slice_z - 0.01, 0.01)

    update(viewer)

    return False


if __name__ == "__main__":
    keys = {"space": "toggle showing input mesh, output mesh or slice through tet-mesh of convex hull",
            ". / ,": "push back/pull forward slicing plane"}

    print_usage(keys)

    V = igl.eigen.MatrixXd()
    BC = igl.eigen.MatrixXd()
    W = igl.eigen.MatrixXd()
    T = igl.eigen.MatrixXi()
    F = igl.eigen.MatrixXi()
    G = igl.eigen.MatrixXi()

    slice_z = 0.5
    overlay = 0

    # Load mesh: (V,T) tet-mesh of convex hull, F contains facets of input
    # surface mesh _after_ self-intersection resolution
    igl.readMESH(TUTORIAL_SHARED_PATH + "big-sigcat.mesh", V, T, F)

    # Compute barycenters of all tets
    igl.barycenter(V, T, BC)

    # Compute generalized winding number at all barycenters
    print("Computing winding number over all %i tets..." % T.rows())
    igl.winding_number(V, F, BC, W)

    # Extract interior tets
    Wt = sum(W > 0.5)
    CT = igl.eigen.MatrixXi(Wt, 4)
    k = 0
    for t in range(T.rows()):
        if W[t] > 0.5:
            CT.setRow(k, T.row(t))
            k += 1

    # find bounary facets of interior tets
    igl.boundary_facets(CT, G)

    # boundary_facets seem to be reversed...
    G = G.rowwiseReverse()

    # normalize
    W = (W - W.minCoeff()) / (W.maxCoeff() - W.minCoeff())

    # Plot the generated mesh
    viewer = igl.glfw.Viewer()
    update(viewer)
    viewer.callback_key_down = key_down
    viewer.launch()
