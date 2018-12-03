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
import math

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["glfw"]
check_dependencies(dependencies)


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
T = igl.eigen.MatrixXi()
tree = igl.AABB()
FN = igl.eigen.MatrixXd()
VN = igl.eigen.MatrixXd()
EN = igl.eigen.MatrixXd()
E = igl.eigen.MatrixXi()
EMAP = igl.eigen.MatrixXi()

max_distance = 1
slice_z = 0.5
overlay = False

viewer = igl.glfw.Viewer()


def append_mesh(C_vis, F_vis, V_vis, V, F, color):
    F_vis.conservativeResize(F_vis.rows() + F.rows(), 3)
    F_vis.setBottomRows(F.rows(), F + V_vis.rows())
    V_vis.conservativeResize(V_vis.rows() + V.rows(), 3)
    V_vis.setBottomRows(V.rows(), V)
    C_vis.conservativeResize(C_vis.rows() + V.rows(), 3)
    colorM = igl.eigen.MatrixXd(V.rows(), C_vis.cols())
    colorM.rowwiseSet(color)
    C_vis.setBottomRows(V.rows(), colorM)


def update_visualization(viewer):
    global V, F, T, tree, FN, VN, EN, E, EMAP, max_distance, slice_z, overlay
    plane = igl.eigen.MatrixXd([0.0, 0.0, 1.0, -((1 - slice_z) * V.col(2).minCoeff() + slice_z * V.col(2).maxCoeff())])
    V_vis = igl.eigen.MatrixXd()
    F_vis = igl.eigen.MatrixXi()

    # Extract triangle mesh slice through volume mesh and subdivide nasty triangles
    J = igl.eigen.MatrixXi()
    bary = igl.eigen.SparseMatrixd()
    igl.marching_tets(V, T, plane, V_vis, F_vis, J, bary)
    max_l = 0.03
    while True:
        l = igl.eigen.MatrixXd()
        igl.edge_lengths(V_vis, F_vis, l)
        l /= (V_vis.colwiseMaxCoeff() - V_vis.colwiseMinCoeff()).norm()

        if l.maxCoeff() < max_l:
            break

        bad = l.rowwiseMaxCoeff() > max_l
        notbad = l.rowwiseMaxCoeff() <= max_l  # TODO replace by ~ operator
        F_vis_bad = igl.eigen.MatrixXi()
        F_vis_good = igl.eigen.MatrixXi()
        igl.slice_mask(F_vis, bad, 1, F_vis_bad)
        igl.slice_mask(F_vis, notbad, 1, F_vis_good)
        igl.upsample(V_vis, F_vis_bad)
        F_vis = igl.cat(1, F_vis_bad, F_vis_good)

    # Compute signed distance
    S_vis = igl.eigen.MatrixXd()
    I = igl.eigen.MatrixXi()
    N = igl.eigen.MatrixXd()
    C = igl.eigen.MatrixXd()

    # Bunny is a watertight mesh so use pseudonormal for signing
    igl.signed_distance_pseudonormal(V_vis, V, F, tree, FN, VN, EN, EMAP, S_vis, I, C, N)

    # push to [0,1] range
    S_vis = 0.5 * (S_vis / max_distance) + 0.5
    C_vis = igl.eigen.MatrixXd()
    # color without normalizing
    igl.parula(S_vis, False, C_vis)

    if overlay:
        append_mesh(C_vis, F_vis, V_vis, V, F, igl.eigen.MatrixXd([[0.8, 0.8, 0.8]]))

    viewer.data().clear()
    viewer.data().set_mesh(V_vis, F_vis)
    viewer.data().set_colors(C_vis)
    viewer.core.lighting_factor = overlay


def key_down(viewer, key, modifier):
    global slice_z, overlay

    if key == ord(' '):
        overlay = not overlay
    elif key == ord('.'):
        slice_z = min(slice_z + 0.01, 0.99)
    elif key == ord(','):
        slice_z = max(slice_z - 0.01, 0.01)
    else:
        return False

    update_visualization(viewer)
    return True


print("Press [space] to toggle showing surface.")
print("Press '.'/',' to push back/pull forward slicing plane.")

# Load mesh: (V,T) tet-mesh of convex hull, F contains original surface triangles
igl.readMESH(TUTORIAL_SHARED_PATH + "bunny.mesh", V, T, F)

# Call to point_mesh_squared_distance to determine bounds
sqrD = igl.eigen.MatrixXd()
I = igl.eigen.MatrixXi()
C = igl.eigen.MatrixXd()
igl.point_mesh_squared_distance(V, V, F, sqrD, I, C)
max_distance = math.sqrt(sqrD.maxCoeff())

# Precompute signed distance AABB tree
tree.init(V, F)

# Precompute vertex, edge and face normals
igl.per_face_normals(V, F, FN)
igl.per_vertex_normals(V, F, igl.PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, FN, VN)
igl.per_edge_normals(V, F, igl.PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM, FN, EN, E, EMAP)

# Plot the generated mesh
update_visualization(viewer)
viewer.callback_key_down = key_down
viewer.data().show_lines = False
viewer.launch()
