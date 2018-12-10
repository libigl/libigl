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


def pre_draw(viewer):
    global pose, anim_t, C, BE, P, U, M, anim_t_dir

    if viewer.core.is_animating:
        # Interpolate pose and identity
        anim_pose = igl.RotationList(len(pose))

        for e in range(len(pose)):
            anim_pose[e] = pose[e].slerp(anim_t, igl.eigen.Quaterniond.Identity())

        # Propagate relative rotations via FK to retrieve absolute transformations
        vQ = igl.RotationList()
        vT = []
        igl.forward_kinematics(C, BE, P, anim_pose, vQ, vT)
        dim = C.cols()
        T = igl.eigen.MatrixXd(BE.rows() * (dim + 1), dim)
        for e in range(BE.rows()):
            a = igl.eigen.Affine3d.Identity()
            a.translate(vT[e])
            a.rotate(vQ[e])
            T.setBlock(e * (dim + 1), 0, dim + 1, dim, a.matrix().transpose().block(0, 0, dim + 1, dim))

        # Compute deformation via LBS as matrix multiplication
        U = M * T

        # Also deform skeleton edges
        CT = igl.eigen.MatrixXd()
        BET = igl.eigen.MatrixXi()
        igl.deform_skeleton(C, BE, T, CT, BET)

        viewer.data().set_vertices(U)
        viewer.data().set_edges(CT, BET, sea_green)
        viewer.data().compute_normals()
        anim_t += anim_t_dir
        anim_t_dir *= -1.0 if (0.0 >= anim_t or anim_t >= 1.0) else 1.0

    return False


def key_down(viewer, key, mods):
    global selected, W
    if key == ord('.'):
        selected += 1
        selected = min(max(selected, 0), W.cols()-1)
        set_color(viewer)
    elif key == ord(','):
        selected -= 1
        selected = min(max(selected, 0), W.cols()-1)
        set_color(viewer)
    elif key == ord(' '):
        viewer.core.is_animating = not viewer.core.is_animating

    return True


def set_color(viewer):
    global selected, W
    C = igl.eigen.MatrixXd()
    igl.jet(W.col(selected), True, C)
    viewer.data().set_colors(C)


if __name__ == "__main__":
    keys = {".": "show next weight function",
            ",": "show previous weight function",
            "space": "toggle animation"}

    print_usage(keys)

    V = igl.eigen.MatrixXd()
    W = igl.eigen.MatrixXd()
    U = igl.eigen.MatrixXd()
    C = igl.eigen.MatrixXd()
    M = igl.eigen.MatrixXd()
    Q = igl.eigen.MatrixXd()
    T = igl.eigen.MatrixXi()
    F = igl.eigen.MatrixXi()
    BE = igl.eigen.MatrixXi()
    P = igl.eigen.MatrixXi()

    sea_green = igl.eigen.MatrixXd([[70. / 255., 252. / 255., 167. / 255.]])

    selected = 0
    pose = igl.RotationList()
    anim_t = 1.0
    anim_t_dir = -0.03

    igl.readMESH(TUTORIAL_SHARED_PATH + "hand.mesh", V, T, F)
    U = igl.eigen.MatrixXd(V)
    igl.readTGF(TUTORIAL_SHARED_PATH + "hand.tgf", C, BE)

    # Retrieve parents for forward kinematics
    igl.directed_edge_parents(BE, P)

    # Read pose as matrix of quaternions per row
    igl.readDMAT(TUTORIAL_SHARED_PATH + "hand-pose.dmat", Q)
    igl.column_to_quats(Q, pose)
    assert (len(pose) == BE.rows())

    # List of boundary indices (aka fixed value indices into VV)
    b = igl.eigen.MatrixXi()
    # List of boundary conditions of each weight function
    bc = igl.eigen.MatrixXd()

    igl.boundary_conditions(V, T, C, igl.eigen.MatrixXi(), BE, igl.eigen.MatrixXi(), b, bc)

    # compute BBW weights matrix
    bbw_data = igl.BBWData()
    # only a few iterations for sake of demo
    bbw_data.active_set_params.max_iter = 8
    bbw_data.verbosity = 2
    if not igl.bbw(V, T, b, bc, bbw_data, W):
        exit(-1)

    # Normalize weights to sum to one
    igl.normalize_row_sums(W, W)
    # precompute linear blend skinning matrix
    igl.lbs_matrix(V, W, M)

    # Plot the mesh with pseudocolors
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(U, F)
    set_color(viewer)
    viewer.data().set_edges(C, BE, sea_green)
    viewer.data().show_lines = False
    viewer.data().show_overlay_depth = False
    viewer.data().line_width = 1
    viewer.core.trackball_angle.normalize()
    viewer.callback_pre_draw = pre_draw
    viewer.callback_key_down = key_down
    viewer.core.is_animating = False
    viewer.core.animation_max_fps = 30.0
    viewer.launch()
