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
import math

sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies, print_usage

dependencies = ["glfw"]
check_dependencies(dependencies)


def pre_draw(viewer):
    global recompute, anim_t, poses, C, BE, P, U, M, anim_t_dir

    if recompute:
        # Find pose interval
        begin = int(math.floor(anim_t)) % len(poses)
        end = int(math.floor(anim_t) + 1) % len(poses)
        t = anim_t - math.floor(anim_t)

        # Interpolate pose and identity
        anim_pose = igl.RotationList()
        for e in range(len(poses[begin])):
            anim_pose.append(poses[begin][e].slerp(t, poses[end][e]))

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
        if use_dqs:
            igl.dqs(V, W, vQ, vT, U)
        else:
            U = M * T

        # Also deform skeleton edges
        CT = igl.eigen.MatrixXd()
        BET = igl.eigen.MatrixXi()
        igl.deform_skeleton(C, BE, T, CT, BET)

        viewer.data().set_vertices(U)
        viewer.data().set_edges(CT, BET, sea_green)
        viewer.data().compute_normals()
        if viewer.core.is_animating:
            anim_t += anim_t_dir
        else:
            recompute = False

    return False


def key_down(viewer, key, mods):
    global recompute, use_dqs, animation
    recompute = True
    if key == ord('D') or key == ord('d'):
        use_dqs = not use_dqs
        viewer.core.is_animating = False
        animation = False
        if use_dqs:
            print("Switched to Dual Quaternion Skinning")
        else:
            print("Switched to Linear Blend Skinning")
    elif key == ord(' '):
        if animation:
            viewer.core.is_animating = False
            animation = False
        else:
            viewer.core.is_animating = True
            animation = True
    return False


if __name__ == "__main__":
    keys = {"d": "toggle between LBS and DQS",
            "space": "toggle animation"}

    print_usage(keys)

    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()
    C = igl.eigen.MatrixXd()
    BE = igl.eigen.MatrixXi()
    P = igl.eigen.MatrixXi()
    W = igl.eigen.MatrixXd()
    M = igl.eigen.MatrixXd()

    sea_green = igl.eigen.MatrixXd([[70. / 255., 252. / 255., 167. / 255.]])

    anim_t = 0.0
    anim_t_dir = 0.015
    use_dqs = False
    recompute = True
    animation = False  # Flag needed as there is some synchronization problem with viewer.core.is_animating

    poses = [[]]

    igl.readOBJ(TUTORIAL_SHARED_PATH + "arm.obj", V, F)
    U = igl.eigen.MatrixXd(V)
    igl.readTGF(TUTORIAL_SHARED_PATH + "arm.tgf", C, BE)

    # retrieve parents for forward kinematics
    igl.directed_edge_parents(BE, P)
    rest_pose = igl.RotationList()
    igl.directed_edge_orientations(C, BE, rest_pose)
    poses = [[igl.eigen.Quaterniond.Identity() for i in range(4)] for j in range(4)]

    twist = igl.eigen.Quaterniond(pi, igl.eigen.MatrixXd([1, 0, 0]))
    poses[1][2] = rest_pose[2] * twist * rest_pose[2].conjugate()
    bend = igl.eigen.Quaterniond(-pi * 0.7, igl.eigen.MatrixXd([0, 0, 1]))
    poses[3][2] = rest_pose[2] * bend * rest_pose[2].conjugate()

    igl.readDMAT(TUTORIAL_SHARED_PATH + "arm-weights.dmat", W)
    igl.lbs_matrix(V, W, M)

    # Plot the mesh with pseudocolors
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(U, F)
    viewer.data().set_edges(C, BE, sea_green)
    viewer.data().show_lines = False
    viewer.data().show_overlay_depth = False
    viewer.data().line_width = 1
    viewer.core.trackball_angle.normalize()
    viewer.callback_pre_draw = pre_draw
    viewer.callback_key_down = key_down
    viewer.core.is_animating = False
    viewer.core.camera_zoom = 2.5
    viewer.core.animation_max_fps = 30.0
    viewer.launch()
