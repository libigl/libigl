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

dependencies = ["copyleft", "glfw"]
check_dependencies(dependencies)


def key_down(viewer, key, modifier):
    if key == ord('1'):
        viewer.data().clear()
        viewer.data().set_mesh(V, F)
    elif key == ord('2'):
        viewer.data().clear()
        viewer.data().set_mesh(SV, SF)
    elif key == ord('3'):
        viewer.data().clear()
        viewer.data().set_mesh(BV, BF)

    return True


if __name__ == "__main__":
    keys = {"1": "show original mesh",
            "2": "show marching cubes contour of signed distance",
            "3": "show marching cubes contour of indicator function"}

    print_usage(keys)

    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()

    # Read in inputs as double precision floating point meshes
    igl.read_triangle_mesh(TUTORIAL_SHARED_PATH + "armadillo.obj", V, F)

    # number of vertices on the largest side
    s = 50
    Vmin = V.colwiseMinCoeff()
    Vmax = V.colwiseMaxCoeff()
    h = (Vmax - Vmin).maxCoeff() / s
    res = (s * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).castint()

    def lerp(res, Vmin, Vmax, di, d):
        return Vmin[d] + float(di) / (res[d] - 1) * (Vmax[d] - Vmin[d])

    # create grid
    print("Creating grid...")
    GV = igl.eigen.MatrixXd(res[0] * res[1] * res[2], 3)
    for zi in range(res[2]):
        z = lerp(res, Vmin, Vmax, zi, 2)
        for yi in range(res[1]):
            y = lerp(res, Vmin, Vmax, yi, 1)
            for xi in range(res[0]):
                x = lerp(res, Vmin, Vmax, xi, 0)
                GV.setRow(xi + res[0] * (yi + res[1] * zi), igl.eigen.MatrixXd([[x, y, z]]))

    # compute values
    print("Computing distances...")
    S = igl.eigen.MatrixXd()
    B = igl.eigen.MatrixXd()
    I = igl.eigen.MatrixXi()
    C = igl.eigen.MatrixXd()
    N = igl.eigen.MatrixXd()

    igl.signed_distance(GV, V, F, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N)
    # Convert distances to binary inside-outside data --> aliasing artifacts
    B = S.copy()
    for e in range(B.rows()):
        if B[e] > 0:
            B[e] = 1
        else:
            if B[e] < 0:
                B[e] = -1
            else:
                B[e] = 0

    print("Marching cubes...")
    SV = igl.eigen.MatrixXd()
    BV = igl.eigen.MatrixXd()
    SF = igl.eigen.MatrixXi()
    BF = igl.eigen.MatrixXi()

    igl.copyleft.marching_cubes(S, GV, res[0], res[1], res[2], SV, SF)
    igl.copyleft.marching_cubes(B, GV, res[0], res[1], res[2], BV, BF)

    # Plot the generated mesh
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(SV, SF)
    viewer.callback_key_down = key_down
    viewer.launch()
