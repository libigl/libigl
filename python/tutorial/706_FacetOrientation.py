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

dependencies = ["embree", "glfw"]
check_dependencies(dependencies)


def key_down(viewer, key, modifier):
    global facetwise, is_showing_reoriented, FF
    if key == ord('F') or key == ord('f'):
        facetwise = (facetwise + 1) % 2
    elif key == ord('S') or key == ord('s'):
        scramble_colors()
    elif key == ord(' '):
        is_showing_reoriented = ~is_showing_reoriented

    viewer.data().clear()
    viewer.data().set_mesh(V, FF[facetwise] if is_showing_reoriented else F)
    viewer.data().set_colors(RGBcolors[facetwise])

    return True

def scramble_colors():
    global C, viewer, RGBcolors
    for p in range(2):
        R = igl.eigen.MatrixXi()
        igl.randperm(C[p].maxCoeff() + 1, R)
        C[p] = igl.slice(R, igl.eigen.MatrixXi(C[p]))
        HSV = igl.eigen.MatrixXd(C[p].rows(), 3)
        HSV.setCol(0, 360.0 * C[p].castdouble() / C[p].maxCoeff())
        HSVright = igl.eigen.MatrixXd(HSV.rows(), 2)
        HSVright.setConstant(1.0)
        HSV.setRightCols(2, HSVright)
        igl.hsv_to_rgb(HSV, RGBcolors[p])
    viewer.data().set_colors(RGBcolors[facetwise])



if __name__ == "__main__":
    keys = {"space": "toggle between original and reoriented faces",
            "F,f": "toggle between patchwise and facetwise reorientation",
            "S,s": "scramble colors"}
    print_usage(keys)

    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()
    C = [igl.eigen.MatrixXi(), igl.eigen.MatrixXi()]
    RGBcolors = [igl.eigen.MatrixXd(), igl.eigen.MatrixXd()]
    FF = [igl.eigen.MatrixXi(), igl.eigen.MatrixXi()]
    is_showing_reoriented = False
    facetwise = 0

    igl.read_triangle_mesh(TUTORIAL_SHARED_PATH + "truck.obj", V, F)

    # Compute patches
    for p in range(2):
        I = igl.eigen.MatrixXi()
        igl.embree.reorient_facets_raycast(V, F, F.rows() * 100, 10, p == 1, False, False, I, C[p])
        # apply reorientation
        FF[p].conservativeResize(F.rows(), F.cols())
        for i in range(I.rows()):
            if I[i]:
                FF[p].setRow(i, F.row(i).rowwiseReverse())
            else:
                FF[p].setRow(i, F.row(i))

    # Plot the generated mesh
    viewer = igl.glfw.Viewer()
    viewer.data().set_mesh(V, FF[facetwise] if is_showing_reoriented else F)
    viewer.data().set_face_based(True)
    scramble_colors()
    viewer.callback_key_down = key_down
    viewer.launch()
