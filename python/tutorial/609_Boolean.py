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

dependencies = ["cgal", "glfw"]
check_dependencies(dependencies)

boolean_type_names = {igl.MESH_BOOLEAN_TYPE_UNION: "Union", igl.MESH_BOOLEAN_TYPE_INTERSECT: "Intersect", igl.MESH_BOOLEAN_TYPE_MINUS: "Minus", igl.MESH_BOOLEAN_TYPE_XOR: "XOR", igl.MESH_BOOLEAN_TYPE_RESOLVE: "Resolve"}

boolean_types = list(boolean_type_names.keys())




def update(viewer):
    print("Calculating A %s B..." % boolean_type_names[boolean_type])
    igl.cgal.mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J)
    C = igl.eigen.MatrixXd(FC.rows(), 3)

    for f in range(C.rows()):
        if J[f] < FA.rows():
            C.setRow(f, Red)
        else:
            C.setRow(f, Green)

    viewer.data().clear()
    viewer.data().set_mesh(VC, FC)
    viewer.data().set_colors(C)
    print("Done.")


def key_down(viewer, key, modifier):
    global boolean_type

    if key == ord('.'):
        boolean_type = boolean_types[(boolean_types.index(boolean_type) + 1) % (len(boolean_types))]
    elif key == ord(','):
        boolean_type = boolean_types[(boolean_types.index(boolean_type) + len(boolean_types) - 1) % len(boolean_types)]
    elif key == ord('['):
        viewer.core.camera_dnear -= 0.1
    elif key == ord(']'):
        viewer.core.camera_dnear += 0.1
    else:
        return False

    update(viewer)

    return False


if __name__ == "__main__":

    VA = igl.eigen.MatrixXd()
    FA = igl.eigen.MatrixXi()
    VB = igl.eigen.MatrixXd()
    FB = igl.eigen.MatrixXi()
    VC = igl.eigen.MatrixXd()
    FC = igl.eigen.MatrixXi()
    J = igl.eigen.MatrixXi()

    Red = igl.eigen.MatrixXd([[1, 0, 0]])
    Green = igl.eigen.MatrixXd([[0, 1, 0]])

    # Load meshes in OFF format
    igl.readOFF(TUTORIAL_SHARED_PATH + "cheburashka.off", VA, FA)
    igl.readOFF(TUTORIAL_SHARED_PATH + "decimated-knight.off", VB, FB)

    boolean_type = igl.MESH_BOOLEAN_TYPE_UNION

    viewer = igl.glfw.Viewer()
    update(viewer)

    print(
        "Usage: Press '.' to switch to next boolean operation type. \nPress ',' to switch to previous boolean operation type. \nPress ']' to push near cutting plane away from camera. \nPress '[' to pull near cutting plane closer to camera. \nHint: investigate _inside_ the model to see orientation changes. \n")

    viewer.data().show_lines = True
    viewer.callback_key_down = key_down
    viewer.core.camera_dnear = 3.9
    viewer.launch()
