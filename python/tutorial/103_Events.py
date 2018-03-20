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


V1 = igl.eigen.MatrixXd()
F1 = igl.eigen.MatrixXi()

V2 = igl.eigen.MatrixXd()
F2 = igl.eigen.MatrixXi()

def key_pressed(viewer, key, modifier):
    print("Key: ", chr(key))

    if key == ord('1'):
        # # Clear should be called before drawing the mesh
        viewer.data().clear()
        # # Draw_mesh creates or updates the vertices and faces of the displayed mesh.
        # # If a mesh is already displayed, draw_mesh returns an error if the given V and
        # # F have size different than the current ones
        viewer.data().set_mesh(V1, F1)
        viewer.core.align_camera_center(V1,F1)
    elif key == ord('2'):
        viewer.data().clear()
        viewer.data().set_mesh(V2, F2)
        viewer.core.align_camera_center(V2,F2)
    return False


#  Load two meshes
igl.readOFF(TUTORIAL_SHARED_PATH + "bumpy.off", V1, F1)
igl.readOFF(TUTORIAL_SHARED_PATH + "fertility.off", V2, F2)

print("1 Switch to bump mesh")
print("2 Switch to fertility mesh")

viewer = igl.glfw.Viewer()

# Register a keyboard callback that allows to switch between
# the two loaded meshes
viewer.callback_key_pressed = key_pressed
viewer.data().set_mesh(V1, F1)
viewer.launch()
