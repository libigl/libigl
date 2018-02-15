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
F = igl.eigen.MatrixXi()

N_vertices = igl.eigen.MatrixXd()
N_faces = igl.eigen.MatrixXd()
N_corners = igl.eigen.MatrixXd()


# This function is called every time a keyboard button is pressed
def key_pressed(viewer, key, modifier):
    if key == ord('1'):
        viewer.data().set_normals(N_faces)
        return True
    elif key == ord('2'):
        viewer.data().set_normals(N_vertices)
        return True
    elif key == ord('3'):
        viewer.data().set_normals(N_corners)
        return True
    return False


# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "fandisk.off", V, F)

# Compute per-face normals
N_faces = igl.eigen.MatrixXd()
igl.per_face_normals(V, F, N_faces)

# Compute per-vertex normals
N_vertices = igl.eigen.MatrixXd()
igl.per_vertex_normals(V, F, igl.PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, N_vertices)

# Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
N_corners = igl.eigen.MatrixXd()
igl.per_corner_normals(V, F, 20, N_corners)

# Plot the mesh
viewer = igl.glfw.Viewer()
viewer.callback_key_pressed = key_pressed
viewer.data().show_lines = False
viewer.data().set_mesh(V, F)
viewer.data().set_normals(N_faces)
print("Press '1' for per-face normals.")
print("Press '2' for per-vertex normals.")
print("Press '3' for per-corner normals.")
viewer.launch()
