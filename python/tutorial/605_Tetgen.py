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

dependencies = ["tetgen", "glfw"]
check_dependencies(dependencies)


# Input polygon
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
B = igl.eigen.MatrixXd()

# Tetrahedralized interior
TV = igl.eigen.MatrixXd()
TT = igl.eigen.MatrixXi()
TF = igl.eigen.MatrixXi()

viewer = igl.glfw.Viewer()


def key_down(viewer, key, modifier):
    if key >= ord('1') and key <= ord('9'):
        t = float((key - ord('1')) + 1) / 9.0
        v = igl.eigen.MatrixXd()
        v = B.col(2) - B.col(2).minCoeff()
        v /= v.col(0).maxCoeff()

        s = []
        for i in range(v.size()):
            if v[i, 0] < t:
                s.append(i)

        V_temp = igl.eigen.MatrixXd(len(s) * 4, 3)
        F_temp = igl.eigen.MatrixXd(len(s) * 4, 3).castint()

        for i in range(len(s)):
            V_temp.setRow(i * 4 + 0, TV.row(TT[s[i], 0]))
            V_temp.setRow(i * 4 + 1, TV.row(TT[s[i], 1]))
            V_temp.setRow(i * 4 + 2, TV.row(TT[s[i], 2]))
            V_temp.setRow(i * 4 + 3, TV.row(TT[s[i], 3]))

            F_temp.setRow(i * 4 + 0, igl.eigen.MatrixXd([[(i*4)+0, (i*4)+1, (i*4)+3]]).castint())
            F_temp.setRow(i * 4 + 1, igl.eigen.MatrixXd([[(i*4)+0, (i*4)+2, (i*4)+1]]).castint())
            F_temp.setRow(i * 4 + 2, igl.eigen.MatrixXd([[(i*4)+3, (i*4)+2, (i*4)+0]]).castint())
            F_temp.setRow(i * 4 + 3, igl.eigen.MatrixXd([[(i*4)+1, (i*4)+2, (i*4)+3]]).castint())

        viewer.data().clear()
        viewer.data().set_mesh(V_temp, F_temp)
        viewer.data().set_face_based(True)

    else:
        return False

    return True


# Load a surface mesh
igl.readOFF(TUTORIAL_SHARED_PATH + "fertility.off", V, F)

# Tetrahedralize the interior
igl.tetgen.tetrahedralize(V, F, "pq1.414Y", TV, TT, TF)

# Compute barycenters
igl.barycenter(TV, TT, B)

# Plot the generated mesh
key_down(viewer, ord('5'), 0)
viewer.callback_key_down = key_down
viewer.launch()
