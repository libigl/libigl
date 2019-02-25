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
from math import pi

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["comiso", "glfw"]
check_dependencies(dependencies)


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Face barycenters
B = igl.eigen.MatrixXd()

# Scale for visualizing the fields
global_scale = 1
extend_arrows = False

# Cross field
X1 = igl.eigen.MatrixXd()
X2 = igl.eigen.MatrixXd()

# Bisector field
BIS1 = igl.eigen.MatrixXd()
BIS2 = igl.eigen.MatrixXd()

# Combed bisector
BIS1_combed = igl.eigen.MatrixXd()
BIS2_combed = igl.eigen.MatrixXd()

# Per-corner, integer mismatches
MMatch = igl.eigen.MatrixXi()

# Field singularities
isSingularity = igl.eigen.MatrixXi()
singularityIndex = igl.eigen.MatrixXi()

# Per corner seams
Seams = igl.eigen.MatrixXi()

# Combed field
X1_combed = igl.eigen.MatrixXd()
X2_combed = igl.eigen.MatrixXd()

# Global parametrization (with seams)
UV_seams = igl.eigen.MatrixXd()
FUV_seams = igl.eigen.MatrixXi()

# Global parametrization
UV = igl.eigen.MatrixXd()
FUV = igl.eigen.MatrixXi()

# Texture
texture_R = igl.eigen.MatrixXuc()
texture_G = igl.eigen.MatrixXuc()
texture_B = igl.eigen.MatrixXuc()


# Create a texture that hides the integer translation in the parametrization
def line_texture():
    size = 128
    size2 = int(size / 2)
    lineWidth = 3
    texture_R.setConstant(size, size, 255)

    for i in range(0, size):
        for j in range(size2 - lineWidth, size2 + lineWidth + 1):
            texture_R[i, j] = 0

    for i in range(size2 - lineWidth, size2 + lineWidth + 1):
        for j in range(0, size):
            texture_R[i, j] = 0

    texture_G = texture_R.copy()
    texture_B = texture_R.copy()
    return (texture_R, texture_G, texture_B)


def key_down(viewer, key, modifier):
    global extend_arrows, texture_R, texture_G, texture_B

    if key == ord('E'):
        extend_arrows = not extend_arrows

    if key < ord('1') or key > ord('8'):
        return False

    viewer.data().clear()
    viewer.data().show_lines = False
    viewer.data().show_texture = False

    if key == ord('1'):
        # Cross field
        viewer.data().set_mesh(V, F)
        viewer.data().add_edges(B - global_scale * X1 if extend_arrows else B, B + global_scale * X1,
                              igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data().add_edges(B - global_scale * X2 if extend_arrows else B, B + global_scale * X2,
                              igl.eigen.MatrixXd([[0, 0, 1]]))

    if key == ord('2'):
        # Bisector field
        viewer.data().set_mesh(V, F)
        viewer.data().add_edges(B - global_scale * BIS1 if extend_arrows else B, B + global_scale * BIS1,
                              igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data().add_edges(B - global_scale * BIS2 if extend_arrows else B, B + global_scale * BIS2,
                              igl.eigen.MatrixXd([[0, 0, 1]]))

    if key == ord('3'):
        # Bisector field combed
        viewer.data().set_mesh(V, F)
        viewer.data().add_edges(B - global_scale * BIS1_combed if extend_arrows else B, B + global_scale * BIS1_combed,
                              igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data().add_edges(B - global_scale * BIS2_combed if extend_arrows else B, B + global_scale * BIS2_combed,
                              igl.eigen.MatrixXd([[0, 0, 1]]))

    if key == ord('4'):
        # Singularities and cuts
        viewer.data().set_mesh(V, F)

        # Plot cuts
        l_count = Seams.sum()
        P1 = igl.eigen.MatrixXd(l_count, 3)
        P2 = igl.eigen.MatrixXd(l_count, 3)

        for i in range(0, Seams.rows()):
            for j in range(0, Seams.cols()):
                if Seams[i, j] != 0:
                    P1.setRow(l_count - 1, V.row(F[i, j]))
                    P2.setRow(l_count - 1, V.row(F[i, (j + 1) % 3]))
                    l_count -= 1

        viewer.data().add_edges(P1, P2, igl.eigen.MatrixXd([[1, 0, 0]]))

        # Plot the singularities as colored dots (red for negative, blue for positive)
        for i in range(0, singularityIndex.size()):
            if 2 > singularityIndex[i] > 0:
                viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[1, 0, 0]]))
            elif singularityIndex[i] > 2:
                viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[1, 0, 0]]))

    if key == ord('5'):
        # Singularities and cuts, original field
        # Singularities and cuts
        viewer.data().set_mesh(V, F)
        viewer.data().add_edges(B - global_scale * X1_combed if extend_arrows else B, B + global_scale * X1_combed,
                              igl.eigen.MatrixXd([[1, 0, 0]]))
        viewer.data().add_edges(B - global_scale * X2_combed if extend_arrows else B, B + global_scale * X2_combed,
                              igl.eigen.MatrixXd([[0, 0, 1]]))

        # Plot cuts
        l_count = Seams.sum()

        P1 = igl.eigen.MatrixXd(l_count, 3)
        P2 = igl.eigen.MatrixXd(l_count, 3)

        for i in range(0, Seams.rows()):
            for j in range(0, Seams.cols()):
                if Seams[i, j] != 0:
                    P1.setRow(l_count - 1, V.row(F[i, j]))
                    P2.setRow(l_count - 1, V.row(F[i, (j + 1) % 3]))
                    l_count -= 1

        viewer.data().add_edges(P1, P2, igl.eigen.MatrixXd([[1, 0, 0]]))

        # Plot the singularities as colored dots (red for negative, blue for positive)
        for i in range(0, singularityIndex.size()):
            if 2 > singularityIndex[i] > 0:
                viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[1, 0, 0]]))
            elif singularityIndex[i] > 2:
                viewer.data().add_points(V.row(i), igl.eigen.MatrixXd([[0, 1, 0]]))

    if key == ord('6'):
        # Global parametrization UV
        viewer.data().set_mesh(UV, FUV)
        viewer.data().set_uv(UV)
        viewer.data().show_lines = True

    if key == ord('7'):
        # Global parametrization in 3D
        viewer.data().set_mesh(V, F)
        viewer.data().set_uv(UV, FUV)
        viewer.data().show_texture = True

    if key == ord('8'):
        # Global parametrization in 3D with seams
        viewer.data().set_mesh(V, F)
        viewer.data().set_uv(UV_seams, FUV_seams)
        viewer.data().show_texture = True

    viewer.data().set_colors(igl.eigen.MatrixXd([[1, 1, 1]]))

    viewer.data().set_texture(texture_R, texture_B, texture_G)

    viewer.core.align_camera_center(viewer.data().V, viewer.data().F)

    return False


# Load a mesh in OFF format
igl.readOFF(TUTORIAL_SHARED_PATH + "3holes.off", V, F)

# Compute face barycenters
igl.barycenter(V, F, B)

# Compute scale for visualizing fields
global_scale = .5 * igl.avg_edge_length(V, F)

# Contrain one face
b = igl.eigen.MatrixXd([[0]]).castint()
bc = igl.eigen.MatrixXd([[1, 0, 0]])

# Create a smooth 4-RoSy field
S = igl.eigen.MatrixXd()

igl.comiso.nrosy(V, F, b, bc, igl.eigen.MatrixXi(), igl.eigen.MatrixXd(), igl.eigen.MatrixXd(), 4, 0.5, X1, S)

# Find the orthogonal vector
B1 = igl.eigen.MatrixXd()
B2 = igl.eigen.MatrixXd()
B3 = igl.eigen.MatrixXd()

igl.local_basis(V, F, B1, B2, B3)

X2 = igl.rotate_vectors(X1, igl.eigen.MatrixXd.Constant(1, 1, pi / 2), B1, B2)

gradient_size = 50
iterations = 0
stiffness = 5.0
direct_round = False

# Always work on the bisectors, it is more general
igl.compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2)

# Comb the field, implicitly defining the seams
igl.comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed)

# Find the integer mismatches
igl.cross_field_missmatch(V, F, BIS1_combed, BIS2_combed, True, MMatch)

# Find the singularities
igl.find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex)

# Cut the mesh, duplicating all vertices on the seams
igl.cut_mesh_from_singularities(V, F, MMatch, Seams)

# Comb the frame-field accordingly
igl.comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed)

# Global parametrization
igl.comiso.miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams, UV, FUV, gradient_size, stiffness,
               direct_round, iterations, 5, True, True)

# Global parametrization (with seams, only for demonstration)
igl.comiso.miq(V, F, X1_combed, X2_combed, MMatch, isSingularity, Seams, UV_seams, FUV_seams, gradient_size,
               stiffness, direct_round, iterations, 5, False)

# Plot the mesh
viewer = igl.glfw.Viewer()

# Replace the standard texture with an integer shift invariant texture
(texture_R, texture_G, texture_B) = line_texture()

# Plot the original mesh with a texture parametrization
key_down(viewer, ord('7'), 0)

# Launch the viewer
viewer.callback_key_down = key_down
viewer.launch()
