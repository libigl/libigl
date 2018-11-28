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
import math

sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = []
check_dependencies(dependencies)

if __name__ == "__main__":
    V = igl.eigen.MatrixXd()
    F = igl.eigen.MatrixXi()

    # Load meshes in OFF format
    igl.readOBJ(TUTORIAL_SHARED_PATH + "horse_quad.obj", V, F)

    # Count the number of irregular vertices, the border is ignored
    irregular = igl.is_irregular_vertex(V, F)
    vertex_count = V.rows()
    irregular_vertex_count = sum(irregular)
    irregular_ratio = irregular_vertex_count / vertex_count

    print("Irregular vertices: \n%d/%d (%.2f%%)\n" % (
        irregular_vertex_count, vertex_count, irregular_ratio * 100))

    # Compute areas, min, max and standard deviation
    area = igl.eigen.MatrixXd()
    igl.doublearea(V, F, area)
    area /= 2.0

    area_avg = area.mean()
    area_min = area.minCoeff() / area_avg
    area_max = area.maxCoeff() / area_avg
    area_ns = (area - area_avg) / area_avg
    area_sigma = math.sqrt(area_ns.squaredMean())

    print("Areas (Min/Max)/Avg_Area Sigma: \n%.2f/%.2f (%.2f)\n" % (
        area_min, area_max, area_sigma))

    # Compute per face angles, min, max and standard deviation
    angles = igl.eigen.MatrixXd()
    igl.internal_angles(V, F, angles)
    angles = 360.0 * (angles / (2 * math.pi))

    angle_avg = angles.mean()
    angle_min = angles.minCoeff()
    angle_max = angles.maxCoeff()
    angle_ns = angles - angle_avg
    angle_sigma = math.sqrt(angle_ns.squaredMean())

    print("Angles in degrees (Min/Max) Sigma: \n%.2f/%.2f (%.2f)\n" % (
        angle_min, angle_max, angle_sigma))
