# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import sys, os
import time
# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

import tcpviewer

from shared import TUTORIAL_SHARED_PATH


## This is a test application for the TCPViewer
# Make sure to launch the tcpviewer.py first

# Read a mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "beetle.off", V, F)

# Send it to the viewer
viewer = tcpviewer.TCPViewer()
viewer.data().set_mesh(V, F)
viewer.launch()
