## This is a test application for the TCPViewer

# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import os
import time

# Launch the tcp viewer
os.system("python ../tcpviewer.py&")

# Wait for it to set up the socket
time.sleep(1)

import igl
import tcpviewer

# Read a mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF('../../tutorial/shared/beetle.off', V, F)

# Send it to the viewer
viewer = tcpviewer.TCPViewer()
viewer.data.set_mesh(V, F)
viewer.launch()
