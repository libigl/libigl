import sys, os
import time
# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

import tcpviewer

from shared import TUTORIAL_SHARED_PATH


## This is a test application for the TCPViewer

# Launch the tcp viewer
os.system("python ../tcpviewer.py&")

# Wait for it to set up the socket
time.sleep(1)

# Read a mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "beetle.off", V, F)

# Send it to the viewer
viewer = tcpviewer.TCPViewer()
viewer.data.set_mesh(V, F)
viewer.launch()
