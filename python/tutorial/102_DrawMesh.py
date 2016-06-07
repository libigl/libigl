import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "beetle.off", V, F)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.launch()
