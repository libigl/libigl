import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


# Load mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "bumpy.off", V, F)

# Compute Gaussian curvature
K = igl.eigen.MatrixXd()
igl.gaussian_curvature(V, F, K)

# Compute pseudocolor
C = igl.eigen.MatrixXd()
igl.jet(K, True, C)

# Plot the mesh with pseudocolors
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.data.set_colors(C)
viewer.launch()
