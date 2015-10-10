# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
C = igl.eigen.MatrixXd()

# Load a mesh in OFF format
igl.readOFF("../../tutorial/shared/screwdriver.off", V, F)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)

# Use the z coordinate as a scalar field over the surface
Z = V.col(2);

# Compute per-vertex colors
igl.jet(Z,True,C)

# Add per-vertex colors
viewer.data.set_colors(C)

# Launch the viewer
viewer.launch()
