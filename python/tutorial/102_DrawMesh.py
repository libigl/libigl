# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../../tutorial/shared/beetle.off", V, F)

# Plot the mesh
viewer = igl.viewer.Viewer();
viewer.data.set_mesh(V, F);
viewer.launch();
