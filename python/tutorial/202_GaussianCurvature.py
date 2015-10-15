# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

# Load mesh
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../../tutorial/shared/bumpy.off",V,F);

# Compute Gaussian curvature
K = igl.eigen.MatrixXd();
igl.gaussian_curvature(V,F,K);

# Compute pseudocolor
C = igl.eigen.MatrixXd();
igl.jet(K,True,C);

# Plot the mesh with pseudocolors
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.data.set_colors(C)
viewer.launch()
