import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

igl.readOFF(TUTORIAL_SHARED_PATH + "decimated-knight.off", V, F)

# Sort barycenters lexicographically
BC = igl.eigen.MatrixXd()
sorted_BC = igl.eigen.MatrixXd()

igl.barycenter(V, F, BC)

I = igl.eigen.MatrixXi()
J = igl.eigen.MatrixXi()

# sorted_BC = BC(I,:)
igl.sortrows(BC, True, sorted_BC, I)

# Get sorted "place" from sorted indices
J.resize(I.rows(), 1)
# J(I) = 1:numel(I)

igl.slice_into(igl.coloni(0, I.size() - 1), I, J)

# Pseudo-color based on sorted place
C = igl.eigen.MatrixXd()
igl.jet(J.castdouble(), True, C)

# Plot the mesh with pseudocolors
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.data.set_colors(C)
viewer.launch()
