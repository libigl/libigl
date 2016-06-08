import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import check_dependencies

dependencies = ["triangle", "viewer"]
check_dependencies(dependencies)


# Input polygon
V = igl.eigen.MatrixXd([[-1, -1], [1, -1], [1, 1], [-1, 1], [-2, -2], [2, -2], [2, 2], [-2, 2]])
E = igl.eigen.MatrixXi([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6,7], [7,4]])
H = igl.eigen.MatrixXd([[0, 0]])

# Triangulated Interior
V2 = igl.eigen.MatrixXd()
F2 = igl.eigen.MatrixXi()

igl.triangle.triangulate(V, E, H, "a0.005q", V2, F2)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V2, F2)
viewer.launch()
