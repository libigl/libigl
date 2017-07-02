import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH


# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF(TUTORIAL_SHARED_PATH + "cube.off", V, F)

# Print the vertices and faces matrices (commented out to make this file compatible with python 2.x and 3.x)
# print("Vertices: \n", V, sep='')
# print("Faces: \n", F, sep='')

# Save the mesh in OBJ format
igl.writeOBJ("cube.obj",V,F)
