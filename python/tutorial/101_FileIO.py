from __future__ import print_function
# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../../tutorial/shared/cube.off", V, F)

# Print the vertices and faces matrices
print("Vertices: \n", V, sep='')
print("Faces: \n", F, sep='')

# Save the mesh in OBJ format
igl.writeOBJ("cube.obj",V,F)
