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

# 100 random indicies into rows of F
I = igl.eigen.MatrixXi()
igl.floor((0.5 * (igl.eigen.MatrixXd.Random(100, 1) + 1.) * F.rows()), I)

# 50 random indicies into rows of I
J = igl.eigen.MatrixXi()
igl.floor((0.5 * (igl.eigen.MatrixXd.Random(50, 1) + 1.) * I.rows()), J)

# K = I(J);
K = igl.eigen.MatrixXi()
igl.slice(I, J, K)

# default green for all faces
# C = p2e(np.array([[0.4,0.8,0.3]])).replicate(F.rows(),1)
C = igl.eigen.MatrixXd([[0.4, 0.8, 0.3]]).replicate(F.rows(), 1)

# Red for each in K
R = igl.eigen.MatrixXd([[1.0, 0.3, 0.3]]).replicate(K.rows(), 1)
# C(K,:) = R
igl.slice_into(R, K, 1, C)

# Plot the mesh with pseudocolors
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.data.set_colors(C)
viewer.launch()
