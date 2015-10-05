# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Load a mesh in OFF format
igl.readOFF("../../tutorial/shared/bunny.off", V, F)

# Find the bounding box
m = V.colwiseMinCoeff()
M = V.colwiseMaxCoeff()


# Corners of the bounding box
V_box = igl.eigen.MatrixXd(
[
[m[0,0], m[0,1], m[0,2]],
[M[0,0], m[0,1], m[0,2]],
[M[0,0], M[0,1], m[0,2]],
[m[0,0], M[0,1], m[0,2]],
[m[0,0], m[0,1], M[0,2]],
[M[0,0], m[0,1], M[0,2]],
[M[0,0], M[0,1], M[0,2]],
[m[0,0], M[0,1], M[0,2]]
]
)

E_box = igl.eigen.MatrixXi(
[
[0, 1],
[1, 2],
[2, 3],
[3, 0],
[4, 5],
[5, 6],
[6, 7],
[7, 4],
[0, 4],
[1, 5],
[2, 6],
[7 ,3]
]
)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)

# Plot the corners of the bounding box as points
viewer.data.add_points(V_box,igl.eigen.MatrixXd([[1,0,0]]))

# Plot the edges of the bounding box
for i in range(0,E_box.rows()):
    viewer.data.add_edges(
        V_box.row(E_box[i,0]),
        V_box.row(E_box[i,1]),
        igl.eigen.MatrixXd([[1,0,0]]))

# Plot labels with the coordinates of bounding box vertices
l1 = 'x: ' + str(m[0,0]) + ' y: ' + str(m[0,1]) + ' z: ' + str(m[0,2])
viewer.data.add_label(m.transpose(),l1)

l2 = 'x: ' + str(M[0,0]) + ' y: ' + str(M[0,1]) + ' z: ' + str(M[0,2])
viewer.data.add_label(M.transpose(),l2);

# Launch the viewer
viewer.launch();
