import sys, os

# Add the igl library to the modules search path
sys.path.insert(0, os.getcwd() + "/../")
import pyigl as igl

from shared import TUTORIAL_SHARED_PATH, check_dependencies

dependencies = ["viewer"]
check_dependencies(dependencies)


b = igl.eigen.MatrixXi()
B = igl.eigen.MatrixXd()
bc = igl.eigen.MatrixXd()
lx = igl.eigen.MatrixXd()
ux = igl.eigen.MatrixXd()
Beq = igl.eigen.MatrixXd()
Bieq = igl.eigen.MatrixXd()
Z = igl.eigen.MatrixXd()

Q = igl.eigen.SparseMatrixd()
Aeq = igl.eigen.SparseMatrixd()
Aieq = igl.eigen.SparseMatrixd()


def solve(viewer):
    global Q, B, b, bc, Aeq, Beq, Aieq, Bieq, lx, ux, Z
    params = igl.active_set_params()
    params.max_iter = 8

    igl.active_set(Q, B, b, bc, Aeq, Beq, Aieq, Bieq, lx, ux, params, Z)

    C = igl.eigen.MatrixXd()
    igl.jet(Z, 0, 1, C)
    viewer.data.set_colors(C)


def key_down(viewer, key, mod):
    global Beq, solve
    if key == ord('.'):
        Beq[0, 0] = Beq[0, 0] * 2.0
        solve(viewer)
        return True
    elif key == ord(','):
        Beq[0, 0] = Beq[0, 0] / 2.0
        solve(viewer)
        return True
    elif key == ord(' '):
        solve(viewer)
        return True
    return False


V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

igl.readOFF(TUTORIAL_SHARED_PATH + "cheburashka.off", V, F)

# Plot the mesh
viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)
viewer.core.show_lines = False
viewer.callback_key_down = key_down

# One fixed point on belly
b = igl.eigen.MatrixXi([[2556]])
bc = igl.eigen.MatrixXd([[1]])

# Construct Laplacian and mass matrix
L = igl.eigen.SparseMatrixd()
M = igl.eigen.SparseMatrixd()
Minv = igl.eigen.SparseMatrixd()

igl.cotmatrix(V, F, L)
igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_VORONOI, M)
igl.invert_diag(M, Minv)

# Bi-Laplacian
Q = L.transpose() * (Minv * L)

# Zero linear term
B = igl.eigen.MatrixXd.Zero(V.rows(), 1)

# Lower and upper bound
lx = igl.eigen.MatrixXd.Zero(V.rows(), 1)
ux = igl.eigen.MatrixXd.Ones(V.rows(), 1)

# Equality constraint constrain solution to sum to 1
Beq = igl.eigen.MatrixXd([[0.08]])
Aeq = M.diagonal().sparseView().transpose()

# (Empty inequality constraints)
solve(viewer)
print("Press '.' to increase scale and resolve.")
print("Press ',' to decrease scale and resolve.")

viewer.launch()
