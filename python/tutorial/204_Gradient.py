# Add the igl library to the modules search path
import sys, os
sys.path.insert(0, os.getcwd() + "/../")

import igl

V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()

# Load a mesh in OFF format
igl.readOFF("../../tutorial/shared/cheburashka.off", V, F)

# Read scalar function values from a file, U: #V by 1
U = igl.eigen.MatrixXd()
igl.readDMAT("../../tutorial/shared/cheburashka-scalar.dmat",U)
U = U.col(0)

# Compute gradient operator: #F*3 by #V
G = igl.eigen.SparseMatrixd()
igl.grad(V,F,G)

# Compute gradient of U
GU = (G*U).MapMatrix(F.rows(),3)

# Compute gradient magnitude
GU_mag = GU.rowwiseNorm()

viewer = igl.viewer.Viewer()
viewer.data.set_mesh(V, F)

# Compute pseudocolor for original function
C = igl.eigen.MatrixXd()

igl.jet(U,True,C)

# Or for gradient magnitude
# igl.jet(GU_mag,True,C)

viewer.data.set_colors(C);

# Average edge length divided by average gradient (for scaling)
max_size = igl.avg_edge_length(V,F) / GU_mag.mean()

# Draw a black segment in direction of gradient at face barycenters
BC = igl.eigen.MatrixXd()
igl.barycenter(V,F,BC)

black = igl.eigen.MatrixXd([[0.0,0.0,0.0]])
viewer.data.add_edges(BC,BC+max_size*GU, black)

# Hide wireframe
viewer.core.show_lines = False

viewer.launch()
