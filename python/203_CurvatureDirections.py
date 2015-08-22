import igl

V = igl.eigen.MatrixXd();
F = igl.eigen.MatrixXi();
igl.read_triangle_mesh("../tutorial/shared/fertility.off", V, F);

# Alternative discrete mean curvature
HN = igl.eigen.MatrixXd();
L = igl.eigen.SparseMatrixd();
M = igl.eigen.SparseMatrixd();
Minv = igl.eigen.SparseMatrixd();


igl.cotmatrix(V,F,L);
igl.massmatrix(V,F,igl.MASSMATRIX_TYPE_VORONOI,M);

igl.invert_diag(M,Minv);

# Laplace-Beltrami of position
temp = igl.eigen.SparseMatrixd();
temp = L*V
HN = -Minv*(L*V);

# Extract magnitude as mean curvature
#   VectorXd H = HN.rowwise().norm();
#
#   // Compute curvature directions via quadric fitting
#   MatrixXd PD1,PD2;
#   VectorXd PV1,PV2;
#   igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
#   // mean curvature
#   H = 0.5*(PV1+PV2);
#
#   igl::viewer::Viewer viewer;
#   viewer.data.set_mesh(V, F);
#
#
#   // Compute pseudocolor
#   MatrixXd C;
#   igl::parula(H,true,C);
#   viewer.data.set_colors(C);
#
#   // Average edge length for sizing
#   const double avg = igl::avg_edge_length(V,F);
#
#   // Draw a blue segment parallel to the minimal curvature direction
#   const RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
#   viewer.data.add_edges(V + PD1*avg, V - PD1*avg, blue);
#
#   // Draw a red segment parallel to the maximal curvature direction
#   viewer.data.add_edges(V + PD2*avg, V - PD2*avg, red);
#
#   // Hide wireframe
#   viewer.core.show_lines = false;
#
#   viewer.launch();
# }
