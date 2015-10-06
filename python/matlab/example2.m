% Launch the external viewer
launch_viewer;

V = py.igl.eigen.MatrixXd();
F = py.igl.eigen.MatrixXi();
py.igl.read_triangle_mesh('../tutorial/shared/fertility.off', V, F);

% Alternative discrete mean curvature
HN = py.igl.eigen.MatrixXd();
L = py.igl.eigen.SparseMatrixd();
M = py.igl.eigen.SparseMatrixd();
Minv = py.igl.eigen.SparseMatrixd();


py.igl.cotmatrix(V,F,L);
py.igl.massmatrix(V,F,py.igl.MASSMATRIX_TYPE_VORONOI,M);

py.igl.invert_diag(M,Minv);

% Laplace-Beltrami of position
HN = -Minv*(L*V);

% Extract magnitude as mean curvature
H = HN.rowwiseNorm();

% Compute curvature directions via quadric fitting
PD1 = py.igl.eigen.MatrixXd();
PD2 = py.igl.eigen.MatrixXd();

PV1 = py.igl.eigen.MatrixXd();
PV2 = py.igl.eigen.MatrixXd();

py.igl.principal_curvature(V,F,PD1,PD2,PV1,PV2);

% Mean curvature
H = 0.5*(PV1+PV2);

viewer = py.tcpviewer_single.TCPViewer();
viewer.data.set_mesh(V, F);

% Compute pseudocolor
C = py.igl.eigen.MatrixXd();
py.igl.parula(H,true,C);

viewer.data.set_colors(C);

% Average edge length for sizing
avg = py.igl.avg_edge_length(V,F);

% Draw a blue segment parallel to the minimal curvature direction
red  = m2p([0.8,0.2,0.2]);
blue = m2p([0.2,0.2,0.8]);

viewer.data.add_edges(V + PD1*avg, V - PD1*avg, blue);

% Draw a red segment parallel to the maximal curvature direction
viewer.data.add_edges(V + PD2*avg, V - PD2*avg, red);

% Plot
viewer.launch()
