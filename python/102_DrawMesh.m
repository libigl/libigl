% Load a mesh in OFF format
V = py.igl.eigen.MatrixXd();
F = py.igl.eigen.MatrixXi();
py.igl.readOFF('../tutorial/shared/beetle.off', V, F);

V

% Plot the mesh
viewer = py.igl.viewer.Viewer();
viewer.data.set_mesh(V, F);
viewer.launch();

