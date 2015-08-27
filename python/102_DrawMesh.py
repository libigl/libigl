import igl

# Load a mesh in OFF format
V = igl.eigen.MatrixXd()
F = igl.eigen.MatrixXi()
igl.readOFF("../tutorial/shared/beetle.off", V, F)

# Plot the mesh
viewer = igl.viewer.Viewer();
viewer.data.set_mesh(V, F);
viewer.launch();
