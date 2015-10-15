%% Launch the external viewer
launch_viewer;

%% Load a mesh in OFF format
V = py.igl.eigen.MatrixXd();
F = py.igl.eigen.MatrixXi();
py.igl.readOFF('../tutorial/shared/beetle.off', V, F);

%% Scale the x coordinate in matlab
V = p2m(V);
V(:,1) = V(:,1) * 2;
V = m2p(V);

%% Plot the mesh
viewer = py.tcpviewer_single.TCPViewer();
viewer.data.set_mesh(V, F);
viewer.launch();

